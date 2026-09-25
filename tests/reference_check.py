#!/usr/bin/env python3
"""Independent CBCL reference decoder used to validate cuda-demux output.

Decodes an Illumina CBCL run straight from the on-disk format (no shared code
with the C++ tool), assigns barcodes with the same matching rule the tool uses,
then streams the tool's FASTQ output and compares a sample of reads position by
position.  Requires only numpy.

    tests/reference_check.py --run RUN_DIR --samplesheet SS.csv --output OUT_DIR
                             [--stride 997] [--edge 32] [--mismatches 1]

Exit status is non-zero when any sampled read differs or per-sample counts
disagree, so the script doubles as an end-to-end regression test.
"""
import argparse
import glob
import gzip
import os
import re
import struct
import sys
import xml.etree.ElementTree as ET
import zlib

import numpy as np

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def revcomp(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


# --------------------------------------------------------------------------- run structure
def parse_run_info(run_dir):
    """Returns (segments, i5_is_rc) where segments is a list of (kind, ncycles)
    in cycle order with kind in {R1, I1, I2, R2}."""
    root = ET.parse(os.path.join(run_dir, "RunInfo.xml")).getroot()
    reads = root.find("Run").find("Reads")
    segs, seen_r, seen_i, i5_rc = [], 0, 0, False
    for rd in reads.findall("Read"):
        n = int(rd.get("NumCycles"))
        if rd.get("IsIndexedRead") == "Y":
            seen_i += 1
            kind = "I1" if seen_i == 1 else "I2"
            if seen_i == 2 and rd.get("IsReverseComplement") == "Y":
                i5_rc = True
        else:
            seen_r += 1
            kind = "R1" if seen_r == 1 else "R2"
        segs.append((kind, n))
    return segs, i5_rc


def parse_samplesheet(path):
    samples, in_data, cols = [], False, None
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("["):
                sec = line[1:line.index("]")]
                in_data = sec in ("BCLConvert_Data", "Data")
                cols = None
                continue
            if not in_data:
                continue
            f = [x.strip() for x in line.split(",")]
            if cols is None:
                low = [x.lower() for x in f]
                cols = {}
                for i, h in enumerate(low):
                    if h in ("sample_id", "sampleid", "sample"):
                        cols["id"] = i
                    elif h in ("index", "index1", "i7_index_id"):
                        cols["i1"] = i
                    elif h in ("index2", "i5_index_id"):
                        cols["i2"] = i
                    elif h == "lane":
                        cols["lane"] = i
                continue
            sid = f[cols["id"]] if cols.get("id", 99) < len(f) else ""
            if not sid:
                continue
            i1 = f[cols["i1"]] if "i1" in cols and cols["i1"] < len(f) else ""
            i2 = f[cols["i2"]] if "i2" in cols and cols["i2"] < len(f) else ""
            lane = 0
            if "lane" in cols and cols["lane"] < len(f) and f[cols["lane"]]:
                lane = int(f[cols["lane"]])
            samples.append((sid, i1.upper(), i2.upper(), lane))
    return samples


# --------------------------------------------------------------------------- filters / tiles
def parse_filters(lane_dir, lane_no):
    """Returns ordered list of (tile_id, pf_mask) sorted by tile id."""
    out = []
    for path in glob.glob(os.path.join(lane_dir, f"s_{lane_no}_*.filter")):
        tile = int(re.search(r"_(\d+)\.filter$", path).group(1))
        with open(path, "rb") as fh:
            data = fh.read()
        version = struct.unpack_from("<I", data, 0)[0]
        if version == 0:
            n = struct.unpack_from("<I", data, 8)[0]
            body = data[12:12 + n]
        else:
            n = struct.unpack_from("<I", data, 4)[0]
            body = data[8:8 + n]
        mask = (np.frombuffer(body, dtype=np.uint8) & 1).astype(bool)
        out.append((tile, mask))
    out.sort(key=lambda t: t[0])
    return out


class CbclFile:
    def __init__(self, path):
        self.path = path
        with open(path, "rb") as fh:
            head = fh.read(1 << 20)
        ver, hs, bpb, bpq, nbins = struct.unpack_from("<HIBBI", head, 0)
        off = 12
        self.bins = [struct.unpack_from("<II", head, off + 8 * i)[1] for i in range(nbins)]
        off += 8 * nbins
        nt = struct.unpack_from("<I", head, off)[0]
        off += 4
        self.tiles = []
        for i in range(nt):
            self.tiles.append(struct.unpack_from("<IIII", head, off + 16 * i))
        off += 16 * nt
        self.nonpf_excluded = head[off] != 0
        assert off + 1 == hs, f"{path}: header size mismatch {off + 1} vs {hs}"
        assert bpb == 2 and bpq == 2, f"{path}: unsupported bit widths {bpb}/{bpq}"
        self.header_size = hs
        self.qmap = np.zeros(4, dtype=np.uint8)
        for i, q in enumerate(self.bins[:4]):
            self.qmap[i] = q

    def blocks(self):
        """Yields (tile_id, num_clusters, base_idx[uint8], qbin[uint8])."""
        with open(self.path, "rb") as fh:
            fh.seek(self.header_size)
            for tile_id, n, unc, comp in self.tiles:
                raw = zlib.decompress(fh.read(comp), 16 + zlib.MAX_WBITS)
                assert len(raw) == unc, f"{self.path} tile {tile_id}: bad block size"
                b = np.frombuffer(raw, dtype=np.uint8)
                nib = np.empty(len(b) * 2, dtype=np.uint8)
                nib[0::2] = b & 0x0F
                nib[1::2] = b >> 4
                nib = nib[:n]
                yield tile_id, n, nib & 3, nib >> 2


def encode_read(base_idx, qbin, qmap):
    """Applies the tool's convention: q-bin 0 is a no-call -> 'N' with Q2 ('#')."""
    seq = BASES[base_idx].copy()
    qual = qmap[qbin] + 33
    nocall = qbin == 0
    seq[nocall] = ord("N")
    qual[nocall] = ord("#")
    return seq, qual


# --------------------------------------------------------------------------- decode lane
def decode_lane(run_dir, lane_dir, lane_no, segs, sample_idx):
    """Decodes every cycle for the clusters in `sample_idx` (compact PF indices)
    and every index cycle for all PF clusters.

    Returns dict with 'seq', 'qual' [len(sample_idx) x total_cycles] and
    'idx_seq' [n_pf x n_index_cycles] (uint8 ASCII)."""
    filters = parse_filters(lane_dir, lane_no)
    tile_info = {}
    raw_off = pf_off = 0
    for tile, mask in filters:
        tile_info[tile] = (raw_off, pf_off, mask)
        raw_off += len(mask)
        pf_off += int(mask.sum())
    n_pf = pf_off
    total_cycles = sum(n for _, n in segs)

    kinds = []
    for kind, n in segs:
        kinds += [kind] * n
    index_cycles = [c for c, k in enumerate(kinds) if k in ("I1", "I2")]
    index_pos = {c: i for i, c in enumerate(index_cycles)}

    seq = np.full((len(sample_idx), total_cycles), ord("N"), dtype=np.uint8)
    qual = np.full((len(sample_idx), total_cycles), ord("#"), dtype=np.uint8)
    idx_seq = np.full((n_pf, len(index_cycles)), ord("N"), dtype=np.uint8)

    # For the sampled clusters we need a scatter from compact index -> row.
    row_of = np.full(n_pf, -1, dtype=np.int64)
    row_of[sample_idx] = np.arange(len(sample_idx))

    for c in range(total_cycles):
        cdir = os.path.join(lane_dir, f"C{c + 1}.1")
        files = sorted(glob.glob(os.path.join(cdir, "*.cbcl")))
        assert files, f"no cbcl in {cdir}"
        for path in files:
            cb = CbclFile(path)
            for tile_id, n, bidx, qbin in cb.blocks():
                r0, p0, mask = tile_info[tile_id]
                if cb.nonpf_excluded:
                    assert n == int(mask.sum()), f"{path} tile {tile_id}: PF count mismatch"
                    compact = np.arange(p0, p0 + n)
                else:
                    assert n == len(mask), f"{path} tile {tile_id}: raw count mismatch"
                    bidx, qbin = bidx[mask], qbin[mask]
                    compact = np.arange(p0, p0 + int(mask.sum()))
                s, q = encode_read(bidx, qbin, cb.qmap)
                if c in index_pos:
                    idx_seq[compact, index_pos[c]] = s
                rows = row_of[compact]
                hit = rows >= 0
                seq[rows[hit], c] = s[hit]
                qual[rows[hit], c] = q[hit]
        if (c + 1) % 50 == 0 or c + 1 == total_cycles:
            print(f"  decoded cycle {c + 1}/{total_cycles}", file=sys.stderr)
    return dict(seq=seq, qual=qual, idx_seq=idx_seq, n_pf=n_pf, kinds=kinds)


# --------------------------------------------------------------------------- barcode assignment
def assign_samples(idx_seq, samples, lane_no, i5_rc, max_mm, i1_len):
    """Returns int32 array of sample index per cluster (-1 = undetermined).

    Rule (same as the tool / bcl-convert): the candidate with the fewest total
    mismatches wins if it is unique and has at most `max_mm` mismatches in each
    index separately. 'N' on either side counts as a mismatch."""
    barcodes, owner = [], []
    seen = {}
    for si, (sid, i1, i2, lane) in enumerate(samples):
        if lane not in (0, lane_no):
            continue
        bc = i1 + (revcomp(i2) if i5_rc else i2)
        if bc in seen:
            continue
        seen[bc] = len(barcodes)
        barcodes.append(bc)
        owner.append(si)
    L = idx_seq.shape[1]
    assert all(len(b) == L for b in barcodes), "barcode length != index cycles"
    bc_arr = np.frombuffer("".join(barcodes).encode(), dtype=np.uint8).reshape(len(barcodes), L)
    n = idx_seq.shape[0]
    best = np.full(n, 1000, dtype=np.int16)
    second = np.full(n, 1000, dtype=np.int16)
    best_b = np.full(n, -1, dtype=np.int32)
    best_ok = np.zeros(n, dtype=bool)
    read_n = idx_seq == ord("N")
    for b in range(len(barcodes)):
        mmv = (idx_seq != bc_arr[b]) | read_n | (bc_arr[b] == ord("N"))
        mm1 = mmv[:, :i1_len].sum(axis=1).astype(np.int16)
        mm2 = mmv[:, i1_len:].sum(axis=1).astype(np.int16)
        mm = mm1 + mm2
        better = mm < best
        second = np.where(better, best, np.minimum(second, mm))
        best = np.where(better, mm, best)
        best_b = np.where(better, b, best_b)
        best_ok = np.where(better, (mm1 <= max_mm) & (mm2 <= max_mm), best_ok)
    ok = best_ok & ((second - best) >= 1)
    out = np.full(n, -1, dtype=np.int32)
    out[ok] = np.asarray(owner, dtype=np.int32)[best_b[ok]]
    return out


# --------------------------------------------------------------------------- fastq streaming
def open_any(path):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")


def check_file(path, wanted, label):
    """wanted: dict ordinal -> (seq bytes, qual bytes). Returns (#checked, #bad, #records)."""
    if not wanted and not os.path.exists(path):
        return 0, 0, 0
    checked = bad = n = 0
    pending = dict(wanted)
    with open_any(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().rstrip(b"\n")
            fh.readline()
            q = fh.readline().rstrip(b"\n")
            if n in pending:
                es, eq = pending.pop(n)
                checked += 1
                if s != es or q != eq:
                    bad += 1
                    if bad <= 5:
                        print(f"MISMATCH {label} record {n} ({h.strip().decode()})")
                        print(f"   got  {s.decode()[:80]}\n   want {es.decode()[:80]}")
                        print(f"   got  {q.decode()[:80]}\n   want {eq.decode()[:80]}")
            n += 1
    if pending:
        bad += len(pending)
        print(f"MISSING {label}: {len(pending)} expected records beyond EOF (file has {n})")
    return checked, bad, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--stride", type=int, default=997, help="check every Nth PF cluster")
    ap.add_argument("--edge", type=int, default=32, help="also check first/last N PF clusters")
    ap.add_argument("--mismatches", type=int, default=1, help="allowed per index (tool's --barcode-mismatches)")
    ap.add_argument("--lane", type=int, default=0, help="only this lane (0 = all)")
    args = ap.parse_args()

    segs, i5_rc = parse_run_info(args.run)
    samples = parse_samplesheet(args.samplesheet)
    ids = [s[0] for s in samples]
    print(f"run structure: {segs}  i5 reverse-complement: {i5_rc}  samples: {len(samples)}")

    base = os.path.join(args.run, "Data", "Intensities", "BaseCalls")
    lane_dirs = sorted(glob.glob(os.path.join(base, "L00*")))
    total_bad = total_checked = 0
    for lane_dir in lane_dirs:
        lane_no = int(os.path.basename(lane_dir)[1:])
        if args.lane and lane_no != args.lane:
            continue
        # cluster sample: every stride-th cluster plus the edges
        filters = parse_filters(lane_dir, lane_no)
        n_pf = sum(int(m.sum()) for _, m in filters)
        pick = set(range(0, n_pf, args.stride))
        pick.update(range(min(args.edge, n_pf)))
        pick.update(range(max(0, n_pf - args.edge), n_pf))
        sample_idx = np.array(sorted(pick), dtype=np.int64)
        print(f"lane {lane_no}: {n_pf} PF clusters, checking {len(sample_idx)}", file=sys.stderr)

        dec = decode_lane(args.run, lane_dir, lane_no, segs, sample_idx)
        i1_len = sum(n for k, n in segs if k == "I1")
        assign = assign_samples(dec["idx_seq"], samples, lane_no, i5_rc, args.mismatches, i1_len)

        # ordinal of each cluster within its sample's output file
        ordinal = np.zeros(n_pf, dtype=np.int64)
        for si in list(range(len(samples))) + [-1]:
            m = assign == si
            ordinal[m] = np.arange(int(m.sum()))
        counts = {sid: int((assign == si).sum()) for si, sid in enumerate(ids)}
        counts["undetermined"] = int((assign == -1).sum())

        kinds = np.array(dec["kinds"])
        r1 = np.where(kinds == "R1")[0]
        r2 = np.where(kinds == "R2")[0]
        wanted = {}  # (sample, read) -> {ordinal: (seq, qual)}
        for row, cidx in enumerate(sample_idx):
            sid = ids[assign[cidx]] if assign[cidx] >= 0 else "undetermined"
            o = int(ordinal[cidx])
            wanted.setdefault((sid, 1), {})[o] = (dec["seq"][row, r1].tobytes(), dec["qual"][row, r1].tobytes())
            if len(r2):
                wanted.setdefault((sid, 2), {})[o] = (dec["seq"][row, r2].tobytes(), dec["qual"][row, r2].tobytes())

        for sid in ids + ["undetermined"]:
            for rd in (1, 2) if len(r2) else (1,):
                stem = os.path.join(args.output, f"{sid}_L{lane_no:03d}_R{rd}_001.fastq")
                path = stem + ".gz" if os.path.exists(stem + ".gz") else stem
                w = wanted.get((sid, rd), {})
                if not w and not os.path.exists(path):
                    if counts[sid]:
                        print(f"MISSING file {path} (expected {counts[sid]} reads)")
                        total_bad += 1
                    continue
                checked, bad, nrec = check_file(path, w, f"{sid} R{rd}")
                total_checked += checked
                total_bad += bad
                if nrec != counts[sid]:
                    print(f"COUNT {sid} R{rd}: file has {nrec}, reference expects {counts[sid]}")
                    total_bad += 1
        print(f"lane {lane_no}: matched {n_pf - counts['undetermined']} / {n_pf} "
              f"({100.0 * (n_pf - counts['undetermined']) / n_pf:.2f}%)")

    print(f"checked {total_checked} reads, {total_bad} problems")
    sys.exit(1 if total_bad else 0)


if __name__ == "__main__":
    main()
