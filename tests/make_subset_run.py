#!/usr/bin/env python3
"""Cuts a small, self-consistent run folder out of a real CBCL run.

Keeps only the requested tiles (default: the first tile of every surface) in
every cycle's CBCL files and the matching .filter files, and copies
RunInfo.xml / RunParameters.xml / SampleSheet.csv. The result is a valid run
that cuda-demux and tests/reference_check.py can process in seconds, which
makes it a handy regression fixture for a given instrument layout.

    tests/make_subset_run.py --run RUN_DIR --out SUBSET_DIR [--tiles 1101,2101]
                             [--lanes 1] [--per-surface 1]
"""
import argparse
import glob
import os
import re
import shutil
import struct
import sys


def read_header(fh):
    head = fh.read(12)
    ver, hs, bpb, bpq, nbins = struct.unpack("<HIBBI", head)
    bins = fh.read(8 * nbins)
    (nt,) = struct.unpack("<I", fh.read(4))
    tiles = [struct.unpack("<IIII", fh.read(16)) for _ in range(nt)]
    flag = fh.read(1)
    assert fh.tell() == hs, "header size mismatch"
    return dict(ver=ver, bpb=bpb, bpq=bpq, bins=bins, tiles=tiles, flag=flag, header_size=hs)


def subset_cbcl(src, dst, keep):
    with open(src, "rb") as fh:
        h = read_header(fh)
        kept, blocks = [], []
        offset = h["header_size"]
        for t in h["tiles"]:
            tile_id, n, unc, comp = t
            if tile_id in keep:
                fh.seek(offset)
                blocks.append(fh.read(comp))
                kept.append(t)
            offset += comp
    header_size = 12 + len(h["bins"]) + 4 + 16 * len(kept) + 1
    with open(dst, "wb") as out:
        out.write(struct.pack("<HIBBI", h["ver"], header_size, h["bpb"], h["bpq"], len(h["bins"]) // 8))
        out.write(h["bins"])
        out.write(struct.pack("<I", len(kept)))
        for t in kept:
            out.write(struct.pack("<IIII", *t))
        out.write(h["flag"])
        for b in blocks:
            out.write(b)
    return len(kept)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tiles", help="comma-separated tile ids to keep (default: first per surface)")
    ap.add_argument("--per-surface", type=int, default=1, help="tiles to keep per surface when --tiles is not given")
    ap.add_argument("--lanes", help="comma-separated lane numbers (default: all)")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    for name in ("RunInfo.xml", "RunParameters.xml", "SampleSheet.csv"):
        src = os.path.join(args.run, name)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(args.out, name))

    base = os.path.join(args.run, "Data", "Intensities", "BaseCalls")
    lanes = [int(x) for x in args.lanes.split(",")] if args.lanes else None
    for lane_dir in sorted(glob.glob(os.path.join(base, "L00*"))):
        lane_no = int(os.path.basename(lane_dir)[1:])
        if lanes and lane_no not in lanes:
            continue
        filters = {}
        for f in glob.glob(os.path.join(lane_dir, f"s_{lane_no}_*.filter")):
            filters[int(re.search(r"_(\d+)\.filter$", f).group(1))] = f
        if args.tiles:
            keep = {int(x) for x in args.tiles.split(",")}
        else:
            keep = set()
            by_surface = {}
            for t in sorted(filters):
                by_surface.setdefault(t // 1000, []).append(t)
            for ts in by_surface.values():
                keep.update(ts[:args.per_surface])
        missing = keep - set(filters)
        if missing:
            sys.exit(f"lane {lane_no}: no filter file for tiles {sorted(missing)}")

        out_lane = os.path.join(args.out, "Data", "Intensities", "BaseCalls", f"L{lane_no:03d}")
        os.makedirs(out_lane, exist_ok=True)
        for t in sorted(keep):
            shutil.copy(filters[t], os.path.join(out_lane, os.path.basename(filters[t])))

        cycles = sorted(glob.glob(os.path.join(lane_dir, "C*.1")),
                        key=lambda p: int(os.path.basename(p)[1:].split(".")[0]))
        for cdir in cycles:
            out_c = os.path.join(out_lane, os.path.basename(cdir))
            os.makedirs(out_c, exist_ok=True)
            for f in sorted(glob.glob(os.path.join(cdir, "*.cbcl"))):
                subset_cbcl(f, os.path.join(out_c, os.path.basename(f)), keep)
        print(f"lane {lane_no}: kept tiles {sorted(keep)} across {len(cycles)} cycles")


if __name__ == "__main__":
    main()
