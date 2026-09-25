#!/usr/bin/env python3
"""
Benchmark runner for cuda-demux

Discovers Illumina run folders under ~/Desktop/illumina, sorts by folder size
(ascending), runs cuda-demux sequentially (next starts only when previous finishes),
and records benchmark telemetry: CPU, GPU, memory, I/O, durations, and environment
details. Produces a machine-readable JSON per run and a summary CSV.

Requirements (optional fallbacks apply when unavailable):
  - psutil (process CPU/memory and IO; otherwise uses /proc and ps)
  - pynvml (GPU stats; otherwise falls back to nvidia-smi if present)

Example:
  python3 scripts/benchmark_cuda_demux.py \
    --base "~/Desktop/illumina" \
    --output-root ./benchmarks \
    --binary ./build/cuda-demux \
    --device 0

Notes:
  - Sample sheets: if a run folder contains a top-level SampleSheet.csv, it is used.
    Otherwise, pass --default-samplesheet to apply the same sheet to all runs lacking one.
  - Binary: if --binary is omitted, the script tries `cuda-demux` from PATH, then
    `./build/cuda-demux` relative to repo root.
  - GPU device selection: pass --device to query/stick to a specific GPU.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# Optional deps
try:
    import psutil  # type: ignore
except Exception:  # pragma: no cover
    psutil = None  # type: ignore

try:
    import pynvml  # type: ignore
except Exception:  # pragma: no cover
    pynvml = None  # type: ignore


def expand(p: str) -> Path:
    return Path(os.path.expanduser(os.path.expandvars(p))).resolve()


def find_cuda_demux(binary_arg: Optional[str]) -> str:
    if binary_arg:
        return str(expand(binary_arg))
    # Try PATH
    path_bin = shutil.which("cuda-demux")
    if path_bin:
        return path_bin
    # Try local build default
    local = Path(__file__).resolve().parents[1] / "build" / "cuda-demux"
    return str(local)


def list_run_dirs(base_dir: Path) -> List[Path]:
    """Return immediate subdirectories of base_dir that look like Illumina runs."""
    if not base_dir.is_dir():
        return []
    candidates = []
    for child in sorted(base_dir.iterdir()):
        if child.is_dir():
            # Heuristic: has Data/Intensities/BaseCalls or contains CBCLs
            basecalls = child / "Data" / "Intensities" / "BaseCalls"
            if basecalls.exists():
                candidates.append(child)
                continue
            # Fallback: has RunInfo.xml + any *.cbcl files under tree
            if (child / "RunInfo.xml").exists():
                for root, _dirs, files in os.walk(child):
                    if any(f.lower().endswith(".cbcl") for f in files):
                        candidates.append(child)
                        break
    return candidates


def folder_size_bytes(path: Path) -> int:
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += (Path(root) / f).stat().st_size
            except FileNotFoundError:
                pass
    return total


def detect_samplesheet(run_dir: Path, default: Optional[Path]) -> Optional[Path]:
    direct = run_dir / "SampleSheet.csv"
    if direct.exists():
        return direct
    # Search variations commonly used
    for pat in ["SampleSheet*.csv", "*sample*sheet*.csv", "*SampleSheet*.csv"]:
        matches = list(run_dir.glob(pat))
        if matches:
            return matches[0]
    return default


def cpu_info() -> Dict[str, str]:
    info: Dict[str, str] = {}
    # CPU model
    try:
        with open("/proc/cpuinfo", "r", encoding="utf-8", errors="ignore") as f:
            text = f.read()
        m = re.search(r"^model name\s*:\s*(.+)$", text, re.MULTILINE)
        if m:
            info["cpu_model"] = m.group(1).strip()
        cores = len([l for l in text.splitlines() if l.startswith("processor\t:")])
        if cores:
            info["cpu_cores_logical"] = str(cores)
    except Exception:
        pass
    # RAM total
    try:
        if psutil:
            vm = psutil.virtual_memory()
            info["mem_total_gb"] = f"{vm.total / (1024**3):.2f}"
        else:
            with open("/proc/meminfo", "r", encoding="utf-8", errors="ignore") as f:
                text = f.read()
            m = re.search(r"MemTotal:\s*(\d+) kB", text)
            if m:
                info["mem_total_gb"] = f"{int(m.group(1)) / (1024**2):.2f}"
    except Exception:
        pass
    info["platform"] = sys.platform
    info["python"] = sys.version.split()[0]
    return info


def init_nvml() -> bool:
    try:
        if pynvml is None:
            return False
        pynvml.nvmlInit()
        return True
    except Exception:
        return False


def gpu_env_info(device_index: Optional[int]) -> Dict[str, str]:
    info: Dict[str, str] = {}
    try:
        if pynvml is not None:
            if not init_nvml():
                raise RuntimeError
            info["nvidia_driver"] = pynvml.nvmlSystemGetDriverVersion().decode()
            info["cuda_demux_device"] = str(device_index) if device_index is not None else "default"
            if device_index is not None:
                h = pynvml.nvmlDeviceGetHandleByIndex(device_index)
                info["gpu_name"] = pynvml.nvmlDeviceGetName(h).decode()
                mem = pynvml.nvmlDeviceGetMemoryInfo(h)
                info["gpu_mem_total_gb"] = f"{mem.total / (1024**3):.2f}"
            else:
                # Report first GPU
                count = pynvml.nvmlDeviceGetCount()
                if count > 0:
                    h = pynvml.nvmlDeviceGetHandleByIndex(0)
                    info["gpu_name"] = pynvml.nvmlDeviceGetName(h).decode()
                    mem = pynvml.nvmlDeviceGetMemoryInfo(h)
                    info["gpu_mem_total_gb"] = f"{mem.total / (1024**3):.2f}"
            pynvml.nvmlShutdown()
            return info
    except Exception:
        pass
    # Fallback: nvidia-smi
    try:
        q = [
            "nvidia-smi",
            "--query-gpu=name,memory.total,driver_version",
            "--format=csv,noheader,nounits",
        ]
        if device_index is not None:
            q.extend(["-i", str(device_index)])
        out = subprocess.check_output(q, text=True).strip().splitlines()
        if out:
            name, mem_total, drv = [x.strip() for x in out[0].split(",")]
            info["gpu_name"] = name
            info["gpu_mem_total_gb"] = f"{float(mem_total) / 1024:.2f}"
            info["nvidia_driver"] = drv
    except Exception:
        pass
    return info


def get_cuda_demux_version(binary: str) -> Optional[str]:
    try:
        out = subprocess.check_output([binary, "--version"], text=True, stderr=subprocess.STDOUT, timeout=10)
        return out.strip()
    except Exception:
        return None


class Monitor:
    def __init__(self, pid: int, device_index: Optional[int], interval: float = 1.0):
        self.pid = pid
        self.interval = max(0.25, float(interval))
        self.device_index = device_index
        self._stop = threading.Event()
        self.thread = threading.Thread(target=self._run, daemon=True)
        self.samples: List[Dict[str, float]] = []
        self._ps_proc = None
        if psutil:
            try:
                self._ps_proc = psutil.Process(pid)
                # Prime CPU percent
                self._ps_proc.cpu_percent(interval=None)
            except Exception:
                self._ps_proc = None

        # NVML setup
        self._nvml_inited = init_nvml()
        self._nvml_handle = None
        if self._nvml_inited and self.device_index is not None:
            try:
                self._nvml_handle = pynvml.nvmlDeviceGetHandleByIndex(self.device_index)
            except Exception:
                self._nvml_handle = None

    def start(self):
        self.thread.start()

    def stop(self):
        self._stop.set()
        self.thread.join(timeout=self.interval * 2)
        if self._nvml_inited:
            try:
                pynvml.nvmlShutdown()
            except Exception:
                pass

    def _gpu_sample(self) -> Tuple[Optional[float], Optional[float], Optional[float]]:
        """Return (gpu_util_pct, gpu_mem_used_gb, gpu_proc_mem_gb)."""
        # Prefer NVML
        try:
            if self._nvml_inited:
                handle = self._nvml_handle
                if handle is None:
                    # default to device 0
                    handle = pynvml.nvmlDeviceGetHandleByIndex(0)
                util = pynvml.nvmlDeviceGetUtilizationRates(handle)
                mem = pynvml.nvmlDeviceGetMemoryInfo(handle)
                proc_mem_gb = None
                try:
                    procs = pynvml.nvmlDeviceGetComputeRunningProcesses_v3(handle)
                    for p in procs:
                        if int(p.pid) == self.pid:
                            proc_mem_gb = float(p.usedGpuMemory) / (1024 ** 3)
                            break
                except Exception:
                    pass
                return float(util.gpu), float(mem.used) / (1024 ** 3), proc_mem_gb
        except Exception:
            pass

        # Fallback to nvidia-smi
        try:
            q = [
                "nvidia-smi",
                "--query-gpu=utilization.gpu,memory.used",
                "--format=csv,noheader,nounits",
            ]
            if self.device_index is not None:
                q.extend(["-i", str(self.device_index)])
            out = subprocess.check_output(q, text=True).strip().splitlines()
            util_s, mem_s = [x.strip() for x in out[0].split(",")]
            gpu_util = float(util_s)
            gpu_mem_gb = float(mem_s) / 1024.0
        except Exception:
            gpu_util = None
            gpu_mem_gb = None

        # Per-process GPU mem via nvidia-smi
        proc_mem_gb = None
        try:
            out = subprocess.check_output(
                [
                    "nvidia-smi",
                    "--query-compute-apps=pid,used_memory",
                    "--format=csv,noheader,nounits",
                ],
                text=True,
            )
            for line in out.strip().splitlines():
                if not line.strip():
                    continue
                pid_s, used_s = [x.strip() for x in line.split(",")]
                if int(pid_s) == self.pid:
                    proc_mem_gb = float(used_s) / 1024.0
                    break
        except Exception:
            pass

        return gpu_util, gpu_mem_gb, proc_mem_gb

    def _proc_sample(self) -> Tuple[Optional[float], Optional[float], Optional[int], Optional[int]]:
        """Return (cpu_pct, rss_gb, read_bytes, write_bytes)."""
        if self._ps_proc is not None:
            try:
                cpu = float(self._ps_proc.cpu_percent(interval=None))  # percent of a single CPU
                rss_gb = float(self._ps_proc.memory_info().rss) / (1024 ** 3)
                io = None
                try:
                    io = self._ps_proc.io_counters()
                except Exception:
                    pass
                rbytes = int(io.read_bytes) if io else None
                wbytes = int(io.write_bytes) if io else None
                return cpu, rss_gb, rbytes, wbytes
            except Exception:
                return None, None, None, None
        # Fallback: ps for %CPU and RSS (kB)
        try:
            out = subprocess.check_output(["ps", "-p", str(self.pid), "-o", "%cpu=,rss="], text=True)
            parts = out.strip().split()
            if len(parts) >= 2:
                cpu = float(parts[0])
                rss_gb = float(int(parts[1])) / (1024 ** 2)
            else:
                cpu = None
                rss_gb = None
        except Exception:
            cpu = None
            rss_gb = None
        # Read/Write bytes via /proc/PID/io
        rbytes = None
        wbytes = None
        try:
            with open(f"/proc/{self.pid}/io", "r", encoding="utf-8", errors="ignore") as f:
                txt = f.read()
            m = re.search(r"read_bytes:\s*(\d+)", txt)
            if m:
                rbytes = int(m.group(1))
            m = re.search(r"write_bytes:\s*(\d+)", txt)
            if m:
                wbytes = int(m.group(1))
        except Exception:
            pass
        return cpu, rss_gb, rbytes, wbytes

    def _run(self):
        while not self._stop.is_set():
            cpu, rss_gb, rbytes, wbytes = self._proc_sample()
            gpu_util, gpu_mem_gb, proc_gpu_mem_gb = self._gpu_sample()
            ts = time.time()
            self.samples.append(
                {
                    "ts": ts,
                    "cpu_pct": cpu if cpu is not None else float("nan"),
                    "rss_gb": rss_gb if rss_gb is not None else float("nan"),
                    "gpu_util_pct": gpu_util if gpu_util is not None else float("nan"),
                    "gpu_mem_gb": gpu_mem_gb if gpu_mem_gb is not None else float("nan"),
                    "proc_gpu_mem_gb": proc_gpu_mem_gb if proc_gpu_mem_gb is not None else float("nan"),
                    "read_bytes": rbytes if rbytes is not None else -1,
                    "write_bytes": wbytes if wbytes is not None else -1,
                }
            )
            time.sleep(self.interval)


def agg(values: List[float]) -> Dict[str, Optional[float]]:
    vals = [v for v in values if v == v]  # filter NaN
    if not vals:
        return {"avg": None, "min": None, "max": None}
    return {
        "avg": sum(vals) / len(vals),
        "min": min(vals),
        "max": max(vals),
    }


def write_json(path: Path, data: Dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def write_csv(path: Path, rows: List[Dict[str, object]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    ap = argparse.ArgumentParser(description="Benchmark cuda-demux across Illumina run folders")
    ap.add_argument("--base", default="~/Desktop/illumina", help="Base directory containing run folders (default: ~/Desktop/illumina)")
    ap.add_argument("--output-root", default="./benchmarks", help="Directory to store benchmark outputs (JSON/CSV and per-run logs)")
    ap.add_argument("--binary", default=None, help="Path to cuda-demux binary; default: PATH or ./build/cuda-demux")
    ap.add_argument("--device", type=int, default=None, help="CUDA device index to pass to cuda-demux and monitor")
    ap.add_argument("--interval", type=float, default=1.0, help="Sampling interval in seconds (default: 1.0)")
    ap.add_argument("--default-samplesheet", default=None, help="Fallback SampleSheet.csv if a run lacks one")
    ap.add_argument("--extra-args", default="", help="Extra arguments to pass to cuda-demux (quoted string)")
    ap.add_argument("--gzip", action="store_true", help="Add --gzip to cuda-demux for compressed FASTQs")
    ap.add_argument("--dry-run", action="store_true", help="Only list discovered runs and the planned order; do not execute")

    args = ap.parse_args()

    base_dir = expand(args.base)
    out_root = expand(args.output_root)
    binary = find_cuda_demux(args.binary)
    if not Path(binary).exists():
        print(f"[WARN] cuda-demux binary not found at {binary}. The script will fail when executed.")

    default_ss = expand(args.default_samplesheet) if args.default_samplesheet else None

    runs = list_run_dirs(base_dir)
    if not runs:
        print(f"[ERROR] No run folders detected under {base_dir}")
        sys.exit(1)

    sized_runs: List[Tuple[Path, int]] = []
    print(f"Discovered {len(runs)} candidate run(s); computing sizes for ordering...")
    for r in runs:
        sz = folder_size_bytes(r)
        sized_runs.append((r, sz))
    sized_runs.sort(key=lambda x: x[1])

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    bench_dir = out_root / f"cuda_demux_benchmark_{timestamp}"
    bench_dir.mkdir(parents=True, exist_ok=True)

    # Summary info
    env = {**cpu_info(), **gpu_env_info(args.device)}
    env["cuda_demux_binary"] = binary
    env["cuda_demux_version"] = get_cuda_demux_version(binary) or "unknown"
    write_json(bench_dir / "environment.json", env)

    order_list = [
        {
            "run_path": str(p),
            "run_name": p.name,
            "size_bytes": sz,
            "size_gb": round(sz / (1024 ** 3), 3),
        }
        for p, sz in sized_runs
    ]
    write_json(bench_dir / "planned_order.json", {"runs": order_list})

    print("Planned order (smallest to largest):")
    for i, entry in enumerate(order_list, 1):
        print(f"  {i:2d}. {entry['run_name']}  ({entry['size_gb']} GB)")

    if args.dry_run:
        print("Dry-run enabled; exiting before execution.")
        return

    # Catch Ctrl-C gracefully
    stop_flag = {"stop": False}

    def on_sigint(_sig, _frm):
        stop_flag["stop"] = True
        print("\n[INFO] Interrupt received. Will stop after current run.")

    signal.signal(signal.SIGINT, on_sigint)

    summary_rows: List[Dict[str, object]] = []

    for idx, (run_dir, run_sz) in enumerate(sized_runs, 1):
        if stop_flag["stop"]:
            break
        run_name = run_dir.name
        print(f"\n[{idx}/{len(sized_runs)}] Starting run: {run_name}")

        # Resolve sample sheet
        ss_path = detect_samplesheet(run_dir, default_ss)
        if ss_path is None or not ss_path.exists():
            print(f"[WARN] No SampleSheet.csv for {run_name}. Skipping run.")
            continue

        # Prepare output directories
        run_out_root = bench_dir / "cuda_demux_outputs" / run_name
        run_out_root.mkdir(parents=True, exist_ok=True)
        log_path = bench_dir / "logs" / f"{run_name}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        metrics_json = bench_dir / "metrics" / f"{run_name}.json"

        # Build command
        cmd = [binary, "--input", str(run_dir), "--samplesheet", str(ss_path), "--output", str(run_out_root)]
        if args.device is not None:
            cmd += ["--device", str(args.device)]
        if args.gzip:
            cmd += ["--gzip"]
        if args.extra_args:
            cmd += shlex.split(args.extra_args)

        # Launch process
        start_time = time.time()
        with open(log_path, "w", encoding="utf-8") as log_f:
            log_f.write(f"Command: {' '.join(shlex.quote(c) for c in cmd)}\n")
            log_f.flush()
            proc = subprocess.Popen(cmd, stdout=log_f, stderr=subprocess.STDOUT, cwd=str(run_dir))
        mon = Monitor(proc.pid, args.device, interval=args.interval)
        mon.start()

        rc = None
        try:
            rc = proc.wait()
        finally:
            mon.stop()
        end_time = time.time()
        elapsed = end_time - start_time

        # Aggregate metrics
        samples = mon.samples
        cpu_stats = agg([s.get("cpu_pct", float("nan")) for s in samples])
        rss_stats = agg([s.get("rss_gb", float("nan")) for s in samples])
        gpu_stats = agg([s.get("gpu_util_pct", float("nan")) for s in samples])
        gpumem_stats = agg([s.get("gpu_mem_gb", float("nan")) for s in samples])
        proc_gpumem_stats = agg([s.get("proc_gpu_mem_gb", float("nan")) for s in samples])

        # IO delta
        read_b = [s.get("read_bytes", -1) for s in samples if s.get("read_bytes", -1) >= 0]
        write_b = [s.get("write_bytes", -1) for s in samples if s.get("write_bytes", -1) >= 0]
        read_bytes_delta = (max(read_b) - min(read_b)) if len(read_b) >= 2 else None
        write_bytes_delta = (max(write_b) - min(write_b)) if len(write_b) >= 2 else None

        run_record = {
            "run_name": run_name,
            "run_path": str(run_dir),
            "size_bytes": run_sz,
            "size_gb": round(run_sz / (1024 ** 3), 3),
            "command": cmd,
            "return_code": rc,
            "start_time": datetime.fromtimestamp(start_time).isoformat(),
            "end_time": datetime.fromtimestamp(end_time).isoformat(),
            "elapsed_seconds": round(elapsed, 3),
            "metrics": {
                "cpu_pct": cpu_stats,
                "rss_gb": rss_stats,
                "gpu_util_pct": gpu_stats,
                "gpu_mem_gb": gpumem_stats,
                "proc_gpu_mem_gb": proc_gpumem_stats,
                "io_bytes": {
                    "read_delta": read_bytes_delta,
                    "write_delta": write_bytes_delta,
                },
                "sample_interval_s": args.interval,
                "samples": len(samples),
            },
        }

        write_json(metrics_json, run_record)

        summary_rows.append(
            {
                "run_name": run_name,
                "size_gb": round(run_sz / (1024 ** 3), 3),
                "rc": rc,
                "elapsed_s": round(elapsed, 3),
                "cpu_avg_pct": None if cpu_stats["avg"] is None else round(cpu_stats["avg"], 2),
                "cpu_max_pct": None if cpu_stats["max"] is None else round(cpu_stats["max"], 2),
                "rss_avg_gb": None if rss_stats["avg"] is None else round(rss_stats["avg"], 3),
                "rss_max_gb": None if rss_stats["max"] is None else round(rss_stats["max"], 3),
                "gpu_avg_pct": None if gpu_stats["avg"] is None else round(gpu_stats["avg"], 2),
                "gpu_max_pct": None if gpu_stats["max"] is None else round(gpu_stats["max"], 2),
                "gpu_mem_avg_gb": None if gpumem_stats["avg"] is None else round(gpumem_stats["avg"], 3),
                "gpu_mem_max_gb": None if gpumem_stats["max"] is None else round(gpumem_stats["max"], 3),
                "proc_gpu_mem_max_gb": None if proc_gpumem_stats["max"] is None else round(proc_gpumem_stats["max"], 3),
                "io_read_delta_bytes": read_bytes_delta,
                "io_write_delta_bytes": write_bytes_delta,
            }
        )

        print(
            f"Completed {run_name}: rc={rc}, elapsed={elapsed:.1f}s, "
            f"CPU avg={summary_rows[-1]['cpu_avg_pct']}%, GPU avg={summary_rows[-1]['gpu_avg_pct']}%"
        )

    # Write summary
    if summary_rows:
        write_csv(bench_dir / "summary.csv", summary_rows)
        write_json(bench_dir / "summary.json", {"rows": summary_rows, "environment": env, "planned": order_list})
        print(f"\nBenchmark complete. Results in: {bench_dir}")
    else:
        print("\nNo runs executed. Nothing to summarize.")


if __name__ == "__main__":
    main()

