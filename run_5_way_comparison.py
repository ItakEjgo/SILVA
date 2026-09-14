import subprocess
import re
import sys
import argparse
import os

def run_cmd(algo, ratio, is_single_version, base_file, update_file, threads, use_numa):
    cmd = [
        "build/main",
        "-t", "batch-insert",
        "-a", algo,
        "-i", base_file,
        "-u", update_file,
        "-br", str(ratio),
        "-p", "0.2"
    ]
    if use_numa:
        cmd = ["numactl", "-i", "all"] + cmd
        
    if threads > 0:
        cmd = ["taskset", "-c", f"0-{threads-1}"] + cmd
        
    if is_single_version:
        cmd.append("-sv")
        
    env = dict(os.environ)
    if threads > 0:
        env["PARLAY_NUM_THREADS"] = str(threads)
        env["OMP_NUM_THREADS"] = str(threads)
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=env)
        if res.returncode != 0:
            return "Error (Crash)"
        
        match = re.search(r'batch insert time \(avg\):\s+([\d\.]+)', res.stdout)
        if match:
            return f"{float(match.group(1)):.3f} s"
        else:
            return "Parse Error"
    except subprocess.TimeoutExpired:
        return "> 300.0 s (Timeout)"

def main():
    parser = argparse.ArgumentParser(description="Run 5-way algorithm comparison benchmark")
    parser.add_argument("--dataset-base", default="dataset/uniform/1M_2_1.in", help="Base dataset path")
    parser.add_argument("--dataset-update", default="dataset/uniform/1M_2_2.in", help="Update dataset path")
    parser.add_argument("--ratios", type=str, default="0.001", help="Comma-separated ratios, e.g. 0.001,0.01,0.1")
    parser.add_argument("--numa", action="store_true", help="Use numactl -i all")
    parser.add_argument("--threads", type=int, default=0, help="Threads (default 0 means all available cores)")
    args = parser.parse_args()

    if not os.path.exists("build/main"):
        print("Error: build/main not found! Please compile the program first: cd build && make -j8")
        return

    configs = [
        ("SPaCtree", "spac", True, "Single Version (In-place)"),
        ("PaCZUtree", "paczu", True, "Single Version (In-place)"),
        ("PaCZtree", "pacz", True, "Single Version (In-place)"),
        ("PaCZUtree", "paczu", False, "Multi Version (COW)"),
        ("PaCZtree", "pacz", False, "Multi Version (COW)"),
    ]

    print(f"=== Running 5-Way Benchmark Comparison ===")
    print(f"Batch Ratios: {args.ratios}")
    print(f"Threads: {args.threads}")
    print("-" * 60)

    ratios_list = [float(r.strip()) for r in args.ratios.split(",")]

    for r in ratios_list:
        results = []
        print(f"\n>>> Starting Test Ratio = {r} <<<")
        for name, algo, sv, desc in configs:
            print(f"Running: {name:10} | Mode: {desc:<25} ... ", end="", flush=True)
            t = run_cmd(algo, r, sv, args.dataset_base, args.dataset_update, args.threads, args.numa)
            print(t)
            results.append((name, desc, t))

        print("\n" + "=" * 70)
        print(f" 📊 Ultimate Benchmark Comparison (Ratio = {r}) ")
        print("=" * 70)
        print(f"| {'Data Structure':<14} | {'Version Mode':<27} | {'Total Time':<20} |")
        print(f"|{'-'*16}|{'-'*29}|{'-'*22}|")
        for name, desc, t in results:
            print(f"| {name:<14} | {desc:<27} | {t:<20} |")


if __name__ == "__main__":
    main()
