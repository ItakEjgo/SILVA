import subprocess
import re
import sys
import argparse
import os

def run_cmd(algo, ratio, is_single_version, base_file, update_file, threads):
    cmd = [
        "taskset", f"-c", f"0-{threads-1}",
        "build/main",
        "-t", "batch-insert",
        "-a", algo,
        "-i", base_file,
        "-u", update_file,
        "-br", str(ratio),
        "-p", "0.2"
    ]
    if is_single_version:
        cmd.append("-sv")
        
    env = dict(os.environ)
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
    parser = argparse.ArgumentParser(description="运行 5 种算法组合对比实验")
    parser.add_argument("--dataset-base", default="dataset/uniform/1M_2_1.in", help="基础数据集路径")
    parser.add_argument("--dataset-update", default="dataset/uniform/1M_2_2.in", help="更新数据集路径")
    parser.add_argument("--ratio", type=float, default=0.001, help="批量更新比例")
    parser.add_argument("--threads", type=int, default=8, help="线程数")
    args = parser.parse_args()

    if not os.path.exists("build/main"):
        print("错误: 找不到 build/main，请先执行 cd build && make -j8 编译程序！")
        return

    configs = [
        ("SPaCtree", "spac", True, "单版本 (原地修改)"),
        ("PaCZUtree", "paczu", True, "单版本 (原地修改)"),
        ("PaCZtree", "pacz", True, "单版本 (原地修改)"),
        ("PaCZUtree", "paczu", False, "多版本 (写时复制)"),
        ("PaCZtree", "pacz", False, "多版本 (写时复制)"),
    ]

    print(f"=== 运行 5 种组合的基准实验 ===")
    print(f"数据集 Base: {args.dataset_base}")
    print(f"数据集 Update: {args.dataset_update}")
    print(f"Batch Ratio: {args.ratio}")
    print(f"Threads: {args.threads}")
    print("-" * 60)

    results = []
    for name, algo, sv, desc in configs:
        print(f"正在运行: {name:10} | 模式: {desc} ... ", end="", flush=True)
        t = run_cmd(algo, args.ratio, sv, args.dataset_base, args.dataset_update, args.threads)
        print(t)
        results.append((name, desc, t))

    print("\n" + "=" * 60)
    print(" 📊 终极基准对比表 ")
    print("=" * 60)
    print(f"| {'数据结构':<12} | {'运行场景 (版本模式)':<20} | {'实际总耗时':<15} |")
    print(f"|{'-'*14}|{'-'*22}|{'-'*17}|")
    for name, desc, t in results:
        print(f"| {name:<12} | {desc:<20} | {t:<15} |")

if __name__ == "__main__":
    main()
