import os
import subprocess
import re
import csv
import argparse
from datetime import datetime

class SilvaExperimentRunner:
    def __init__(self, threads=None, ratios=None, dataset_base=None, dataset_update=None, snapshot_ratio="0.2"):
        self.threads = threads
        self.ratios = ratios
        self.dataset_base_file = dataset_base
        self.dataset_update_file = dataset_update
        self.snapshot_ratio = snapshot_ratio
        self.base_dir = "."
        self.binary_path = "build/main"
        self.datasets_dir = "dataset"
        self.distributions = ["uniform"]
        self.sizes = ["10M"]
        self.algorithms = ["mvq", "pacz", "rlog", "pkdlog"]
        self.results_dir = "results_experiments"
        
        os.makedirs(self.results_dir, exist_ok=True)

    def get_datasets(self):
        if self.dataset_base_file:
            dist_name = "custom"
            size_name = "custom"
            yield (dist_name, size_name, self.dataset_base_file, self.dataset_update_file)
        else:
            for dist in self.distributions:
                for size in self.sizes:
                    db = f"{self.datasets_dir}/{dist}/{size}_2_1.in"
                    du = f"{self.datasets_dir}/{dist}/{size}_2_2.in"
                    yield (dist, size, db, du)

    def run_command(self, cmd, timeout=1800):
        # Apply thread restrictions
        env = os.environ.copy()
        if self.threads is not None:
            env["PARLAY_NUM_THREADS"] = str(self.threads)
            env["OMP_NUM_THREADS"] = str(self.threads)
            cmd = ["taskset", "-c", "0"] + cmd

        print(f"  [RUN] {' '.join(cmd)}")
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
            if res.returncode != 0:
                print(f"  [ERROR] Process exited with code {res.returncode}")
                return res.stdout + f"\n[PROCESS_CRASHED_WITH_CODE_{res.returncode}]"
            return res.stdout
        except subprocess.TimeoutExpired:
            print("  [TIMEOUT] Process took too long.")
            return "TIMEOUT"
        except Exception as e:
            print(f"  [CRASH] {e}")
            return "CRASH"

    def parse_time(self, output, keyword):
        if output in ["TIMEOUT", "CRASH", None]:
            return output
        pattern = rf'\[{keyword}\].*?time \(avg\):\s*([\d\.]+)|\[.*{keyword}.*\].*?time \(avg\):\s*([\d\.]+)'
        # Some outputs might be [MVQ] or [Zorder-PACZ] or [PkdTree]
        # Actually our C++ output is like: [PkdTree]: build time (avg): 0.123
        pattern2 = r'(?:\]:?.*?time \(avg\):?)\s*([\d\.]+)'
        match = re.search(pattern2, output)
        if match:
            return match.group(1)
        return "PARSE_ERROR"
    def parse_memory(self, output):
        if output in ["TIMEOUT", "CRASH", None]:
            return "N/A"
        pattern = r'\[memory_MB\]:\s*([\d\.]+)'
        match = re.search(pattern, output)
        if match:
            return match.group(1)
        return "N/A"

    def parse_batch_time(self, output, keyword):
        if output in ["TIMEOUT", "CRASH", None]:
            return []
        
        results = []
        pattern = r'(?:\[memory_MB\]:\s*([\d\.]+)\n)?.*?\[batch_(?:size|ratio)\]:\s*([\d\.]+)\n.*?(?:time \(avg\)):\s*([\d\.]+)'
        matches = re.findall(pattern, output)
        for mem, batch_size, time_s in matches:
            mem_val = mem if mem else "N/A"
            results.append((batch_size, time_s, mem_val))
        return results

    def run_build_experiment(self):
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        suffix = f"_threads_{self.threads}" if self.threads else "_threads_all"
        csv_path = os.path.join(self.results_dir, f"exp_build{suffix}_{timestamp}.csv")
        
        print(f"=== Starting BUILD Experiment ===")
        print(f"Results will be saved to: {csv_path}")
        
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Distribution", "Size", "Algorithm", "Threads", "Build_Time_Seconds", "Memory_MB"])
            
            for dist, size, dataset_base, _ in self.get_datasets():
                if not dataset_base or not os.path.exists(dataset_base):
                    print(f"File not found: {dataset_base}")
                    continue
                        
                print(f"\n--- Testing Dataset: {dist} - {size} ---")
                for algo in self.algorithms:
                    cmd = [self.binary_path, "-t", "build", "-a", algo, "-i", dataset_base]
                    output = self.run_command(cmd, timeout=300)
                    audit_log_path = os.path.join(self.results_dir, f"audit_log_build{suffix}_{timestamp}.txt")
                    with open(audit_log_path, 'a') as af:
                        af.write(f"\n{'='*60}\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] CMD: {' '.join(cmd)}\n{'='*60}\n")
                        if output: af.write(output + "\n")
                     # Increased timeout for single core 50M
                    
                    time_s = self.parse_time(output, "build")
                    mem_mb = self.parse_memory(output)
                    threads_str = str(self.threads) if self.threads else "ALL"
                    writer.writerow([dist, size, algo, threads_str, time_s, mem_mb])
                    f.flush()
                    print(f"  -> Result: {time_s}s")
                    
        self.plot_build_results(csv_path)
        print(f"=== BUILD Experiment Completed ===")


    def plot_build_results(self, csv_path):
        try:
            import matplotlib.pyplot as plt
            import numpy as np
            from collections import defaultdict
        except ImportError:
            print("[Warning] matplotlib or numpy not installed. Skipping plot generation.")
            return

        time_data = defaultdict(lambda: defaultdict(list))
        mem_data = defaultdict(lambda: defaultdict(list))
        
        with open(csv_path, 'r') as f:
            import csv
            reader = csv.DictReader(f)
            for row in reader:
                if row["Build_Time_Seconds"] in ["N/A", "PARSE_ERROR", "TIMEOUT", "CRASH"]: continue
                dist = row["Distribution"]
                size = row["Size"]
                algo = row["Algorithm"]
                time_s = float(row["Build_Time_Seconds"])
                mem_mb = row["Memory_MB"]
                
                time_data[dist][size].append((algo, time_s))
                if mem_mb not in ["N/A", "PARSE_ERROR"]:
                    mem_data[dist][size].append((algo, float(mem_mb)))

        import os
        plot_dir = os.path.join(self.results_dir, "build_plots")
        os.makedirs(plot_dir, exist_ok=True)

        for dist, sizes_dict in time_data.items():
            sizes = list(sizes_dict.keys())
            # Ensure sizes are sorted properly (1M, 10M, 20M...)
            sizes.sort(key=lambda x: int(x.replace('M', '')) if 'M' in x else 0)
            
            algos = self.algorithms
            
            # Plot Time
            plt.figure(figsize=(12, 6))
            x = np.arange(len(sizes))
            width = 0.15
            
            for i, algo in enumerate(algos):
                y = []
                for sz in sizes:
                    val = next((v for a, v in sizes_dict[sz] if a == algo), 0)
                    y.append(val)
                plt.bar(x + i*width - width*len(algos)/2, y, width, label=algo)
            
            plt.title(f'Build Time (Seconds) - {dist}')
            plt.xlabel('Dataset Size')
            plt.ylabel('Time (Seconds)')
            plt.xticks(x, sizes)
            plt.legend()
            plt.grid(axis='y')
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, f'build_time_{dist}.png'))
            plt.close()

        for dist, sizes_dict in mem_data.items():
            sizes = list(sizes_dict.keys())
            sizes.sort(key=lambda x: int(x.replace('M', '')) if 'M' in x else 0)
            algos = self.algorithms
            
            plt.figure(figsize=(12, 6))
            x = np.arange(len(sizes))
            width = 0.15
            
            for i, algo in enumerate(algos):
                y = []
                for sz in sizes:
                    val = next((v for a, v in sizes_dict[sz] if a == algo), 0)
                    y.append(val)
                plt.bar(x + i*width - width*len(algos)/2, y, width, label=algo)
            
            plt.title(f'Build Memory (MB) - {dist}')
            plt.xlabel('Dataset Size')
            plt.ylabel('Memory (MB)')
            plt.xticks(x, sizes)
            plt.legend()
            plt.grid(axis='y')
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, f'build_mem_{dist}.png'))
            plt.close()
        print(f"Generated build plots in {plot_dir}")

    def run_batch_experiment(self, task_name, action_keyword):
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        suffix = f"_threads_{self.threads}" if self.threads else "_threads_all"
        csv_path = os.path.join(self.results_dir, f"exp_{task_name}{suffix}_{timestamp}.csv")
        
        print(f"=== Starting {task_name.upper()} Experiment ===")
        print(f"Results will be saved to: {csv_path}")
        
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Distribution", "Size", "Algorithm", "Threads", "Batch_Size", f"{task_name.capitalize()}_Time_Seconds", "Memory_MB"])
            
            for dist, size, dataset_base, dataset_update in self.get_datasets():
                if not dataset_base or not dataset_update or not os.path.exists(dataset_base) or not os.path.exists(dataset_update):
                    print(f"Dataset pairs not found: {dataset_base}, {dataset_update}")
                    continue
                        
                print(f"\n--- Testing Dataset: {dist} - {size} ---")
                for algo in self.algorithms:
                    cmd = [self.binary_path, "-t", task_name, "-a", algo, "-i", dataset_base, "-u", dataset_update]
                    if self.ratios:
                        cmd.extend(["-br", self.ratios])
                    if self.snapshot_ratio:
                        cmd.extend(["-p", str(self.snapshot_ratio)])
                    output = self.run_command(cmd, timeout=1800)
                    audit_log_path = os.path.join(self.results_dir, f"audit_log_{task_name}{suffix}_{timestamp}.txt")
                    with open(audit_log_path, 'a') as af:
                        af.write(f"\n{'='*60}\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] CMD: {' '.join(cmd)}\n{'='*60}\n")
                        if output: af.write(output + "\n")
                    
                    
                    threads_str = str(self.threads) if self.threads else "ALL"
                    
                    res_list = self.parse_batch_time(output, action_keyword)
                    # Capture and plot per-batch data!
                    batch_times_matches = re.findall(r'\[per_batch_time\]:\s*(.*)', output) if output else []
                    batch_mems_matches = re.findall(r'\[per_batch_mem\]:\s*(.*)', output) if output else []
                    
                    if not batch_times_matches and output and 'step_time' in output:
                        step_times = re.findall(r'\[step_time\]:\s*([\d\.]+)', output)
                        step_mems = re.findall(r'\[step_mem\]:\s*([\d\.]+)', output)
                        if step_times and step_mems:
                            batch_times_matches = [','.join(step_times)]
                            batch_mems_matches = [','.join(step_mems)]

                    if not batch_times_matches and output and '[per_batch_time_val]:' in output:
                        times_vals = re.findall(r'\[per_batch_time_val\]:\s*([\d\.]+)', output)
                        mems_vals = re.findall(r'\[per_batch_mem_val\]:\s*([\d\.]+)', output)
                        if times_vals: batch_times_matches = [','.join(times_vals)]
                        if mems_vals: batch_mems_matches = [','.join(mems_vals)]

                    if output in ["TIMEOUT", "CRASH", None] or (output and "[PROCESS_CRASHED" in output):
                        writer.writerow([dist, size, algo, threads_str, "N/A", "CRASH", "OOM"])
                        f.flush()
                        # Do not continue, let it plot the partial data if any!

                    if batch_times_matches and batch_mems_matches:
                        ratios_matches = re.findall(r'\[Testing Ratio\]:\s*([\d\.]+)', output) if output else []
                        # Fallback if testing ratio not found
                        ratios_to_plot = [r[0] for r in res_list] if res_list else ratios_matches
                        if not ratios_to_plot and self.ratios:
                            ratios_to_plot = self.ratios.split(',')
                        
                        for idx, batch_size in enumerate(ratios_to_plot):
                            if idx < len(batch_times_matches) and idx < len(batch_mems_matches):
                                t_str = batch_times_matches[idx].strip()
                                m_str = batch_mems_matches[idx].strip()
                                if t_str and m_str:
                                    t_arr = [float(x) for x in t_str.split(',') if x]
                                    m_arr = [float(x) for x in m_str.split(',') if x]
                                    self.plot_single_ratio(task_name, dist, size, algo, batch_size, t_arr, m_arr)

                    if not res_list:
                        writer.writerow([dist, size, algo, threads_str, "N/A", "PARSE_ERROR", "N/A"])
                    else:
                        for batch_size, time_s, mem_val in res_list:
                            writer.writerow([dist, size, algo, threads_str, batch_size, time_s, mem_val])
                    f.flush()
                    f.flush()
                    print(f"  -> Parsed {len(res_list)} batch sizes.")
                    
        self.plot_results(csv_path, task_name)
        print(f"=== {task_name.upper()} Experiment Completed ===")

    def plot_single_ratio(self, task_name, dist, size, algo, ratio, times, mems):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return
        import os
        ratio_str = f"{float(ratio):g}"
        plot_dir = os.path.join(self.results_dir, "batch_plots", f"ratio_{ratio_str}")
        os.makedirs(plot_dir, exist_ok=True)
        
        fig, ax1 = plt.subplots(figsize=(10, 6))
        
        color = 'tab:red'
        ax1.set_xlabel('Batch Index')
        ax1.set_ylabel('Time (ms)', color=color)
        ax1.plot(range(1, len(times) + 1), times, color=color, marker='o', markersize=3, label='Time (ms)')
        ax1.tick_params(axis='y', labelcolor=color)
        
        ax2 = ax1.twinx()  
        color = 'tab:blue'
        ax2.set_ylabel('Memory (MB)', color=color)  
        ax2.plot(range(1, len(mems) + 1), mems, color=color, marker='s', markersize=3, label='Memory (MB)')
        ax2.tick_params(axis='y', labelcolor=color)
        
        plt.title(f'{task_name.capitalize()} Per-Batch - {algo} (Ratio={ratio}, {dist}-{size})')
        fig.tight_layout()
        
        out_png = os.path.join(plot_dir, f'{task_name}_{algo}_ratio{ratio}_{dist}_{size}.png')
        plt.savefig(out_png)
        plt.close()
        print(f"  -> Generated per-batch plot for {algo} (Ratio={ratio})")

    def plot_results(self, csv_path, task_name):
        try:
            import matplotlib.pyplot as plt
            from collections import defaultdict
        except ImportError:
            print("[Warning] matplotlib not installed. Skipping plot generation.")
            return

        time_data = defaultdict(lambda: defaultdict(list))
        mem_data = defaultdict(lambda: defaultdict(list))
        
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            time_col = f"{task_name.capitalize()}_Time_Seconds"
            for row in reader:
                if row["Batch_Size"] == "N/A": continue
                ds = f'{row["Distribution"]}-{row["Size"]}'
                algo = row["Algorithm"]
                batch = float(row["Batch_Size"])
                time_s = float(row[time_col])
                mem_mb = row["Memory_MB"]
                
                time_data[ds][algo].append((batch, time_s))
                if mem_mb != "N/A" and mem_mb != "PARSE_ERROR":
                    mem_data[ds][algo].append((batch, float(mem_mb)))

        import os
        for ds, algos in time_data.items():
            plt.figure(figsize=(10, 6))
            for algo, vals in algos.items():
                vals.sort()
                x = [v[0] for v in vals]
                y = [v[1] for v in vals]
                plt.plot(x, y, marker='o', label=algo)
            plt.title(f'{task_name.capitalize()} Time vs Batch Size ({ds})')
            plt.xlabel('Batch Size (Ratio)')
            plt.ylabel('Time (Seconds)')
            plt.xscale('log')
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(csv_path.replace('.csv', f'_{ds}_time.png'))
            plt.close()

        for ds, algos in mem_data.items():
            plt.figure(figsize=(10, 6))
            for algo, vals in algos.items():
                vals.sort()
                x = [v[0] for v in vals]
                y = [v[1] for v in vals]
                plt.plot(x, y, marker='o', label=algo)
            plt.title(f'{task_name.capitalize()} Memory vs Batch Size ({ds})')
            plt.xlabel('Batch Size (Ratio)')
            plt.ylabel('Memory (MB)')
            plt.xscale('log')
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(csv_path.replace('.csv', f'_{ds}_mem.png'))
            plt.close()
        print(f"Generated plots for {csv_path}")


    def run_batch_insert_experiment(self):
        self.run_batch_experiment("batch-insert", "insert")

    def run_batch_delete_experiment(self):
        self.run_batch_experiment("batch-delete", "delete")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SILVA Automated Experiment Framework")
    parser.add_argument("--task", choices=["build", "insert", "delete", "query", "all"], required=True, help="Experiment task to run")
    parser.add_argument("--threads", type=int, default=None, help="Restrict to N threads (uses taskset -c 0 to N-1 if specified)")
    parser.add_argument("--ratios", type=str, default=None, help="Comma-separated batch ratios (e.g., '0.01,0.1,0.25,0.5,1.0')")
    parser.add_argument("--dataset-base", type=str, default=None, help="Specific base dataset file (e.g. dataset/uniform/10M_2_1.in)")
    parser.add_argument("--dataset-update", type=str, default=None, help="Specific update dataset file for batch ops (e.g. dataset/uniform/10M_2_2.in)")
    parser.add_argument("--snapshot-ratio", type=str, default="0.2", help="Snapshot ratio parameter (-p) passed to binary (default: 0.2)")
    
    args = parser.parse_args()
    
    runner = SilvaExperimentRunner(
        threads=args.threads,
        ratios=args.ratios,
        dataset_base=args.dataset_base,
        dataset_update=args.dataset_update,
        snapshot_ratio=args.snapshot_ratio
    )
    
    if args.task in ["build", "all"]:
        runner.run_build_experiment()
    if args.task in ["insert", "all"]:
        runner.run_batch_insert_experiment()
    if args.task in ["delete", "all"]:
        runner.run_batch_delete_experiment()
    if args.task in ["query", "all"]:
        runner.run_query_experiment()
