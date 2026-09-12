import os
import subprocess
import re
import csv
import argparse
from datetime import datetime

class SilvaExperimentRunner:
    def __init__(self, threads=None, ratios=None):
        self.threads = threads
        self.ratios = ratios
        self.base_dir = "."
        self.binary_path = "build/main"
        self.datasets_dir = "dataset"
        self.results_dir = "results_experiments"
        
        self.distributions = ["uniform", "varden"]
        self.sizes = ["1M", "10M", "20M", "30M", "40M", "50M"]
        self.algorithms = ["mvq", "pacz", "boost", "rlog", "pkdtree", "pkdlog"]
        
        os.makedirs(self.results_dir, exist_ok=True)

    def run_command(self, cmd, timeout=300):
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
                return None
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
            
            for dist in self.distributions:
                for size in self.sizes:
                    dataset_base = f"{self.datasets_dir}/{dist}/{size}_2_1.in"
                    if not os.path.exists(dataset_base):
                        continue
                        
                    print(f"\n--- Testing Dataset: {dist} - {size} ---")
                    for algo in self.algorithms:
                        cmd = [self.binary_path, "-t", "build", "-a", algo, "-i", dataset_base]
                        output = self.run_command(cmd, timeout=300) # Increased timeout for single core 50M
                        
                        time_s = self.parse_time(output, "build")
                        mem_mb = self.parse_memory(output)
                        threads_str = str(self.threads) if self.threads else "ALL"
                        writer.writerow([dist, size, algo, threads_str, time_s, mem_mb])
                        f.flush()
                        print(f"  -> Result: {time_s}s")
                        
        print(f"=== BUILD Experiment Completed ===")

    def run_batch_experiment(self, task_name, action_keyword):
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        suffix = f"_threads_{self.threads}" if self.threads else "_threads_all"
        csv_path = os.path.join(self.results_dir, f"exp_{task_name}{suffix}_{timestamp}.csv")
        
        print(f"=== Starting {task_name.upper()} Experiment ===")
        print(f"Results will be saved to: {csv_path}")
        
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Distribution", "Size", "Algorithm", "Threads", "Batch_Size", f"{task_name.capitalize()}_Time_Seconds", "Memory_MB"])
            
            for dist in self.distributions:
                for size in self.sizes:
                    dataset_base = f"{self.datasets_dir}/{dist}/{size}_2_1.in"
                    dataset_update = f"{self.datasets_dir}/{dist}/{size}_2_2.in"
                    if not os.path.exists(dataset_base) or not os.path.exists(dataset_update):
                        continue
                        
                    print(f"\n--- Testing Dataset: {dist} - {size} ---")
                    for algo in self.algorithms:
                        cmd = [self.binary_path, "-t", task_name, "-a", algo, "-i", dataset_base, "-u", dataset_update]
                        if self.ratios:
                            cmd.extend(["-br", self.ratios])
                        output = self.run_command(cmd, timeout=300)
                        
                        threads_str = str(self.threads) if self.threads else "ALL"
                        
                        res_list = self.parse_batch_time(output, action_keyword)

                        # Capture and plot per-chunk data!
                        chunk_times_matches = re.findall(r'\[per_chunk_time\]:\s*(.*)', output)
                        chunk_mems_matches = re.findall(r'\[per_chunk_mem\]:\s*(.*)', output)
                        if not chunk_times_matches and '[per_chunk_time_val]:' in output:
                            times_vals = re.findall(r'\[per_chunk_time_val\]:\s*([\d\.]+)', output)
                            mems_vals = re.findall(r'\[per_chunk_mem_val\]:\s*([\d\.]+)', output)
                            if times_vals: chunk_times_matches = [','.join(times_vals)]
                            if mems_vals: chunk_mems_matches = [','.join(mems_vals)]

                        if chunk_times_matches and chunk_mems_matches and res_list:
                            for idx, (batch_size, _, _) in enumerate(res_list):
                                if idx < len(chunk_times_matches) and idx < len(chunk_mems_matches):
                                    t_str = chunk_times_matches[idx].strip()
                                    m_str = chunk_mems_matches[idx].strip()
                                    if t_str and m_str:
                                        t_arr = [float(x) for x in t_str.split(',') if x]
                                        m_arr = [float(x) for x in m_str.split(',') if x]
                                        self.plot_single_batch(task_name, dist, size, algo, batch_size, t_arr, m_arr)

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

    def plot_single_batch(self, task_name, dist, size, algo, ratio, times, mems):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return
            
        import os
        plot_dir = os.path.join(self.results_dir, "chunk_plots")
        os.makedirs(plot_dir, exist_ok=True)
        
        fig, ax1 = plt.subplots(figsize=(10, 6))
        
        color = 'tab:red'
        ax1.set_xlabel('Chunk Index')
        ax1.set_ylabel('Time (ms)', color=color)
        ax1.plot(range(1, len(times) + 1), times, color=color, marker='o', markersize=3, label='Time (ms)')
        ax1.tick_params(axis='y', labelcolor=color)
        
        ax2 = ax1.twinx()  
        color = 'tab:blue'
        ax2.set_ylabel('Memory (MB)', color=color)  
        ax2.plot(range(1, len(mems) + 1), mems, color=color, marker='s', markersize=3, label='Memory (MB)')
        ax2.tick_params(axis='y', labelcolor=color)
        
        plt.title(f'{task_name.capitalize()} Per-Chunk - {algo} (Ratio={ratio}, {dist}-{size})')
        fig.tight_layout()
        
        out_png = os.path.join(plot_dir, f'{task_name}_{algo}_ratio{ratio}_{dist}_{size}.png')
        plt.savefig(out_png)
        plt.close()
        print(f"  -> Generated per-chunk plot for {algo} (Ratio={ratio})")

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
    
    args = parser.parse_args()
    
    runner = SilvaExperimentRunner(threads=args.threads, ratios=args.ratios)
    
    if args.task in ["build", "all"]:
        runner.run_build_experiment()
    if args.task in ["insert", "all"]:
        runner.run_batch_insert_experiment()
    if args.task in ["delete", "all"]:
        runner.run_batch_delete_experiment()
    if args.task in ["query", "all"]:
        runner.run_query_experiment()
