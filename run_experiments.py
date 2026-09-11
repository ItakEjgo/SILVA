import os
import subprocess
import re
import csv
import argparse
from datetime import datetime

class SilvaExperimentRunner:
    def __init__(self, threads=None):
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

    def parse_batch_time(self, output, keyword):
        if output in ["TIMEOUT", "CRASH", None]:
            return []
        
        results = []
        pattern = r'\[batch_size\]:\s*(\d+)\n.*?(?:time \(avg\)):\s*([\d\.]+)'
        matches = re.findall(pattern, output)
        for batch_size, time_s in matches:
            results.append((batch_size, time_s))
        return results

    def run_build_experiment(self):
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        suffix = f"_threads_{self.threads}" if self.threads else "_threads_all"
        csv_path = os.path.join(self.results_dir, f"exp_build{suffix}_{timestamp}.csv")
        
        print(f"=== Starting BUILD Experiment ===")
        print(f"Results will be saved to: {csv_path}")
        
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Distribution", "Size", "Algorithm", "Threads", "Build_Time_Seconds"])
            
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
                        threads_str = str(self.threads) if self.threads else "ALL"
                        writer.writerow([dist, size, algo, threads_str, time_s])
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
            writer.writerow(["Distribution", "Size", "Algorithm", "Threads", "Batch_Size", f"{task_name.capitalize()}_Time_Seconds"])
            
            for dist in self.distributions:
                for size in self.sizes:
                    dataset_base = f"{self.datasets_dir}/{dist}/{size}_2_1.in"
                    dataset_update = f"{self.datasets_dir}/{dist}/{size}_2_2.in"
                    if not os.path.exists(dataset_base) or not os.path.exists(dataset_update):
                        continue
                        
                    print(f"\n--- Testing Dataset: {dist} - {size} ---")
                    for algo in self.algorithms:
                        cmd = [self.binary_path, "-t", task_name, "-a", algo, "-i", dataset_base, "-u", dataset_update]
                        output = self.run_command(cmd, timeout=300)
                        
                        threads_str = str(self.threads) if self.threads else "ALL"
                        
                        if output in ["TIMEOUT", "CRASH", None]:
                            writer.writerow([dist, size, algo, threads_str, "N/A", output])
                            print(f"  -> Result: {output}")
                            continue
                            
                        res_list = self.parse_batch_time(output, action_keyword)
                        if not res_list:
                            writer.writerow([dist, size, algo, threads_str, "N/A", "PARSE_ERROR"])
                            print("  -> Result: PARSE_ERROR")
                        for b_size, t in res_list:
                            writer.writerow([dist, size, algo, threads_str, b_size, t])
                        f.flush()
                        print(f"  -> Parsed {len(res_list)} batch sizes.")
                        
        print(f"=== {task_name.upper()} Experiment Completed ===")

    def run_batch_insert_experiment(self):
        self.run_batch_experiment("batch-insert", "insert")

    def run_batch_delete_experiment(self):
        self.run_batch_experiment("batch-delete", "delete")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SILVA Automated Experiment Framework")
    parser.add_argument("--task", choices=["build", "insert", "delete", "query", "all"], required=True, help="Experiment task to run")
    parser.add_argument("--threads", type=int, default=None, help="Restrict to N threads (uses taskset -c 0 to N-1 if specified)")
    
    args = parser.parse_args()
    
    runner = SilvaExperimentRunner(threads=args.threads)
    
    if args.task in ["build", "all"]:
        runner.run_build_experiment()
    if args.task in ["insert", "all"]:
        runner.run_batch_insert_experiment()
    if args.task in ["delete", "all"]:
        runner.run_batch_delete_experiment()
    if args.task in ["query", "all"]:
        runner.run_query_experiment()
