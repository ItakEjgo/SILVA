import os
import subprocess

sizes = [
    (1000000, "1M"),
    (10000000, "10M"),
    (20000000, "20M"),
    (30000000, "30M"),
    (40000000, "40M"),
    (50000000, "50M")
]
dim = 2
varden_dir = "dataset/varden"

print("Generating Varden Datasets...")
dbscan_path = "baselines/pkdtree/script/DBSCAN"
wash_script = "baselines/pkdtree/script/wash_varden.py"

for size_num, size_name in sizes:
    outfile = f"{varden_dir}/{size_name}_2.in"
    if os.path.exists(outfile):
        print(f"Skipping {outfile} (already exists)")
        continue
    
    temp_raw = f"{varden_dir}/raw_{size_name}.in"
    
    print(f"Running DBSCAN for {size_name} Varden...")
    subprocess.run([dbscan_path, "-algo", "0", "-ds", temp_raw, "-n", str(size_num), "-d", str(dim), "-vd", "1"], check=True)
    
    print(f"Washing data for {size_name} Varden...")
    subprocess.run(["python3", wash_script, temp_raw, outfile], check=True)
    
    os.remove(temp_raw)
    print(f"Generated {outfile}")

print("All Varden datasets generated successfully!")
