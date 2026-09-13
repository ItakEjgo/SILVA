import os
import subprocess
import shutil

sizes = [
    (1000000, "1M"),
    (10000000, "10M"),
    (20000000, "20M"),
    (30000000, "30M"),
    (40000000, "40M"),
    (50000000, "50M")
]
dim = 2

uniform_dir = "dataset/uniform"
varden_dir = "dataset/varden"

# 1. Generate Uniform
print("Resuming Uniform Datasets...")
generator_path = "baselines/spacetreelib/build/data_generator"
for size_num, size_name in sizes:
    outfile = f"{uniform_dir}/{size_name}_2.in"
    if os.path.exists(outfile):
        print(f"Skipping {outfile} (already exists)")
        continue
    
    temp_path = "dataset/temp_uniform"
    os.makedirs(temp_path, exist_ok=True)
    print(f"Generating {size_name} Uniform...")
    subprocess.run([generator_path, temp_path, str(size_num), str(dim), "1", "0"], check=True)
    
    gen_file = f"{temp_path}/{size_num}_{dim}/1.in"
    if os.path.exists(gen_file):
        shutil.move(gen_file, outfile)
        print(f"Generated {outfile}")
    
    shutil.rmtree(temp_path)

# 2. Generate Varden
print("Generating Varden Datasets...")
dbscan_path = "baselines/spacetreelib/script/DBSCAN"
wash_script = "baselines/spacetreelib/script/wash_varden.py"
varDensity = "5"

for size_num, size_name in sizes:
    outfile = f"{varden_dir}/{size_name}_2.in"
    if os.path.exists(outfile):
        print(f"Skipping {outfile} (already exists)")
        continue
    
    temp_raw = f"{varden_dir}/raw_{size_name}.in"
    
    print(f"Running DBSCAN for {size_name} Varden...")
    subprocess.run([dbscan_path, "-algo", "0", "-ds", temp_raw, "-n", str(size_num), "-d", str(dim), "-vd", varDensity], check=True)
    
    print(f"Washing data for {size_name} Varden...")
    subprocess.run(["python3", wash_script, temp_raw, outfile], check=True)
    
    os.remove(temp_raw)
    print(f"Generated {outfile}")

print("All datasets generated successfully!")
