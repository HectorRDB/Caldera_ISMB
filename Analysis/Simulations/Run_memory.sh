# Add path
dir="/Output/Simulations/Speed"
mkdir -p $dir

python3 /Analysis/Simulations/Memory.py > ${dir}/out_mem 2>&1
