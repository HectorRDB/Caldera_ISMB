out="/Output/Simulations/Imbalance"
mkdir -p $out
# Add path

python3 /Analysis/Simulations/Imbalance.py > ${out}/out 2>&1
