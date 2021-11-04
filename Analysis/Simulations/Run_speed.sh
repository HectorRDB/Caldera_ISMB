# Add path
dir="/Output/Simulations/Speed"
mkdir -p $dir

# /Analysis/Simulations/Speed.py -i 1  > ${dir}/out_1 2>&1
# /Analysis/Simulations/Speed.py -i 2 -n 50 > ${dir}/out_2 2>&1
# /Analysis/Simulations/Speed.py -i 3 -p 0.2 > ${dir}/out_3 2>&1
# /Analysis/Simulations/Speed.py -i 4 -p 0.2 --alpha=0.0001 > ${dir}/out_4 2>&1
# /Analysis/Simulations/Speed.py -i 5 -C 2 -n 50 > ${dir}/out_5 2>&1
/Analysis/Simulations/Speed.py -i 2 -n 50 --long > ${dir}/out_6 2>&1