# Add path
raw="/Raw/Pseudomonas"
mkdir -p ${raw}
data="/Data/Pseudomonas/"
out="/Output/Pseudomonas/"
mkdir -p ${out}

# Download data
mkdir -p ${data}
cd ${data}
wget https://plmbox.math.cnrs.fr/seafhttp/files/0cbc10b6-1748-414e-a10e-4e6c06003008/step1.tgz
tar -xvf step1.tgz
rm step1.tgz
cd /

# Running the actual COIN script
date >> ${out}/all_steps
python3 /src/COIN/bin/coin-script -l ${data} -o ${out} --verbose True  \
    -P /Raw/Pseudomonas/PA_pop_kmeans2.txt  --Lmax 2000 -s 7 --alpha 0.000001 \
    >> ${out}/COIN.out 2>&1
date >> ${out}/all_steps
echo "Finished COIN" >> ${out}/all_steps
