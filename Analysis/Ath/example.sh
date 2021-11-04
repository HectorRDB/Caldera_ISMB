data="Data/Ath/"
mkdir -p ${data}
out="Output/Ath/"
mkdir -p ${out}

# downloads kegg pathways, extracts a list of edges
# Running the actual CALDERA script
R CMD BATCH Analysis/Ath/cleanFormat.R ${out}/ath.out
python3 CALDERA/Pre-process/toMajor.py -l ${data} >> ${out}/ath.out 2>&1

python3 CALDERA/Scripts/Core.py -l ${data} -o ${out} --verbose True -t 1\
    -C 1 --Lmax 50000 --kmax 10000000 --alpha 0.05>> ${out}/ath.out 2>&1
