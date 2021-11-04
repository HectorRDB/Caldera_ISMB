out=/Example/Output
mkdir -p $out
raw=/Example/Raw
data=/Example/Data/
mkdir -p $raw
mkdir -p $raw/Genomes
# Create data
R CMD BATCH Analysis/Example/generate.R ${out}/generate.out

# Current process:
DBGWAS -strains ${raw}/strains -only1 -nb-cores 2 -output ${data} 
python3 src/CALDERA/bin/Pre-Process/toMajor.py -l ${data}
python3 src/CALDERA/bin/caldera-script -l ${data} -o ${out}/caldera/ -v -t 2 -C 1
DBGWAS -strains ${raw}/strains -skip1 -skip2  -nb-cores 2 \
    -output ${out}/caldera/ -nh 200 -caldera 




