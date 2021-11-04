echo "Getting data" >> ${out}/all_steps.out
# Add path ----
raw="/Raw/Pseudomonas"
mkdir -p ${raw}
data="/Data/Pseudomonas/"
out="/Output/Pseudomonas/"
mkdir -p ${out}
# Download data ----
wget https://www.dropbox.com/s/0g1llvdbfv1jys6/pseudomonas_aeruginosa_full_dataset.zip?dl=1
mv 'pseudomonas_aeruginosa_full_dataset.zip?dl=1' \
    ${raw}/pseudomonas_aeruginosa_full_dataset.zip
cd ${raw}
pwd >> ${out}/all_steps.out
unzip ${raw}/pseudomonas_aeruginosa_full_dataset.zip
rm pseudomonas_aeruginosa_full_dataset.zip
mv ${raw}/pseudomonas_aeruginosa_full_dataset/strains.newick ${raw}/
echo "##########################################################" >> ${out}/all_steps.out

# Run DBGWAS ------
echo "Running DBGWAS" >> ${out}/all_steps.out
DBGWAS -strains ${raw}/pseudomonas_aeruginosa_full_dataset/strains \
    -nb-cores 6 -output ${data} \
    -no-preview >> ${out}/all_steps.out 2>&1
echo "Finished DBGWAS" >> ${out}/all_steps.out
echo "##########################################################" >> ${out}/all_steps.out

# Preprocess ------
python3 /src/CALDERA/bin/Pre-Process/toMajor.py -l ${data} > ${out}/caldera.out 2>&1
R CMD BATCH /Analysis/Pseudomonas/makePop.R ${out}/makePop.out
echo "Finished pre-processing" >> ${out}/all_steps.out
echo "" >> ${out}/all_steps.out
echo "##########################################################" >> ${out}/all_steps.out

# Caldera ------
echo "Running CALDERA" >> ${out}/all_steps.out
python3 /src/CALDERA/bin/caldera-script -l ${data} -o ${out} -v -t 6 \
    -P /Raw/Pseudomonas/PA_pop_kmeans2.txt -s 7 --alpha=0 --Lmax 2000 >> ${out}/caldera.out 2>&1
DBGWAS -strains ${raw}/pseudomonas_aeruginosa_full_dataset/strains -skip1 -skip2 \
    -nb-cores 2 -output ${out}/caldera/ -caldera >> ${out}/caldera.out 2>&1
echo "Finished CALDERA" >> ${out}/all_steps.out
echo "##########################################################" >> ${out}/all_steps.out

echo "Running CALDERA with many alpha values" >> ${out}/all_steps.out
alpha=0.000000000001
for i in {1..10}
do 
    echo ${alpha} >> ${out}/all_steps.out
    python3 /src/CALDERA/bin/caldera-script -l ${data} -o ${out}/Caldera_S7_${alpha} -v -t 6 \
        -P /Raw/Pseudomonas/PA_pop_kmeans2.txt -s 7 --alpha=${alpha} \
        --Lmax 2000 >> ${out}/caldera.out 2>&1
    DBGWAS -strains ${raw}/pseudomonas_aeruginosa_full_dataset/strains -skip1 -skip2 \
        -nb-cores 2 -output ${out}/Caldera_S7_${alpha} -caldera \
        -SFF p1.0 -nh 2 >> ${out}/caldera.out 2>&1
    echo "" >> ${out}/all_steps.out
    rm -rf ${out}/Caldera_S7_${alpha}/step*
    rm -rf ${out}/Caldera_S7_${alpha}/visualisations
    rm -rf ${out}/Caldera_S7_${alpha}/textualOutput/components
    alpha=$(echo $alpha | sed 's/^0.0/0./g')
done

Rscript --no-save --verbose /Analysis/Pseudomonas/Comp_Results.R \
    -a ${alpha} -l Caldera_S7_ >> ${out}/analyze.out 2>&1