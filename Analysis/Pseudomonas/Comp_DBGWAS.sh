# Add path
raw="/Raw/Pseudomonas"
mkdir -p ${raw}
data="/Data/Pseudomonas/"
out="/Output/Pseudomonas/"
mkdir -p ${out}

# Download data
wget https://www.dropbox.com/s/0g1llvdbfv1jys6/pseudomonas_aeruginosa_full_dataset.zip?dl=1
mv 'pseudomonas_aeruginosa_full_dataset.zip?dl=1' \
    ${raw}/pseudomonas_aeruginosa_full_dataset.zip
cd ${raw}
pwd >> ${out}/all_steps.out
unzip ${raw}/pseudomonas_aeruginosa_full_dataset.zip
rm pseudomonas_aeruginosa_full_dataset.zip
mv ${raw}/pseudomonas_aeruginosa_full_dataset/strains.newick ${raw}/

# Run DBGWAS ------
echo "Running DBGWAS" >> ${out}/all_steps.out
DBGWAS -strains ${raw}/pseudomonas_aeruginosa_full_dataset/strains \
    -nb-cores 6 -output ${data} \
    -no-preview >> ${out}/all_steps.out 2>&1
alpha=0.000000000001
n=$(wc -l ${data}/step2/patterns.txt | awk '{print $1}')
calc(){ awk "BEGIN { print "$*" }"; }
for i in {1..12}
do 
    echo ${alpha} >> ${out}/all_steps.out
    cut=$(calc $alpha/$n)
    mkdir ${out}/DBGWAS_p_${alpha}/
    DBGWAS -strains ${raw}/pseudomonas_aeruginosa_full_dataset/strains \
        -nb-cores 6 -output ${data} -skip1 -skip2 \
        -SFF p${cut} -no-preview -nh 2 \
        >> ${out}/all_steps.out 2>&1
    mv ${data}/textualOutput ${out}/DBGWAS_p_${alpha}/
    mkdir ${out}/DBGWAS_q_${alpha}/
    DBGWAS -strains ${raw}/pseudomonas_aeruginosa_full_dataset/strains \
        -nb-cores 6 -output ${data} -skip1 -skip2 \
        -SFF q${alpha} -no-preview -nh 2  \
       >> ${out}/all_steps.out 2>&1
    echo "" >> ${out}/all_steps.out
    mv ${data}/textualOutput ${out}/DBGWAS_q_${alpha}/
    alpha=$(echo $alpha | sed 's/^0.0/0./g')
done

Rscript --no-save --verbose /Analysis/Pseudomonas/Comp_Results.R \
    -a ${alpha} -l DBGWAS_ -d TRUE >> ${out}/analyze.out 2>&1

DBGWAS -strains ${raw}/pseudomonas_aeruginosa_full_dataset/strains \
    -nb-cores 6 -output ${data} -skip1 skip2 -nh 2\
    -no-preview >> ${out}/all_steps.out 2>&1
mkdir ${out}/DBGWAS/
cp -r ${data}/visualisations ${out}/DBGWAS/
cp -r ${data}/textualOutput ${out}/DBGWAS/
