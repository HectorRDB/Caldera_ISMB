out=/Output/Akkermansia/
mkdir -p $out
raw=/Raw/Akkermansia/
mkdir -p ${raw}
data=/Data/Akkermansia/
# Get metadata
wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-021-02427-7/MediaObjects/13059_2021_2427_MOESM1_ESM.xlsx
mv 13059_2021_2427_MOESM1_ESM.xlsx ${raw}/metadata.xlsx
R CMD BATCH /Analysis/Akkermansia/clean_metadata.R ${out}/clean_metadata.out
# Get genomes
url=https://zenodo.org/record/5018705/files/genomes.tar.gz?download=0
mkdir -p ${raw}/Genomes/
wget -O ${raw}/Genomes/genomes.tar.gz ${url}
cd ${raw}/Genomes/ && tar -xvf genomes.tar.gz && rm genomes.tar.gz
# Remove non-useful files
for file in $(ls)
do
    if ! grep -qxFe "$file" ${raw}/strains.tsv; then
        echo "Deleting: $file"
        rm "$file"
    fi
done

# Getting the reference database
sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
PATH=${PATH}:${HOME}/edirect
esearch -db protein -query "refseq[filter] AND txid239935[Organism:exp]" | efetch -format fasta > ${raw}/Akkermansia.refseq.sequences
cp ${raw}/Akkermansia.refseq.sequences ${out}

# Running CALDERA
echo "Running DBGWAS" >> ${out}/all_steps.out
cd /
DBGWAS -strains ${raw}/strains -SFF q0.05 -nh 2\
    -nb-cores 4 -output ${data} -pt-db ${raw}/Akkermansia.refseq.sequences \
    >> ${out}/all_steps.out 2>&1
cp -r ${data}/visualisations ${out}/DBGWAS/visualisations
cp -r ${data}/textualOutput ${out}/DBGWAS/textualOutput
echo "Finished DBGWAS" >> ${out}/all_steps.out
python3 /src/CALDERA/bin/Pre-Process/toMajor.py -l ${data} > ${out}/caldera.out 2>&1
echo "Finished pre-processing" >> ${out}/all_steps.out
echo "" >> ${out}/all_steps.out
echo "##########################################################" >> ${out}/all_steps.out
echo "Running CALDERA" >> ${out}/all_steps.out
python3 /src/CALDERA/bin/caldera-script -l ${data} -o ${out}/CALDERA -v -t 4 \
    -s 5 --alpha=0 --Lmax 2000 --batch_size=50000\
    >> ${out}/caldera.out 2>&1
DBGWAS -strains ${raw}/strains -skip1 -skip2 -caldera \
    -pt-db ${raw}/Akkermansia.refseq.sequences -output ${out}/CALDERA > ${out}/viz.out 2>&1
echo "Finished CALDERA" >> ${out}/all_steps.out
echo "Finished All Steps" >> ${out}/all_steps.out

