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

start=`date +%s`
mkdir $out/pyseer
awk -v OFS='\t' '{print $1, $2}' $raw//pseudomonas_aeruginosa_full_dataset/strains > $raw/phenotypes.tsv
awk '{print $1, $3}' $raw/pseudomonas_aeruginosa_full_dataset/strains | tail -n +2 > $raw/fsm_file_list.txt
cd $raw
/scr/fsm-lite/fsm-lite -l $raw/fsm_file_list.txt -m 30 -M 31 -s 20 -S 260 -v -t kmers | gzip -c - > ${raw}/kmers.txt.gz
awk '{print $1}' $raw/phenotypes.tsv | tail -n +2 > $raw/samples.txt
python /scr/pyseer/pyseer-runner.py --phenotypes $raw/phenotypes.tsv --kmers $raw/kmers.txt.gz \
    --no-distances --cpu 4 --wg enet > $out/pyseer/pyseer_wg.assoc
end=`date +%s`
runtime=$((end-start))
echo $runtime