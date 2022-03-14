out=/Output/Explo/
mkdir -p $out
raw=/Raw/Explo/
data=/Data/Explo/
mkdir -p $raw
mkdir -p $raw/Genomes
# Create data
R CMD BATCH /Analysis/Explo/generate.R ${out}/generate.out

# Run DBGWAS
echo "Running DBGWAS" >> ${out}/all_steps.out
DBGWAS -strains ${raw}/strains -nh 200 \
    -nb-cores 2 -output ${data} -SFF q1.0 \
    >> ${out}/all_steps.out 2>&1

cp -r ${data} ${out}/DBGWAS/
echo "Finished DBGWAS" >> ${out}/all_steps.out

python3 /src/CALDERA/bin/Pre-Process/toMajor.py -l ${data} >> ${out}/all_steps.out 2>&1
echo "Finished pre-processing" >> ${out}/all_steps.out

# Running the actual CALDERA script with one stage
alpha=0.1
for i in {1..15}
do  
    echo ${alpha} >> ${out}/all_steps.out
    echo "... Running CALDERA" >> ${out}/all_steps.out
    # All unitigs
    python3 /src/AllUnitigs/OneStage_script.py -l ${data} -o ${out}/AllUnitigs_${alpha}/ \
        -v -C 1 --alpha=${alpha} >> ${out}/caldera.out 2>&1
    DBGWAS -strains ${raw}/strains -skip1 -skip2 -caldera \
        -nb-cores 2 -output ${out}/AllUnitigs_${alpha}/ > \
        ${out}/AllUnitigs_${alpha}/step3.out 2>&1
    # Stages
    for stage in $(echo '1 2 3 5 10 15 20')
    do
        echo "... Stage "${stage} >> ${out}/all_steps.out
        python3 /src/CALDERA/bin/caldera-script -l ${data} -o ${out}/Stage${stage}_${alpha}/ -v -t 2 \
            -C 1 -s ${stage} --Lmax 20000 --alpha=${alpha} >> ${out}/caldera.out 2>&1
        DBGWAS -strains ${raw}/strains -skip1 -skip2 -caldera\
            -nb-cores 2 -output ${out}/Stage${stage}_${alpha}/ > \
            ${out}/Stage${stage}_${alpha}/step3.out 2>&1
    done 
    # All stages
    python3 /src/CALDERA/bin/caldera-script -l ${data} -o ${out}/AllStages_${alpha}/ -v -t 2 \
        -C 1 --batch_size=20000 --Lmax 20000 --alpha=${alpha} >> ${out}/caldera.out 2>&1
    DBGWAS -strains ${raw}/strains -skip1 -skip2 -caldera\
        -nb-cores 2 -output ${out}/AllStages_${alpha}/ -nh 200 > \
        ${out}/AllStages_${alpha}/step3.out 2>&1
    echo "... Finished CALDERA" >> ${out}/all_steps.out
    # Post-processing
    Rscript --no-save --verbose /Analysis/Explo/analyze_dbg.R \
        -a ${alpha} >> ${out}/analyze_dbg.out 2>&1 
    echo "... Finished All Steps" >> ${out}/all_steps.out
    echo "" >> ${out}/all_steps.out
        alpha=$(echo $alpha | sed 's/^0./0.0/g')
done

# Run pyseer
mkdir $out/pyseer
awk -v OFS='\t' '{print $1, $2}' $raw/strains > $raw/phenotypes.tsv
awk '{print $1, $3}' $raw/strains | tail -n +2 > $raw/fsm_file_list.txt
/scr/fsm-lite/fsm-lite -l $raw/fsm_file_list.txt -m 30 -M 31 -s 2 -S 610 -v -t kmers | gzip -c - > ${raw}/kmers.txt.gz
/scr/mash/mash sketch -s 10000 -o ${raw}/samples ${raw}/Genomes/*
/scr/mash/mash dist ${raw}/samples.msh ${raw}/samples.msh | python /scr/pyseer/square_mash-runner.py > ${raw}/mash.tsv
python /scr/pyseer/pyseer-runner.py --phenotypes $raw/phenotypes.tsv --kmers $raw/kmers.txt.gz \
    --no-distances --cpu 4 > $out/pyseer/pyseer.assoc
python /scr/pyseer/pyseer-runner.py --phenotypes $raw/phenotypes.tsv --kmers $raw/kmers.txt.gz \
    --distances ${raw}/mash.tsv --cpu 4 > $out/pyseer/pyseer_dist.assoc
python /scr/pyseer/pyseer-runner.py --phenotypes $raw/phenotypes.tsv --kmers $raw/kmers.txt.gz \
    --no-distances --cpu 4 --wg enet > $out/pyseer/pyseer_wg.assoc
awk '{print $1}' $raw/phenotypes.tsv | tail -n +2 > $raw/samples.txt
python /scr/pyseer/similarity-runner.py --kmers $raw/kmers.txt.gz $raw/samples.txt > ${raw}/sim.tsv
python /scr/pyseer/pyseer-runner.py --phenotypes $raw/phenotypes.tsv --kmers $raw/kmers.txt.gz \
    --similarity ${raw}/sim.tsv --cpu 4 --lmm > $out/pyseer/pyseer_lmm.assoc
cp ${raw}/kmers.txt.gz $out/pyseer/kmers.txt.gz

# # Klover 
# pyenv global 2.7.18
# awk -v OFS='\t' '{print $1, $3}' $raw/strains | tail -n +2> $raw/genomic-data.tsv
# awk -v OFS='\t' '{print $1, $2}' $raw/strains | tail -n +2> $raw/phenotypes.tsv
# kover dataset create from-contigs --phenotype-description "Resistance example" \
#     --genomic-data $raw/genomic-data.tsv --phenotype-metadata $raw/phenotypes.tsv \
#     --output example.kover
# kover dataset split --dataset example.kover --id example_split \
#     --train-size 0.666 --folds 5 --random-seed 72 --progress
# kover learn scm --dataset example.kover --split example_split \
#     --model-type conjunction disjunction --p 0.1 1.0 10.0 \
#     --hp-choice cv --n-cpu 2 --progress

# Final merge of results
R CMD BATCH /Analysis/Explo/analyze_kmer.R ${out}/analyze_pyseer.out

# Clean
cp -r ${out}/AllStages_0.000000000000001 ${out}/Setup
R CMD BATCH /Analysis/Explo/Viz.R ${out}/viz.out
DBGWAS -strains ${raw}/strains -skip1 -skip2 -caldera\
    -nb-cores 2 -output ${out}/Setup -SFF p2 > \
    ${out}/all_steps.out 2>&1

cp -r ${out}/AllUnitigs_0.00000001 ${out}/AllUnitigs
cp -r ${out}/AllStages_0.00000001 ${out}/AllStages
rm -rf ${out}/Stage*
rm -rf ${out}/AllUnitigs_*
rm -rf ${out}/AllStages_*
echo "Done" >> ${out}/all_steps.out