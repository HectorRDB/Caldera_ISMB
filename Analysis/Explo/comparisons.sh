out=/Output/Explo/
mkdir -p $out
raw=/Raw/Explo/
data=/Data/Explo/
mkdir -p $raw
mkdir -p $raw/Genomes
# Create data
R CMD BATCH /Analysis/Explo/generate.R ${out}/generate.out

echo "Running DBGWAS" >> ${out}/all_steps.out
DBGWAS -strains ${raw}/strains -nh 200\
    -nb-cores 2 -output ${data} -SFF q1.0\
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
    Rscript --no-save --verbose /Analysis/Explo/analyze.R \
        -a ${alpha} >> ${out}/analyze.out 2>&1 
    echo "... Finished All Steps" >> ${out}/all_steps.out
    echo "" >> ${out}/all_steps.out
        alpha=$(echo $alpha | sed 's/^0./0.0/g')
done

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