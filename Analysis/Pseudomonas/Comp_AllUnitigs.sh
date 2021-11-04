# Add path
raw="/Raw/Pseudomonas"
mkdir -p ${raw}
data="/Data/Pseudomonas/"
out="/Output/Pseudomonas/"
mkdir -p ${out}

# Download data
mkdir -p ${data}
cd ${data}
wget https://plmbox.math.cnrs.fr/seafhttp/files/3350a8aa-ff2e-435a-8aa8-8798682ea6cf/step1.tgz
tar -xvf step1.tgz
rm step1.tgz
cd /
echo "##########################################################" >> ${out}/all_unitigs.out
echo "Running unitigs_analysis" >> ${out}/all_unitigs.out
alpha=0.000000000001
for i in {12..1}
do 
    echo ${alpha} >> ${out}/all_unitigs.out
    echo "...Running unitigs_analysis" >> ${out}/all_unitigs.out
    python3 /src/AllUnitigs/OneStage_script.py --alpha=${alpha} -l $data \
        -o ${out}/AllUnitigs_${alpha} -P /Raw/Pseudomonas/PA_pop_kmeans2.txt \
        -v >> ${out}/testing.out 2>&1
    DBGWAS -strains ${raw}/pseudomonas_aeruginosa_full_dataset/strains -skip1 -skip2 \
        -nb-cores 2 -output ${out}/AllUnitigs_${alpha} -caldera -SFF p1.0 -nh 2 \
        >> ${out}/testing.out 2>&1
    echo "" >> ${out}/all_unitigs.out
    rm -rf ${out}/AllUnitigs_${alpha}/step*
    rm -rf ${out}/AllUnitigs_${alpha}/visualisations
    rm -rf ${out}/AllUnitigs_${alpha}/textualOutput/components
    alpha=$(echo $alpha | sed 's/^0.0/0./g')
done

Rscript --no-save --verbose /Analysis/Pseudomonas/Comp_Results.R \
    -l AllUnitigs_ >> ${out}/analyze.out 2>&1