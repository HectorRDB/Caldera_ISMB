raw="/Raw/Ath/"
mkdir -p ${raw}
data="/Data/Ath/"
mkdir -p ${data}
out="/Output/Ath/"
mkdir -p ${out}

# downloads kegg pathways, extracts a list of edges
# and a list of genes that are present in at least one pathway.
cd $raw
mkdir -p ${raw}/kgmlAth
mkdir -p ${data}/step1/edges
R CMD BATCH /Analysis/Ath/getAthPathways.R ${out}/getAthPathways.out
# Get the genomes
wget https://easygwas.ethz.ch/down/dataset/download/7/
unzip ${raw}/index.html
rm ${raw}/index.html

plink --ped ${raw}/genotype.ped --map ${raw}/genotype.map --recode vcf --out ${raw}/athV
java -Xmx4g -jar snpEff.jar closest Arabidopsis_thaliana athV.vcf > athV.ann.vcf
rm genotype.ped genotype.map athV.vcf

# neurips-getPatterns.sh extracts the genotypes of all SNPs corresponding to
# each gene, creating one file per gene in a geneToPat directory.
mkdir -p ${raw}/geneToPat
/Analysis/Ath/getPatterns.sh genes
head -n 573 /Raw/Ath/geneToPat/AT5G46570 > /Raw/Ath/geneToPat/AT5G46570
head -n 14 ${raw}/athV.ann.vcf | tail -n 1 > ${data}/ath.names

# uilds an ath.patterns file with one row per gene and one column per sample,
# with a 1 when the gene has at least one mutation (even heterozygous) for a
# linked SNP. The word "ps" should be manually added at the beginning of the
# first row to match the DBGWAS output format.
R CMD BATCH /Analysis/Ath/genoFromSNPs.R ${out}/genoFromSNPs.out
rm athV.ann.vcf

# Pre-process to make compatible with data
R CMD BATCH /Analysis/Ath/cleanFormat.R ${out}/finish-pre-process.out
python3 /CALDERA/Pre-process/toMajor.py -l ${data}

# Running the actual CALDERA script
python3 /CALDERA/Scripts/Core.py -l ${data} -o ${out} --verbose True -t 4 \
    -C 1 --Lmax 50000 --kmax 10000000 --alpha 0.05 > ${out}/ath.out 2>&1

R CMD BATCH /Analysis/Ath/post-process.R ${out}/post-process.out
