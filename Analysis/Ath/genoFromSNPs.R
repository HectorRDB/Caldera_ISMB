## 1001 Genomes data obtained from https://easygwas.ethz.ch/down/1/

dir <- "/Raw/Ath/geneToPat/"
all.files <- list.files(dir)
pat <- matrix(0, nrow = length(all.files), ncol = 1135)
rownames(pat) <- all.files
info <- file.info(sapply(all.files, FUN = function(ss) file.path(dir, ss)))

## Empty files are genes with no mutation, pattern should stay all 0
all.files <- all.files[info$size != 0]

for (ff in all.files) {
  snp.tab <- read.table(file.path(dir, ff), stringsAsFactors = FALSE)

  ## sanity check
  missing <- sapply(snp.tab[[8]], FUN = function(ss) (length(grep(ff, ss)) == 0))
  if (any(missing)) {
    stop(sprintf("Gene %s missing from some of its assigned SNPs", ff))
  }

  m <- as.matrix(snp.tab[-(1:9)])

  pat[ff, ] <- apply(m, 2, FUN = function(v) 1 * any(v != "0/0"))
}

write.table(rownames(pat),  row.names = FALSE, col.names = FALSE, quote = FALSE,
            file = "/Data/Ath/step1/ath.genes")

ath.names <- strsplit(readLines("/Data/Ath/ath.names", 1), "\t")[[1]][-(1:9)]
pat <- cbind(1:nrow(pat), pat)
colnames(pat) <- c("ps", ath.names)

write.table(pat, row.names = FALSE, col.names = TRUE, quote = FALSE, 
            file = "/Data/Ath/step1/ath.patterns")

## phenotype
pheno <- read.table("/Raw/Ath/phenotypes.pheno", header = TRUE)
pheno.names <- sapply(pheno[["FID"]], FUN = function(ss) paste(ss, ss, sep = "_"))
DTF.pheno <- cbind(pheno.names, 1 * (pheno[["DTF1"]] > median(pheno[["DTF1"]], na.rm = TRUE)))
colnames(DTF.pheno) <- c("ID", "pheno")
write.table(DTF.pheno, row.names = FALSE, col.names = TRUE, quote = FALSE,
            file = "/Data/Ath/step1/ath.pheno")
