library(here)

X <- read.table(here("Data", "Pseudomonas", "step1", "bugwas_input.unique_rows.binary"), header=TRUE, sep=" ")
rownames(X) <- as.character(X[, 1])
X <- X[, -1]

## PCA: scale matrix, compute SVD
scX  <-  t(scale(t(X), scale=FALSE))
svdscX <- svd(scX)

## Gap after 1st or 4 first PCs
pdf(file = here("Figures", "Pseudomonas", 'PA-pseudomonas-PC-variance.pdf'))
plot(svdscX$d^2/(sum(svdscX$d^2)), pch=19, xlab='Principal components', ylab = 'Percentage of explained variance')
dev.off()

km.res <- list()
klist <- c(2, 5)
for(ii in 1:length(klist)){
    print(klist[ii])
    km.res[[ii]] <- kmeans(t(X), klist[ii], nstart=10)
}
    
## Cluster 1 (red, 9 samples) could be merged with cluster 4 (yellow)
pdf(file = here("Figures", "Pseudomonas", 'clusters-5.pdf'))
colVec <- c('red', 'green', 'blue', 'yellow', 'orange')
plot(svdscX$v[, 1], svdscX$v[, 2], pch=19, col=colVec[km.res[[3]]$cluster], xlab='PC1', ylab='PC2')
#legend('topright', legend=as.character(klist), col=colVec, pch=19)
dev.off()

## 2-cluster pop
write.table(file = here("Raw", "Pseudomonas", 'PA_pop_kmeans2.txt'),
    km.res[[1]]$cluster - 1, quote = FALSE, row.names = FALSE, col.names = FALSE)

## 4-cluster version (5-cluster k-means, then merge #4 with #1)
pop4 <- km.res[[3]]$cluster
pop4[pop4 == 4]  <- 1 # Now contains 1, 2, 3, 5
pop4[pop4 == 5]  <- 4 # Now contains 1, 2, 3, 4

## Visual sanity check
## plot(svdscX$v[, 1], svdscX$v[, 2], pch=19, col=colVec[pop4], xlab='PC1', ylab='PC2')
write.table(file = here("Raw", "Pseudomonas", 'PA_pop_kmeans4.txt'), 
    pop4 - 1, quote = FALSE, row.names = FALSE, col.names = FALSE)
