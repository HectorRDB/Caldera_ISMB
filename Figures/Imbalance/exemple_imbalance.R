library(phylogram)
library(here)
library(dendextend)
library(magrittr)
pheno <- read.table( here("Data", "Amikacin", "step1", "bugwas_input.id_phenotype"), header = TRUE)
tree <- read.dendrogram(file = here("data", "strains.newick"))
dist <- as.matrix(cophenetic(tree))
pop <- cutree(tree, k = 3)
pop <- pop[pheno$ID]
write.table(pop, here("data", "PA_pop.txt"), row.names = FALSE, col.names = FALSE)
table(pop)
coords <- cmdscale(dist, k = 2)
wss <- sapply(2:20, function(k){
  kmeans(coords, centers = k, nstart = 10)$tot.withinss
})
plot(2:20, wss)
# We choose k = 3
cl <- kmeans(coords, centers = 3, nstart = 10)
(table(cl$cluster) / length(cl$cluster)) %>%
  round(2) %>%
  min() %>%
  `*`(100) %>%
  paste0("%")
  
