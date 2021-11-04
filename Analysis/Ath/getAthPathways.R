library(KEGGREST)
library(DEGraph)
library(here)

## This script downloads kegg pathways for Arabidopsis thaliana (ath),
## and extracts a list of edges for each pathway.

## It requires a ./kgmlAth directory

## Get ath pathway list from kegg
ath.list <- keggList("pathway", "ath")
pathway.list <- sapply(ath.list, FUN = function(ss) strsplit(ss, " - ")[[1]][1])

## Save pathway names
write.table(pathway.list, row.name = TRUE, col.name = FALSE, quote = FALSE,
            file = "/Data/Ath/step1/ath.pathways")

## Download kgml files for all pathways
pathway.codes <- sub("path:", "", names(ath.list))
for (cc in pathway.codes) {
  write(keggGet(cc, "kgml"), file = file.path("kgmlAth", paste(cc, "kgml", sep = ".")))
}

## Build graph objects from kgml files
grList <- getKEGGPathways(path = "kgmlAth", organism = "ath", 
                          pattern = "^ath([0-9]+).kgml", verbose = TRUE)

## Extract edge list from a graphNEL object (returned by
## getKEGGPathways), as required by CALDERA
listEdges <- function(gnel, return.table = FALSE) {
  nodesv <- gnel@nodes
  edge.list <- gnel@edgeL
  edge.sources <- sapply(names(edge.list), 
                         FUN = function(ss) strsplit(ss, "ath:")[[1]][2])
  edge.table <- c()
  for (e in 1:length(edge.list)) {
    if (length(edge.list[[e]]$edges) != 0) {
      add.edges <- sapply(nodesv[edge.list[[e]]$edges], 
                          FUN = function(ss) {
                            paste(edge.sources[e], 
              sapply(ss, FUN = function(sss) strsplit(sss, "ath:")[[1]][2]))})
      names(add.edges) <- NULL
      edge.table <- c(edge.table, add.edges)
    }
  }
  write.table(edge.table, 
              file = paste0("/Data/Ath/step1/edges/", gnel@graphData$info@number), 
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  if (return.table) {
    return(edge.table)
  }
}
full.table <- unlist(lapply(grList, listEdges, TRUE))

## Build a list of all genes involved in at least one of the downloaded pathways
full.list <- as.vector(sapply(full.table, FUN = function(ss) strsplit(ss, " ")[[1]]))
unique.list <- unique(full.list)
write.table(unique.list, file = "/Raw/Ath/genes", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)
