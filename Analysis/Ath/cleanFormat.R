library(here)
library(dplyr)
Patterns <- read.table(here("Data", "Ath", "step1", "ath.patterns"), sep = "",
                       header = TRUE, stringsAsFactors = FALSE)
colnames(Patterns) <- sub("X", "", colnames(Patterns))
Pheno <- read.table(here("Data", "Ath", "step1", "ath.pheno"), header = TRUE,
                    stringsAsFactors = FALSE)
Pheno <- Pheno[!is.na(Pheno$pheno), ]
Patterns <- Patterns[, c("ps", Pheno$ID)]
Patterns$ps <- factor(Patterns$ps, levels = Patterns$ps)

edges <- list.files(here("Data", "Ath", "step1", "edges"))
edges <- lapply(edges, function(edge){
  df <- read.table(here("Data", "Ath", "step1", "edges", edge), stringsAsFactors = F,
                   col.names = c("from", "to")) %>%
    mutate(from = factor(from, levels = levels(Patterns$ps)),
           to = factor(to, levels = levels(Patterns$ps))) %>%
    mutate(from = as.numeric(from) - 1, to = as.numeric(to) - 1)
  return(df)
}) %>% 
  bind_rows() %>%
  distinct()
  
write.table(edges, file = here("Data", "Ath", "step1", "graph.edges.dbg"), 
            row.names = F, col.names = F, sep = "\t")
write.table(rep(1, nrow(Patterns)), here("Data", "Ath", "step1", "graph.nodes"),
            col.names = F, sep = "  ")
write.table(rep(1, nrow(Patterns)), col.names = F, row.names = F,
            file = here("Data", "Ath", "step1", "weight_correction"))
write.table(Pheno, here("Data", "Ath", "step1", "bugwas_input.id_phenotype"), 
            row.names = F)

Numbers <- Patterns %>%
  mutate(number = as.numeric(ps) - 1,
         ps = as.character(ps)) %>%
  select(ps, number)
write.table(Numbers, here("Output", "Ath", "numbers.matching"),
            row.names = F,  sep = "\t")
Patterns <- Patterns %>%
  mutate(ps = as.numeric(ps) - 1)
write.table(Patterns, row.names = F,  sep = "\t",
            file = here("Data", "Ath", "step1", "bugwas_input.all_rows.binary"))

