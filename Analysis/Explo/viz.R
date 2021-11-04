library(dplyr)
library(ggplot2)
library(readr)
library(here)
library(stringr)
library(tidyr)
options(scipen = 999)
pheno_nodes <- read_delim(here("Output", "Explo", "DBGWAS", "step1", "bugwas_input.all_rows.binary"),
                          delim = " ", col_names = c("NodeId", paste0("Sample", 1:100)),
                          skip = 1) %>%
  pivot_longer(-NodeId, names_to = "Sample", values_to = "Present") %>%
  mutate(pheno = if_else(Sample %in% paste0("Sample", 1:50), 1, 0)) %>%
  filter(Present == 1) %>%
  group_by(NodeId) %>%
  summarise(res_node = mean(pheno == 0) == 0) 

bugwas <- c(
    paste0(pheno_nodes$NodeId, collapse = " "),
    paste0(pheno_nodes$NodeId[pheno_nodes$res_node], collapse = " ")
)
write_lines(bugwas, here("Output", "Explo", "Setup", "step1", "bugwas_input.unique_rows_to_all_rows.binary"))
gemma <- pheno_nodes %>%
    mutate(Pheno = if_else(res_node, 1, 0)) %>%
    select(NodeId, Pheno)
write_tsv(gemma, here("Output", "Explo", "Setup", "step1", "gemma_input.unitig_to_pattern.binary"),
    col_names = FALSE)
nodes_to_css <- pheno_nodes %>%
    mutate(ccs = if_else(res_node, "0 1", "0")) %>%
    select(ccs)
write_tsv(gemma, here("Output", "Explo", "Setup", "step2", "nodes_to_css"),
    col_names = FALSE)
write_tsv(data.frame("weight" = rep(1, nrow(gemma))), here("Output", "Explo", "Setup", "step1", "weight_correction"),
    col_names = FALSE)
pats <- data.frame("ID" = c(0, 1), "pvals" = c("0.001", "0.0000001"),
                   "statistics" = c(70, 70), "Pheno0" = c(50, 0), "Pheno1" = c(50, 50),
                   "PhenoNA" = c(0, 0))
write_tsv(pats, here("Output", "Explo", "Setup", "step2", "patterns.txt"))

