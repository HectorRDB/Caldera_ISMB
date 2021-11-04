#!/usr/bin/env Rscript
suppressWarnings(library(here))
suppressWarnings(library(dplyr))
suppressWarnings(library(readr))
suppressWarnings(library(stringr))
suppressWarnings(library(tidyr))
suppressWarnings(library(optparse))
# Arguments for R Script ----
option_list <- list(
  make_option(c("-a", "--alpha"),
              action = "store", default = NA, type = "double",
              help = "The value of alpha"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$a)) {
  alpha <- opt$a
  cat("The value of alpha is ", alpha, "\n")
} else {
  stop("Missing alpha argument")
}

options(scipen = 999)
# Load input
pheno_nodes <- read_delim(here("Output", "Explo", "DBGWAS", "step1", "bugwas_input.all_rows.binary"),
                          delim = " ", col_names = c("NodeId", paste0("Sample", 1:100)),
                          skip = 1) %>%
  pivot_longer(-NodeId, names_to = "Sample", values_to = "Present") %>%
  mutate(pheno = if_else(Sample %in% paste0("Sample", 1:50), 1, 0)) %>%
  filter(Present == 1) %>%
  group_by(NodeId) %>%
  summarise(res_node = mean(pheno == 0) == 0) %>%
  mutate(NodeId = paste0("n", NodeId))

# Metrics
metrics <- function(result, pheno_nodes) {
  if (is.null(result)) {
    return(data.frame("TN" = sum(!pheno_nodes$res_node),
                      "FN" = sum(pheno_nodes$res_node), "TP" = 0, "FP" = 0))
  }
  df <- pheno_nodes %>%
    mutate(positive = NodeId %in% result$NodeId)
  return(data.frame("TN" = sum(!df$positive & !df$res_node),
                    "FN" = sum(!df$positive & df$res_node),
                    "TP" = sum(df$positive & df$res_node),
                    "FP" = sum(df$positive & !df$res_node)))
}

read_sig_nodes <- function(folder, alpha) {
  file <- here("Output", "Explo", paste0(folder, alpha), "step3", "significant_unitigs.txt")
  if (file.exists(file)) {
    df <- read_tsv(file, col_names = "NodeId") %>%
      mutate(NodeId = paste0("n", NodeId))
  } else {
    df <- NULL
  }
  return(df)
}

res_methods <- bind_rows(
  "All Unitigs" = metrics(read_sig_nodes("AllUnitigs_", alpha), pheno_nodes),
  "All Stages" = metrics(read_sig_nodes("AllStages_", alpha), pheno_nodes),
  "DBGWAS" = metrics(
    read_tsv(here("Output", "Explo", "DBGWAS", "textualOutput", "all_comps_nodes_info.tsv")) %>%
      filter(`p-value`<= alpha / nrow(.)) %>%
      select(NodeId),
    pheno_nodes),
  .id = "Analysis"
)

stages <- lapply(c(1, 2, 3, 5, 10, 15, 20), function(stage) {
  metrics(read_sig_nodes(paste0("Stage", stage, "_"), alpha), pheno_nodes) %>%
    mutate(Analysis = paste0("Stage", stage)) %>%
    return()
}) %>% bind_rows()

write_csv(bind_rows(res_methods, stages),
          here("Output", "Explo", "Results", paste0("results_", -log10(alpha), ".csv")))
