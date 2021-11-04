#!/usr/bin/env Rscript
suppressWarnings(library(here))
suppressWarnings(library(dplyr))
suppressWarnings(library(readr))
suppressWarnings(library(stringr))
suppressWarnings(library(purrr))
suppressWarnings(library(tidyr))
suppressWarnings(library(optparse))
# Arguments for R Script ----
option_list <- list(
  make_option(c("-l", "--loc"),
              action = "store", default = NA, type = "character",
              help = "Default floder name"
  ),
  make_option(c("-d", "--dbgwas"),
              action = "store", default = TRUE, type = "logical",
              help = "Are we running on DBGWAS"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
} else {
  stop("Missing loc argument")
}

# Load input
aac6 <- read_lines(here("Raw", "Amikacin", "aac-significant-unitigs.txt"))
plasmid <- read_lines(here("Raw", "Amikacin", "plasmid-all-unitigs.txt"))

# Analyze function
metrics <- function(loc, alpha) {
  print(loc)
  file <- here("Output", "Amikacin", loc, "textualOutput", "all_comps_nodes_info.tsv")
  if (!file.exists(file)) {
    res <- data.frame("AAC6" = 0, "Plasmid" = 0, "Not_in_either" = 0,
                      "Rank_AAC6" = -1, "Rank_plasmid" = -1,
                      "n_comp_plasmid" = -1)
    return(res)
  }
  nodes <- read_tsv(file)
  if (nrow(nodes) == 0) {
    res <- data.frame("AAC6" = 0, "Plasmid" = 0, "Not_in_either" = 0,
                      "Rank_AAC6" = -1, "Rank_plasmid" = -1,
                      "n_comp_plasmid" = -1)
    return(res)
  }
  if (opt$d) {
    n <- min(nodes$`q-Value`)/ min(nodes$`p-value`) 
    nodes <- nodes %>%
      filter(`Significant?` == "Yes") %>%
      rename("CCSPvalue" = `p-value`)
  } else {
    nodes <- nodes %>%
      filter(CCSId != -1)
  }
  if (nrow(nodes) == 0) {
    res <- data.frame("AAC6" = 0, "Plasmid" = 0, "Not_in_either" = 0,
                      "Rank_AAC6" = -1, "Rank_plasmid" = -1,
                      "n_comp_plasmid" = 0)
  } else {
    sig_nodes <- nodes$NodeId
    comps <- nodes %>%
      group_by(CompId) %>%
      slice_min(CCSPvalue, with_ties = FALSE) %>%
      arrange(CCSPvalue)
    aac6_comp <- nodes %>% filter(NodeId == aac6)
    aac6_comp <- max(-1, which(comps$CompId == aac6_comp$CompId))
    plasmid_comp <- nodes %>% filter(NodeId %in% plasmid)
    n_plasmid_comp <- n_distinct(plasmid_comp$CompId)
    plasmid_comp <- min(which(comps$CompId %in% plasmid_comp$CompId))
    if (is.infinite(plasmid_comp)) plasmid_comp <- -1
    res <- data.frame("AAC6" = aac6 %in% sig_nodes %>% as.numeric(),
                      "Plasmid" = mean(plasmid %in% sig_nodes),
                      "Not_in_either" = mean(!sig_nodes %in% c(plasmid, aac6)),
                      "Rank_AAC6" = aac6_comp,
                      "Rank_plasmid" = plasmid_comp,
                      "n_comp_plasmid" = n_plasmid_comp
                      )
  }
  return(res)
}

locs <- list.dirs(here("Output", "Amikacin"), recursive = FALSE, full.names = FALSE) %>%
  str_subset(loc)
alphas <- word(locs, -1, sep = "_") %>% as.numeric()
res <- map2(locs, alphas, metrics)
names(res) <- alphas
res <- bind_rows(res, .id = "alpha")

write_csv(res,  here("Output", "Amikacin", paste0(loc, "results.csv")))
