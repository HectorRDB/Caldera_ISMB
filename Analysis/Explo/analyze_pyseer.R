suppressWarnings(library(here))
suppressWarnings(library(dplyr))
suppressWarnings(library(readr))
suppressWarnings(library(stringr))
suppressWarnings(library(tidyr))
suppressWarnings(library(optparse))
# Load kmer input
res <- paste0("Sample_", 0:49)
sens <- paste0("Sample_", 50:99)
kmers <- read_lines(here("Output", "Explo", "pyseer", "kmers.txt.gz")) %>%
  data.frame("line" = .) %>%
  mutate(variant = word(line, 1, sep = " "))

kmers$g_res <- lapply(kmers$line, function(i) {
  any(lapply(res, str_detect, string = i) %>% unlist()) %>% return()
}) %>% unlist()

kmers$g_sens <- lapply(kmers$line, function(i) {
  any(lapply(sens, str_detect, string = i) %>% unlist()) %>% return()
}) %>% unlist()

kmers <- kmers %>%
  mutate(res_node = g_res & !g_sens) %>%
  select(variant, res_node)

metrics <- function(df, alpha) {
  df %>% 
    mutate(Sig = if_else(is.na(pval), FALSE, pval < alpha / n_distinct(variant))) %>%
    group_by(Analysis) %>%
    summarise("TN" = sum(Sig & !res_node),
              "FN" = sum(!Sig & res_node),
              "TP" = sum(Sig & res_node),
              "FP" = sum(Sig & !res_node)) %>%
    return()
}

pyseer <- lapply(
  c("pyseer.assoc", "pyseer_wg.assoc", "pyseer_lmm.assoc"),
  function(file) {
    df <- read_table(here("Output", "Explo", "pyseer", file)) %>%
      select(variant, `lrt-pvalue`) %>%
      rename(pval = `lrt-pvalue`) %>%
      right_join(kmers)
}) %>%
  `names<-`(c("kmer fixed", "kmer wg", "kmer lmm")) %>%
  bind_rows(.id = "Analysis")


# Analysis
for (i in 1:15) {
  alpha <- 10^(-i)
  res_pyseer <- metrics(pyseer, alpha)
  res <- read_csv(here("Output", "Explo", "Results", paste0("results_", i, ".csv")))
  write_csv(bind_rows(res, res_pyseer),
            here("Output", "Explo", "Results", paste0("results_pyseer_", -log10(alpha), ".csv")))
}

