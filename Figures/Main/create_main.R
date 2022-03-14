library(dplyr)
library(stringr)
library(here)
library(ggplot2)
library(cowplot)
library(tidyr)
library(readr)
library(network)
library(ggnetwork)

# Figure 1a  ----
files <- list.files(here("Output", "Simulations", "Speed")) %>%
  str_subset("_res_2.txt")
res <- lapply(files, function(file){
  df <- read.table(here("Output", "Simulations", "Speed", file)) %>%
    as.matrix() %>% t()
}) %>%
  do.call(what = 'rbind') %>%
  as.data.frame() %>%
  rename("CALDERA_1_core_BFS" = "V1",
         "CALDERA_1_core_DFS" = "V2", 
         "COIN+LAMP2" = "V3",
         "CALDERA_5_core_BFS" = "V4") %>%
  mutate(file = files,
         N = word(file, 1, sep = "_") %>% as.numeric()) %>%
  select(-file) %>%
  pivot_longer(-N, names_to = "Method", values_to = "time") %>%
  mutate(Method = factor(Method,
                         levels = c("CALDERA_1_core_BFS", "CALDERA_1_core_DFS", "COIN+LAMP2",
                                    "CALDERA_5_core_BFS"))) %>%
  mutate(main = if_else(str_detect(Method, "BFS"), "BFS", "DFS")) %>%
  mutate(Method = case_when(
    str_detect(Method, "CALDERA_1_core") ~ "CALDERA (1 core)",
    str_detect(Method, "CALDERA_5_core") ~ "CALDERA (5 cores)",
    TRUE ~ "COIN + LAMP2"))
p1 <- ggplot(res, aes(x = N, y = time, col = Method, linetype = main)) +
  geom_line(size = 2) +
  scale_y_log10(breaks = 10^(c(1, 3, 5))) +
  scale_x_log10() +
  labs(x = "Number of nodes in the graph",
       y = "Time (seconds), log scale",
       col = "Method",
       title = "Time needed to explore a graph with different number of nodes",
       linetype = "Exploration") +
  theme_bw() +
  scale_linetype_manual(values = c("BFS" = "dashed", "DFS"  = 'solid')) +
  scale_color_manual(values = c("COIN + LAMP2" = "#377EB8",
                                "CALDERA (1 core)" = "#E41A1C",
                                "CALDERA (5 cores)" = "#67000D"),
                     labels = c( "COIN+\nLAMP2", "CALDERA\n1 core", "CALDERA\n5 cores")) +
  theme(legend.justification = 'center') +
  guides(linetype = guide_legend(keywidth = 3, override.aes = list(size = 1.3)))
p1
# Figure 1b ----
res <- list.files(here("Output", "Explo", "Results")) %>%
  str_subset("pyseer") %>%
  lapply(., function(name) {
    alpha <- word(name, -1, sep = "_") %>% word(1, sep = "\\.")
    alpha <- 10^(-as.numeric(alpha))
    res_alpha <- read_csv(here("Output", "Explo", "Results", name)) %>%
      mutate(alpha = alpha)
  }) %>%
  bind_rows() %>%
  mutate(coverage = 100 * TP / (TP + FN)) %>%
  filter(!str_detect(Analysis, "Stage") | Analysis == "All Stages") %>%
  mutate(Analysis = factor(Analysis, 
                           levels = c("All Stages", "All Unitigs",
                                      "kmer fixed", "kmer lmm", "kmer wg",
                                      "DBGWAS")))

p2 <- ggplot(res, aes(x = alpha, y = coverage, col = Analysis)) +
  geom_line(size = 2) +
  scale_x_log10(breaks = 10^(7:1 * -2)) +
  scale_y_continuous(breaks = 0:4 * 25, labels = paste0(0:4 * 25, "%")) +
  theme_bw() +
  labs(x = "Value of alpha",
       col = "Method",
       y = "Coverage of the accessory genome",
       title = "Coverage of the accessory genome by different methods,\nwhen changing the alpha threshold") +
  scale_color_manual(values = c("All Stages" = "#67000D",
                                "All Unitigs" = "#E41A1C",
                                "kmer fixed" = "#C7E9C0",
                                "kmer lmm" = "#74C476",
                                "kmer wg" = "#006D2C",
                                "DBGWAS" = "#984EA3"),
                     labels = c("CALDERA", "All unitigs",
                                "Fixed effect model on kmers",
                                "Mixed model effect on kmers",
                                "Elastic net model on kmers",
                                "DBGWAS"))

# Main figure ----
theme <- theme(
  plot.title = element_text(hjust = .5),
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.box = "vertical",
  plot.background = element_blank()
)


plotList <- list(
  p1 + theme,
  p2 + theme)

p <- plot_grid(plotlist = plotList, ncol = 1, labels = c('a)', "b)"),
               label_size = 15, scale = .95)
save_plot(p, filename = here("Figures", "Main", "fig1.pdf"),
          base_height = 12, base_width = 7, limitsize = F)

