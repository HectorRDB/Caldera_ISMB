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
  NULL +
  labs(x = "Number of nodes in the graph",
       y = "Time (seconds)",
       col = "Method",
       linetype = "Exploration") +
  theme_bw() +
  scale_linetype_manual(values = c("BFS" = "dashed", "DFS"  = 'solid')) +
  scale_color_manual(values = c("COIN + LAMP2" = "#377EB8",
                                "CALDERA (1 core)" = "#E41A1C",
                                "CALDERA (5 cores)" = "#67000D"),
                     labels = c( "COIN+\nLAMP2", "CALDERA\n1 core", "CALDERA\n5 cores")) +
  theme(legend.justification = 'left')

# Figure 1b ----
res <- list.files(here("Output", "Explo")) %>% str_subset("results_") %>%
  lapply(., function(name) {
    alpha <- word(name, 2, sep = "_") %>% word(1, sep = "\\.")
    alpha <- 10^(-as.numeric(alpha))
    res_alpha <- read_csv(here("Output", "Explo", name)) %>%
      mutate(alpha = alpha)
  }) %>%
  bind_rows() %>%
  mutate(coverage = 100 * TP / (TP + FN))

p2 <- ggplot(res %>%
              filter(!str_detect(Analysis, "Stage") | Analysis == "All Stages"), 
            aes(x = alpha, y = coverage, col = Analysis)) +
  geom_line(size = 2) +
  scale_x_log10() +
  scale_y_continuous(breaks = 0:4 * 25, labels = paste0(0:4 * 25, "%")) +
  theme_bw() +
  labs(x = "Value of alpha", y = "Coverage of the accessory genome", col = "Method") +
  scale_color_manual(values = c("All Stages" = "#67000D",
                                "All Unitigs" = "#984EA3",
                                "DBGWAS" = "#4DAF4A"),
                     labels = c("CALDERA", "All unitigs", "DBGWAS")) +
  theme(legend.justification = 'right')
# Main figure ----
theme <- theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.box = "horizontal",
  plot.background = element_blank())


plotList <- list(
  p1 + theme + guides(linetype = guide_legend(keywidth = 2)),
  p2 + theme)

p <- plot_grid(plotlist = plotList, ncol = 2, labels = c('a)', "b)"),
               label_size = 15, scale = .95)
save_plot(p, filename = here("Figures", "Main", "fig1.pdf"),
          base_height = 5, base_width = 12, limitsize = F)

