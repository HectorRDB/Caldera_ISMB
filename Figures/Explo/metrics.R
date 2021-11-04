library(dplyr)
library(ggplot2)
library(readr)
library(here)
library(stringr)
library(cowplot)
res <- list.files(here("Output", "Explo", "Results")) %>%
  lapply(., function(name) {
    alpha <- word(name, 2, sep = "_") %>% word(1, sep = "\\.")
    alpha <- 10^(-as.numeric(alpha))
    res_alpha <- read_csv(here("Output", "Explo", "Results", name)) %>%
      mutate(alpha = alpha)
  }) %>%
  bind_rows() %>%
  mutate(coverage = 100 * TP / (TP + FN))

p <- ggplot(res %>%
              filter(!str_detect(Analysis, "Stage") | Analysis == "All Stages"), 
            aes(x = alpha, y = coverage, col = Analysis)) +
  geom_line(size = 2) +
  scale_x_log10() +
  scale_y_continuous(breaks = 0:4 * 25, labels = paste0(0:4 * 25, "%")) +
  theme_bw() +
  labs(x = "Value of \u03b1", y = "Coverage of the accessory genome") +
  scale_color_manual(values = c("All Stages" = "#67000D",
                                "All Unitigs" = "#E41A1C",
                                "DBGWAS" = "#C6DBEF"),
                     labels = c("CALDERA", "All unitigs", "DBGWAS"))
ggsave(filename = here("Figures", "Explo", "Coverage_different_methods.png"),
       height = 6, width = 10)

p <- ggplot(res %>%
              filter(str_detect(Analysis, "Stage")) %>%
              mutate(Stage = str_remove(Analysis, "Stage"),
                     Stage = if_else(Stage == "All s", "All", Stage),
                     Stage = factor(Stage, levels = c(1, 2, 3, 5, 10, 15, 20, "All"))), 
            aes(x = alpha, y = coverage, col = Stage)) +
  geom_line(size = 2) +
  scale_x_log10() +
  scale_y_continuous(breaks = 0:4 * 25, labels = paste0(0:4 * 25, "%")) +
  theme_bw() +
  labs(x = "Value of \u03b1", y = "Coverage of the accessory 
       genome") +
  scale_color_brewer(palette = "Spectral", direction = -1)
ggsave(filename = here("Figures", "Explo", "Coverage_different_stages.png"),
       height = 6, width = 10)

p1 <- ggplot(res %>%
              filter(!str_detect(Analysis, "Stage") | Analysis == "All Stages"),
            aes(x = alpha, y = 100 * FP / (FP + TN), col = Analysis)) +
  geom_line(size = 2, alpha = .8) +
  scale_x_log10() +
  scale_y_continuous(breaks = 0:5 * 2, labels = paste0(0:5 * 2, "%")) +
  expand_limits(y = 10) +
  theme_bw() +
  labs(x = "Value of \u03b1", y = "False Positive Rate at the unitig level") +
  scale_color_manual(values = c("All Stages" = "#67000D",
                                "All Unitigs" = "#E41A1C",
                                "DBGWAS" = "#C6DBEF"),
                     labels = c("CALDERA", "All unitigs", "DBGWAS")) +
  theme(legend.position = "bottom") +
  guides(col = guide_legend(ncol = 2))


p2 <- ggplot(res %>%
         filter(str_detect(Analysis, "Stage")) %>%
         mutate(Stage = str_remove(Analysis, "Stage"),
                Stage = if_else(Stage == "All s", "All", Stage),
                Stage = factor(Stage, levels = c(1, 2, 3, 5, 10, 15, 20, "All"))),
       aes(x = alpha, y = 100 * FP / (FP + TN), col = Stage)) +
  geom_line(size = 2) +
  scale_x_log10() +
  scale_y_continuous(breaks = 0:5 * 2, labels = paste0(0:5 * 2, "%")) +
  expand_limits(y = 10) +
  theme_bw() +
  labs(x = "Value of \u03b1", y = "False Positive Rate at the unitig level") +
  scale_color_brewer(palette = "Spectral", direction = -1) +
  theme(legend.position = "bottom")

p <- plot_grid(p1, p2, labels = c("a)", "b)"), scale = .95)

save_plot(filename = here("Figures", "Explo", "FPR.png"), plot = p,
          base_height = 5, base_width = 10)
