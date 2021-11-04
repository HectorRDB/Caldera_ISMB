library(here)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(cowplot)
files <- list.files(here("Output", "Simulations", "Imbalance")) %>%
  str_subset("txt") %>%
  str_subset("times", negate = TRUE)

res <- lapply(files, function(file){
  df <- readLines(here("Output", "Simulations", "Imbalance", file))[-1]
  df <- str_remove(df, "We continue with \\(")
  df <- str_remove(df, ",\\)")
  type <- word(file, 1, sep = "_")
  if (type == "full") {
    method <- "New Bound"
  } else {
    method <- "FACS Bound"
  }
  imbalance <- word(file, 2, sep = "_") %>%
    str_remove("\\.txt") %>% as.numeric()
  return(c("N" = sum(as.numeric((df))),
           "Imb" = imbalance,
           "method" = method))  
}) %>%
  do.call(what = 'rbind') %>%
  as.data.frame(stringsAsFactors = F) %>%
  mutate(Imb = as.numeric(Imb), N = as.numeric(N)) %>%
  filter(Imb <= .25)
p1 <- ggplot(res, aes(x = Imb, y = N, col = method)) +
  geom_line(size = 2) +
  scale_y_log10() +
  NULL +
  labs(x = "Ratio of the smallest\nover the largest population",
       y = "Number of subgraphs\nexplored",
       col = "CALDERA with") +
  theme_bw() +
  scale_color_brewer(palette = "Set1", direction = -1) +
  theme(legend.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 22),
        legend.title = element_text(size = 24),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = 'left',
        legend.box = "vertical") +
  guides(linetype = guide_legend(override.aes = list(size = .5)))

# Times
df <- read.table(here("Output", "Simulations", "Imbalance", "times.txt"), 
                 col.names = c("FACS Bound", "New Bound")) %>%
  mutate(Imb = unique(res$Imb)) %>%
  pivot_longer(-Imb, names_to = "method", values_to = "time")
p2 <- ggplot(df, aes(x = Imb, y = time, col = method)) +
  geom_line(size = 2) +
  scale_y_log10() +
  NULL +
  labs(x = "Ratio of the smallest\nover the largest population",
       y = "Time (s)",
       col = "CALDERA with") +
  theme_bw() +
  scale_color_brewer(palette = "Set1", direction = -1) +
  theme(legend.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 22),
        legend.title = element_text(size = 24),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = 'left',
        legend.box = "vertical") +
  guides(linetype = guide_legend(override.aes = list(size = .5)))

leg <- get_legend(p1)
p <- plot_grid(p1 + guides(col = 'none'), p2 + guides(col = 'none'),
               labels = c("a)", "b)"))
p <- plot_grid(p, plot_grid(NULL, leg, rel_widths = c(.5, 1)),
               ncol = 1, rel_heights = c(.8, .1))
save_plot(here("Figures", "Imbalance", "Imb.png"), plot = p,
          base_height = 7, base_width = 15)
