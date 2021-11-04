library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(cowplot)

df2 <- read.table(here("Figures", "Imbalance", "env_1.txt"),
                 col.names = c("x1", "x2", "CALDERA", "FACS"))
df1 <- read.table(here("Figures", "Imbalance", "env_2.txt"),
                  col.names = c("x1", "x2", "CALDERA", "FACS"))
TH <- qchisq(10^(-8), df = 1, lower.tail = FALSE) 

plots <- lapply(list(df1, df2), function(df){
  df <- df %>%
    mutate(x1 = x1, x2 = x2) %>%
    mutate(Prunable = (CALDERA < TH) + (FACS < TH)) %>%
    mutate(Prunable = case_when(Prunable == 0 ~ "By neither",
                                Prunable == 1 ~ "By CALDERA only",
                                TRUE ~ "By CALDERA and FACS"))
  
  
  p <- ggplot(df, aes(x = x1, y = x2, fill = Prunable)) +
    geom_tile(height = 1, width = 1, col = "transparent") +
    theme_minimal() +
    geom_vline(xintercept = max(df$x1 / 2)) +
    geom_hline(yintercept = max(df$x2 / 2) - .5) +
    labs(x = expression(x[1]), y = expression(x[2]), fill = "Is the subgraph\nprunable?") +
    scale_x_continuous(breaks = c(0, max(df$x1) / 2, max(df$x1)),
                       labels = c(0, expression(n[1] / 2), expression(n[1]))) +
    scale_y_continuous(breaks = c(-.5, (max(df$x2) / 2) - .5, max(df$x2) + .5),
                       labels = c(0, expression(n[2]/2), expression(n[2]))) +
    scale_fill_manual(
      values = c("By neither" = "#FEE0D2",
                 "By CALDERA only" = "#FC9272",
                 "By CALDERA and FACS" = "#DE2D26")
    )
  
  return(p)
})
plots[[1]] <- plots[[1]] + guides(fill = FALSE)
legend <- get_legend(plots[[2]] + theme(legend.position = "bottom"))
plots[[2]] <- plots[[2]] + guides(fill = FALSE)
first_row <- plot_grid(plotlist = plots, ncol = 2, rel_widths  = c(1, 1),
                       labels = c("a)", "b)"))
second_row <- plot_grid(NULL, legend, NULL, ncol = 3, rel_widths = c(1, 4, 1))
p <- plot_grid(first_row, second_row, rel_heights = c(4, 1), ncol = 1)

save_plot(filename = here("Figures", "Imbalance", "intuition_imbalance.png"),
          plot = p)
