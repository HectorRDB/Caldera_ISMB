library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(here)
library(stringr)
library(cowplot)

res_all_unitigs <- read_csv(here("Output", "Pseudomonas", "AllUnitigs_results.csv"))
res_caldera <- read_csv(here("Output", "Pseudomonas", "Caldera_S7_results.csv"))
res_dbgwas <- read_csv(here("Output", "Pseudomonas", "DBGWAS_q_results.csv"))

cov <- bind_rows("All Unitigs" = res_all_unitigs, "CALDERA" = res_caldera,
                 "DBGWAS" = res_dbgwas, .id = "Method") %>%
  select(alpha, Method, Plasmid, AAC6, Not_in_either) %>%
  pivot_longer(c("Plasmid", "AAC6", "Not_in_either"), values_to = "Coverage", names_to = "Object")

p <- ggplot(cov %>% filter(Object != "Not_in_either"),
            aes(x = alpha, y = Coverage, col = Method)) +
  geom_line(size = 2, alpha = .8) +
  facet_wrap(~Object) +
  scale_x_log10(breaks = 10^(-12:3)) +
  scale_y_continuous(breaks = 0:4 * 25/100, labels = paste0(0:4 * 25, "%")) +
  theme_bw() +
  labs(x = "Value of \u03b1", y = "Coverage of the plasmid and the AAC6' mutation") +
  scale_color_manual(values = c("CALDERA" = "#67000D",
                                "All Unitigs" = "#E41A1C",
                                "DBGWAS" = "#C6DBEF"),
                     labels = c("CALDERA", "All unitigs", "DBGWAS"))
p
ggsave(filename = here("Figures", "Pseudomonas", "Coverage_different_methods.png"),
       height = 4, width = 10)

cov <- cov %>%
  filter(Object != "AAC6") %>%
  pivot_wider(names_from = "Object", values_from = "Coverage") %>%
  arrange(alpha)
ggplot(cov, aes(x = Plasmid, y = Not_in_either, col = Method)) +
  geom_path(size = 1) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = c("CALDERA" = "#67000D",
                                "All Unitigs" = "#E41A1C",
                                "DBGWAS" = "#C6DBEF")) +
  labs(x = "Coverage of the plasmid", y = "Percentage of unitigs not in the plasmid or AAC6'")


ranks <- bind_rows("All Unitigs" = res_all_unitigs, "CALDERA" = res_caldera,
                   "DBGWAS" = res_dbgwas, .id = "Method") %>%
  select(alpha, Method, starts_with("rank"), n_comp_plasmid) %>%
  pivot_longer(-c("alpha", "Method"), values_to = "Rank", names_to = "Object") %>%
  filter(Rank != -1) %>%
  mutate(Type = if_else(str_detect(Object, "Rank"), "Ranking", "Number of components"),
         Object = if_else(str_detect(Object, "AAC6"), "AAC6", "Plasmid")) %>%
  pivot_wider(names_from = Type, values_from = Rank, values_fill = 1) %>%
  pivot_longer(c("Ranking", "Number of components"), names_to = "Type", values_to = "Values") %>%
  filter(Method != "All Unitigs")

ggplot(ranks, aes(x = alpha, y = Values, col = Method)) +
  geom_line(size = 2, alpha = .8) +
  facet_grid(Type~Object) +
  scale_x_log10(breaks = 10^(-12:3)) +
  theme_bw() +
  labs(x = "Value of \u03b1", y = "") +
  scale_color_manual(values = c("CALDERA" = "#67000D",
                                # "All Unitigs" = "#E41A1C",
                                "DBGWAS" = "#C6DBEF"),
                     labels = c("CALDERA", 
                                # "All unitigs",
                                "DBGWAS"))
 