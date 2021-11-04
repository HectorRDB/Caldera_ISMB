library(here)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# Figure 1 ----
files <- list.files(here("Output", "Simulations", "Speed")) %>%
  str_subset("_res_2.txt")

res <- lapply(files, function(file){
  df <- read.table(here("Output", "Simulations", "Speed", file),
                   stringsAsFactors = F) %>%
    as.matrix() %>% t()
}) %>%
  do.call(what = 'rbind') %>%
  as.data.frame(stringsAsFactors = F) %>%
  rename("CALDERA_1_core_BFS" = "V1",
         "CALDERA_1_core_DFS" = "V2", 
         "COIN+LAMP2" = "V3",
         "CALDERA_5_core_BFS" = "V4") %>%
  mutate(file = files,
         N = word(file, 1, sep = "_"),
         N = as.numeric(N)) %>%
  select(-file) %>%
  pivot_longer(-N, names_to = "Method", values_to = "time") %>%
  mutate(Method = factor(Method,
    levels = c("CALDERA_1_core_BFS", "CALDERA_1_core_DFS", "COIN+LAMP2",
               "CALDERA_5_core_BFS"))) %>%
  mutate(main = if_else(str_detect(Method, "BFS"), "BFS", "DFS")) %>%
  mutate(Method = case_when(
    str_detect(Method, "CALDERA_1_core") ~ "CALDERA (1 core)",
    str_detect(Method, "CALDERA_5_core") ~ "CALDERA (5 cores)",
    TRUE ~ "COIN + LAMP2")) %>%
  filter(time < 2 * 24 * 60^2)
p <- ggplot(res, aes(x = N, y = time, col = Method, linetype = main)) +
  geom_line(size = 2) +
  scale_y_log10() +
  scale_x_log10() +
  NULL +
  labs(x = "Number of nodes in the graph",
       y = "Time (seconds)",
       col = "Method",
       linetype = "Exploration") +
  theme_bw() +
  scale_linetype_manual(values = c("BFS" = "dashed", "DFS"  = 'solid')) +
  scale_color_manual(values = c("COIN + LAMP2" = "#C6DBEF",
                                "CALDERA (1 core)" = "#E41A1C",
                                "CALDERA (5 cores)" = "#67000D"),
                     labels = c("CALDERA\n1 core", "CALDERA\n5 cores", "COIN+\nLAMP2")) +
  theme(legend.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 22),
        legend.title = element_text(size = 24),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = 'left',
        legend.box = "vertical") +
  guides(linetype = guide_legend(override.aes = list(size = .5))) 
p
ggsave(filename = here("Figures", "Speed", "results_speed.pdf"), plot = p,
       width = 9, height = 9.5)

# Mean ----
files <- list.files(here("Output", "Simulations", "Speed"))%>%
  str_subset("_res_") %>%
  # str_subset("_res_5.txt", negate = TRUE) %>%
  identity()
res <- lapply(files, function(file){
  df <- read.table(here("Output", "Simulations", "Speed", file),
                   stringsAsFactors = F) %>%
    as.matrix() %>% t()
}) %>%
  do.call(what = 'rbind') %>%
  as.data.frame(stringsAsFactors = F) %>%
  rename("CALDERA_1_core_BFS" = "V1",
         "CALDERA_1_core_DFS" = "V2", 
         "COIN+LAMP2" = "V3",
         "CALDERA_5_core_BFS" = "V4") %>%
  mutate(file = files,
         N = word(file, 1, sep = "_"),
         N = as.numeric(N),
         scenario = word(file, 3, sep = "_"),
         scenario = as.numeric(str_remove(scenario, "\\.txt"))) %>%
  select(-file) %>%
  pivot_longer(-c(N, scenario), names_to = "Method", values_to = "time")

res %>%
  group_by(scenario, N) %>%
  mutate(ratio = time / min(time)) %>%
  ungroup() %>%
  group_by(Method) %>%
  summarise(sd = sd(ratio),
            ratio = mean(ratio))


# Supp figures ----
res <- res %>%
  mutate(timeout = if_else(scenario == 1, .95 * 2 * 24 * 60^2, .95 * 24 * 60^2),
           Method = factor(Method, 
                         levels = c("CALDERA_1_core_BFS", "CALDERA_1_core_DFS", "COIN+LAMP2",
                                                  "CALDERA_5_core_BFS")),
         main = if_else(str_detect(Method, "BFS"), "BFS", "DFS"),
         Method = case_when(str_detect(Method, "CALDERA_1_core") ~ "CALDERA (1 core)",
                            str_detect(Method, "CALDERA_5_core") ~ "CALDERA (5 cores)",
                            TRUE ~ "COIN + LAMP2"))
p <- ggplot(res, aes(x = N, y = time, col = Method, linetype = main)) +
  geom_line(size = 2) +
  scale_y_log10() +
  scale_x_log10(breaks = c(10^(2:4))) +
  geom_text(aes(y = .6 * timeout), label = "timeout", col = "black", size = 8,
            data = res %>% group_by(scenario) %>%
              dplyr::summarise(N = if_else(max(N) >= 10^4, 300, 150), timeout = min(timeout),
                               main = "BFS")) +
  geom_hline(aes(yintercept = timeout), col = "black", size = 2) +
  NULL +
  labs(x = "Number of nodes in the graph",
       y = "Time (seconds)",
       col = "Method",
       linetype = "Exploration") +
  theme_bw() +
  scale_linetype_manual(values = c("BFS" = "dashed", "DFS"  = 'solid')) +
  scale_color_manual(values = c("COIN + LAMP2" = "#C6DBEF",
                                "CALDERA (1 core)" = "#E41A1C",
                                "CALDERA (5 cores)" = "#67000D"),
                     labels = c("CALDERA\n1 core", "CALDERA\n5 cores", "COIN+\nLAMP2")) +
  theme(legend.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 22),
        legend.title = element_text(size = 24),
        strip.text = element_text(size = 22),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = 'left',
        legend.box = "horizontal") +
  guides(linetype = guide_legend(override.aes = list(size = .5),
                                 title.position="top", title.hjust = 0.5),
         color = guide_legend(title.position = "top", title.hjust = 0.5)) +
  facet_wrap(~scenario, scales = "free")
p
ggsave(filename = here("Figures", "Speed", "supp_speed.pdf"), plot = p,
       width = 15, height = 12)

# Mem ----
files <- list.files(here("Output", "Simulations", "Speed"))%>%
  str_subset("_mem")
res <- lapply(files, function(file){
  df <- read.table(here("Output", "Simulations", "Speed", file),
                   stringsAsFactors = F) %>%
    as.matrix() %>% t()
}) %>%
  do.call(what = 'rbind') %>%
  as.data.frame(stringsAsFactors = F) %>%
  rename("CALDERA_1_core_BFS" = "V1",
         "CALDERA_1_core_DFS" = "V2", 
         "COIN+LAMP2" = "V3") %>%
  mutate(file = files,
         N = word(file, 1, sep = "_"),
         N = as.numeric(N)) %>%
  select(-file) %>%
  pivot_longer(-N, names_to = "Method", values_to = "Mem") %>%
  mutate(main = if_else(str_detect(Method, "BFS"), "BFS", "DFS"),
         Method = case_when(str_detect(Method, "CALDERA_1_core") ~ "CALDERA (1 core)",
                            TRUE ~ "COIN + LAMP2"))

p <- ggplot(res, aes(x = N, y = Mem, col = Method, linetype = main)) +
  geom_line(size = 2) +
  labs(x = "Number of nodes in the graph",
       y = "Peak memory usage (Gb)",
       col = "Method",
       linetype = "Exploration") +
  theme_bw() +
  scale_linetype_manual(values = c("BFS" = "dashed", "DFS"  = 'solid')) +
  scale_color_manual(values = c("COIN + LAMP2" = "#C6DBEF",
                                "CALDERA (1 core)" = "#E41A1C"),
                     labels = c("CALDERA\n1 core", "COIN+\nLAMP2")) +
  theme(legend.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 22),
        legend.title = element_text(size = 24),
        legend.position = "bottom",
        legend.direction = "horizontal",
        # legend.justification = 'left',
        legend.box = "vertical") +
  guides(linetype = guide_legend(override.aes = list(size = .5)))
p
ggsave(filename = here("Figures", "Speed", "results_mem.pdf"), plot = p,
       width = 9, height = 7)
