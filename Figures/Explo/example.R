library(dplyr)
library(ggplot2)
library(readr)
library(here)
library(tidyr)
library(stringr)
library(cowplot)
library(network)
library(GGally)
library(magick)
# Small node

p1 <- image_read(here("Figures", "Explo", "schema.png"))

# full graph
pheno_nodes <- read_tsv(here("Output", "Explo", "Setup", "textualOutput", "all_comps_nodes_info.tsv")) %>%
  select(NodeId, CCSId) %>%
  mutate(pheno = if_else(CCSId == 0, "gene A", "gene B"),
         NodeId = str_remove(NodeId, "n") %>% as.numeric())


edges <- read_tsv(here("Output", "Explo", "DBGWAS", "step1", "graph.edges.dbg"),
                  col_names = FALSE) %>%
  group_by(X1, X2, X3) %>%
  mutate(from = min(X1, X2), to = max(X1, X2)) %>%
  ungroup() %>%
  select(from, to) %>%
  distinct()

net <- network::network(edges %>% select(from, to), directed = FALSE,
                        vertices = pheno_nodes)
set.seed(14)
p2 <- ggnet2(net, mode = "kamadakawai", node.size = 1, color = "pheno",
       palette = c("gene A" = "steelblue", "gene B" = "tomato"), label = "") +
  labs(color = '') +
  theme_void() + 
  theme(legend.position = c(.8, .8),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 20)
        )
p2
p <- ggdraw() +
  draw_image(p1, height = .65, width = .65, x = .05, y = -.05) +
  draw_plot(p2, x = .13, width = .9, height = 1) +
  draw_text(x = .05, y = .9, text = "a)", size = 17) +
  draw_text(x = .05, y = .5, text = "b)", size = 17)
p
save_plot(here("Figures", "Explo", "Example.png"), plot = p2,
          base_height = 6, base_width = 10)
