library(here)
library(dplyr)
library(stringr)

# CLean nodes ----
nodes <- readLines(here("Output", "Ath", "Sig_subgraphs_nodes"))

# Find pathways ----
node_matching <- read.table(here("Output", "Ath", "numbers.matching"),
                            header = TRUE) %>%
  select(ps, number) %>%
  arrange(number) %>%
  mutate(ps = factor(ps))
edges <- list.files(here("Data", "Ath", "step1", "edges"))
names(edges) <- edges
edges <- lapply(edges, function(edge){
  df <- read.table(here("Data", "Ath", "step1", "edges", edge), stringsAsFactors = F,
                   col.names = c("from", "to")) %>%
    mutate(from = factor(from, levels = levels(node_matching$ps)),
           to = factor(to, levels = levels(node_matching$ps))) %>%
    mutate(from = as.numeric(from) - 1, to = as.numeric(to) - 1)
  return(df)
}) %>% 
  bind_rows(.id = "pathway") %>%
  distinct()
  
pathways <- list()
for (i in 1:length(nodes)) {
  node <- nodes[i] %>% str_split(",") %>% unlist() %>% as.numeric()
  paths <- edges %>%
    filter(from %in% node | to %in% node) %>%
    select(pathway) %>%
    distinct()
  pathways[[i]] <- paths
}
pathways <- unlist(pathways) %>% table()


genes <- nodes %>% str_split(",") %>% unlist() %>% as.numeric()
genes <- node_matching %>%
  filter(number %in% genes)
