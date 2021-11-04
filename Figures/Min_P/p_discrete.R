library(ggplot2)
library(dplyr)
library(here)
n <- 100
n1 <- 50
x <- 64
as <- 14:50
names(as) <- as
df <- lapply(as, function(a) {
  mat <- matrix(c(a, n1 - a, x - a, n - n1 - x + a), ncol = 2, byrow = TRUE)
  return(data.frame(p_value = fisher.test(mat)$p.value))
}) %>%
  bind_rows(.id = "as") %>%
  mutate(as = as.numeric(as))
p <- ggplot(df, aes(x = as, y = log10(p_value))) +
  geom_point() +
  theme_bw() +
  labs(y = expression(log[10](p(S))),
       x = expression(a[S]))

ggsave(plot = p, filename = here("Figures", "Min_P", "p_discrete.pdf"), height = 5)
