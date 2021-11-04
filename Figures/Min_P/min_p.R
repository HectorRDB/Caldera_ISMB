library(ggplot2)
library(dplyr)
library(here)
n <- 100
n1 <- 25
xs <- 0:100
names(xs) <- xs

minimal_p <- function(x, n1, n) {
  # Compute marginals
  n2 <- n - n1
  
  # Compute the two possible p-values
  al <- max(0, x - n2)
  ar <- min(x, n1)
  mat1 <- matrix(c(al, x - al, n1 - al, n2 - x + al), ncol = 2, byrow = TRUE)
  mat2 <- matrix(c(ar, x - ar, n1 - ar, n2 - x + ar), ncol = 2, byrow = TRUE)
  pval1 <- fisher.test(mat1)$p.value
  pval2 <- fisher.test(mat2)$p.value
  
  # Return the right matrice.
  return(min(pval1, pval2))
}

df <- lapply(xs, FUN = minimal_p, n1 = n1, n = n) %>%
  unlist() %>%
  data.frame("min_p" = .) %>%
  mutate(xs = xs)
p <- ggplot(df, aes(x = xs, y = log10(min_p))) +
  geom_point() +
  theme_bw() +
  labs(y = expression(log[10](p^{"*"}~(S))),
       x = expression(x[S])) +
  scale_x_continuous(breaks = c(0, n1, n/2, n - n1, n),
                     labels = c("0", expression(n[1]), expression(n/2), expression(n[2]), "n"))
p

ggsave(plot = p, filename = here("Figures", "Min_P", "min_p.pdf"), height = 5)