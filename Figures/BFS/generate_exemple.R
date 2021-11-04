library(dplyr)
library(ggplot2)
library(here)

thresh <- function(alpha, n){
  n1 <- ceiling(n / 2)
  pvals <- lapply(n1:n, function(x) {
    mat <- matrix(c(n1, 0, x - n1, n - x), byrow = TRUE, ncol = 2)
    return(fisher.test(x = mat, )$p.value)
  }) %>% unlist()
  sigma <- max((n1:n)[pvals <= alpha])
}

prob_binom <- function(s, p) {
  if (s == 1) {
    return(p)
  } else {
    return(p + prob_binom(s - 1, p) * (1 - p))
  }
}

Prob_0 <- function(n, s, p, alpha, N) {
  sigma <- thresh(alpha, n)
  not_test <- pbinom(n - sigma - 1, n, prob_binom(s, p)) +
    pbinom(sigma, n, prob_binom(s, p), lower.tail = FALSE)
  if (not_test == 1) {
    not_test <- max(
      pbinom(n - sigma - 1, n, prob_binom(s, p), log.p = TRUE),
      pbinom(sigma, n, prob_binom(s, p), lower.tail = FALSE, log.p = TRUE))
    not_test <- 
    return(exp(N * not_test))
  } else {
    return(not_test^N)
  }
}

Prob_test <- function(n, s, p, alpha, N) {
  sigma <- thresh(alpha, n)
  not_test <- pbinom(n - sigma - 1, n, prob_binom(s, p)) +
    pbinom(sigma, n, prob_binom(s, p), lower.tail = FALSE)
  return(1 - not_test)
}

Prune <- function(n, s, p, alpha, N) {
  sigma <- thresh(alpha, n)
  not_test <- pbinom(sigma, n, prob_binom(s, p), lower.tail = FALSE)
  return(not_test)
}

df1 <- data.frame(s = rep(c(1, 2, 3), 3),
                 N = rep(c(100, 1000, 10000), each = 3))
df2 <- data.frame(s = rep(c(1, 2, 3), 99),
                  p = rep(1:99/100, each = 3))
df <- full_join(df1, df2) %>%
  as.matrix()
ps <- rep(0, nrow(df))
for (i in 1:nrow(df)) {
  cond <- df[i, ]
  ps[i] <- (Prob_0(n = 100, s = cond[1], p = cond[3], alpha = 10^(-4), N = cond[2]))
}

p_test <- rep(0, nrow(df))
for (i in 1:nrow(df)) {
  cond <- df[i, ]
  p_test[i] <- (Prob_test(n = 100, s = cond[1], p = cond[3], alpha = 10^(-4), N = cond[2]))
}

prune <- rep(0, nrow(df))
for (i in 1:nrow(df)) {
  cond <- df[i, ]
  prune[i] <- (Prune(n = 100, s = cond[1], p = cond[3], alpha = 10^(-4), N = cond[2]))
}

df <- as.data.frame(df)
df$ps <- ps
df$p_test <-  p_test
df$prune <-  prune

p1 <- ggplot(df %>% filter(N == 100),
            aes(x = p, y = p_test, col = as.character(s))) +
  geom_line(lwd = 2) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(col = "Stage", x = "Probabilty of a one",
       y = "Probability that a subgraph is 1-testable") +
  theme()
p2 <- ggplot(df %>% filter(N == 100),
            aes(x = p, y = prune, col = as.character(s))) +
  geom_line(lwd = 2) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(col = "Stage", x = "Probabilty of a one",
       y = "Probability that a subgraph is prunable")
legend <- get_legend(p1 + theme(legend.position = "bottom"))

first_row <- plot_grid(p1 + guides(col = FALSE), p2 + guides(col = FALSE),
                       ncol = 2)
second_row <- plot_grid(NULL, legend, NULL, ncol = 3, rel_widths = c(1, 4, 1))
p <- plot_grid(first_row, second_row, ncol = 1, rel_heights = c(5, 1))
save_plot(filename = here("Figures", "BFS", "simu.pdf"), plot = p,
          base_width = 9)
