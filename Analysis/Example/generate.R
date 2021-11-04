library(here)
library(dplyr)
library(readr)
# Helper functions ----
bases <- c("A", "T", "G", "C")
.modif <- function(seq, mu, insert, delete = insert) {
  new_seq <- lapply(seq, function(base) {
    # Mutation
    if (rbinom(1, 1, mu) == 1) {
      return(sample(bases[bases != base], 1))
    }
    # Insertion
    if (rbinom(1, 1, insert) == 1) {
      return(paste0(base, sample(bases, 1)))
    }
    # Deletion
    if (rbinom(1, 1, delete) == 1) {
      return("")
    }
    return(base)
  })
  return(paste0(new_seq, collapse = ""))
}

.generate_buffer <- function(buffer, n, mu, insert, n_distinct) {
  b <- sample(bases, buffer, replace = TRUE)
  buffer <- lapply(seq_len(n_distinct), function(seq) {
    return(.modif(b, mu = mu, insert = insert))
  })
  return(lapply(seq_len(n), sample, x = buffer, size = 1))
}

.generate_modif <- function(L, m, mu_res, insert) {
  accessory <- sample(bases, L, replace = TRUE)
  accessories <- lapply(seq_len(m), function(seq) {
    return(.modif(accessory, mu = mu_res, insert = insert))
  })
  return(accessories)
}

# Parameters ----
set.seed(21)
n <- 100
p <- .5
L <- 10^3
mu_res <- .01
n_distinct <- 10
mu <- .01
insert <- .01
buffer <- L
m <- round(n * p)

# Generate simulation data ----
buffer <- .generate_buffer(buffer, n, mu, insert, n_distinct)
accessories <- .generate_modif(L, n_distinct, mu_res, insert) %>%
  Reduce(f = `c`)
accessories <- lapply(seq_len(m), function(i) {
  return(sample(accessories, 1))
})
seq_res <- lapply(seq_len(m), function(i) {
  paste0(buffer[[i]], accessories[[i]])
}) %>% unlist()
seq_sens <- lapply(seq_len(n - m), function(i) {
  buffer[[i + m]]
}) %>% unlist()
seqs <- c(seq_res, seq_sens)
seqs <- lapply(seq_len(n), function(i) {
  header <- paste0(">Sample_", i, "Generated Strain")
  write_lines(c(header, seqs[i]),
              here("Example", "Raw", "Genomes", paste0("Sample_", i, ".fna")))
  return()
})

df <- data.frame(
  "GenomeID" = paste0("Sample_", seq_len(n)),
  "Pheno" = c(rep(1, m), rep(0, n - m)),
  "path" = paste0(here("Example", "Raw", "Genomes"), "/Sample_", seq_len(n), ".fna")
)

write_tsv(df, here("Example", "Raw", "strains"))
