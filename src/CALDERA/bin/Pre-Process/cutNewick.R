#!/usr/bin/env Rscript

suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the data"
  ),
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the output"
  ),
  make_option(c("-k", "--kuts"),
              action = "store", default = NA, type = "integer",
              help = "How many groups to create"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at ", loc, "\n")
} else {
  stop("Missing l argument")
}

if (!is.na(opt$o)) {
  output <- opt$o
} else {
  stop("Missing o argument")
  cat("The output will be stored at ", output, "\n")
}

if (!is.na(opt$k)) {
  k <- opt$k
  cat("Cutting into ", k, " groups.\n")
} else {
  stop("Missing k argument")
}


library(phylogram)
library(dendextend)

tree <- read.dendrogram(file = loc)
pop <- cutree(tree, k = k)
write.table(pop - 1, output, row.names = FALSE, col.names = FALSE)

