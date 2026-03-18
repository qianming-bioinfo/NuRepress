#!/usr/bin/env Rscript
cran_pkgs <- c(
  "readr", "dplyr", "tidyr", "stringr", "ggplot2",
  "data.table", "purrr", "cluster", "mclust"
)
bioc_pkgs <- c(
  "GenomicRanges", "IRanges", "S4Vectors",
  "rtracklayer", "Rsamtools", "Biostrings"
)

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(NULL))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

cat("Installed / verified CRAN packages:\n")
cat(paste(cran_pkgs, collapse = ", "), "\n")
cat("Installed / verified Bioconductor packages:\n")
cat(paste(bioc_pkgs, collapse = ", "), "\n")
