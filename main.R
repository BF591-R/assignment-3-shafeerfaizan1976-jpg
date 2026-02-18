#!/usr/bin/Rscript
## Author: Faizan Shafeer
## faizi225@bu.edu
## BU BF591
## Assignment Bioinformatics Basics

#### Bioconductor ####
# it is standard among R packages to define libraries and packages at the 
# beginning of a script. Also note that a package should NOT be installed every 
# time a script runs.
# The bioconductor repository has installation instructions for biomaRt: 
# https://bioconductor.org/install/

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("biomaRt", quietly = TRUE)){
  BiocManager::install("biomaRt")
}
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(tidyverse))

#### Loading and processing data ####
#' Load Expression Data
#'
#' @param filepath A text string of the full filepath to the file to load.
#'
#' @return A tibble containing the data loaded from the CSV in `filepath`. 
#' 
#' @details Note that not all CSVs are created equal, and there are often cases where 
#' the data will not load in correctly on the first try. You may want to write this functon to 
#' adjust the CSV file being loaded into this assignment so that it can be formed into a 
#' tibble correctly.
#'
#' @examples 
#' `data <- load_expression('/project/bf528/project_1/data/example_intensity_data.csv')`
load_expression <- function(filepath) {
  data <- read.csv(filepath, check.names = FALSE)
  colnames(data)[1] <- "probe"
    return(as_tibble(data))
}

#' Filter 15% of the gene expression values.
#'
#' @param tibble A tibble of expression values, rows by probe and columns by sample.
#'
#' @return A tibble of affymetrix probe names from the input expression data tibble. 
#' These names match the rows with 15% or more of the expression values about log2(15).
#' 
#'
#' @examples `samples <- filter_15(data_tib)`
#' `> str(samples)`
#' `tibble [40,158 × 1] (S3: tbl_df/tbl/data.frame)`
#' `$ probe: chr [1:40158] "1007_s_at" "1053_at" "117_at" "121_at" ...`
filter_15 <- function(tibble){
  expr_cols <- tibble[, -1]
  threshold <- log2(15)
  n_samples  <- ncol(expr_cols)
  
  # apply() with MARGIN=1 runs the function on each row
  row_passes <- apply(expr_cols, 1, function(row) {
    sum(row > threshold) / n_samples >= 0.15
  })
  
  return(select(tibble[row_passes, ], 1))
}

#### Gene name conversion ####

#' Convert affymetrix array names into hgnc_symbol IDs using biomaRt. Inputs and 
#' outputs will likely not be the same size.
#'
#' @param affy_tib A single column tibble of strings containing array names.
#'
#' @return A 2 column tibble that contains affy IDs in the first column,
#' and their corresponding HGNC gene ID in the second column. Note that not all affy IDs 
#' will necessarily correspond with a gene ID, and one gene may have multiple affy IDs.
#' 
#' @details Connecting to ensembl via biomaRt can be...hit or miss...so you may 
#' want to check if data was correctly returned (or if it was just empty). The 
#' `getBM()` function may not accept a tibble, so you might need to convert your 
#' input into a flat vector.
#'
#' @examples 
#' `> affy_to_hgnc(tibble(c('202860_at', '1553551_s_at')))`
#' `affy_hg_u133_plus_2 hgnc_symbol`
#' `1        1553551_s_at      MT-ND1`
#' `2        1553551_s_at       MT-TI`
#' `3        1553551_s_at       MT-TM`
#' `4        1553551_s_at      MT-ND2`
#' `5           202860_at     DENND4B`
affy_to_hgnc <- function(affy_vector) {
  # Convert tibble to plain character vector
  if (is.data.frame(affy_vector) || is_tibble(affy_vector)) {
    affy_vector <- pull(affy_vector, 1)
  }
  
  mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                     dataset = "hsapiens_gene_ensembl")
  
  results <- getBM(
    attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"),
    filters    = "affy_hg_u133_plus_2",
    values     = affy_vector,
    mart       = mart
  )
  
  return(as_tibble(results))
}

#' Reduce a tibble of expression data to only the rows in good_genes or bad_genes.
#'
#' @param expr_tibble A tibble holding the expression data, each row corresponding to
#' one affymetrix probe ID and each column to a sample.
#' @param names_ids A two column tibble that associates affy IDs with HGNC gene IDs. 
#' Generated `with affy_to_hgnc()`.
#' @param good_genes A list of gene names stored as a vector of strings.
#' @param bad_genes A list of gene names stored as a vector of strings.
#'
#' @return A tibble with two additional columns added:
#' 1. HGNC gene IDs 
#' 2. Does the gene is this row fall into "good" or "bad" genes?
#' This tibble should be reduced to only rows present in good or bad genes. All
#' other rows can be discarded.
#' 
#' @details In order to plot only our genes of interest, we need to rearrange our 
#' data to include only the elements we want to see. We also want to add to columns, 
#' one that associates the probeids with the HGNC gene name, and one that says if 
#' that gene is in the good or bad sets of genes.
#'
#' @examples 
#' `plot_tibble <- reduce_data(expr_tibble = expr, names_ids = sample_names,`
#' `                           goodGenes, badGenes)`
#' `> head(plot_tibble)`
#' `A tibble: 6 × 38`
#' `  probeids    hgnc    gene_set    GSM972389 ...`
#' `  <chr>       <chr>   <chr>       <dbl>     ...`
#' `1 202860_at   DENND4B good        7.16      ...`
#' `2 204340_at   TMEM187 good        6.40      ...`
reduce_data <- function(expr_tibble, names_ids, good_genes, bad_genes){
  # Match probe IDs in the expression tibble to the affy IDs in names_ids
  probe_ids  <- pull(expr_tibble, 1)
  match_idx  <- match(probe_ids, names_ids$affy_hg_u133_plus_2)
  hgnc_names <- names_ids$hgnc_symbol[match_idx]
  
  # Add hgnc_symbol right after the probe column
  result <- add_column(expr_tibble, hgnc_symbol = hgnc_names, .after = 1)
  
  # Keep only rows that belong to good or bad gene sets
  keep   <- which(result$hgnc_symbol %in% c(good_genes, bad_genes))
  result <- result[keep, ]
  
  # Label each row as "good" or "bad"
  gene_set_labels <- ifelse(result$hgnc_symbol %in% good_genes, "good", "bad")
  
  # Add gene_set right after hgnc_symbol
  result <- add_column(result, gene_set = gene_set_labels, .after = "hgnc_symbol")
  
  return(result)
}

#' Convert a wide format tibble to long for easy plotting
#'
#' @param tibble A tibble of data in wide format. Specifically, it should be 
#' operating on the reduced tibble created by the previous function
#'
#' @return A tibble properly converted from wide to long format, the old sample 
#' columns should now be contained within a single column called "sample"
#' @export
#'
#' @examples
convert_to_long <- function(tibble) {
  result <- pivot_longer(
    tibble,
    cols      = starts_with("GSM"),
    names_to  = "sample",
    values_to = "value"
  )
  return(result)
}

