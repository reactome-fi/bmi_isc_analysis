usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE, repos = "http://cran.us.r-project.org")
  require(p, character.only = TRUE)
}
usePackage("pacman")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("mgsa")
# --------------------------------------------------------
p_load(dplyr)
p_load(Matrix)
p_load(patchwork)
p_load(tidyverse)
p_load(plotly)
p_load(devtools)
p_load(reticulate)
# Required by iCytotrace
p_load(HiClimR)
p_load(ccaPP)
p_load(ggpubr)
p_load(nnls)
p_load(egg)
p_load(sva)
p_load(Rtsne)
p_load(plyr)
# p_load(mgsa)

set.seed(1234)

library(reticulate)
# Change to use the local miniconda env
use_python("path/to/env/for/cytotrace/bin/python")
reticulate::import('scanoramaCT')
reticulate::import('numpy')
# The following is needed for Python3.7
Sys.setenv(RETICULATE_PYTHON="path/to/env/for/cytotrace/bin/python")
# Running modified iCytotrace source - refer to manuscript
source('path/to/modified/iCytoTRACE.R')
source('path/to/unmodified/plotCytoTRACE.R')
# The following python package should be configured first after set the following
have_scanoramaCT <- T
have_numpy <- T

# Change these to local seq data directories
in.path <- "path/to/repo/"
E12.dir <- paste0(in.path, "raw/Adult_Irr_gfp/filtered_feature_bc_matrix//filtered_feature_bc_matrix/")
E17.dir <- paste0(in.path, "raw/17_5_gfp/filtered_feature_bc_matrix/")
adult.dir <- paste0(in.path, "raw/Adult_Hom_gfp/filtered_feature_bc_matrix/")
adult.dir <- paste0(in.path, "raw/Adult_Irr_gfp/filtered_feature_bc_matrix/")


outs2df <- function(path, case){
  
  barcode.path <- paste0(path, "barcodes.tsv.gz")
  features.path <- paste0(path, "features.tsv.gz")
  matrix.path <- paste0(path, "matrix.mtx.gz")
  
  mat <- readMM(file = matrix.path)
  feature.names <- read.delim(features.path, 
                              header = FALSE,
                              stringsAsFactors = FALSE)
  
  barcode.names <- read.delim(barcode.path, 
                              header = FALSE,
                              stringsAsFactors = FALSE)
  
  barcode.names$V2 <-  gsub("-1", "", barcode.names$V1)
  colnames(mat) <- barcode.names$V2
  rownames(mat) <- feature.names$V1
  
  mat <- as.matrix(mat)
  mat.tibb <- as_tibble(mat)
  colnames(mat.tibb) <- paste(case, colnames(mat.tibb), sep = "_")

  
  mat.tibb <- mat.tibb %>% add_column(gene_name = feature.names$V2, .before = 1)
  mat.tibb <- mat.tibb[!duplicated(mat.tibb$gene_name), ]
  
  mat.df <- as.data.frame(mat.tibb)
  rownames(mat.df) <- mat.df$gene_name
  mat.df <- mat.df[, -1]
  mat <-  as.matrix(mat.df)
  return(mat)
}

E12 <- outs2df(path = E12.dir, case = "12_5")
E17 <- outs2df(path = E17.dir, case = "17_5")
Hom <- outs2df(path = Hom.dir, case = "Hom")
Irr <- outs2df(path = Irr.dir, case = "Irr")

# Results for Figure 5E
#print('working on E12, E17, Hom')
#datasets.E12E17Hom <- list(E12, E17, Hom)
#results.E12E17Hom <- iCytoTRACE(datasets.E12E17Hom)
#save(results.E12E17Hom, file = "results_E12E17Hom.RData")
#plotCytoTRACE(results.E12E17Hom)

# Results for Figure 6G
# print('working on E12, E17, Hom, Irr')
# datasets.E12E17IrrHom <- list(E12, E17, Irr, Hom)
# rm(E12, E17, Irr, Hom)
# results.E12E17IrrHom <- iCytoTRACE(datasets.E12E17IrrHom)
# save(results.E12E17IrrHom, file = "results_E12E17IrrHom.RData")
# plotCytoTRACE(results.E12E17IrrHom)





