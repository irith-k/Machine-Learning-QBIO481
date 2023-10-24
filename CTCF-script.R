## Install packages
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
# DNAshapeR
BiocManager::install("DNAshapeR")
# Caret
install.packages("caret", repos = "http://cran.us.r-project.org")
install.packages("glmnet", repos = "http://cran.us.r-project.org")

## Initialization
library(DNAshapeR)
library(caret)
library(glmnet)
workingPath <- "/Users/irithkatiyar/Desktop/qbio481/Machine-Learning-QBIO481/CTCF/"