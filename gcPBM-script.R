######################################
# 10.02.2019
# Multiple Linear Regression (MLR) example
# BISC 481
######################################

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
install.packages("ggplot2", repos = "http://cran.us.r-project.org")

## Initialization
library(DNAshapeR)
library(caret)
library(glmnet)
library(ggplot2)
workingPath <- "/Users/irithkatiyar/Desktop/qbio481/Machine-Learning-QBIO481/gcPBM/"

## Iterate through datasets
results <- list()
datafiles <- list("Mad.txt", "Max.txt", "Myc.txt")
for (datafile in datafiles) {
    ## Predict DNA shapes
    fn_fasta <- paste0(workingPath, paste0(datafile, ".fa"))
    pred <- getShape(fn_fasta)

    ## Iterate through 1-mer and 1-mer+shape feature types
    res <- list()
    featureTypes <- list(c("1-mer"), c("1-mer", "1-shape"))
    for(featureType in featureTypes) {
        ## Encode feature vectors    
        featureVector <- encodeSeqShape(fn_fasta, pred, featureType)
        head(featureVector)

        ## Build MLR model by using Caret
        # Data preparation
        fn_exp <- paste0(workingPath, datafile)
        exp_data <- read.table(fn_exp)
        df <- data.frame(affinity=exp_data$V2, featureVector)

        # Arguments setting for Caret
        trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

        # Prediction with L2-regularized
        model <- train(affinity~., data = df, trControl=trainControl, 
                    method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
        model
        res <- append(res, model$results$Rsquared[1])
    }
    results <- append(results, res)
}
print(results)

## Print average R-squared values for each dataset
avg_results <- list()
for(i in seq_along(datafiles)) {
  res <- c(results[[2*i-1]], results[[2*i]])
  datafile <- datafiles[[i]]
  avg <- mean(unlist(res))
  avg_results <- append(avg_results, avg)
  print(datafile)
  print(avg)
}

## Create plot
# Create a data frame with R-squared values and dataset names
data <- data.frame(
  Mad = c(results[[1]], results[[2]]),
  Mat = c(results[[3]], results[[4]]),
  Myc = c(results[[5]], results[[6]])
)

# Convert data frame to matrix
data_matrix <- as.matrix(data)

# Define the dataset names
dataset_names <- data$Dataset

# Create a bar graph
barplot(data_matrix, beside = TRUE, col = c("blue", "red"),
        names.arg = dataset_names, 
        main = "R-squared Values for Datasets",
        xlab = "Dataset", ylab = "R-squared")
legend("bottomright", legend = c("1-mer", "1-mer+shape"), fill = c("blue", "red"))
