# R script calling the main function and running on the demo datasets (Data1 from GBM)

library(glmnet)
library(grpreg)
library(survival)
library(MCMCglmm)

source("MainFunction.R")
source("UtilFunctions.R")

options(warn=-1)
# import the demo datasets -- data from each biological level
Gene <- scale(as.matrix(read.csv("Gene.csv", sep=";")))
CNV <- scale(as.matrix(read.csv("CNV.csv", sep=";")))
methy <- scale(as.matrix(read.csv("methy.csv", sep=";")))

# import information about markers, the gene they map to, and their functional networks
Gene_Methy <- read.csv("Gene_Methy.csv", sep=";")
Gene_CNV <- read.csv("Gene_CNV.csv",sep=";")
Gene_Pathway <- read.csv("Gene_Pathway.csv", sep=";")

# import survival data (survival time and censoring indicator) and other clinical covariates
y <- read.csv("Y.csv", header=TRUE, sep=";")

# The examples below run the main function that implements the integrative genomic model 
# allowing for association between CNVs and methylations

# Run the Integrative-gene scenario using
# univariate models for the association between methylation sites and CNVs
ana = Integrated_Original(Gene = Gene, methy = methy, CNV = CNV, y = y,
                          Gene_CNV = Gene_CNV, Gene_Methy = Gene_Methy, 
                          Gene_Pathway = Gene_Pathway, multi_methy = FALSE,
                          intra = TRUE, pathway = FALSE,nfolds = 50,
                          alpha = 1)

ana$R2.adj

# Run the Integrative-network scenario using
# univariate models for the association between methylation sites and CNVs
ana.network = Integrated_Original(Gene = Gene, methy = methy, CNV = CNV, y = y,
                          Gene_CNV = Gene_CNV, Gene_Methy = Gene_Methy, 
                          Gene_Pathway = Gene_Pathway, multi_methy = FALSE,
                          intra = FALSE, pathway = TRUE,nfolds = 50,
                          alpha = 1)
ana.network$R2.adj

