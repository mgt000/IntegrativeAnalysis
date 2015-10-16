This repository contains the source code, sample R scripts, and example datasets to implement the integrative genomic models presented in the Bioinformatics paper "Evaluation of hierarchical models for integrative genomic analyses" by M. Denis and M.G. Tadesse.

R code:
- MainFunction.R -- main function for integrative models allowing for CNV-methylation association
- MainFunction-noCNVmethyl.R -- main function for integrative models with no CNV-methylation association
- UtilFunctions.R -- defines "regression" function called in main function
- call.R -- sample script to run integrative models on demo datasets

Example datasets:
- CNV.csv -- matrix of CNV data
- methy.csv -- matrix of methylation data 
- Gene.csv -- matrix of gene expression data
- Gene_CNV.csv -- CNV markers with the identifier of the gene they map to and their functional network
- Gene_Methy.csv -- methylation markers with the identifier of the gene they map to and their functional network
- Gene_Pathway.csv -- gene identifier and associated functional network
- Y.csv -- matrix with survival time, censoring indicator, and additional covariates
