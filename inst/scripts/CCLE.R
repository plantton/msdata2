library("MSnbase")
library("openxlsx")
library("dplyr")
library("tidyr")

## Set workdir as the current path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Table of Contents
## Part 1: Normalized Protein Data
## Part 2: Non-Normalized Protein Data
## Part 3: Biological triplicate on the first two multiplexed experiments

## Part 1:  Normalized Protein Data

## Read Sample Information
sampleInfoPath <- "../extdata/Table_S1_Sample_Information.xlsx"
## Read Normalized Protein Data (Excel sheet from CCLE)
normDataPath <- "../extdata/Table_S2_Protein_Quant_Normalized.xlsx"
## Load these 2 files
dataInfo <- read.xlsx(sampleInfoPath, sheet = 2, startRow = 1, colNames = T)
normalizedData <- read.xlsx(normDataPath,
                            sheet = 2,
                            startRow = 1,
                            colNames = T,
                            cols = 1:426)

## Check duplicate cell lines
## 10 cell lines are quantified twice, not all of which are "Bridge" samples.
duplicateCL <- dataInfo %>%
                group_by(CCLE.Code) %>%
                summarise(no_rows = length(CCLE.Code))

## Format the PhenoData of MSnSet

## Read the normalized expression data
