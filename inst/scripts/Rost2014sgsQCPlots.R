library("MSnbase")
library("readr")
library("dplyr")
library("tidyr")

## Firstly we need the manual Skyline analysis result
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sgsManualFile <- paste("../../inst/extdata/",
                      "OpenSWATH_SM3_GoldStandardManualResults.csv",
                      sep = "")

sgsManualData <- read.csv(sgsManualFile)

## Select Replicate Names from Human background data
sgsHumanReplicateName <- gsub("split_(.*)_combined.*",
                           "\\1", colnames(exprs(Rost2014sgsHuman)))
## Select Skylike manual data which has corresponding Human background data
sgsManualDataHuman <- sgsManualData %>%
                        filter(ReplicateName %in% sgsHumanReplicateName)

