library(MSnbase)
library(BiocFileCache)
library(HDF5Array)
library(mzR)
library(beachmat)
library(Matrix)

## Set workdir as the current path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## Sed to load the MSnbase package and
setMSnbaseVerbose(FALSE)

## Set Raw Data Path
rawMSPath <- "../../../RawMS/napedro_L120420_010_SW.mzXML.gz"

## Read MS raw datasets as "MSnExp" instance (One Raw dataset)
## https://lgatto.github.io/MSnbase/articles/v04-benchmarking.html
rawMSdata <- readMSData(rawMSPath, msLevel. = 1L,
                        mode = "onDisk")

