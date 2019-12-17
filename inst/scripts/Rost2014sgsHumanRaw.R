library(MSnbase)
## library(BiocFileCache)
library(HDF5Array)
library(mzR)
library(beachmat)
library(Matrix)

## Set workdir as the current path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## Sed to load the MSnbase package and
setMSnbaseVerbose(FALSE)

## Set Raw Data Path
rawMSPath <- "../../../RawMSData/napedro_L120420_010_SW.mzXML.gz"

## Read MS raw datasets as "MSnExp" instance (One Raw dataset)
## https://lgatto.github.io/MSnbase/articles/v04-benchmarking.html
rawMSdata <- readMSData(rawMSPath,
                        msLevel. = 1L,
                        mode = "onDisk")
rawMSdata
length(rawMSdata)
rawMSdata[1:3]
## Chromatograms
chr <- chromatogram(rawMSdata)
plot(chr)
##
dim(fData(rawMSdata))
colnames(fData(rawMSdata))

## Plotting MS spectra
spc1000 <- rawMSdata[[1000]]
plot(spc1000)

## Read MS2 raw data
rawMSdata2 <- readMSData(rawMSPath,
                         msLevel. = 2L,
                         mode = "onDisk")
rawMSdata2
## Chromatograms
spcMS2 <- rawMSdata2[[1]]
plot(spcMS2)
## Filter and Visualize
chr2 <- filterRt(rawMSdata[1:5], rt = c(6, 7))
chr2 <- filterMz(chr2, mz = c(1005, 1006))
plot(chr2, type = "XIC")
