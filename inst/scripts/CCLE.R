## Ref: MSnbase IO capabilities; check MSnbase vignette for more info.
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

## Replace "NA"s with "Biological Samples" on the column of "Notes"
dataInfo$Notes[is.na(dataInfo$Notes)] <- "biological samples"

## Subset "dataInfo" - the meta data for all biological lines (remove bridges)
dataInfoSamples <- subset(dataInfo, Notes %in% "biological samples")

## Extract the quantitation data, from the original data
## Protein_Id: Composite unique protein identifier. It concatenates the
##      UniProt subdatabase and two different UniProt accession types.
exprsData <- normalizedData[, c(1, 49:426)]

## Format the pData(sample meta-data) of MSnSet
exprsTenPx <- as.integer(gsub(".*_TenPx", "", names(exprsData)[-1]))
exprsCellLines <- gsub("_TenPx.*$", "", names(exprsData)[-1])  ## CCLE Code
pdataCCLEleft <- data.frame(exprsCellLines, exprsTenPx)
names(pdataCCLEleft)[1] <- "CCLE.Code"
## Join pdataCCLE and dataInfoSamples as pData(exprsData)
## warning is compatible with the metadata
pdataCCLE <- pdataCCLEleft %>%
                 dplyr::left_join(dataInfoSamples,
                                  by=c("CCLE.Code" = "CCLE.Code",
                                       "exprsTenPx" = "Protein.10-Plex.ID"))
## Integration verified:
##    identical(exprsCellLines, pdataCCLE$CCLE.Code) = TRUE
##    identical(exprsTenPx, as.integer(pdataCCLE$exprsTenPx)) = TRUE

## Rename column of pdataCCLE
pdataCCLE <- pdataCCLE %>% dplyr::rename(ccleCode = CCLE.Code,
                                         Protein10PlexID = exprsTenPx,
                                         CellLine = Cell.Line,
                                         TissueOfOrigin = Tissue.of.Origin,
                                         ProteinTMTLabel=Protein.TMT.Label)
rownames(pdataCCLE) <- names(exprsData)[-1]


## Experiment info format as "AnnotatedDataFrame"
## nrow(phenoData) == dim(e)[2]
pd <- new("AnnotatedDataFrame", pdataCCLE)

## Extract the feature meta-data, from the original data
featuresDT <- normalizedData[, 1:48]
rownames(featuresDT) <- featuresDT$Protein_Id
## Feature data format as "AnnotatedDataFrame"
fd <- new("AnnotatedDataFrame", featuresDT)


## Read the normalized expression data
normalizedCCLE <- readMSnSet2(file = exprsData, ecol = 2:379, skip = 0, fnames = 1)
## Expression data
e <- exprs(normalizedCCLE)

## Experimental data to add
experiment <- new("MIAPE",
                  lab = "Gygi Lab, Harvard",
                  name = "Steve Gygi",
                  contact = "David Nusinow",
                  email = "david_nusinow@hms.harvard.edu",
                  samples = list(
                    species = "375 cell lines of the Cancer Cell Line Encyclopedia",
                    operator = "Gygi Lab"
                  ),
                  title = "Quantitative Proteomics of the Cancer Cell Line Encyclopedia",
                  abstract = "Proteins are essential agents of biological processes. To date, large-scale profiling of cell line collections including the Cancer Cell Line Encyclopedia (CCLE) has focused primarily on genetic information whereas deep interrogation of the proteome has remained out of reach. Here, we expand the CCLE through quantitative profiling of thousands of proteins by mass spectrometry across 375 cell lines from diverse lineages to reveal information undiscovered by DNA and RNA methods. We observe unexpected correlations within and between pathways that are largely absent from RNA. An analysis of microsatellite instable (MSI) cell lines reveals the dysregulation of specific protein complexes associated with surveillance of mutation and translation. These and other protein complexes were associated with sensitivity to knockdown of several different genes. These data in conjunction with the wider CCLE are a broad resource to explore cellular behavior and facilitate cancer research.",
                  pubMedIds = "31978347",
                  url = "https://doi.org/10.1016/j.cell.2019.12.023",
                  instrumentModel = "Thermo-Fisher Orbitrap Fusion or Orbitrap Fusion Lumos",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "ESI",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = "Sequest/R script",
                  collisionEnergy = "",
                  dateStamp = "Jan 2020"
)

## The "MSnProcess" Class
process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("Various Normalisation Approaches")),
               normalised=TRUE)

## CCLE Normalized Data
CCLENormalized <- new("MSnSet",
                        exprs = e,
                        phenoData = pd,
                        experimentData = experiment,
                        featureData = fd,
                        processingData = process)

## checks
stopifnot(dim(pData(CCLENormalized))[1] == ncol(e),
          dim(fData(CCLENormalized))[1] == nrow(e),
          validObject(CCLENormalized))

save(CCLENormalized, file="../../data/CCLENormalized.rda",
     compress = "xz", compression_level = 9)

