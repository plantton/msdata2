library("MSnbase")
library("readr")
library("SWATH2stats")


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sgshumantxtfile <- "../../inst/extdata/OpenSWATH_SM3_GoldStandardAutomatedResults_human_peakgroups.txt"
sgshumantxt <- read.delim(sgshumantxtfile)
## import 422 stable isotope-labeled standard (SIS) peptides data
sispeptidefile <- "../extdata/NBT-L30171D-OpenSWATH_SM2_GoldStandardPeptides.csv"
sis.peptides <- read.csv(sispeptidefile,sep=";")


rownames(sgshumantxt) <- gsub("AQUA4SWATH_(.*)/2_run0_split_napedro_(.*)_SW_combined.*",
                              "\\1_\\2_", sgshumantxt$transition_group_id)


## Import the entire data matrix as MsnSet - assayData
## Note: "transition_group_id", "id" is unique for each row
##       "filename" contains 30 different strings. - Check SGS data description.
##       "Sequence", "FullPeptideName" ,"aggr_Fragment_Annotation" contain 345 different levels/numbers.
##       "ProteinName" contains 16 different levels.
## Attention: As suggested by "SWATH2stats":
##             For some R packages such as  imsbInfer,
##             All the columns should be preserved.
##             Hence I will keep all the columns from the OpenSWATH results.
e.sgs.human <- readMSnSet2(file = sgshumantxt, ecol = 1:54)

## Create a "Study design" table for experimental design description
Study_design.sgs.human <- data.frame(Filename = unique(sgshumantxt$filename))
## Add "Condition" column: Control vs Disease
## "Condition" represents different dilution steps
Study_design.sgs.human$Condition <- gsub(".*_(.*)_SW.*", "\\1", Study_design.sgs.human$Filename)
## "BioReplicate" column
Study_design.sgs.human$BioReplicate <- factor(gsub(".*_(.*)_0.*", "\\1", Study_design.sgs.human$Filename))
levels(Study_design.sgs.human$BioReplicate) <- seq(nlevels(Study_design.sgs.human$BioReplicate))
## "Run" column:
Study_design.sgs.human$Run <- seq_along(Study_design.sgs.human$Filename)


## Experimental data to add
experiment <- new("MIAPE",
                  lab = "Aebersold Group, ETH Zürich",
                  name = "Ruedi Aebersold",
                  contact = "Ruedi Aebersold",
                  email = "aebersold@imsb.biol.ethz.ch",
                  samples = list(
                    species = c("Human","Yeast","Streptococcus pyogenes")
                  ),
                    title = "OpenSWATH enables automated, targeted analysis of data-independent acquisition MS data:SWATH-MS Gold Standard (SGS) Dataset, S. pyogenes Dataset",
                  abstract = "The SWATH-MS Gold Standard (SGS) dataset consists of 90 SWATH-MS runs of 422 synthetic stable isotope-labeled standard (SIS) peptides in ten different dilution steps, spiked into three protein backgrounds of varying complexity (water, yeast and human), acquired in three technical replicates. The SGS dataset was manually annotated, resulting in 342 identified and quantified peptides with three or four transitions each. In total, 30,780 chromatograms were inspected and 18,785 were annotated with one true peak group, whereas in 11,995 cases no peak was detected. See also http://www.openswath.org/openswath_data.html for details.",
                  pubMedIds = "24727770",
                  url = "http://compms.org/resources/reference-data/50",
                  instrumentModel = "AB SCIEX TripleTOF 5600 System",
                  instrumentManufacturer = "AB SCIEX",
                  ionSource = "ESI/IonDrive™ Turbo V, DuoSpray™ source",
                  analyser = "QTOF/TripleTOF",
                  detectorType = "microchannel plate/4-channel MCP detector",
                  softwareName = "OpenSWATH",
                  collisionEnergy = "",
                  dateStamp = "2014-03-07"
)

## Expression data
e <- exprs(e.sgs.human)

## Experiment info
## nrow(phenoData) == 54 == dim(e)[2]
toName <- paste0(colnames(e))
colnames(e) <- toName
pd <- data.frame(toName,
                 row.names=colnames(e))
pd <- new("AnnotatedDataFrame", pd)


## feature data
fd <- sgshumantxt$Sequence
fd <- as.data.frame(fd)
colnames(fd) <- "sequence"
fd$Filename <- sgshumantxt$filename
fd$sortslot <- seq_along(fd$sequence)
# left outer join "fd" with "sis.peptides"
fd <- merge(x = fd, y = sis.peptides, by = "sequence", all.x = TRUE)
# Sort "fd" by "sort slot" column

## Annotate featureData with Study_Design table: Study_design.sgs.human
fd <- merge(x = fd, y = Study_design.sgs.human, by = "Filename", all.x = TRUE)
fd <- fd[order(fd$sortslot),]
rownames(fd) <- rownames(sgshumantxt)
fd$sortslot <- NULL
fd <- new("AnnotatedDataFrame", fd)

## The "MSnProcess" Class
process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("No Normalisation")),
               normalised=FALSE)


Rost2014sgs <- new("MSnSet",
                      exprs = e,
                      phenoData = pd,
                      experimentData = experiment,
                      featureData = fd,
                      processingData = process)



## checks
stopifnot(dim(pData(Rost2014sgs))[1] == ncol(e),
          dim(fData(Rost2014sgs))[1] == nrow(e),
          validObject(Rost2014sgs))

save(Rost2014sgs, file="../../data/Rost2014Humansgs.rda",
     compress = "xz", compression_level = 9)

