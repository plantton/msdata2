library("MSnbase")
library("readr")
library("dplyr")
library("tidyr")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sgsyeasttxtfile <- "../../inst/extdata/OpenSWATH_SM3_GoldStandardAutomatedResults_yeast_peakgroups.txt"
sgsyeasttxt <- read.delim(sgsyeasttxtfile)
## import 422 stable isotope-labeled standard (SIS) peptides data
sispeptidefile <- "../extdata/NBT-L30171D-OpenSWATH_SM2_GoldStandardPeptides.csv"
sis.peptides <- read.csv(sispeptidefile,sep=";")
names(sis.peptides) <- c("Sequence", names(sis.peptides)[2:3])


## Create a "Study design" table for experimental design description
Study_design.sgs.yeast <- data.frame(filename = unique(sgsyeasttxt$filename))
# Add "Condition" column: Control vs Disease
# "Condition" represents different dilution steps
Study_design.sgs.yeast$Condition <- gsub(".*_(.*)_SW.*", "\\1", Study_design.sgs.yeast$filename)
# Manually modify the table, according to the published reference **
Study_design.sgs.yeast$Condition[10] <- "010"
# "BioReplicate" column
Study_design.sgs.yeast$BioReplicate <- gsub(".*_(.*)_0.*", "\\1", Study_design.sgs.yeast$filename)
# Manually modify the table, according to the published reference **
Study_design.sgs.yeast$BioReplicate[10] <- "L120228"
Study_design.sgs.yeast$BioReplicate <- factor(Study_design.sgs.yeast$BioReplicate)
levels(Study_design.sgs.yeast$BioReplicate) <- seq(nlevels(Study_design.sgs.yeast$BioReplicate))
# "Run" column:
Study_design.sgs.yeast$Run <- seq_along(Study_design.sgs.yeast$filename)

## Create Peptide level aggregation data frame as Expression assay data
## left outer join "sgsyeasttxt" with "sis.peptides"
sgsyeasttxt$sortslot <- seq_along(sgsyeasttxt$Sequence)
sgsyeasttxt <- merge(x = sgsyeasttxt, y = sis.peptides, by = "Sequence", all.x = TRUE)
## Annotate sgs data frame with Study_Design table:
sgsyeasttxt <- merge(x = sgsyeasttxt, y = Study_design.sgs.yeast, by = "filename", all.x = TRUE)
sgsyeasttxt <- sgsyeasttxt[order(sgsyeasttxt$sortslot),]
## assign rownames to the data frame
rownames(sgsyeasttxt) <- gsub("AQUA4SWATH_(.*)/2_run0_split_napedro_(.*)_SW_combined.*",
                              "\\1_\\2_", sgsyeasttxt$transition_group_id)
sgsyeasttxt$sortslot <- NULL


## Column - "total_xic" as assay data, reshape the original data frame
## Peptide sequences existed in all 30 elution runs?
sequence.RunOverview <- sgsyeasttxt %>% group_by(Sequence) %>% summarise(n())
## Sequence: No missing values. **
## Arrange the data frame: firstly group all rows on "Sequence";
## then arrange all rows within one "Sequence" category on "Run"
sgsyeasttxt <- sgsyeasttxt %>%
  group_by(Sequence) %>%
  arrange(Run) %>%
  arrange(Sequence)

## Reshape "total_xic" to new matrix: each row represents one peptide sequence
exprs.xic <- sgsyeasttxt %>%
  select(Sequence, Run, total_xic) %>%
  spread(Run, total_xic, fill = NA) %>%
  as.data.frame()
rownames(exprs.xic) <- exprs.xic$Sequence
colnames(exprs.xic)[!names(exprs.xic) %in% "Sequence"] <- as.character(Study_design.sgs.yeast$filename[match(colnames(exprs.xic)[-1],
                                                                                                             Study_design.sgs.yeast$Run)])


## Reshape the dataframe as Sequence/Peptide level data
## With columns as 30 "run"s

## Import the entire data matrix as MsnSet - assayData
## Note: "transition_group_id", "id" is unique for each row
##       "filename" contains 30 different strings. - Check SGS data description.
##       "Sequence", "FullPeptideName" ,"aggr_Fragment_Annotation" contain 345 different levels/numbers.
##       "ProteinName" contains 16 different levels.
## Attention: As suggested by "SWATH2stats":
##             For some R packages such as  imsbInfer,
##             All the columns should be preserved.
##             Hence I will try to keep as much info from OpenSWATH results as possible.
e.sgs.yeast <- readMSnSet2(file = exprs.xic, ecol = 2:31)

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
e <- exprs(e.sgs.yeast)

## Experiment info
## nrow(phenoData) == dim(e)[2]
pd <- data.frame(Study_design.sgs.yeast,
                 row.names=colnames(e))
rownames(pd) <- pd$filename
pd <- new("AnnotatedDataFrame", pd)

## feature data
## Drop redundant columns from original OpenSWATH output
sgsyeasttxt <- sgsyeasttxt %>%
  select(-c(filename,  # Put columns' names here to remove them
            transition_group_record,
            decoy,
            transition_group_id,
            run_id,
            FullPeptideName,
            xx_swath_prelim_score,
            aggr_Peak_Apex,
            aggr_Fragment_Annotation,
            log10_total_xic,
            peak_group_rank,
            Condition,
            total_xic,
            id,
            BioReplicate))

## Function to reshape column into data frame (fData)
## Create an empty assay
toName <- colnames(sgsyeasttxt[,!(names(sgsyeasttxt) %in% c("Sequence", "Run"))])
fd <- data.frame(toName,
                 row.names = colnames(sgsyeasttxt[,!(names(sgsyeasttxt) %in% c("Sequence", "Run"))]))
fd <- t(fd)

## Columns contains redundant info
# Charge: Always 2
#  ProteinName: n_distinct() == 1
#  nr_peaks: n_distinct() == 1
#  organism: n_distinct() == 1
#  Picked: n_distinct() == 1
fd <- fd[, !colnames(fd) %in% c("Charge",
                                "ProteinName",
                                "nr_peaks",
                                "organism",
                                "picked"),  drop = FALSE]

## fData
f.df <- array(0, dim = c(dim(e)[1], dim(fd)[2]))
f.df <- data.frame(f.df, row.names = rownames(exprs.xic))

# colnames(f.df) <- colnames(fd)
for (i in 1:dim(fd)[2]) {
  i.mat <- sgsyeasttxt %>%
    select(Sequence, Run, colnames(fd)[i]) %>%
    spread(Run, colnames(fd)[i], fill = NA)
  f.df[, i] <- as.matrix(i.mat[, -1])
}
colnames(f.df) <- colnames(fd)

## Add redundant/repeated values
df.1 <- sgsyeasttxt %>% select(Sequence, Charge, ProteinName, nr_peaks, organism, picked) %>%
  group_by(Sequence) %>%
  summarise(Charge = unique(Charge),
            ProteinName = unique(ProteinName),
            nr_peaks = unique(nr_peaks),
            organism = unique(organism),
            picked = unique(picked)) %>%
  arrange(Sequence) %>%
  as.data.frame()
rownames(df.1) <- df.1$Sequence
# rownames(df.1) <- df.1$Sequence
f.df <- cbind(f.df, df.1[,c("Charge",
                            "ProteinName",
                            "nr_peaks",
                            "organism",
                            "picked")][match(rownames(f.df), rownames(df.1)),])
f.df <- new("AnnotatedDataFrame", f.df)


## The "MSnProcess" Class
process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("No Normalisation")),
               normalised=FALSE)


Rost2014sgsYeast <- new("MSnSet",
                   exprs = e,
                   phenoData = pd,
                   experimentData = experiment,
                   featureData = f.df,
                   processingData = process)
## Normalise
#  Rost2014sgsYeast <- normalise(Rost2014sgsYeast, method = "sum")

## checks
stopifnot(dim(pData(Rost2014sgsYeast))[1] == ncol(e),
          dim(fData(Rost2014sgsYeast))[1] == nrow(e),
          validObject(Rost2014sgsYeast))

save(Rost2014sgsYeast, file="../../data/Rost2014Yeastsgs.rda",
     compress = "xz", compression_level = 9)
