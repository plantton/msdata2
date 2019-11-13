library("MSnbase")
library("readxl")
library("dplyr")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  data.ramus <- "../extdata/1-s2.0-S187439191530186X-mmc1.xlsx"


exdata.1 <- read_xlsx(data.ramus, sheet = 1)
exdata.2 <- read_xlsx(data.ramus, sheet = 2)
exdata.3 <- read_xlsx(data.ramus, range = cell_rows(1:2976), sheet = 3)
exdata.4 <- read_xlsx(data.ramus, range = cell_rows(1:2381), sheet = 4)
exdata.5 <- read_xlsx(data.ramus, range = cell_rows(1:2722), sheet = 5)
exdata.6 <- read_xlsx(data.ramus, range = cell_rows(1:2645), sheet = 6)
exdata.7 <- read_xlsx(data.ramus, range = cell_rows(1:2626), sheet = 7)
exdata.8 <- read_xlsx(data.ramus, range = cell_rows(1:2620), sheet = 8)

## Data 1: SC MFPaQ
workflow1 <- readMSnSet2(file = exdata.1, ecol = 5:10)

## Experimental data
## Experiment info
experiment1 <- new("MIAPE",
                  samples = list(
                    species = c("Human","Yeast")
                    ),
                  title = "Benchmarking quantitative label-free LC-MS data processing workflows using a complex spiked proteomic standard dataset.",
                  abstract = 'Proteomic workflows based on nanoLC-MS/MS data-dependent-acquisition analysis have progressed tremendously in recent years. High-resolution and fast sequencing instruments have enabled the use of label-free quantitative methods, based either on spectral counting or on MS signal analysis, which appear as an attractive way to analyze differential protein expression in complex biological samples. However, the computational processing of the data for label-free quantification still remains a challenge. Here, we used a proteomic standard composed of an equimolar mixture of 48 human proteins (Sigma UPS1) spiked at different concentrations into a background of yeast cell lysate to benchmark several label-free quantitative workflows, involving different software packages developed in recent years. This experimental design allowed to finely assess their performances in terms of sensitivity and false discovery rate, by measuring the number of true and false-positive (respectively UPS1 or yeast background proteins found as differential). The spiked standard dataset has been deposited to the ProteomeXchange repository with the identifier PXD001819 and can be used to benchmark other label-free workflows, adjust software parameter settings, improve algorithms for extraction of the quantitative metrics from raw MS data, or evaluate downstream statistical methods.
                              BIOLOGICAL SIGNIFICANCE: Bioinformatic pipelines for label-free quantitative analysis must be objectively evaluated in their ability to detect variant proteins with good sensitivity and low false discovery rate in large-scale proteomic studies. This can be done through the use of complex spiked samples, for which the "ground truth" of variant proteins is known, allowing a statistical evaluation of the performances of the data processing workflow. We provide here such a controlled standard dataset and used it to evaluate the performances of several label-free bioinformatics tools (including MaxQuant, Skyline, MFPaQ, IRMa-hEIDI and Scaffold) in different workflows, for detection of variant proteins with different absolute expression levels and fold change values. The dataset presented here can be useful for tuning software tool parameters, and also testing new algorithms for label-free quantitative analysis, or for evaluation of downstream statistical methods.',
                  pubMedIds = "26585461",
                  instrumentModel = "LTQ Velos-Orbitrap",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "ESI",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = c("ExtractMSn", "Mascot Search Engine", "MFPaQ")
                  )

## combine {BiocGenerics} **
## Expression data
e1 <- exprs(workflow1)

## Experiment info
pd1 <- data.frame(Replicate = c(seq(1,3),seq(1,3)),
                 Concentration.Level = c(rep('A',3),rep('B', 3)) ,
                 Quantification.Method = "Spectral counting",
                 row.names=colnames(e1))
pd1 <- new("AnnotatedDataFrame", pd1)

## feature data
fd1 <- fData(workflow1)
names(fd1) <- c(names(fd1)[1:8],"Filtering 1", "Filtering 2",
                "Filtering 3", "BH rank")
fd1 <- new("AnnotatedDataFrame", fd1)

# The "MSnProcess" Class
process1 <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".", sep=""),
                 paste("No Normalisation")),
               normalised=FALSE)

##
ramus2015SCMfPaQ <- new("MSnSet",
                          exprs = e1,
                          phenoData = pd1,
                          experimentData = experiment1,
                          featureData = fd1,
                          processingData = process1)

stopifnot(validObject(ramus2015SCMfPaQ))

## Data 2: SC MaxQuant
workflow2 <- readMSnSet2(file = exdata.2, ecol = 5:10)

## Experiment info
experiment2 <- new("MIAPE",
                   samples = list(
                     species = c("Human","Yeast")
                   ),
                   title = "Benchmarking quantitative label-free LC-MS data processing workflows using a complex spiked proteomic standard dataset.",
                   abstract = 'Proteomic workflows based on nanoLC-MS/MS data-dependent-acquisition analysis have progressed tremendously in recent years. High-resolution and fast sequencing instruments have enabled the use of label-free quantitative methods, based either on spectral counting or on MS signal analysis, which appear as an attractive way to analyze differential protein expression in complex biological samples. However, the computational processing of the data for label-free quantification still remains a challenge. Here, we used a proteomic standard composed of an equimolar mixture of 48 human proteins (Sigma UPS1) spiked at different concentrations into a background of yeast cell lysate to benchmark several label-free quantitative workflows, involving different software packages developed in recent years. This experimental design allowed to finely assess their performances in terms of sensitivity and false discovery rate, by measuring the number of true and false-positive (respectively UPS1 or yeast background proteins found as differential). The spiked standard dataset has been deposited to the ProteomeXchange repository with the identifier PXD001819 and can be used to benchmark other label-free workflows, adjust software parameter settings, improve algorithms for extraction of the quantitative metrics from raw MS data, or evaluate downstream statistical methods.
                              BIOLOGICAL SIGNIFICANCE: Bioinformatic pipelines for label-free quantitative analysis must be objectively evaluated in their ability to detect variant proteins with good sensitivity and low false discovery rate in large-scale proteomic studies. This can be done through the use of complex spiked samples, for which the "ground truth" of variant proteins is known, allowing a statistical evaluation of the performances of the data processing workflow. We provide here such a controlled standard dataset and used it to evaluate the performances of several label-free bioinformatics tools (including MaxQuant, Skyline, MFPaQ, IRMa-hEIDI and Scaffold) in different workflows, for detection of variant proteins with different absolute expression levels and fold change values. The dataset presented here can be useful for tuning software tool parameters, and also testing new algorithms for label-free quantitative analysis, or for evaluation of downstream statistical methods.',
                   pubMedIds = "26585461",
                   instrumentModel = "LTQ Velos-Orbitrap",
                   instrumentManufacturer = "ThermoScientific",
                   ionSource = "ESI",
                   analyser = "Orbitrap",
                   detectorType = "Orbitrap",
                   softwareName = c("Andromeda", "Andromeda", "MaxQuant")
)

## combine {BiocGenerics} **
## Expression data
e2 <- exprs(workflow2)

## Experiment info
pd2 <- data.frame(Replicate = c(seq(1,3),seq(1,3)),
                  Concentration.Level = c(rep('A',3),rep('B', 3)) ,
                  Quantification.Method = "Spectral counting",
                  row.names=colnames(e2))
pd2 <- new("AnnotatedDataFrame", pd2)

## feature data
fd2 <- fData(workflow2)
names(fd2) <- c(names(fd2)[1:8],"Filtering 1", "Filtering 2",
                "Filtering 3", "BH rank")
fd2 <- new("AnnotatedDataFrame", fd2)

# The "MSnProcess" Class
process2 <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".", sep=""),
                  paste("No Normalisation")),
                normalised=FALSE)

##
ramus2015SCMaxQuant <- new("MSnSet",
                        exprs = e2,
                        phenoData = pd2,
                        experimentData = experiment2,
                        featureData = fd2,
                        processingData = process2)

stopifnot(validObject(ramus2015SCMaxQuant))

## Data 3: SC Irma-Heidi
workflow3 <- readMSnSet2(file = exdata.3, ecol = 5:10)

## Experiment info
experiment3 <- new("MIAPE",
                   samples = list(
                     species = c("Human","Yeast")
                   ),
                   title = "Benchmarking quantitative label-free LC-MS data processing workflows using a complex spiked proteomic standard dataset.",
                   abstract = 'Proteomic workflows based on nanoLC-MS/MS data-dependent-acquisition analysis have progressed tremendously in recent years. High-resolution and fast sequencing instruments have enabled the use of label-free quantitative methods, based either on spectral counting or on MS signal analysis, which appear as an attractive way to analyze differential protein expression in complex biological samples. However, the computational processing of the data for label-free quantification still remains a challenge. Here, we used a proteomic standard composed of an equimolar mixture of 48 human proteins (Sigma UPS1) spiked at different concentrations into a background of yeast cell lysate to benchmark several label-free quantitative workflows, involving different software packages developed in recent years. This experimental design allowed to finely assess their performances in terms of sensitivity and false discovery rate, by measuring the number of true and false-positive (respectively UPS1 or yeast background proteins found as differential). The spiked standard dataset has been deposited to the ProteomeXchange repository with the identifier PXD001819 and can be used to benchmark other label-free workflows, adjust software parameter settings, improve algorithms for extraction of the quantitative metrics from raw MS data, or evaluate downstream statistical methods.
                              BIOLOGICAL SIGNIFICANCE: Bioinformatic pipelines for label-free quantitative analysis must be objectively evaluated in their ability to detect variant proteins with good sensitivity and low false discovery rate in large-scale proteomic studies. This can be done through the use of complex spiked samples, for which the "ground truth" of variant proteins is known, allowing a statistical evaluation of the performances of the data processing workflow. We provide here such a controlled standard dataset and used it to evaluate the performances of several label-free bioinformatics tools (including MaxQuant, Skyline, MFPaQ, IRMa-hEIDI and Scaffold) in different workflows, for detection of variant proteins with different absolute expression levels and fold change values. The dataset presented here can be useful for tuning software tool parameters, and also testing new algorithms for label-free quantitative analysis, or for evaluation of downstream statistical methods.',
                   pubMedIds = "26585461",
                   instrumentModel = "LTQ Velos-Orbitrap",
                   instrumentManufacturer = "ThermoScientific",
                   ionSource = "ESI",
                   analyser = "Orbitrap",
                   detectorType = "Orbitrap",
                   softwareName = c("Mascot Distiller", "Mascot", "IRMa/hEIDI")
)

## combine {BiocGenerics} **
## Expression data
e3 <- exprs(workflow3)

## Experiment info
pd3 <- data.frame(Replicate = c(seq(1,3),seq(1,3)),
                  Concentration.Level = c(rep('A',3),rep('B', 3)) ,
                  Quantification.Method = "Spectral counting",
                  row.names=colnames(e3))
pd3 <- new("AnnotatedDataFrame", pd3)

## feature data
fd3 <- fData(workflow3)
names(fd3) <- c(names(fd3)[1:8],"Filtering 1", "Filtering 2",
                "Filtering 3", "BH rank")
fd3 <- new("AnnotatedDataFrame", fd3)

# The "MSnProcess" Class
process3 <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".", sep=""),
                  paste("No Normalisation")),
                normalised=FALSE)

##
ramus2015SCIrmaHeidi <- new("MSnSet",
                           exprs = e3,
                           phenoData = pd3,
                           experimentData = experiment3,
                           featureData = fd3,
                           processingData = process3)

stopifnot(validObject(ramus2015SCIrmaHeidi))

## Data 4: SC Irma-Heidi
workflow3 <- readMSnSet2(file = exdata.3, ecol = 5:10)

## Experiment info
experiment3 <- new("MIAPE",
                   samples = list(
                     species = c("Human","Yeast")
                   ),
                   title = "Benchmarking quantitative label-free LC-MS data processing workflows using a complex spiked proteomic standard dataset.",
                   abstract = 'Proteomic workflows based on nanoLC-MS/MS data-dependent-acquisition analysis have progressed tremendously in recent years. High-resolution and fast sequencing instruments have enabled the use of label-free quantitative methods, based either on spectral counting or on MS signal analysis, which appear as an attractive way to analyze differential protein expression in complex biological samples. However, the computational processing of the data for label-free quantification still remains a challenge. Here, we used a proteomic standard composed of an equimolar mixture of 48 human proteins (Sigma UPS1) spiked at different concentrations into a background of yeast cell lysate to benchmark several label-free quantitative workflows, involving different software packages developed in recent years. This experimental design allowed to finely assess their performances in terms of sensitivity and false discovery rate, by measuring the number of true and false-positive (respectively UPS1 or yeast background proteins found as differential). The spiked standard dataset has been deposited to the ProteomeXchange repository with the identifier PXD001819 and can be used to benchmark other label-free workflows, adjust software parameter settings, improve algorithms for extraction of the quantitative metrics from raw MS data, or evaluate downstream statistical methods.
                              BIOLOGICAL SIGNIFICANCE: Bioinformatic pipelines for label-free quantitative analysis must be objectively evaluated in their ability to detect variant proteins with good sensitivity and low false discovery rate in large-scale proteomic studies. This can be done through the use of complex spiked samples, for which the "ground truth" of variant proteins is known, allowing a statistical evaluation of the performances of the data processing workflow. We provide here such a controlled standard dataset and used it to evaluate the performances of several label-free bioinformatics tools (including MaxQuant, Skyline, MFPaQ, IRMa-hEIDI and Scaffold) in different workflows, for detection of variant proteins with different absolute expression levels and fold change values. The dataset presented here can be useful for tuning software tool parameters, and also testing new algorithms for label-free quantitative analysis, or for evaluation of downstream statistical methods.',
                   pubMedIds = "26585461",
                   instrumentModel = "LTQ Velos-Orbitrap",
                   instrumentManufacturer = "ThermoScientific",
                   ionSource = "ESI",
                   analyser = "Orbitrap",
                   detectorType = "Orbitrap",
                   softwareName = c("Mascot Distiller", "Mascot", "IRMa/hEIDI")
)

## combine {BiocGenerics} **
## Expression data
e3 <- exprs(workflow3)

## Experiment info
pd3 <- data.frame(Replicate = c(seq(1,3),seq(1,3)),
                  Concentration.Level = c(rep('A',3),rep('B', 3)) ,
                  Quantification.Method = "Spectral counting",
                  row.names=colnames(e3))
pd3 <- new("AnnotatedDataFrame", pd3)

## feature data
fd3 <- fData(workflow3)
names(fd3) <- c(names(fd3)[1:8],"Filtering 1", "Filtering 2",
                "Filtering 3", "BH rank")
fd3 <- new("AnnotatedDataFrame", fd3)

# The "MSnProcess" Class
process3 <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".", sep=""),
                  paste("No Normalisation")),
                normalised=FALSE)

##
ramus2015SCIrmaHeidi <- new("MSnSet",
                            exprs = e3,
                            phenoData = pd3,
                            experimentData = experiment3,
                            featureData = fd3,
                            processingData = process3)

stopifnot(validObject(ramus2015SCIrmaHeidi))

