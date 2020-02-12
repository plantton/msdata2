library("MSnbase")
library("readxl")
library("dplyr")

## Set workdir as the current path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## read quantitative results
dataRamus <- "../extdata/1-s2.0-S187439191530186X-mmc1.xlsx"

## There's no value missing for spectral counting workflows: 1 - 4
exdata1 <- read_xlsx(dataRamus, range = cell_rows(1:2728), sheet = 1)
exdata2 <- read_xlsx(dataRamus, range = cell_rows(1:2803), sheet = 2)
exdata3 <- read_xlsx(dataRamus, range = cell_rows(1:2976), sheet = 3)
exdata4 <- read_xlsx(dataRamus, range = cell_rows(1:2381), sheet = 4)
exdata5 <- read_xlsx(dataRamus, range = cell_rows(1:2722), sheet = 5)
exdata6 <- read_xlsx(dataRamus, range = cell_rows(1:2645), sheet = 6)
exdata7 <- read_xlsx(dataRamus, range = cell_rows(1:2626), sheet = 7)
exdata8 <- read_xlsx(dataRamus, range = cell_rows(1:2620), sheet = 8)

## Workflow 1: Spectral Counting - ExtractMSn, Mascot, MFPaQ
workflow1 <- readMSnSet2(file = exdata1, ecol = 5:10)

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
                 ConcentrationLevel = c(rep('A',3),rep('B', 3)) ,
                 QuantificationMethod = "Spectral counting",
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
ramus2016SCMfPaQ <- new("MSnSet",
                          exprs = e1,
                          phenoData = pd1,
                          experimentData = experiment1,
                          featureData = fd1,
                          processingData = process1)

stopifnot(dim(pData(ramus2016SCMfPaQ))[1] == ncol(e1),
          dim(fData(ramus2016SCMfPaQ))[1] == nrow(e1),
          validObject(ramus2016SCMfPaQ))

## Data 2: SC MaxQuant
workflow2 <- readMSnSet2(file = exdata2, ecol = 5:10)

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
                  ConcentrationLevel = c(rep('A',3),rep('B', 3)),
                  QuantificationMethod = "Spectral counting",
                  row.names=colnames(e2))
pd2 <- new("AnnotatedDataFrame", pd2)

## feature data
fd2 <- fData(workflow2)
names(fd2) <- c(names(fd2)[1:8],"Filtering 1", "Filtering 2",
                "Filtering 3", "BH rank")
fd2 <- new("AnnotatedDataFrame", fd2)

## The "MSnProcess" Class
process2 <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".", sep=""),
                  paste("No Normalisation")),
                normalised=FALSE)

##
ramus2016SCMaxQuant <- new("MSnSet",
                        exprs = e2,
                        phenoData = pd2,
                        experimentData = experiment2,
                        featureData = fd2,
                        processingData = process2)

stopifnot(dim(pData(ramus2016SCMaxQuant))[1] == ncol(e2),
          dim(fData(ramus2016SCMaxQuant))[1] == nrow(e2),
          validObject(ramus2016SCMaxQuant))


## Data 3: SC Irma-Heidi
workflow3 <- readMSnSet2(file = exdata3, ecol = 5:10)

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
                  ConcentrationLevel = c(rep('A',3),rep('B', 3)),
                  QuantificationMethod = "Spectral counting",
                  row.names=colnames(e3))
pd3 <- new("AnnotatedDataFrame", pd3)

## feature data
fd3 <- fData(workflow3)
names(fd3) <- c(names(fd3)[1:8],"Filtering 1", "Filtering 2",
                "Filtering 3", "BH rank")
fd3 <- new("AnnotatedDataFrame", fd3)

## The "MSnProcess" Class
process3 <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".", sep=""),
                  paste("No Normalisation")),
                normalised=FALSE)

##
ramus2016SCIrmaHeidi <- new("MSnSet",
                           exprs = e3,
                           phenoData = pd3,
                           experimentData = experiment3,
                           featureData = fd3,
                           processingData = process3)

stopifnot(dim(pData(ramus2016SCIrmaHeidi))[1] == ncol(e3),
          dim(fData(ramus2016SCIrmaHeidi))[1] == nrow(e3),
          validObject(ramus2016SCIrmaHeidi))


## Data 4: SC Scaffold
workflow4 <- readMSnSet2(file = exdata4, ecol = 5:10)

## Experiment info
experiment4 <- new("MIAPE",
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
                   softwareName = c("ExtractMSn", "Mascot", "Scaffold")
)

## combine {BiocGenerics} **
## Expression data
e4 <- exprs(workflow4)

## Experiment info
pd4 <- data.frame(Replicate = c(seq(1,3),seq(1,3)),
                  ConcentrationLevel = c(rep('A',3),rep('B', 3)) ,
                  QuantificationMethod = "Spectral counting",
                  row.names=colnames(e4))
pd4 <- new("AnnotatedDataFrame", pd4)

## feature data
fd4 <- fData(workflow4)
names(fd4) <- c(names(fd4)[1:8],"Filtering 1", "Filtering 2",
                "Filtering 3", "BH rank")
fd4 <- new("AnnotatedDataFrame", fd4)

## The "MSnProcess" Class
process4 <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".", sep=""),
                  paste("No Normalisation")),
                normalised=FALSE)

##
ramus2016SCScaffold <- new("MSnSet",
                            exprs = e4,
                            phenoData = pd4,
                            experimentData = experiment4,
                            featureData = fd4,
                            processingData = process4)

stopifnot(dim(pData(ramus2016SCScaffold))[1] == ncol(e4),
          dim(fData(ramus2016SCScaffold))[1] == nrow(e4),
          validObject(ramus2016SCScaffold))

## For Spectral Sounting data, only raw data provided.
save(ramus2016SCMfPaQ, file="../../data/ramus2016SCMfPaQ.rda",
     compress = "xz", compression_level = 9)

save(ramus2016SCMaxQuant, file="../../data/ramus2016SCMaxQuant.rda",
     compress = "xz", compression_level = 9)

save(ramus2016SCIrmaHeidi, file="../../data/ramus2016SCIrmaHeidi.rda",
     compress = "xz", compression_level = 9)

save(ramus2016SCScaffold, file="../../data/ramus2016SCScaffold.rda",
     compress = "xz", compression_level = 9)


## Workflow 5 - 8
## Quantification Method: MS signal analysis

## Workflow 5: ExtractMSn/Mascot/MFPaQ/MFPaQ
## MS signal extraction device - MFPaQ
workflow5 <- readMSnSet2(file = exdata5, ecol = 11:16)

## Expression data
e5 <- exprs(workflow5)

## Missingness in protein intensity values
table(is.na(e5))  ## 107 missing values

## Experiment info
experiment5 <- new("MIAPE",
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
                   softwareName = c("ExtractMSn", "Mascot", "MFPaQ")
)

## Experiment info
pd5 <- data.frame(Replicate = c(seq(1,3),seq(1,3)),
                  ConcentrationLevel = c(rep('A',3),rep('B', 3)) ,
                  QuantificationMethod = "MS signal analysis",
                  row.names=colnames(e5))
pd5 <- new("AnnotatedDataFrame", pd5)

## feature data
fd5 <- fData(workflow5)
fd5 <- new("AnnotatedDataFrame", fd5)

## The "MSnProcess" Class
process5 <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".", sep=""),
                  paste("No Normalisation")),
                normalised=FALSE)


##
ramus2016IntMFPaQ <- new("MSnSet",
                           exprs = e5,
                           phenoData = pd5,
                           experimentData = experiment5,
                           featureData = fd5,
                           processingData = process5)

stopifnot(dim(pData(ramus2016IntMFPaQ))[1] == ncol(e5),
          dim(fData(ramus2016IntMFPaQ))[1] == nrow(e5),
          validObject(ramus2016IntMFPaQ))


## Workflow 6: Andromeda/Andromeda/MaxQuant/MaxQuant(Intensity)
## MS signal extraction device - MaxQuant (Intensity)
workflow6 <- readMSnSet2(file = exdata6, ecol = 12:17)

## Expression data
e6 <- exprs(workflow6)

## Missingness in protein intensity values
table(is.na(e6))  ## 579 missing values

## Experiment info
experiment6 <- new("MIAPE",
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

## Experiment info
pd6 <- data.frame(Replicate = c(seq(1,3),seq(1,3)),
                  ConcentrationLevel = c(rep('A',3),rep('B', 3)) ,
                  QuantificationMethod = "MS signal analysis",
                  row.names=colnames(e6))
pd6 <- new("AnnotatedDataFrame", pd6)

## feature data
fd6 <- fData(workflow6)
fd6 <- new("AnnotatedDataFrame", fd6)

## The "MSnProcess" Class
process6 <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".", sep=""),
                  paste("No Normalisation")),
                normalised=FALSE)

##
ramus2016IntMaxQuant<- new("MSnSet",
                         exprs = e6,
                         phenoData = pd6,
                         experimentData = experiment6,
                         featureData = fd6,
                         processingData = process6)

stopifnot(dim(pData(ramus2016IntMaxQuant))[1] == ncol(e6),
          dim(fData(ramus2016IntMaxQuant))[1] == nrow(e6),
          validObject(ramus2016IntMaxQuant))

## Workflow 7: Andromeda/Andromeda/MaxQuant/MaxQuant(LFQ)
## MS signal extraction device - MaxQuant (LFQ)
workflow7 <- readMSnSet2(file = exdata7, ecol = 12:17)

## Expression data
e7 <- exprs(workflow7)

## Missingness in protein intensity values
table(is.na(e7))  ## 611 missing values

## Experiment info
experiment7 <- new("MIAPE",
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

## Experiment info
pd7 <- data.frame(Replicate = c(seq(1,3),seq(1,3)),
                  ConcentrationLevel = c(rep('A',3),rep('B', 3)) ,
                  QuantificationMethod = "MS signal analysis",
                  row.names=colnames(e7))
pd7 <- new("AnnotatedDataFrame", pd7)

## feature data
fd7 <- fData(workflow7)
fd7 <- new("AnnotatedDataFrame", fd7)

## The "MSnProcess" Class
process7 <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".", sep=""),
                  paste("No Normalisation")),
                normalised=FALSE)

##
ramus2016LFQMaxQuant<- new("MSnSet",
                           exprs = e7,
                           phenoData = pd7,
                           experimentData = experiment7,
                           featureData = fd7,
                           processingData = process7)

stopifnot(dim(pData(ramus2016LFQMaxQuant))[1] == ncol(e7),
          dim(fData(ramus2016LFQMaxQuant))[1] == nrow(e7),
          validObject(ramus2016LFQMaxQuant))


## Workflow 8: Mascot Distiller/Mascot/Scaffold/Skyline
## MS signal extraction device - Skyline
workflow8 <- readMSnSet2(file = exdata8, ecol = 21:26)

## Expression data
e8 <- exprs(workflow8)

## Missingness in protein intensity values
table(is.na(e8))  ## 0 missing values

## Experiment info
experiment8 <- new("MIAPE",
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
                   softwareName = c("Mascor Distiller", "Mascot", "Scaffold")
)

## Experiment info
pd8 <- data.frame(Replicate = c(seq(1,3),seq(1,3)),
                  ConcentrationLevel = c(rep('A',3),rep('B', 3)) ,
                  QuantificationMethod = "MS signal analysis",
                  row.names=colnames(e8))
pd8 <- new("AnnotatedDataFrame", pd8)

## feature data
fd8 <- fData(workflow8)
fd8 <- new("AnnotatedDataFrame", fd8)

## The "MSnProcess" Class
process8 <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".", sep=""),
                  paste("No Normalisation")),
                normalised=FALSE)

##
ramus2016IntSkyline<- new("MSnSet",
                           exprs = e8,
                           phenoData = pd8,
                           experimentData = experiment8,
                           featureData = fd8,
                           processingData = process8)

stopifnot(dim(pData(ramus2016IntSkyline))[1] == ncol(e8),
          dim(fData(ramus2016IntSkyline))[1] == nrow(e8),
          validObject(ramus2016IntSkyline))

## For Spectral Sounting data, only raw data provided.
save(ramus2016IntMFPaQ, file="../../data/ramus2016IntMFPaQ.rda",
     compress = "xz", compression_level = 9)

save(ramus2016IntMaxQuant, file="../../data/ramus2016IntMaxQuant.rda",
     compress = "xz", compression_level = 9)

save(ramus2016LFQMaxQuant, file="../../data/ramus2016LFQMaxQuant.rda",
     compress = "xz", compression_level = 9)

save(ramus2016IntSkyline, file="../../data/ramus2016IntSkyline.rda",
     compress = "xz", compression_level = 9)


## Data cleaning for MSnSet generated from MS signal analysis
library(tidyverse)

## clean dataset from workflow 5
e5 <- exprs(ramus2016IntMFPaQ)

##
class(e5)
dim(e5)
pData(ramus2016IntMFPaQ)
glimpse(e5)
