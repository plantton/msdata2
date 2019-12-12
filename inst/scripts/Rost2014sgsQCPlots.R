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

## We also need the Quantification Result from OpenSWATH for Human Background
sgsHumanTxtfile <- paste("../../inst/extdata/OpenSWATH_SM3_GoldStandard",
                         "AutomatedResults_human_peakgroups.txt",
                         sep = "")
sgsHumanTxt <- read.delim(sgsHumanTxtfile)

## Create Replicate Names for Human-background OpenSWATH data
## Make this "ReplicateName" feature is consistent for the 2 datasets
sgsHumanTxt$ReplicateNames <- gsub("split_(.*)_combined.*",
                                  "\\1", sgsHumanTxt$filename)
## Select Human-background rows from the Skylike manual data
sgsManualDataHuman <- sgsManualData %>%
                        filter(ReplicateName %in% unique(sgsHumanTxt$ReplicateNames))

## Modify "Retention Time" column in "second" metric
sgsManualDataHuman$RetentionTime <- ifelse(sgsManualDataHuman$RetentionTime == "#N/A",
                                           NA,
                                           as.numeric(levels(sgsManualDataHuman$RetentionTime))[sgsManualDataHuman$RetentionTime]*60)


## To compute the identification accuracy, all results reported by mProphet
## above a certain cutoff were taken and it was counted how many of these
## results were false positive (no manual annotation present) or mis-identified
## or mis-identified (manual annotation present but at a different) retention
## time, in our case further than 30 seconds away).

## False positive count
falsePotiveHumanOS <- sgsHumanTxt %>%
                                select(ReplicateNames, Sequence, RT) %>%
                                group_split(ReplicateNames)

falsePotiveHumanManual <- sgsManualDataHuman %>%
                                select(ReplicateName, PeptideSequence, RetentionTime) %>%
                                group_split(ReplicateName)
idenAccurayHuman <- data.frame("ReplicateName" = rep(NA, 30),
                               "FalsePoCount" = rep(NA, 30),
                               "MisIden" = rep(NA, 30))
for (i in 1:30) {
        if (unique(falsePotiveHumanOS[[i]]$ReplicateNames) == unique(falsePotiveHumanManual[[i]]$ReplicateName)) {
              token <- setdiff(unique(falsePotiveHumanOS[[i]]$Sequence), unique(falsePotiveHumanManual[[i]]$PeptideSequence))
              idenAccurayHuman[i,] <- c(unique(falsePotiveHumanOS[[i]]$ReplicateNames), length(token), NA)
        }
}

## Remove records from "falsePotiveHumanManual" which has NA in RT
falsePotiveHumanManualRT <-subset(sgsManualDataHuman,
                                  (!is.na(sgsManualDataHuman$RetentionTime)))
## Group by ReplicateName
misIdenHumanManual <- falsePotiveHumanManualRT %>%
                    group_by(ReplicateName, PeptideSequence) %>%
                    summarise(RT = mean(RetentionTime))
misIdenHumanOS <- sgsHumanTxt %>% select(ReplicateNames, Sequence, RT)

misIdenData <- merge(x = misIdenHumanOS, y = misIdenHumanManual,
                     by.x = c("ReplicateNames", "Sequence"),
                     by.y = c("ReplicateName", "PeptideSequence"),
                     all.y = TRUE)
misIdenData <- misIdenData[which(abs(misIdenData$RT.x - misIdenData$RT.y) > 30), ]
##
misCountHuman <- misIdenData %>% group_by(ReplicateNames) %>%
                     summarise(MisIden = n())

## Final count
idenAccurayHuman$MisIden <- NULL
idenAccurayHuman <- merge(idenAccurayHuman, misCountHuman,
                          by.x = "ReplicateName",
                          by.y = "ReplicateNames",
                          all.x = T)
##
sumSgsHumanTxt <- sgsHumanTxt %>% group_by(ReplicateNames) %>%
                      summarise(total = n())

##
sumAllSGSHuman <- merge(idenAccurayHuman, sumSgsHumanTxt,
                        by.x = "ReplicateName",
                        by.y = "ReplicateNames",
                        all.x = T)
##
sumAllSGSHuman$FalsePoCount <- as.numeric(sumAllSGSHuman$FalsePoCount)
sumAllSGSHuman$AccuRate <- 1 - (sumAllSGSHuman$FalsePoCount + sumAllSGSHuman$MisIden)/sumAllSGSHuman$total




## Heatmap
library(tidyverse)
m3 <- as.data.frame(exprs(Rost2014sgsHuman))
##
m3 <- m3 %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(m3)
##
p <- ggplot(m3,aes(x=colname,y=rowname,fill=value))+
  geom_tile()
##
p1 <- p + labs(x="",y="") +
  geom_tile(colour="white",size=0.25) +
  #theme options
  theme(
    #bold font for legend text
    legend.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank())
##
test <- exprs(object = Rost2014sgsHuman)
colnames(test) <- rep(1:10,3)
test.1 <- colnames(test)
gplots::heatmap.2(test, Colv = F, col = rainbow, dendrogram = 'none')



##
library(ggplot2)
library(GGally)
pairDataFrame <- as.data.frame(exprs(Rost2014sgsHuman)[, c(1, 10, 11, 20, 21, 30)],
                               col.names = c(1, 10, 11, 20, 21, 30))
ggpairs(pairDataFrame,
        upper = list(continuous = "density", combo = "box_no_facet"))
ggpairs(pairDataFrame)


##
library(msdata2)
exDFHuman <- exprs(Rost2014sgsHuman)
exDFYeast <- exprs(Rost2014sgsYeast)
exDFWater <- exprs(Rost2014sgsWater)

## Imputation
exDFHuman <- knnImputation(exDFHuman, k = 10, scale = T, meth = "weighAvg",
                              distData = NULL)
exDFYeast <- knnImputation(exDFYeast, k = 10, scale = T, meth = "weighAvg",
                              distData = NULL)  ## No missingness
exDFWater <- knnImputation(exDFWater, k = 10, scale = T, meth = "weighAvg",
                              distData = NULL)

#
exDFHumanAgg <- exDFHuman[,1:10] + exDFHuman[,11:20] + exDFHuman[,21:30]
exDFYeastAgg <- exDFYeast[,1:10] + exDFYeast[,11:20] + exDFYeast[,21:30]
exDFWaterAgg <- exDFWater[,1:10] + exDFWater[,11:20] + exDFWater[,21:30]


## log
exDFHumanAggLog <- log(exDFHumanAgg)
exDFYeastAggLog <- log(exDFYeastAgg)
exDFWaterAggLog <- log(exDFWaterAgg)

## Standardize to the most indense concentration within each run
exDFHumanNorm <- BBmisc::normalize(exDFHumanAggLog,
                                   method = "standardize",
                                   range = c(0, 1), margin = 1L)
exDFYeastNorm <- BBmisc::normalize(exDFYeastAggLog,
                                   method = "standardize",
                                   range = c(0, 1), margin = 1L)
exDFWaterNorm <- BBmisc::normalize(exDFWaterAggLog,
                                   method = "standardize",
                                   range = c(0, 1), margin = 1L)

##
library(scales)
library(tidyverse)

##
df <-
  tibble(y = (-10:10),
         x = (y^4)*sign(y))


##
log_both <- function(x){ifelse(x == 0, 0, log(abs(x)) * sign(x))}
exp_both <- function(x){exp(abs(x)) * sign(x)} # this is the inverse of log_both

log_both_trans <-
  function(){
    trans_new(name = 'log_both',
              transform = log_both,
              inverse = exp_both)
  }


#transformed
ggplot(df) +
  #no transformation
  geom_point(aes(factor(x), y = 1, fill = x), shape = 21, size = 10) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                       guide = guide_colorbar(order = 1)) +
  ylim(.75, 1.25) +
    labs(colour = "transformed", fill = "default", x = "", y = "")
##
ggplot(df) +
  #no transformation
  geom_tile(aes(factor(x), y = 1, fill = x)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                       guide = guide_colorbar(order = 1)) +
  scale_x_discrete(labels = factor(df$x)) +
  labs(colour = "transformed", fill = "default", x = "", y = "")
