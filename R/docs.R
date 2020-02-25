
#' @title SWATH-MS Gold Standard Dataset
#'
#' @name Rost2014sgs
#'
#' @description
#' The SWATH-MS Gold Standard (SGS) dataset consists of 90 SWATH-MS runs of 422
#' synthetic stable isotope-labeled standard (SIS) peptides in ten different
#' dilution steps (1, 2, 4, 8, ..., 512 times), spiked into three protein
#' backgrounds of varying complexity (water, yeast and human), acquired in
#' three technical replicates. The SGS dataset was manually annotated,
#' resulting in 342 identified and quantified peptides with three or four
#' transitions each. In total, 30,780 chromatograms were inspected and 18,785
#' were annotated with one true peak group, whereas in 11,995 cases no peak was
#' detected.
#'
#' See also http://www.openswath.org/openswath_data.html for details.
#'
#' @aliases Rost2014sgs
#'          Rost2014sgsHuman
#'          Rost2014sgsYeast
#'          Rost2014sgsWater
#'
#' @docType data
#'
#' @usage
#' data(Rost2014sgsHuman)
#' data(Rost2014sgsYeast)
#' data(Rost2014sgsWater)
#'
#' @format The data is an instance of class \code{MSnSet} from package \code{MSnbase}.
#'
#' @details
#'
#' @keywords datasets
#'
#' @references   \emph{OpenSWATH enables automated, targeted analysis of data-independent
#' acquisition MS data.} Röst HL, Rosenberger G, Navarro P, Gillet L, Miladinović SM,
#' Schubert OT, Wolski W, Collins BC, Malmström J, Malmström L, Aebersold R.
#' Nat Biotechnol. 2014 Mar;32   (3):219-23. doi: 10.1038/nbt.2841.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/23979570}{PubMed})
#'
#'
#' @source \href{http://compms.org/resources/reference-data/50}{SWATH-MS Gold Standard Dataset}
#'
#' @examples
#' library(msdata2)
#' data(Rost2014sgsHuman)
#' Rost2014sgsHuman
#' pData(Rost2014sgsHuman)
#' head(exprs(Rost2014sgsHuman))
NULL


#' @title Spiked Proteomic Standard Datasets For Benchmark 8 Different
#' label-free Quantitative Workflows
#'
#' @name ramus2016
#'
#' @description
#' This data consists of quantitative results from 8 different Label-Free
#' quantitative workflows for a controlled, spiked proteomic standard dataset.
#' The proteomics standard data was prepared using a yeast cell lysate
#' accompanied by a serial dilution (respectively 0.05-0.125-0.250-0.5-2.5-5-
#' 12.5-25-50 fmol of UPS1 \eqn{/\mu}g of yeast lysate) of the UPS1 standard mixture
#' (Sigma). Samples were then analyzed in triplicate by nanoLC-MS/MS coupled to
#' an LTQ-Orbitrap Velos mass spectrometer, operated in data-dependent
#' acquisition mode.
#'
#' The dataset was then processed by 8 different workflows for benchmarking,
#' consisting in the following steps: peaklist generation, database search,
#' validation of the identified proteins and extraction of quantitative metric
#' (spectral count or MS signal).
#'
#' @aliases ramus2016
#'          ramus2016SCMfPaQ
#'          ramus2016SCMaxQuant
#'          ramus2016SCIrmaHeidi
#'          ramus2016SCScaffold
#'          ramus2016IntMFPaQ
#'          ramus2016IntMaxQuant
#'          ramus2016LFQMaxQuant
#'          ramus2016IntSkyline
#'
#' @docType data
#'
#' @usage
#' data(ramus2016SCMfPaQ)
#' data(ramus2016SCMaxQuant)
#' data(ramus2016SCIrmaHeidi)
#' data(ramus2016SCScaffold)
#' data(ramus2016IntMFPaQ)
#' data(ramus2016IntMaxQuant)
#' data(ramus2016LFQMaxQuant)
#' data(ramus2016IntSkyline)
#'
#' @format The data is an instance of class \code{MSnSet} from package \code{MSnbase}.
#'
#' @details
#' Full details of the experimental design can be found in the reference.\cr
#'
#' The 8 different datasets correspond to 8 different LFQ workflows:\cr
#' @section
#' ramus2016SCMfPaQ:
#'              \verb{}Peaklist creation device - ExtractMSn\cr
#'                     Database search engine - Mascot\cr
#'                     Validation/Spectral counting device - MFPaQ\cr
#'                     Quantification method - Spectral counting\cr
#' @section
#' ramus2016SCMaxQuant:
#'                     Peaklist creation device - Andromeda\cr
#'                     Database search engine - Andromeda\cr
#'                     Validation/Spectral counting device - MaxQuant\cr
#'                     Quantification method - Spectral counting\cr
#' @section
#' ramus2016SCIrmaHeidi:
#'                     Peaklist creation device - Mascot Distiller\cr
#'                     Database search engine - Mascot\cr
#'                     Validation/Spectral counting device - IRMa/hEIDI\cr
#'                     Quantification method - Spectral counting\cr
#' @section
#' ramus2016SCScaffold:
#'                     Peaklist creation device - ExtractMSn\cr
#'                     Database search engine - Mascot\cr
#'                     Validation/Spectral counting device - Scaffold\cr
#'                     Quantification method - Spectral counting\cr
#' @section
#' ramus2016IntMFPaQ:
#'                     Peaklist creation device - ExtractMSn\cr
#'                     Database search engine - Mascot\cr
#'                     Validation/Spectral counting device - MFPaQ\cr
#'                     MS signal extraction device - MFPaQ\cr
#'                     Quantification method - MS signal analysis\cr
#' @section
#' ramus2016IntMaxQuant:
#'                     Peaklist creation device - Andromeda\cr
#'                     Database search engine - Andromeda\cr
#'                     Validation/Spectral counting device - MaxQuant\cr
#'                     MS signal extraction device - MaxQuant (Intensity)\cr
#'                     Quantification method - MS signal analysis\cr
#' @section
#' ramus2016LFQMaxQuant:
#'                     Peaklist creation device - Andromeda\cr
#'                     Database search engine - Andromeda\cr
#'                     Validation/Spectral counting device - MaxQuant\cr
#'                     MS signal extraction device - MaxQuant (LFQ)\cr
#'                     Quantification method - MS signal analysis\cr
#' @section
#' ramus2016IntSkyline:
#'                     Peaklist creation device - Mascot Distiller\cr
#'                     Database search engine - Mascot\cr
#'                     Validation/Spectral counting device - Scaffold\cr
#'                     MS signal extraction device - Skyline\cr
#'                     Quantification method - MS signal analysis\cr
#'
#'
#'
#' @keywords datasets
#'
#' @references   \emph{Benchmarking quantitative label-free LC-MS data
#' processing workflows using a complex spiked proteomic standard dataset.}
#' Ramus C, Hovasse A, Marcellin M, Hesse AM, Mouton-Barbosa E, Bouyssié D,
#' Vaca S, Carapito C, Chaoui K, Bruley C, Garin J, Cianférani S, Ferro M,
#' Van Dorssaeler A, Burlet-Schiltz O, Schaeffer C, Couté Y, Gonzalez de Peredo A.
#' J Proteomics. 2016 Jan 30;132:51-62. \cr
#' doi: 10.1016/j.jprot.2015.11.011.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/26585461}{PubMed})
#'
#'
#' @source
#' This data was generated from Sup Table 1. Quantitative data obtained from
#' the 8 different workflows. See inst/scripts/ramus2016.R for details.
#'
#' @examples
#' library(msdata2)
#' data(ramus2016IntSkyline)
#' ramus2016IntSkyline
#' pData(ramus2016IntSkyline)
#' head(exprs(ramus2016IntSkyline))
NULL



#' @title Quantitative Proteomics of the Cancer Cell Line Encyclopedia
#'
#' @name CCLENormalized
#'
#' @description
#' In this dataset, we include quantitative proteome profiling by mass
#' spectrometry of 375 cell lines from Cancer Cell Line Encyclopedia. The
#' experiment was performed in multiplex format with 9 biological samples per
#' plex and one common sample to normalize between plexes. The normalized data
#' contains 42 multiplex experiments (TMT) which consists of 504 runs on the
#' mass spectrometer and over 1,500 hours of instrument time.
#'
#' Ref: https://gygi.med.harvard.edu/publications/ccle
#'
#' @aliases CCLENormalized
#'
#' @docType data
#'
#' @usage
#' data(CCLENormalized)
#'
#' @format The data is an instance of class \code{MSnSet} from package \code{MSnbase}.
#'
#' @details
#'
#' @keywords datasets
#'
#' @references   \emph{Quantitative Proteomics of the Cancer Cell Line Encyclopedia.}
#' Nusinow, D., Szpyt, J., Ghandi, M., Rose, C., McDonald, E., Kalocsay, M.,
#' Jané-Valbuena, J., Gelfand, E., Schweppe, D., Jedrychowski, M., Golji, J.,
#' Porter, D., Rejtar, T., Wang, Y., Kryukov, G., Stegmeier, F., Erickson, B.,
#' Garraway, L., Sellers, W. and Gygi, S.
#' Cell. Volume 180, Issue 2, 23 January 2020, Pages 387-402.e16
#' doi: https://doi.org/10.1016/j.cell.2019.12.023.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/31978347}{PubMed})
#'
#'
#' @source \href{https://gygi.med.harvard.edu/publications/ccle}{Quantitative Proteomics of the Cancer Cell Line Encyclopedia}
#'
#' @examples
#' library(msdata2)
#' data(CCLENormalized)
#' CCLENormalized
#' pData(CCLENormalized)
#' head(exprs(CCLENormalized))
NULL











