
#' @title SWATH-MS Gold Standard Dataset
#'
#' @name Rost2014sgs
#'
#' @description
#' The SWATH-MS Gold Standard (SGS) dataset consists of 90 SWATH-MS runs of 422 synthetic
#' stable isotope-labeled standard (SIS) peptides in ten different dilution steps (1, 2, 4,
#' 8, ..., 512 times), spiked into three protein backgrounds of varying complexity (water,
#' yeast and human), acquired in three technical replicates. The SGS dataset was manually
#' annotated, resulting in 342 identified and quantified peptides with three or four
#' transitions each. In total, 30,780 chromatograms were inspected and 18,785 were annotated
#' with one true peak group, whereas in 11,995 cases no peak was detected.
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







