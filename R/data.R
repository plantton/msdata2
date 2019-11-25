msdata2 <- function()
  data(package = "msdata2")


#' Extracts relevant metadata from an \code{MSnSet} instance. See
#' \code{README.md} for a description and explanation of the metadata
#' fields.
#'
#' @title Extract msdata2 metadata
#' @param x A \code{msdata2} data.
#' @return An instance of class \code{msdata2metadata}.
#' @author Chong Tang
#' @aliases print.msdata2metadata
#' @examples
#' library("msdata2")
#' data(Rost2014sgsHuman)
#' data(Rost2014sgsHuman)
##' msdata2metadata(Rost2014sgsHuman)
msdata2metadata <- function(x) {
  ans <- list(Species = experimentData(x)@samples$species,
              Tissue = experimentData(x)@samples$tissue,
              CellLine = ifelse(experimentData(x)@samples$tissue == "Cell",
                                experimentData(x)@samples$cellLine, NA),
              PMID = pubMedIds(x),
              MS = otherInfo(experimentData(x))$MS,
              Experiment = otherInfo(experimentData(x))$spatexp,
              MarkerCol = otherInfo(experimentData(x))$markers.fcol,
              PredictionCol = otherInfo(experimentData(x))$prediction.fcol)
  class(ans) <- c("list", "msdata2metadata")
  ans
}

print.msdata2metadata <- function(x, ...) {
  cat("msdata2 experiment metadata:\n")
  nx <- names(x)
  for (i in nx)
    cat(paste0(" ", i, ": ", x[[i]], "\n"))
}

valid.msdata2metadata <- function(x) {
  stopifnot(inherits(x, "msdata2metadata"))
  !any(sapply(x, is.null))
}

#' @title SWATH-MS Gold Standard Dataset
#'
#' @description
#' The SWATH-MS Gold Standard (SGS) dataset consists of 90 SWATH-MS runs of 422 synthetic
#' stable isotope-labeled standard (SIS) peptides in ten different dilution steps, spiked
#' into three protein backgrounds of varying complexity (water, yeast and human), acquired
#' in three technical replicates. The SGS dataset was manually annotated, resulting in 342
#' identified and quantified peptides with three or four transitions each. In total, 30,780
#' chromatograms were inspected and 18,785 were annotated with one true peak group, whereas
#' in 11,995 cases no peak was detected.
#' See also http://www.openswath.org/openswath_data.html for details.
#'
#' @aliases Rost2014sgsHuman
#'
#' @docType data
#'
#' @usage data(Rost2014sgsHuman)
#'
#' @format The data is an instance of class \code{MSnSet} from package \code{MSnbase}.
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
#' data(Rost2014sgsHuman)
#' Rost2014sgsHuman
#' exprs(Rost2014sgsHuman)
"Rost2014sgsHuman"
