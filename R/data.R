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
#' data(Rost2014Humansgs)
#' data(Rost2014Humansgs)
##' msdata2metadata(Rost2014Humansgs)
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
