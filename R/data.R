#' @title msdata2
#' @aliases msdata2
#' @description  This function lists the data sets available in \code{msdata2}
#'   package by calling \code{data(package = "msdata2")}
#' @usage msdata2()
#' @importFrom utils data packageVersion
#' @import MSnbase
#' @references   See in the respective data sets' manual pages for references to
#' publications.
#' @author Chong Tang <chong.tang@@clouvain.be>
#' @export

msdata2 <- function()
  data(package = "msdata2")


