.onAttach <- function(libname, pkgname) {
  msg <- paste0("\nThis is msdata2 version ",
                packageVersion("msdata2"), ".\n",
                "Use 'msdata2()' to list available data sets.")
  packageStartupMessage(msg)
}
