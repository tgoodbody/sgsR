.onAttach <- function(libname, pkgname) {
  # Runs when attached to search() path such as by library() or require()
  if (!interactive()) {
    return(invisible())
  }

  v <- utils::packageVersion("sgsR")
  packageStartupMessage("sgsR v", v, ". Bug report on <github.com/tgoodbody/sgsR/issues>.")
}
