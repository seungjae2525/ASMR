#' @title Print for \code{asmr} objects
#'
#' @description Print the results for object of class \code{asmr}.
#'
#' @param x An object for class \code{asmr}.
#' @param ... Further arguments (currently not used).
#'
#' @details Print the results for object of class \code{asmr}.
#' From the conditional logistic regression results, \dQuote{Estimate} corresponds to the log odds ratio and \dQuote{SE} corresponds to the standard error of the log odds ratio.
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' ## Load example data
#' data("backPain", package = "gnm")
#' ## Re-express as count data
#' backPainLong <- gnm::expandCategorical(backPain, "pain")
#'
#' @keywords print
#'
#' @seealso
#'  \code{\link[ASMR]{asmr}}, \code{\link[ASMR]{summary.asmr}}, \code{\link[base]{print}}
#'
#' @export

print.asmr <- function(x, ...) {
  if (!(inherits(x, "asmrCLR") | inherits(x, "asmrCPR"))){
    stop("Argument 'x' must be an object of class \"asmrCLR\" or \"asmrCPR\".")
  } else if (inherits(x, "asmrCLR")) {

  } else if (inherits(x, "asmrCLR")) {

  }

}
