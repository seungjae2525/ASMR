#' @title Summary for \code{asmr} objects
#'
#' @description Summary the results for object of class \code{asmr}.
#'
#' @param object An object for class \code{asmr}.
#' @param digits Number of significant digits. Default: \code{max(1L, getOption("digits") - 3L)}.
#' @param ... Further arguments (currently not used).
#'
#' @details Summary the results for object of class \code{asmr}.
#' From the conditional logistic regression results, the odds ratio, confidence interval of the odds ratio, and \emph{P} value are reported.
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' ## Load example data
#' data("backPain", package = "gnm")
#' ## Re-express as count data
#' backPainLong <- gnm::expandCategorical(backPain, "pain")
#'
#' @keywords summary
#'
#' @seealso
#'  \code{\link[ASMR]{asmr}}, \code{\link[ASMR]{print.asmr}}, \code{\link[base]{summary}}
#'
#' @export

summary.asmr <- function(object, digits = max(1L, getOption("digits") - 3L), ...) {
  if (!(inherits(object, "asmrCLR") | inherits(object, "asmrCPR"))){
    stop("Argument 'object' must be an object of class \"asmrCLR\" or \"asmrCPR\".")
  } else if (inherits(object, "asmrCLR")) {

  } else if (inherits(object, "asmrCLR")) {

  }

}
