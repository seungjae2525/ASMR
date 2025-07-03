asmr <- function(...) UseMethod("asmr")

#' @title Analysis of selection-biased matched data
#'
#' @description \code{asmr()} is the main function of \bold{ASMR} package.
#'
#' @param data A dataset that contains a binary or count outcome (\code{Y}) and the binary exposure (\code{X}). If the data is \sQuote{long} format (i.e., for the same subject, measurements before and after the intervention are placed in separate rows), the \code{ID} argument is required. In contrast, if the data is \sQuote{wide} format (i.e., for the same subject, measurements before and after the intervention are placed in separate columns within the same row), the \code{ID} argument is not required.
#' @param Y Variable name(s) of the binary or count outcome. This argument must be a character string for \sQuote{long} format data or a character vector for \sQuote{wide} format data. See Examples.
#' @param X Variable name(s) of the binary exposure. This argument must be a character string for \sQuote{long} format data or a character vector for \sQuote{wide} format data. See Examples.
#' @param ID Variable name(s) which identifies subjects. This argument must be a character string. The length of \code{ID} should be the same as the number of observations (only considered if data is \sQuote{long} format).
#' @param Lambda1 Lower bound of sensitivity parameter (\eqn{\Lambda_{1}}). This value must be positive and less than \code{Lambda2} (i.e., \eqn{0 < \Lambda_{1} < \Lambda_{2}}).
#' @param Lambda2 Upper bound of sensitivity parameter (\eqn{\Lambda_{2}}). This value must be positive and greater than \code{Lambda1} (i.e., \eqn{0 < \Lambda_{1} < \Lambda_{2}}).
#' @param level The confidence level required to obtain the confidence interval for the population sensitivity interval. Default: \code{level = 0.95}.
#' @param method The method to be used. Set \code{method = "CLR"} for binary outcomes and \code{method = "CPR"} for count outcomes. See Details.
#' @param useparallel Logical scalar indicating whether to parallelize our optimization problem for count outcomes. Default: \code{useparallel = FALSE} (only considered if \code{method = "CPR"}).
#' @param n.cores The number of CPU cores to use. Default: \code{n.cores = parallel::detectCores()/2} (ignored if \code{method = "CLR"}).
#'
#'
#' @details
#'   Let \eqn{(Y_{i1}, Y_{i2})} denote the self-matched observations for subject \eqn{i \in \{1, \ldots, N\}} in a finite population, where \eqn{Y_{it} = 1} if the \eqn{i}th individual has the disease of interest at time \eqn{t \in \{1, 2\}} and \eqn{Y_{it} = 0} otherwise.
#'   For example, \eqn{t = 1} may correspond to the pre-intervention period and \eqn{t = 2} to the post-intervention period, or vice versa.
#'   Define \eqn{x_{it}} as a binary indicator where \eqn{x_{it} = 1} if the \eqn{i}th individual received the intervention at time \eqn{t} and \eqn{x_{it} = 0} otherwise.
#'   Assume that \eqn{n} subjects are selected from a finite population of \eqn{N} subjects for the analysis.
#'   For formal notation of the selected sample, define \eqn{S_{i}} as a selection indicator such that \eqn{S_{i} = 1} if the \eqn{i}th individual is included in the sample and \eqn{S_{i} = 0} otherwise, so that \eqn{\sum_{i = 1}^{N} S_{i} = n}.
#'   In \bold{ASMR} package, there are two methods that can be used:
#'
#'   * \code{method = "CLR"}: This is a method for analysis of selection-biased binary matched data.
#'   In this case, the ratio of the selection probabilities for \eqn{i}th individual, conditional on two different outcome patterns observed at time points \eqn{1} and \eqn{2}, is
#'   \deqn{R_{i} = \frac{P(S_{i} = 1 \mid Y_{i1} = 1, Y_{i2} = 0, x_{it}, v_{i})}{P(S_{i} = 1 \mid Y_{i1} = 0, Y_{i2} = 1, x_{it}, v_{i})}.}
#'   Then, for a given \eqn{0 < \Lambda_{1} < \Lambda_{2}}, the set of selection ratio model for the conditional logistic regression is defined by
#'   \deqn{\mathcal{R}(\Lambda_{1}, \Lambda_{2}) = \left\{ R_{i}: \Lambda_{1} \leq R_{i} \leq \Lambda_{2}, \; \forall i \in \mathcal{S} \right\}}
#'   where we assume that the true ratio of the selection probabilities \eqn{R_{i,0}} belongs to \eqn{\mathcal{R}(\Lambda_{1}, \Lambda_{2})} for all \eqn{i \in \mathcal{S} = \{i: S_{i} = 1\}}.
#'   See Section 2.2 in Lee et al. (2025).
#'
#'   * \code{method = "CPR"}: This is a method for analysis of selection-biased count matched data.
#'   In this case, the ratio of the selection probabilities for \eqn{i}th individual, conditional on two different outcome patterns observed at time points \eqn{1} and \eqn{2}, is
#'   \deqn{R_{ik}^{\ast} = \frac{P(S_{i} = 1 \mid Y_{i1} = k, Y_{i2} = y_{i1} + y_{i2} - k, x_{it}, v_{i})}{P(S_{i} = 1 \mid Y_{i1} = y_{i1}, Y_{i2} = y_{i2}, x_{it}, v_{i})}}
#'   for \eqn{k \in \{0, 1, \ldots, y_{i1} + y_{i2}\}}.
#'   Then, for a given \eqn{0 < \Lambda_{1} < \Lambda_{2}}, the set of selection ratio model for the conditional Poisson regression is defined by
#'   \deqn{\mathcal{R}^{\ast}(\Lambda_{1}, \Lambda_{2}) = \left\{ R_{ik}^{\ast}: \Lambda_{1} \leq R_{ik}^{\ast} \leq \Lambda_{2}, \; \forall i \in \mathcal{S}, \; \forall k\right\}}
#'   where we assume that the true ratio of the selection probabilities \eqn{R_{ik,0}^{\ast}} belongs to \eqn{\mathcal{R}^{\ast}(\Lambda_{1}, \Lambda_{2})} for all \eqn{i} and \eqn{k}.
#'   See Section 2.3 in Lee et al. (2025).
#'
#'
#' Here, \eqn{\Lambda_{1}} and \eqn{\Lambda_{2}} are user-specified sensitivity parameters that reflect the selection by \eqn{Y_{i1}}, \eqn{Y_{i2}}, \eqn{x_{it}}, and \eqn{v_{i}}.
#' Based on domain knowledge, we assume that for all \eqn{i} (and \eqn{k}), \eqn{R_{i}} (or \eqn{R_{ik}^{\ast}}) lies between \eqn{\Lambda_{1}} and \eqn{\Lambda_{2}}.
#' Then, we can obtain the sensitivity interval which is the minimum and maximum values of \eqn{\widehat{\beta}} over the provided range of sensitivity parameters (i.e., \eqn{\Lambda_{1}} and \eqn{\Lambda_{2}}).
#' Also, we can obtain the \eqn{100 \times (1 - \alpha)\%} confidence interval for the population sensitivity interval.
#' See Lee et al. (2025) for details.
#'
#'
#' @return An object of class \code{asmr}. The object is a list with the following components:
#' \item{data}{The \code{data} argument.}
#' \item{Y}{The \code{Y} argument.}
#' \item{X}{The \code{X} argument.}
#' \item{ID}{The \code{ID} argument.}
#' \item{Lambda1}{The \code{Lambda1} argument (\eqn{\Lambda_{1}}).}
#' \item{Lambda2}{The \code{Lambda2} argument (\eqn{\Lambda_{2}}).}
#' \item{EST}{It is the maximum likelihood estimate (MLE) \eqn{\widehat{\beta}} from conditional logistic or Poisson regression for selection-biased matched data.}
#' \item{SE}{It is the standard error \eqn{\text{se}(\widehat{\beta})} from conditional logistic regression or Poisson for selection-biased matched data.}
#' \item{SI}{Sensitivity interval.}
#' \item{CI}{Confidence interval for the population sensitivity interval.}
#'
#' In addition, the following additional components exist for \code{method = "CPR"}:
#' \item{Rik_unique_comb_len}{The number of all possible combinations of \eqn{\bm{R}^{\ast}} for optimization problem.}
#' \item{SE_robust}{The robust standard errors for the lower and upper bounds of the sensitivity interval.}
#' \item{Rik}{The set of sensitivity parameters \eqn{\bm{R}^{\ast}} for the lower and upper bounds of the sensitivity interval.}
#'
#' The results for the \code{asmr} are printed with the \code{\link[ASMR]{print.asmr}} function.
#' Also, the results for the \code{asmr} are summarized with the \code{\link[ASMR]{summary.asmr}} function.
#'
#'
#' @examples
#' ## Load example data
#' data("backPain", package = "gnm")
#' ## Re-express as count data
#' backPainLong <- gnm::expandCategorical(backPain, "pain")
#'
#'
#' @seealso
#'  \code{\link[ASMR]{print.asmr}}, \code{\link[ASMR]{summary.asmr}}
#'
#'
#' @references
#' Lee, S., Sim, H., and Lee, W. (2025). Developing statistical methods for selection-biased matched data.
#'
#'
#' @keywords methods
#'
#'
#' @export

##
asmr <- function(data, Y, X, ID = NULL, Lambda1, Lambda2, level = 0.95, method,
                 useparallel = FALSE, n.cores = parallel::detectCores()/2) {
  ## Check if arguments are included in the column name of data!
  if (sum(Y %notin% colnames(x = data)) != 0) {
    stop("Outcome variable name(s) (\"Y\") must be included in data.")
  }
  if (sum(X %notin% colnames(x = data)) != 0) {
    stop("Exposure variable name(s) (\"X\") must be included in data.")
  }

  ## Check outcome, exposure, ID variables
  if (length(x = Y) != length(x = X)) {
    stop("The length of \"Y\" must be equal to the length of \"X\". See Examples")
  }
  if (length(x = Y) > 2) {
    stop("The length of \"Y\" and the length of \"X\" cannot be greater than 2.")
  }
  if (!is.null(x = ID)) {
    if (length(x = ID) > 1) {
      stop("The length of \"ID\" cannot be greater than 1.")
    }
  }

  ## Check sensitivity parameters (Lambda1 and Lambda2)
  if (length(x = Lambda1) != 1 | length(x = Lambda2) != 1) {
    stop("The lengths of both \"Lambda1\" and \"Lambda2\" must be 1.")
  }
  if (Lambda1 <= 0 | Lambda2 <= 0) {
    stop("\"Lambda1\" and \"Lambda2\" must both be positive real numbers.")
  }
  if (Lambda1 > Lambda2) {
    stop("\"Lambda2\" must be larger than \"Lambda1\".")
  }

  ##
  if (method == "CLR") {
    #### Analysis of selection-biased binary matched data
    ## Outcome and exposure variables
    if (length(x = Y) == 1) {
      if (sum(data[, Y] %notin% c(0, 1)) != 0) {
        stop("When method = \"CLR\", all outcome values in \"Y\" must be 0 or 1.")
      }
      if (sum(data[, X] %notin% c(0, 1)) != 0) {
        stop("All exposure values in \"X\" must be 0 or 1.")
      }
      if (is.factor(x = data[, X])) {
        data[, X] <- as.numeric(x = data[, X]) - 1
      }
      if (is.null(x = ID)) {
        stop("If both outcome variable name and exposure variable name have length 1 (i.e., if the data is in wide format), the \"ID\" argument is required.")
      } else if (ID %notin% colnames(x = data)) (
        stop("ID variable name (\"ID\") must be included in data.")
      )
      if (is.factor(x = data[, ID])) {
        data[, ID] <- as.character(x = data[, ID])
      }
    } else {
      if (sum(data[, Y[1]] %notin% c(0, 1)) + sum(data[, Y[2]] %notin% c(0, 1)) != 0) {
        stop("When method = \"CLR\", all outcome values in \"Y\" must be 0 or 1.")
      }
      if (sum(data[, X[1]] %notin% c(0, 1)) + sum(data[, X[2]] %notin% c(0, 1)) != 0) {
        stop("All exposure values must be 0 or 1.")
      }
      if (is.factor(x = data[, X[1]])) {
        data[, X[1]] <- as.numeric(x = data[, X[1]]) - 1
      }
      if (is.factor(x = data[, X[2]])) {
        data[, X[2]] <- as.numeric(x = data[, X[2]]) - 1
      }
    }

    ## Analysis
    result.analysis <- asmrCLR(data = data, Y = Y, X = X, ID = ID,
                               Lambda1 = Lambda1, Lambda2 = Lambda2, level = level)

  } else if (method == "CPR") {
    #### Analysis of selection-biased count matched data
    ## Outcome and exposure variables
    if (length(x = Y) == 1) {
      if (sum(data[, Y] < 0) != 0) {
        stop("When method = \"CPR\", all outcome values in \"Y\" must be greater than or equal to 0.")
      }
      if (sum(data[, X] %notin% c(0, 1)) != 0) {
        stop("All exposure values in \"X\" must be 0 or 1.")
      }
      if (is.factor(x = data[, X])) {
        data[, X] <- as.numeric(x = data[, X]) - 1
      }
      if (is.null(x = ID)) {
        stop("If both outcome variable name and exposure variable name have length 1 (i.e., if the data is in wide format), the \"ID\" argument is required.")
      } else if (ID %notin% colnames(x = data)) (
        stop("ID variable name (\"ID\") must be included in data.")
      )
      if (is.factor(x = data[, ID])) {
        data[, ID] <- as.character(x = data[, ID])
      }
    } else {
      if (sum(data[, Y[1]] < 0) + sum(data[, Y[2]] < 0) != 0) {
        stop("When method = \"CPR\", all outcome values in \"Y\" must be greater than or equal to 0.")
      }
      if (sum(data[, X[1]] %notin% c(0, 1)) + sum(data[, X[2]] %notin% c(0, 1)) != 0) {
        stop("All exposure values in \"X\" must be 0 or 1.")
      }
      if (is.factor(x = data[, X[1]])) {
        data[, X[1]] <- as.numeric(x = data[, X[1]]) - 1
      }
      if (is.factor(x = data[, X[2]])) {
        data[, X[2]] <- as.numeric(x = data[, X[2]]) - 1
      }
    }

    ##
    if (useparallel == TRUE) {
      if (n.cores > parallel::detectCores()) {
        stop("\"n.cores\" is specified to be greater than the number of CPU cores on the current host.")
      }
    }

    ## Analysis
    result.analysis <- asmrCPR(data = data, Y = Y, X = X, ID = ID,
                               Lambda1 = Lambda1, Lambda2 = Lambda2, level = level,
                               useparallel = useparallel, n.cores = n.cores)

  } else {
    stop("\"method\" must be \"CLR\" for binary outcomes or \"CPR\" for count outcomes.")
  }

  ## Return results
  return(result.analysis)
}



