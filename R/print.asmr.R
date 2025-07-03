#' @title Print for \code{asmr} objects
#'
#' @description Print the results for object of class \code{asmr}.
#'
#' @param x An object for class \code{asmr}.
#' @param logscale Whether to report results on a log odds ratio or log risk ratio scale. Default: \code{logscale = FALSE} (i.e., not log scale).
#' @param digits Print digits. Default:  \code{digits = max(1L, getOption("digits") - 3L)}.
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
#'  \code{\link[ASMR]{asmr}}, \code{\link[base]{print}}
#'
#' @export

print.asmr <- function(x, logscale = FALSE,
                       digits = max(1L, getOption("digits") - 3L), ...) {
  if (!(inherits(x, "asmrCLR") | inherits(x, "asmrCPR"))) {
    stop("Argument 'x' must be an object of class \"asmrCLR\" or \"asmrCPR\".")
  }
  xx <- x

  if (is.null(x = xx$ID)) {
    # data_type = "wide"
    data_temp <- xx$data
  } else {
    # data_type = "long"
    data_temp <- make_data_wide(data = xx$data, Y = xx$Y, X = xx$X, ID = xx$ID)
  }

  ##
  cat("##### Summary of data ##### \n")
  cat(paste0("## Number of subjects: ", nrow(x = data_temp), " \n"))

  cat("## Distribution of exposures: \n")
  tmp1 <- as.data.frame(x = expand.grid(x1 = c(0, 1), x2 = c(0, 1)))
  tmp1$count <- c(nrow(x = data_temp[data_temp$x1 == 0 & data_temp$x2 == 0, ]),
                  nrow(x = data_temp[data_temp$x1 == 1 & data_temp$x2 == 0, ]),
                  nrow(x = data_temp[data_temp$x1 == 0 & data_temp$x2 == 1, ]),
                  nrow(x = data_temp[data_temp$x1 == 1 & data_temp$x2 == 1, ]))
  print(x = tmp1)
  cat("where x\u209C (0 or 1) is an intervention incator at time t (1 or 2). \n")
  cat("Here, t = 1 may correspond to the pre-intervention period and \n")
  cat("t = 2 to the post-intervention period, or vice versa. \n")
  cat("\n")

  ##
  data_temp_01 <- data_temp[data_temp$x1 == 0 & data_temp$x2 == 1, c(3, 4)]
  data_temp_10 <- data_temp[data_temp$x1 == 1 & data_temp$x2 == 0, c(3, 4)]
  colnames(x = data_temp_01) <- colnames(x = data_temp_10) <- c("Y1", "Y2")

  cat("## Distribution of outcomes: \n")
  if (nrow(x = data_temp_10) != 0) {
    cat("# For subjects with (x1, x2) = (1, 0), \n")
    print(x = with(data = data_temp_10, expr = table(Y1, Y2)))
    cat("\n")
  }
  if (nrow(x = data_temp_01) != 0) {
    cat("# For subjects with (x1, x2) = (0, 1), \n")
    print(x = with(data = data_temp_01, expr = table(Y1, Y2)))
    cat("\n")
  }

  ##
  cat("##### Results ##### \n")
  cat("## Before considering sample selection bias, \n")
  if (logscale == TRUE) {
    EST <- sprintf(fmt = paste0("%.", digits, "f"), xx$EST)
    SE <- sprintf(fmt = paste0("%.", digits, "f"), xx$SE)
    CI_s1 <- sprintf(fmt = paste0("%.", digits, "f"), xx$CI_s[1])
    CI_s2 <- sprintf(fmt = paste0("%.", digits, "f"), xx$CI_s[2])
  } else {
    EST <- sprintf(fmt = paste0("%.", digits, "f"), exp(xx$EST))
    SE <- sprintf(fmt = paste0("%.", digits, "f"), xx$SE)
    CI_s1 <- sprintf(fmt = paste0("%.", digits, "f"), exp(xx$CI_s[1]))
    CI_s2 <- sprintf(fmt = paste0("%.", digits, "f"), exp(xx$CI_s[2]))
  }

  if (inherits(x, "asmrCLR")) {
    if (logscale == TRUE) {
      cat("Log odds ratio = ")
    } else {
      cat("Odds ratio = ")
    }
    cat(paste0(EST, ", ", xx$level * 100, "% CI = [", CI_s1, ", ", CI_s2, "] \n"))
    cat("(estimated from conditional logistic regression [using survival::clogit()]) \n")
  } else if (inherits(x, "asmrCPR")) {
    if (logscale == TRUE) {
      cat("Log risk ratio = ")
    } else {
      cat("Risk ratio = ")
    }
    cat(paste0(EST, ", ", xx$level * 100, "% CI = [", CI_s1, ", ", CI_s2, "] \n"))
    cat("(estimated from conditional Poisson regression [using gnm::gnm()]) \n")
  }

  ##
  Lambda1.digit <- sapply(X = xx$Lambda1, FUN = function(x) match(TRUE, round(x, 1:20) == x))
  Lambda2.digit <- sapply(X = xx$Lambda2, FUN = function(x) match(TRUE, round(x, 1:20) == x))

  cat("\n")
  cat("## After considering sample selection bias, \n")
  cat("# Sensitivity parameters: ")
  cat(sprintf(fmt = "\u03BB\u2081 = %.*f, \u03BB\u2082 = %.*f\n",
              Lambda1.digit, xx$Lambda1, Lambda2.digit, xx$Lambda2))
  # cat(paste0("Lambda1 = ", Lambda1, ", Lambda2 = ", Lambda2, " \n"))

  # savedig <- options(digits = digits)
  # on.exit(options(savedig))

  ## Analysis results
  if (logscale == TRUE) {
    SI1 <- sprintf(fmt = paste0("%.", digits, "f"), xx$SI[1])
    SI2 <- sprintf(fmt = paste0("%.", digits, "f"), xx$SI[2])
    CI1 <- sprintf(fmt = paste0("%.", digits, "f"), xx$CI[1])
    CI2 <- sprintf(fmt = paste0("%.", digits, "f"), xx$CI[2])
  } else {
    SI1 <- sprintf(fmt = paste0("%.", digits, "f"), exp(xx$SI[1]))
    SI2 <- sprintf(fmt = paste0("%.", digits, "f"), exp(xx$SI[2]))
    CI1 <- sprintf(fmt = paste0("%.", digits, "f"), exp(xx$CI[1]))
    CI2 <- sprintf(fmt = paste0("%.", digits, "f"), exp(xx$CI[2]))
  }
  cat("# Sensitivity interval: ")
  cat(paste0("[Ls, Us] = [", SI1, ", ", SI2,"] \n"))
  cat("# Confidence interval for population sensitivity interval: ")
  cat(paste0("[Lc, Uc] = [", CI1, ", ", CI2,"]"))

  # if (logscale == TRUE) {
  #   tmp4 <- as.data.frame(x = matrix(data = sprintf(paste0("%.", digits, "f"), xx$CI), ncol = 2))
  # } else {
  #   tmp4 <- as.data.frame(x = matrix(data = sprintf(paste0("%.", digits, "f"), exp(xx$CI)), ncol = 2))
  # }
  # rownames(tmp4) <- c("")
  # colnames(tmp4) <- c("Lower bound", "Upper bound")
  # print(tmp4)

  # invisible(x)
  invisible(NULL)
}
