## Analysis of selection-biased binary matched data
asmrCLR <- function(data, Y, X, ID, Lambda1, Lambda2, level = 0.95) {

  ## Convert to a structured dataset
  data_new <- make_data_long(data = data, Y = Y, X = X, ID = ID)

  ## Fit conditional logistic regression for selection-biased binary matched data
  CLR_s <- survival::clogit(formula = Y ~ X + strata(ID),
                            data = data_new,
                            control = coxph.control(iter.max = 1e3))
  EST <- CLR_s$coefficients[[1]]
  SE <- sqrt(x = CLR_s$var)[[1]]
  CI_s <- as.vector(confint.default(CLR_s))

  ## Sensitivity interval
  SI_lower <- EST + log(x = Lambda1)
  SI_upper <- EST + log(x = Lambda2)

  ## Confidence interval for population sensitivity interval
  CI_lower <- EST - qnorm(p = 1 - (1 - level) / 2) * SE + log(x = Lambda1)
  CI_upper <- EST + qnorm(p = 1 - (1 - level) / 2) * SE + log(x = Lambda2)

  ## Return results
  final.result <- list(data = data, Y = Y, X = X, ID = ID, level = level,
                       Lambda1 = Lambda1, Lambda2 = Lambda2,
                       EST = EST, SE = SE, CI_s = CI_s,
                       SI = c(SI_lower = SI_lower, SI_upper = SI_upper),
                       CI = c(CI_lower = CI_lower, CI_upper = CI_upper))

  class(x = final.result) <- c("asmr", "asmrCLR")
  return(final.result)
}
