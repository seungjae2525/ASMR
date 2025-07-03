## Analysis of selection-biased count matched data
asmrCPR <- function(data, Y, X, ID, Lambda1, Lambda2, level = 0.95,
                    useparallel = FALSE, n.cores = parallel::detectCores()/2) {

  ## Fit conditional Poisson regression for selection-biased binary matched data
  data_long <- make_data_long(data = data, Y = Y, X = X, ID = ID)
  CPR_s <- gnm::gnm(formula = Y ~ X, eliminate = factor(x = data_long$ID), data = data_long,
                    family = poisson, tolerance = 1e-10)
  summary_CPR_s <- summary(object = CPR_s)
  EST <- summary_CPR_s$coefficients[[1]]
  SE <- sqrt(x = summary_CPR_s$coefficients[[2]])

  ## Convert to a structured dataset
  data_old <- make_data_wide(data = data, Y = Y, X = X, ID = ID)
  data_old <- base::transform(`_data` = data_old, ss = y1 + y2)

  #### Make dataset for analysis
  ## --> Remove (x_{i1} + x_{i2} = 0), (x_{i1} + x_{i2} = 2), and (y_{i1} + y_{i2} = 0) cases
  ## Data with (x1, x2) = (1, 0) and (y1, y2) != (0, 0) cases
  data_10 <- data_old[(data_old$x1 == 1 & data_old$x2 == 0) &
                        !(data_old$y1 == 0 & data_old$y2 == 0), ]; n10 <- nrow(x = data_10)
  ## Data with (x1, x2) = (0, 1) and (y1, y2) != (0, 0) cases
  data_01 <- data_old[data_old$x1 == 0 & data_old$x2 == 1 &
                        !(data_old$y1 == 0 & data_old$y2 == 0), ]; n01 <- nrow(x = data_01)

  ##
  data_new <- rbind(data_10, data_01); rownames(x = data_new) <- NULL
  y1 <- data_new$y1
  y2 <- data_new$y2
  x1 <- data_new$x1
  x2 <- data_new$x2

  ##
  aa <- get_unique_comb(y1 = y1, y2 = y2,
                        Lambda1 = Lambda1, Lambda2 = Lambda2)
  choices <- lapply(X = aa$rik_comb, FUN = function(z) if (is.list(x = z)) z else list(z))
  rev_choices <- lapply(X = choices, FUN = function(lst) lapply(X = lst, FUN = rev))

  ##
  ss_unique <- sort(x = unique(x = data_new$ss))
  key <- match(x = data_new$ss, table = ss_unique)
  key1 <- if (n10 > 0) key[seq(from = 1, to = n10)] else integer(length = 0)
  key2 <- if (n01 > 0) key[seq(from = n10 + 1, to = n10 + n01)] else integer(length = 0)

  ##
  idx_grid <- expand.grid(lapply(X = choices, FUN = seq_along),
                          KEEP.OUT.ATTRS = FALSE)
  Rik_unique_comb_len <- nrow(x = idx_grid)

  #### For data_10 [(x1, x2) = (1, 0) and (y1, y2) != (0, 0)]
  ## For minimization
  result_10_min <- lapply(X = seq_len(length.out = Rik_unique_comb_len), FUN = function(i) {
    sel <- unlist(x = idx_grid[i, ], use.names = FALSE)
    vals <- lapply(X = key1, FUN = function(k) choices[[k]][[ sel[k] ]])
    unlist(x = vals, use.names = FALSE)
  })

  ## For maximization
  result_10_max <- lapply(X = seq_len(length.out = Rik_unique_comb_len), FUN = function(i) {
    sel <- unlist(x = idx_grid[i, ], use.names = FALSE)
    vals <- lapply(X = key1, FUN = function(k) rev_choices[[k]][[ sel[k] ]])
    unlist(x = vals, use.names = FALSE)
  })

  #### For data_01 [(x1, x2) = (0, 1) and (y1, y2) != (0, 0)]
  ## For minimization
  result_01_min <- lapply(X = seq_len(length.out = Rik_unique_comb_len), FUN = function(i) {
    sel <- unlist(x = idx_grid[i, ], use.names = FALSE)
    vals <- lapply(X = key2, FUN = function(k) rev_choices[[k]][[ sel[k] ]])
    unlist(x = vals, use.names = FALSE)
  })

  ## For maximization
  result_01_max <- lapply(X = seq_len(length.out = Rik_unique_comb_len), FUN = function(i) {
    sel <- unlist(x = idx_grid[i, ], use.names = FALSE)
    vals <- lapply(X = key2, FUN = function(k) choices[[k]][[ sel[k] ]])
    unlist(x = vals, use.names = FALSE)
  })

  #### Set all Rik combinations
  ## Minimum beta
  result_min <- lapply(X = seq_len(length.out = Rik_unique_comb_len), FUN = function(i) {
    c(result_10_min[[i]], result_01_min[[i]])
  })

  ## Maximum beta
  result_max <- lapply(X = seq_len(length.out = Rik_unique_comb_len), FUN = function(i) {
    c(result_10_max[[i]], result_01_max[[i]])
  })

  #### Optimization Setting
  if (useparallel == TRUE) {
    cl <- parallel::makeCluster(spec = n.cores, type = "PSOCK")
    parallel::clusterExport(cl = cl,
                            varlist = c("y1", "y2", "x1", "x2",
                                        "score_sens_beta",
                                        "worker_beta"),
                            envir = environment())
    on.exit(expr = parallel::stopCluster(cl), add = TRUE)
  }

  #### Find optimum beta
  if (useparallel == TRUE) {
    ## Find minimum beta
    betas_min <- parallel::parSapply(cl = cl, X = result_min, FUN = worker_beta,
                                     y1 = y1, y2 = y2, x1 = x1, x2 = x2)
    ## Find maximum beta
    betas_max <- parallel::parSapply(cl = cl, X = result_max, FUN = worker_beta,
                                     y1 = y1, y2 = y2, x1 = x1, x2 = x2)
  } else {
    ## Find minimum beta
    betas_min <- sapply(X = result_min, FUN = worker_beta,
                        y1 = y1, y2 = y2, x1 = x1, x2 = x2)
    ## Find maximum beta
    betas_max <- sapply(X = result_max, FUN = worker_beta,
                        y1 = y1, y2 = y2, x1 = x1, x2 = x2)
  }

  ## Minimum beta and corresponding Rik
  i_min <- which.min(betas_min)
  beta_min <- betas_min[i_min]
  Rik_min <- result_min[[i_min]]

  ## Maximum beta and corresponding Rik
  i_max <- which.max(betas_max)
  beta_max <- betas_max[i_max]
  Rik_max <- result_max[[i_max]]

  ## Robust variances for minimum and maximum betas
  var.lower <- roubst_var_sens_beta(beta = beta_min, rik = Rik_min,
                                    y1 = y1, y2 = y2, x1 = x1, x2 = x2)
  SE_robust_lower <- sqrt(x = var.lower)
  var.upper <- roubst_var_sens_beta(beta = beta_max, rik = Rik_max,
                                    y1 = y1, y2 = y2, x1 = x1, x2 = x2)
  SE_robust_upper <- sqrt(x = var.upper)

  ## Return results
  final.result <- list(data = data, Y = Y, X = X, ID = ID,
                       Lambda1 = Lambda1, Lambda2 = Lambda2,
                       EST = EST, SE = SE,
                       SI = c(SI_lower = beta_min, SI_upper = beta_max),
                       CI = c(CI_lower = beta_min - qnorm(p = 1 - (1 - level) / 2) * SE_robust_lower,
                              CI_upper = beta_max + qnorm(p = 1 - (1 - level) / 2) * SE_robust_upper),
                       Rik_unique_comb_len = Rik_unique_comb_len,
                       SE_robust = c(SE_robust_lower = SE_robust_lower, SE_robust_upper = SE_robust_upper),
                       Rik = list(Rik_min = Rik_min, Rik_max = Rik_max))

  class(x = final.result) <- c("asmr", "asmrCPR")
  return(final.result)
}
