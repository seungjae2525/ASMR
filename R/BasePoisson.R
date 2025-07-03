### Score function for sample selection bias in conditional Poisson regression
score_sens_beta <- function(par, rik, y1, y2, x1, x2) {
  beta <- par
  n <- length(x = y1)
  score <- 0

  ## Pre-calculation for bik
  max_sy <- max(y1 + y2)
  bik_list <- vector(mode = "list", length = max_sy + 1)
  for (jj in 0:max_sy) {
    ks <- 0:jj
    bik_list[[jj + 1]] <- exp(x = lgamma(x = ks + 1) + lgamma(x = jj - ks + 1))
  }

  idx <- 1
  for (i in seq_len(length.out = n)) {
    sy <- y1[i] + y2[i]
    ks <- 0:sy
    xi1 <- x1[i]
    xi2 <- x2[i]

    len_i <- sy + 1
    rik_i <- rik[idx:(idx + len_i - 1)]
    idx <- idx + len_i

    aiks <- ks * xi1 + (sy - ks) * xi2
    biks <- bik_list[[sy + 1]]

    w <- rik_i * exp(x = beta * aiks) / biks

    num <- sum(aiks * w)
    denom <- sum(w)

    score <- score + ((y1[i] * xi1 + y2[i] * xi2) - num / denom)
  }

  return(score)
}

###
get_unique_comb <- function(y1, y2, Lambda1, Lambda2) {
  ##
  ss <- y1 + y2
  ss_unique <- sort(x = unique(x = ss))
  rik_each_len <- ss_unique + 1

  ##
  rik_inx <- 1
  rik_comb_length <- ifelse(test = rik_each_len == 1, yes = 1, no = rik_each_len - 1)
  rik_comb <- vector(mode = "list", length = length(x = ss_unique))

  for (i in 1:length(x = ss_unique)) {
    ss_temp <- ss_unique[i]

    seq_temp <- rik_inx:(rik_inx + rik_comb_length[i] - 1)
    if (ss_temp == 0) {
      rik_comb[[i]] <- rep(x = 1, times = ss_temp + 1)

    } else if (ss_temp == 1) {
      rik_comb[[i]] <- c(Lambda1, Lambda2)

    } else {
      Lambda_mat_Lambda <- lapply(X = 0:(length(x = seq_temp) - 1), FUN = function(ii) {
        c(rep(x = Lambda1, times = ii + 1),
          rep(x = Lambda2, times = length(x = seq_temp) - ii))
      })

      rik_comb[[i]] <- Lambda_mat_Lambda

    }

    rik_inx <- rik_inx + rik_comb_length[i]
  }

  return(list(ss_unique = ss_unique, rik_comb = rik_comb))
}

###
worker_beta <- function(rik_vec, y1, y2, x1, x2) {
  tryCatch(
    stats::uniroot(
      function(b) score_sens_beta(par = b,
                                  rik = rik_vec,
                                  y1 = y1, y2 = y2,
                                  x1 = x1, x2 = x2),
      lower = -10, upper = 10, tol = 1e-8, maxiter = 1000L
    )$root,
    error = function(e) NA_real_
  )
}

### Robust variance estimator for sample selection bias in conditional Poisson regression
roubst_var_sens_beta <- function(beta, rik, y1, y2, x1, x2) {
  n <- length(x = y1)
  score_square <- 0
  hessian <- 0

  ## Pre-calculation for bik
  max_sy <- max(y1 + y2)
  bik_list <- vector(mode = "list", length = max_sy + 1)
  for (jj in 0:max_sy) {
    ks <- 0:jj
    bik_list[[jj + 1]] <- exp(x = lgamma(x = ks + 1) + lgamma(x = jj - ks + 1))
  }

  idx <- 1
  for (i in seq_len(length.out = n)) {
    sy <- y1[i] + y2[i]
    ks <- 0:sy
    xi1 <- x1[i]
    xi2 <- x2[i]

    len_i <- sy + 1
    rik_i <- rik[idx:(idx + len_i - 1)]
    idx <- idx + len_i

    aiks <- ks * xi1 + (sy - ks) * xi2
    biks <- bik_list[[sy + 1]]

    w <- rik_i * exp(x = beta * aiks) / biks

    num <- sum(aiks * w)
    denom <- sum(w)

    score_square <- score_square + ((y1[i] * xi1 + y2[i] * xi2) - num / denom)^{2}


    num2_1 <- sum(aiks^{2} * w)
    denom2_1 <- sum(w)
    num2_2 <- sum(aiks * w)

    hessian <- hessian + ((num2_1 / denom2_1) - (num2_2 / denom2_1)^{2})
  }

  variance <- score_square / (hessian)^{2}

  return(variance)
}

