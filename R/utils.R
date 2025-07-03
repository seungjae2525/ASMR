##
make_data_wide <- function(data, Y, X, ID) {
  ##
  if (length(x = Y) == 2) {
    x1 = data[, X[1]]
    x2 = data[, X[2]]
    y1 = data[, Y[1]]
    y2 = data[, Y[2]]
  } else {
    data <- data.frame(ID = data[, ID], X = data[, X], Y = data[, Y])
    data$pair <- stats::ave(data[, "ID"], data[, "ID"], FUN = seq_along)
    data_wide <- stats::reshape(data = data, direction = "wide",
                                idvar = "ID", timevar = "pair")
    x1 = data_wide$X.1
    x2 = data_wide$X.2
    y1 = data_wide$Y.1
    y2 = data_wide$Y.2
  }

  data <- data.frame(x1 = x1, x2 = x2, y1 = y1, y2 = y2)
  return(data)
}


##
make_data_long <- function(data, Y, X, ID) {
  ##
  if (length(x = Y) == 1) {
    data <- data.frame(ID = data[, ID], X = data[, X], Y = data[, Y])
  } else {
    x1 = data[, X[1]]
    x2 = data[, X[2]]
    y1 = data[, Y[1]]
    y2 = data[, Y[2]]

    X = as.vector(x = rbind(x1, x2))
    Y = as.vector(x = rbind(y1, y2))
    ID = rep(x = seq_len(length.out = length(x = x1)), each = 2)

    data <- data.frame(ID = ID, X = X, Y = Y)
  }

  return(data)
}


##
`%notin%` <- function(x, y) {
  !(x %in% y)
}

##
