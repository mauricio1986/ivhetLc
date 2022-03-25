sech <- function(x) 1 / cosh(x)

## suml function from mlogit (Croissant, 2013)
suml <- function(x){
  n <- length(x)
  if (!is.null(dim(x[[1]]))) {
    d <- dim(x[[1]])
    s <- matrix(0,d[1], d[2])
    for (i in 1:n) {
      x[[i]][is.na(x[[i]])] <- 0
      s <- s + x[[i]]
    }
  }
  else{
    s <- rep(0,length(x[[n]]))
    for (i in 1:n) {
      x[[i]][is.na(x[[i]])] <- 0
      s <- s + x[[i]]
    }
  }
  s
}


repRows <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}

repCols <- function(x, n){
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}


make.shares <- function(theta, obj, ...){
  theta    <- theta
  names.pi <- colnames(model.matrix(obj$formula, data = obj$mf, rhs = 3))
  names.H  <- paste("pi", names.pi, sep = ".")
  H.bar    <- colMeans(model.matrix(obj$formula, data = obj$mf, rhs = 3))
  Q        <- obj$Q
  expsq    <- c(exp(0))
  for (q in 2:Q){
    hd <- c(0)
    for (l in 1:length(H.bar)){
      name.delta <- paste("class", q, names.H[l], sep = ".")
      temp       <- H.bar[l] * theta[name.delta]
      hd         <- hd + temp
    }
    expsq <- c(expsq, exp(hd))
  }
  shares <- expsq / sum(expsq)
  return(shares)
}

make.rhos<- function(theta, obj, ...){
  Q <- obj$Q
  rhos <- c()
  for (q in 1:Q){
    name.rho <- paste("class", q, "athrho", sep = ".")
    rhos <- c(rhos, tanh(theta[name.rho]))
  }
  return(rhos)
}

make.sigmas <- function(theta, obj, ...){
  Q      <- obj$Q
  sigmas <- c()
  for (q in 1:Q){
    name.sigma1 <- paste("class", q, "eq.1.lnsigma", sep = ".")
    name.sigma2 <- paste("class", q, "eq.2.lnsigma", sep = ".")
    sigmas      <- c(sigmas, exp(theta[name.sigma1]), exp(theta[name.sigma2]))
  }
  return(sigmas)
}


make.wlate <- function(theta, obj, ...){
  shares <- make.shares(theta, obj)
  end    <- obj$end.var
  beta.names  <- grep(end, names(theta), value = TRUE)
  wlate <- sum(shares * theta[beta.names])
  return(wlate)
}

#' @export
fit.measures <- function(obj, ...){
  if (!inherits(obj, "lciv")) stop("The model was not estimated using ivhetLc")
  y1 <- obj$y1
  y2 <- obj$y2
  X  <- obj$X
  Z  <- obj$Z
  Hl <- obj$Hl
  Q  <- obj$Q
  res <- ml_lcivd(theta = coef(obj), y1 = y1, y2 = y2, X = X, Z = Z, W = Hl, Q = Q, gradient = FALSE)
  post_pi <- attr(res, 'post_pi')
  pi      <- attr(res, 'pi')
  colnames(post_pi) <- colnames(pi) <- paste("class", 1:Q, sep = "=")

  # Shares
  n <- nrow(pi)
  c.shares <- as.vector(table(apply(post_pi, 1, function(x) which(x == max(x)))) / n)
  u.shares <- as.vector(colMeans(pi))
  # Entropy measure
  E   <- 1 - sum(rowSums(- post_pi * log(post_pi))) / (n * log(Q))
  out <- list(cond.pi  = post_pi,
              unco.pi  = pi,
              entropy  = E,
              c.shares = c.shares,
              u.shares = u.shares)
  return(out)
}
