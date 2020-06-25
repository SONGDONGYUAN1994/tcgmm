#' @title Generate Simulated Data of Two-Condition Gaussian
#'
#' @description
#' This function generates multivariate Gaussian distributed data under case/control conditions.
#' The basic assumption is that for each type, the mean has a shift while the variance keeps the same cross conditions.
#'
#' @param n The sample size of simulated data (cells)
#' @param n.feature Number of features (genes)
#' @param n.group Number of underlying groups (cell types)
#' @param type.prop A vector of same length as \code{n.group}. The proportion of each group.
#' @param con.prop Proportion of condition 1. Default is 0.5.
#' @param mu.mat A matrix which is \code{n.feature} by \code{n.group}. Define the mean parameter of multivariate Gaussian in condition 1.
#' @param sigma.mat A matrix which is \code{n.feature} by \code{n.group}. Define the standard deviation parameter of multivariate Gaussian.
#' @param delta.mat A matrix which is \code{n.feature} by \code{n.group}. Define the shift of mean of multivariate Gaussian.
#'
#' @return
#' A list which contains the data, the group type and the condition.
#' ## y The feature of each data point
#' ## z The group type of each data point
#' ## g The condtion (condition1 / condition 2) of each data point
#'
#' @examples
#' mu1 <- c(5, 7, 9)
#' theta1 <- c(1, 2, 0)
#' sigma1 <- c(1, 2, 3)
#' mu2 <- c(10, 15, 4)
#' theta2 <- c(1, 2, 6)
#' sigma2 <- c(0.4, 0.2, 0.4)
#' mu.mat <- cbind(mu1, mu2)
#' delta.mat <- cbind(theta1, theta2)
#' sigma.mat <- cbind(sigma1, sigma2)
#' dat <- simGen(n = 100, n.feature = 2, n.group = 3, type.prop = c(0.2, 0.3, 0.5),
#' mu.mat = mu.mat, sigma.mat = sigma.mat, delta.mat = delta.mat)
#' @importFrom extraDistr rmnom
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @export
#' @author Dongyuan Song

simGen <- function(n, n.feature, n.group, type.prop, con.prop = 0.5, mu.mat, sigma.mat, delta.mat) {

  ## Check dimension
  stopifnot(length(type.prop) == n.group)
  stopifnot(dim(mu.mat)[1] == n.group)
  stopifnot(dim(mu.mat)[2] == n.feature)
  stopifnot(dim(mu.mat) == dim(sigma.mat))
  stopifnot(dim(mu.mat) == dim(delta.mat))

  ## latend variable
  z <- rmnom(n = n, size = 1, prob = type.prop)

  ## true group variable
  g <- c(rep(0, round(n*con.prop)), rep(1, n - round(n*con.prop)))

  ## multivariate Gaussian
  mu.mat <- data.matrix(mu.mat)
  sigma.mat <- data.matrix(sigma.mat)

  s_control <- (lapply(seq_len(n.group), function(x){
    mu <- mu.mat[x, ]
    sigma <- diag(sigma.mat[x, ])
    rmvnorm(n = n, mean = mu, sigma = sigma)
  }))

  s_case <- (lapply(seq_len(n.group), function(x){
    mu <- mu.mat[x, ] + delta.mat[x, ]
    sigma <- diag(sigma.mat[x, ])
    rmvnorm(n = n, mean = mu, sigma = sigma)
  }))

  y <- lapply(seq_len(n.group), function(x){
    z[, x]*s_control[[x]]*(1-g) + z[, x]*s_case[[x]]*g
  })

  y <- Reduce('+', y)

  res <- list(y = y, z = z, g = g)

  return(res)
}
