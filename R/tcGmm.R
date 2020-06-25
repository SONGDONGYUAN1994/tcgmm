#' @title Fit Two-Condition Gaussian Mixture Model
#'
#' @description
#' This function fits the Two-Condition Gaussian Mixture Model (TCGMM).
#' We assume that the latent groups are consistent between two conditions only with shifts in mean parameters.
#'
#' @param y A matrix with rows as samples (cells) and columns as features (genes)
#' @param g A vector indicating condition 1 (0) and condition 2 (1)
#' @param zInit A matrix indicating the assignment of groups with rows as samples and columns as groups
#' @param maxIter A numeric value of maximum iteration number. Default is 100.
#' @param thresh A numeric value of the converge criteria. Define as the Frobenius norm of the difference of current mean and mean in last iteration. Default is 1e-8.
#' @param verboseN A logical value. Whether to print the iteration number.
#'
#' @return A list with the fitting results
#' @param mu The mean parameter
#' @param sigma The standard deviation parameter
#' @param delta The shift of mean parameter
#' @param z The assignment of groups
#' @param model The fitted regression model of each group
#'
#'
#' @examples
#' library(extraDistr)
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
#' p_int <- c(0.4, 0.3, 0.3)
#' z_int <- rmnom(n = 100, size = 1, prob = p_int)
#' fit <- tcGmm(dat$y, dat$g, zInit = z_int)
#'
#' @importFrom extraDistr rmnom
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom stats coef glm
#' @export
#' @author Dongyuan Song



tcGmm <- function(y, g, zInit, maxIter = 100, thresh = 1e-8, verboseN = TRUE) {

  ## Set dimension
  n <- dim(y)[1]
  n.feature <- dim(y)[2]
  n.group <- dim(zInit)[2]

  for(k in seq_len(maxIter)) {
    if(k == 1) {
      z <- zInit
      gamma_curr <- zInit

      dat_all <- cbind(y, g, z)

      dat_list <- lapply(seq_len(n.group), function(x){
        dat_all[z[, x] == 1, ]
      })

      fit_list <- lapply(seq_len(n.feature), function(i) {
        lapply(seq_len(n.group), function(j){
          dat <- dat_list[[j]]

          fit <- glm(dat[, i] ~ dat[, n.feature + 1], family = "gaussian")
          fit
        })
      })
    }
    else {
      dat_all <- cbind(y, g, gamma_curr, z)
      fit_list <- lapply(seq_len(n.feature), function(i) {
        lapply(seq_len(n.group), function(j){
          dat <- dat_all

          fit <- glm(dat[, i] ~ dat[, n.feature + 1], family = "gaussian", weights = dat[, n.feature + 1 + j])
          fit
        })
      })
    }
    if(k >= 2) {mu_old <- mu_curr; z_old <- z_curr}

    ## Extract fitting values
    mu_curr <- sapply(seq_len(n.feature), function(i) {
      sapply(seq_len(n.group), function(j){
        coef(fit_list[[i]][[j]])[1]
      })
    })

    delta_curr <- sapply(seq_len(n.feature), function(i) {
      sapply(seq_len(n.group), function(j){
        coef(fit_list[[i]][[j]])[2]
      })
    })
    sigma_curr <- sapply(seq_len(n.feature), function(i) {
      sapply(seq_len(n.group), function(j){
        stats::sigma(fit_list[[i]][[j]])
      })
    })

    p_curr <- colMeans(gamma_curr)

    gamma_curr <- apply(dat_all, 1, function(x){

      y_i <- x[seq_len(n.feature)]
      g_i <- x[n.feature + 1]

      ## Calculate density
      if(g_i == 0) {
        d <- sapply(seq_len(n.group), function(i){
          dmvnorm(y_i, mean = mu_curr[i, ], sigma = diag(sigma_curr[i, ]))
        })}
      else {
        d <- sapply(seq_len(n.group), function(i){
          dmvnorm(y_i, mean = mu_curr[i, ] + delta_curr[i, ], sigma = diag(sigma_curr[i, ]))})
      }
      gamma <- sapply(seq_len(n.group), function(i){
        p_curr[i]*d[i]/(sum(p_curr*d))
      })
      gamma
    }
    )
    gamma_curr <- t(gamma_curr)
    group_row <- rep(0, n.group)
    z_curr <- t(apply(gamma_curr, 1, function(x) {
      group_row[which.max(x)] <- 1
      group_row
    }))


    if(verboseN) {
      cat(paste0("Iteration ", k, "\n"))
    }

    if(k >= 2 && norm(mu_curr - mu_old, "F") < thresh && identical(z_old, z_curr)) {
      message(paste0("Iteration ends in ", k, "\n")); break}
  }
  rownames(mu_curr) <- NULL
  rownames(delta_curr) <- NULL
  rownames(sigma_curr) <- NULL

  res <- list(mu = mu_curr, delta = delta_curr, sigma = sigma_curr, z = z_curr, model_fit = fit_list)
  return(res)
}
