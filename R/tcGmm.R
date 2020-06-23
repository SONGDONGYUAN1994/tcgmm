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
#' @param verboseN A llogical value. Whether to print the iteration number.
#'
#' @return
#' A list which contains the fitting result.
#' ## mu The mean parameter
#' ## sigma The standard deviation parameter
#' ## delta The shift of mean parameter
#' ## z The assignment of groups
#' ## model The fitted regression model of each group
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

  for(k in 1:maxIter) {
    if(k == 1) {
      z <- zInit
      gamma_curr <- zInit

      dat_all <- cbind(y, g, z)

      dat_list <- lapply(1:n.group, function(x){
        dat_all[z[, x] == 1, ]
      })

      fit_list <- lapply(1:n.feature, function(i) {
        lapply(1:n.group, function(j){
          dat <- dat_list[[j]]

          fit <- glm(dat[, i] ~ dat[, n.feature + 1], family = "gaussian")
          fit
        })
      })
    }
    else {
      dat_all <- cbind(y, g, gamma_curr, z)
      fit_list <- lapply(1:n.feature, function(i) {
        lapply(1:n.group, function(j){
          dat <- dat_all

          fit <- glm(dat[, i] ~ dat[, n.feature + 1], family = "gaussian", weights = dat[, n.feature + 1 + j])
          fit
        })
      })
    }
    if(k >= 2) mu_old <- mu_curr

    ## Extract fitting values
    mu_curr <- sapply(1:n.feature, function(i) {
      sapply(1:n.group, function(j){
        coef(fit_list[[i]][[j]])[1]
      })
    })

    delta_curr <- sapply(1:n.feature, function(i) {
      sapply(1:n.group, function(j){
        coef(fit_list[[i]][[j]])[2]
      })
    })
    sigma_curr <- sapply(1:n.feature, function(i) {
      sapply(1:n.group, function(j){
        stats::sigma(fit_list[[i]][[j]])
      })
    })

    p_curr <- colMeans(gamma_curr)

    gamma_curr <- apply(dat_all, 1, function(x){

      y_i <- x[1:n.feature]
      g_i <- x[n.feature + 1]

      ## Calculate density
      if(g_i == 0) {
        d <- sapply(1:n.group, function(i){
          dmvnorm(y_i, mean = mu_curr[i, ], sigma = diag(sigma_curr[i, ]))
        })}
      else {
        d <- sapply(1:n.group, function(i){
          dmvnorm(y_i, mean = mu_curr[i, ] + delta_curr[i, ], sigma = diag(sigma_curr[i, ]))})
      }
      gamma <- sapply(1:n.group, function(i){
        p_curr[i]*d[i]/(sum(p_curr*d))
      })
      gamma
    }
    )
    gamma_curr <- t(gamma_curr)

    if(verboseN) {
      cat(paste0("Iteration ", k, "\n"))
    }

    if(k >= 2 &&  norm(mu_curr - mu_old, "F") < thresh) {message(paste0("Iteration ends in ", k, "\n")); break}
  }
  rownames(mu_curr) <- NULL
  rownames(delta_curr) <- NULL
  rownames(sigma_curr) <- NULL

  group_row <- rep(0, n.group)
  z_curr <- t(apply(gamma_curr, 1, function(x) {
    group_row[which.max(x)] <- 1
    group_row
  }))

  res <- list(mu = mu_curr, delta = delta_curr, sigma = sigma_curr, z = z_curr, model_fit = fit_list)
  return(res)
}
