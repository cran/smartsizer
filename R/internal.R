#' @title Internal Functions
#'
#' @description This script contains functions which are implicitly
#' used by other functions in SMARTSizeR. The enduser does not need
#' to work with these functions.
#'
#'
#' @noRd
#'
simulateNormal <- function(V, n_sim = 1000){
  #Simulates 500 sets of n_sim multivariate normal random variables.
  #
  # Args:
  #  V: the covariance matrix
  #  n_sim: number of normal draws for each of the 500 sets.
  #
  # Returns:
  #  A list of 500 sets, each with n_sim normal random draws.
  mc_list <- vector("list", length = 500)
  for (i in 1:500) {
    mc_list[[i]] <- MASS::mvrnorm(n = n_sim, mu = rep.int(0,nrow(V)), Sigma = V)
  }
  return(mc_list)
}


ComputeSVD <- function(V, mc_list){
  # Returns the transformed normal random variables N(0,I) to N(0,V)
  # Uses singular value decomposition
  lapply(mc_list, FUN = function(x) t(eigen(V)$vectors %*%
                                        (diag(sqrt(eigen(V)$values))
                                         %*% t(x))))
}

getOneC <- function(V, Z, alpha = 0.05){
  # Computes a single vector of the ci's given the covariance, a single
  # vector of standard normal random draws, and type I error rate alpha.
  #
  # Args:
  #  V: covariance matrix of \sqrt{n}(\hat{\theta}).
  #  Z: vector of n_sim standard normal random draws.
  #  alpha: type I error rate.
  #
  # Returns:
  #  A single estimate of a nrow(V) length vector the ci's.
  c_alpha <- rep.int(0,nrow(V))

  for (i in 1:nrow(V)) {
    temp <- matrix(0,nrow(Z),nrow(V))
    for (j in setdiff(1:nrow(V),c(i))) {
      temp[,j] <- (Z[,j]-Z[,i])/sqrt(V[i,i]+V[j,j]-2*V[i,j])
    }
    temp <- temp[,-i]
    #temp<-apply(temp,1,sort)
    temp <- apply(temp,1,max)
    c_alpha[i] <- stats::quantile(temp,1-alpha)
  }
  c_alpha
}

computeC <- function(V, mc_list_SVD, alpha = 0.05) {
  # Computes an estimate for ci (quantiles) averaged 500 times.
  #
  # Args:
  #  V: Covariance matrix of \sqrt{n}(\hat{\theta}).
  #  mc_list: A list of 500 sets of n_sim/n_sim standard normal draws.
  #  alpha: Type I Error Probability
  #
  # Depends:
  #  getOneC(V, Z, alpha).
  #
  # Returns:
  #  An nrow(V) length vector of c'i (quantiles) averaged 500 times.
  c_alpha <- rep.int(0,nrow(V))
  for (i in 1:500) {
    c_alpha <- c_alpha + getOneC(V, Z = mc_list_SVD[[i]], alpha)
  }
  c_alpha <- c_alpha/500
  c_alpha
}


getOnePowerForGrid<- function(V, Delta, min_Delta, c_alpha, sample_size, Z){
  power <- 0
  N <- min(which(Delta == 0))
  temp <- matrix(0, nrow(Z), nrow(V))
  for (i in which(Delta >= min_Delta)) {
    temp[,i] <- ((Z[,i] - Z[,N]) + c_alpha[i] * sqrt(V[i,i] + V[N,N] - 2 * V[i,N])) / Delta[i]
  }
  temp <- temp[,-N]
  temp <- apply(temp, 1, max)
  power <- mean(temp < sqrt(sample_size))

  return(power)
}

computePowerForGrid <- function(V, Delta, min_Delta, c_alpha, sample_size, mc_list) {
  ##Averages the power over 500 normal random variable data sets
  ##See documentation for getOnePowerForGrid for more information
  ## Z is a list of N(0,V) simulated random variables matrices
  if (min(Delta) > 0) {
    stop("Delta vector not valid. Difference between best DTR and itself should be 0.")
  }
  if (min(Delta) < 0) {
    stop("Delta vector not valid. Standardized differences in Delta should be non-negative.")
  }
  if (length(Delta) != nrow(V)) {
    stop("The number of DTRs in the Delta vector does not match that of the covariance matrix V.")
  }
  if (det(V) <= 0) {
    stop("Covariance matrix V is not positive definite.")
  }
  if (min_Delta <= 0) {
    stop("Minimum desired detectable effect size min_Delta must be positive.")
  }
  if (sample_size <= 0 | sample_size != round(sample_size)) {
    stop("Sample size n must be a positive integer.")
  }

  avg_power <- 0
  for (i in 1:500) {
    avg_power <- avg_power + getOnePowerForGrid(V, Delta, min_Delta, c_alpha, sample_size, Z = mc_list[[i]])
  }
  avg_power <- avg_power / 500
  return(round(avg_power, 2))
}
