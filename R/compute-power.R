#' @title Compute the Power in a SMART
#'
#' @description Computes the power in an arbitrary SMART design with
#' the goal of identifying optimal embedded dynamic treatment regime (EDTR).
#' The power is the probability of excluding from the set of best EDTRs
#' all EDTRs which are inferior to the best EDTR by min_Delta or more.
#'
#' @param V The covariance matrix of mean EDTR estimators.
#' @param Delta The vector of effect sizes with a zero indicating the best EDTR.
#' @param min_Delta The minimum desired detectable effect size.
#' @param alpha The Type I error rate for not including the true best EDTR.
#' @param sample_size The sample size.
#'
#' @return The power to exclude from the set of best EDTR all EDTR
#' which are inferior to the best EDTR by min_Delta or more.
#'
#' @details  The true best EDTR is included
#' in the set of best with probability at least 1-alpha.
#' Multiple comparisons are adjusted for using the
#' Multiple Comparison with the Best methodology.
#'
#' @seealso
#' \code{\link{computeSampleSize}}
#' @examples
#' #Example 1
#' #V <- rbind(c(1, 0.3, 0.3, 0.3),
#' #           c(0.3, 1, 0.3, 0.3),
#' #           c(0.3, 0.3, 1, 0.3),
#' #           c(0.3, 0.3, 0.3, 1))
#'
#' #Compute power to exclude EDTRs inferior to the best by 0.3 or more
#' #The first DTR is best and the other three are inferior by 0.2, 0.6, and 0.3
#' #The best DTR is included with probability greater than or equal to 95%.
#' #computePower(V,
#' #             Delta = c(0, 0.2, 0.6, 0.3),
#' #             min_Delta = 0.3,
#' #             sample_size = 200)
#' @export
computePower <- function(V, Delta, min_Delta, alpha = 0.05, sample_size) {
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
  if (alpha <= 0 | alpha >= 0.5) {
    stop("Type I error rate alpha must be between 0 and 0.5")
  }
  if (sample_size <= 0 | sample_size != round(sample_size)) {
    stop("Sample size n must be a positive integer.")
  }

  getOnePower <- function(V, Delta, min_Delta, c_alpha, sample_size, Z) {
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

  mc_list <- simulateNormal(V, n_sim = 1000) #from internal.R
  c_alpha <- computeC(V, mc_list, alpha) #from internal.R
  avg_power <- 0
  for (i in 1:500) { #average the power over 500 sets of 1000 individuals for numerical stability
    avg_power <- avg_power + getOnePower(V, Delta, min_Delta, c_alpha, sample_size, Z = mc_list[[i]])
  }
  avg_power <- avg_power / 500
  return(round(avg_power, 2))
}
