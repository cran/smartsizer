#' @title Compute the Sample Size for a SMART.
#'
#' @description Computes the necessary sample size to enroll
#' in an arbitrary SMART design for a specified power with the goal of determining optimal
#' embedded dynamic treatment regime (EDTR). The power is the
#' probability of excluding from the set of best EDTRs all EDTRs inferior to the best
#' by min_Delta or more.
#'
#'
#' @param V The covariance matrix of mean EDTR estimators.
#' @param Delta The vector of effect sizes with the first zero indicating the best EDTR.
#' @param min_Delta The minimum desired detectable effect size.
#' @param alpha The Type I error rate for not including the true best EDTR.
#' @param desired_power The desired power.
#'
#' @return The minimum sample size in order to achieve a power of desired_power to exclude EDTRs
#' from the set of best which are inferior to the optimal EDTR by min_Delta or more.
#'
#' @details The true best EDTR is included in the set of best with probability
#' at least 1-alpha. Multiple comparisons are adjusted for using the
#' Multiple Comparison with the Best methodology.
#'
#' @seealso
#' \code{\link{computePower}}
#'
#' @examples
#' #Example 1
#' #V <- rbind(c(1, 0.3, 0.3, 0.3),
#' #           c(0.3, 1, 0.3, 0.3),
#' #           c(0.3, 0.3, 1, 0.3),
#' #           c(0.3, 0.3, 0.3, 1))
#'
#' #Compute sample size to achieve power of 80% to exclude EDTRs inferior
#' #to the best by 0.3 or more. The first DTR is best and the other
#' #three are inferior by 0.2, 0.6, and 0.3
#' #The best EDTR is included with probability greater than or equal to 95%.
#' #computeSampleSize(V,
#' #                  Delta = c(0, 0.2, 0.6, 0.3),
#' #                  min_Delta = 0.3,
#' #                  alpha = 0.05,
#' #                  desired_power = 0.8)
#' @export
computeSampleSize <- function(V, Delta, min_Delta, alpha = 0.05, desired_power) {
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
  if (desired_power < 0 | desired_power > 1){
    stop("desired_power must be between 0 and 1.")
  }
  getOneSampleSize <- function(V, Delta, min_Delta, c_alpha, desired_power, Z){
    SampleSize <- 0
    N <- min(which(Delta == 0))
    temp <- matrix(0, nrow(Z), nrow(V))
    for (i in which(Delta >= min_Delta)) {
      temp[,i] <- ((Z[,i] - Z[,N]) + c_alpha[i] * sqrt(V[i,i] + V[N,N] - 2 * V[i,N])) / Delta[i]
    }
    temp <- temp[,-N]
    temp <- apply(temp, 1, max)
    SampleSize <- (stats::quantile(temp, desired_power, names = FALSE))^2
    return(SampleSize)
  }
  mc_list <- simulateNormal(V, n_sim = 1000) #from internal.R
  c_alpha <- computeC(V, mc_list, alpha) #from internal.R
  sample_size <- 0
  for (i in 1:500) {  # Averages estimates of necessary sample size over 500 simulated normal data sets
    sample_size <- sample_size + getOneSampleSize(V, Delta, min_Delta, c_alpha, desired_power, Z = mc_list[[i]])
  }
  sample_size <- sample_size / 500
  return(ceiling(sample_size))
}
