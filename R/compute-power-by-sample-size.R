#' @title Compute the Power Over a Grid of Sample Size Values
#'
#' @description Computes the power over
#' a grid of sample size values.
#'
#' @param V The covariance matrix of mean EDTR estimators.
#' @param Delta The vector of effect sizes with a zero indicating the best EDTR.
#' @param min_Delta The minimum desired detectable effect size.
#' @param alpha The Type I error rate for not including the true best EDTR.
#' @param sample_size_grid The vector of sample sizes
#'
#' @return A vector of power for each sample size in the given grid.
#'
#' @details It employs common random variables
#' to reduce the variance. See \code{\link{computePower}} for more details.
#'
#' @seealso \code{\link{computePower}}
#'
#' @examples
#' #V <- rbind(c(1, 0.3, 0.3, 0.3),
#' #           c(0.3, 1, 0.3, 0.3),
#' #           c(0.3, 0.3, 1, 0.3),
#' #           c(0.3, 0.3, 0.3, 1))
#' #computePowerBySampleSize(V,
#' #                         Delta = c(0, 0.2, 0.6, 0.3),
#' #                         min_Delta = 0.3,
#' #                         sample_size_grid = seq(50,300, 50))
#' @export
computePowerBySampleSize <- function(V, Delta, min_Delta, alpha = 0.05, sample_size_grid) {
  mc_list <- simulateNormal(V, n_sim = 1000)
  mc_list_SVD <- ComputeSVD(V, mc_list)
  c_alpha <- computeC(V, mc_list_SVD, alpha) #indep. of n
  sapply(sample_size_grid, function(x) computePowerForGrid(V, Delta, min_Delta, c_alpha, x, mc_list = mc_list_SVD))
}

