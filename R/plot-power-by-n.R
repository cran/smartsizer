#' @title Plot Power by Sample Size
#'
#' @description Plots the power over a grid of sample sizes.
#'
#' @param V The covariance matrix of mean EDTR estimators.
#' @param Delta The vector of effect sizes with a zero indicating the best EDTR.
#' @param min_Delta The minimum desired detectable effect size.
#' @param alpha The Type I error rate for not including the true best EDTR.
#' @param sample_size_grid A vector of sample sizes.
#' @param color The color of the graph.
#'
#' @details It employs common random variables
#' to reduce the variance. See \code{\link{computePower}} for more details.
#'
#'@export
plotPowerByN <- function(V, Delta, min_Delta, alpha = 0.05, sample_size_grid, color = "black") {
  pow <- computePowerBySampleSize(V, Delta, min_Delta, alpha, sample_size_grid)
  graphics::plot(x = sample_size_grid, y = pow,
       type = "o",
       col = color,
       xlab = "Sample size (n)",
       ylab = "Power",
       main = "Power vs. Sample Size")
}
