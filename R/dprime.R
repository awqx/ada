#' Calculate distance d' (d prime)
#'
#' Given two response distributions, calculates d' as described in
#' Physiological Data Analysis - Initial Analysis. Takes advantage of `Rfast`
#' for multiple calculations.
#'
#' @param x1 a numeric matrix
#' @param x2 a numeric matrix
#'
#' @returns numeric d' value

dprime <- function(x1, x2) {
  d_mean <- colMeans(x1) - colMeans(x2)
  sd_comb <- sqrt(0.5 * (Rfast::colVars(x1) + Rfast::colVars(x2)))
  d_mean / sd_comb
}
