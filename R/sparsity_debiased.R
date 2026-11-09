#' Unbiased per-stimulus powers of the true mean
#'
#' For a stimulus with repetitions Y_1, ..., Y_K that are independent given
#' the true mean mu (any noise distribution with finite moments), the mean of
#' products over distinct repetitions is unbiased for mu^p, e.g.,
#' E[Y_r Y_s] = mu^2 for r != s. These are the elementary symmetric
#' polynomials of the repetitions divided by choose(K, p), computed from power
#' sums (Newton's identities).
#'
#' @param Y numeric matrix of stimuli x repetitions (one site)
#' @param p integer power, 1 to 4 (requires K >= p)
#'
#' @returns vector with an unbiased estimate of mu_i^p for each stimulus

rep_power_mean <- function(Y, p) {
  K <- ncol(Y)
  if (p > K) stop("Need at least ", p, " repetitions for power ", p, ".")
  s1 <- rowSums(Y)
  if (p == 1) return(s1 / K)
  s2 <- rowSums(Y^2)
  e <- switch(
    p - 1,
    (s1^2 - s2) / 2,
    {
      s3 <- rowSums(Y^3)
      (s1^3 - 3 * s1 * s2 + 2 * s3) / 6
    },
    {
      s3 <- rowSums(Y^3)
      s4 <- rowSums(Y^4)
      (s1^4 - 6 * s1^2 * s2 + 3 * s2^2 + 8 * s1 * s3 - 6 * s4) / 24
    })
  e / choose(K, p)
}

#' Moments of the tuning curve across stimuli
#'
#' Raw moments m_p = mean over stimuli of mu_i^p, for p = 1 to 4, plus the
#' squared mean m1^2. Plug-in estimates use the trial-averaged responses;
#' debiased estimates use `rep_power_mean()`, which removes the trial-noise
#' bias (of order 1 / K). For m1^2, cross-stimulus products of means are
#' independent, so only the same-stimulus terms need replacing.
#'
#' @param Y numeric matrix of stimuli x repetitions (one site); a vector is
#' treated as noise-free true means
#' @param debias whether to remove the trial-noise bias
#'
#' @returns named vector `m1`, `m1sq`, `m2`, `m3`, `m4`

tuning_moments <- function(Y, debias = T) {
  if (is.null(dim(Y))) {
    mu <- Y
    return(c(m1 = mean(mu), m1sq = mean(mu)^2, m2 = mean(mu^2),
             m3 = mean(mu^3), m4 = mean(mu^4)))
  }
  N <- nrow(Y)
  xbar <- rowMeans(Y)
  if (!debias) {
    return(c(m1 = mean(xbar), m1sq = mean(xbar)^2, m2 = mean(xbar^2),
             m3 = mean(xbar^3), m4 = mean(xbar^4)))
  }
  u2 <- rep_power_mean(Y, 2)
  c(m1 = mean(xbar),
    m1sq = (sum(xbar)^2 - sum(xbar^2) + sum(u2)) / N^2,
    m2 = mean(u2),
    m3 = if (ncol(Y) >= 3) mean(rep_power_mean(Y, 3)) else NA_real_,
    m4 = if (ncol(Y) >= 4) mean(rep_power_mean(Y, 4)) else NA_real_)
}

#' Sparsity measures from tuning moments
#'
#' On non-negative responses, the activity ratio a = m1^2 / m2 (Treves and
#' Rolls), Treves-Rolls sparseness (1 - a) / (1 - 1/N) (Rolls and Tovee, 1995;
#' Vinje and Gallant, 2000), and Hoyer sparsity (sqrt(N) - L1/L2) /
#' (sqrt(N) - 1) = sqrt(N) (1 - sqrt(a)) / (sqrt(N) - 1) are all functions of
#' a. Excess kurtosis is of the distribution of responses across stimuli
#' (lifetime kurtosis). Undefined values (non-positive second moments) are
#' returned as NA.
#'
#' @param m named vector from `tuning_moments()`
#' @param N number of stimuli
#'
#' @returns named vector `activity`, `tr`, `hoyer`, `kurtosis`

moment_sparsity <- function(m, N) {
  a <- if (m[["m2"]] > 0) m[["m1sq"]] / m[["m2"]] else NA_real_
  # central moments from raw moments
  c2 <- m[["m2"]] - m[["m1sq"]]
  c4 <- m[["m4"]] - 4 * m[["m1"]] * m[["m3"]] + 6 * m[["m1sq"]] * m[["m2"]] -
    3 * m[["m1sq"]]^2
  c(activity = a,
    tr = (1 - a) / (1 - 1 / N),
    hoyer = if (!is.na(a) && a >= 0) sqrt(N) * (1 - sqrt(a)) / (sqrt(N) - 1) else NA_real_,
    kurtosis = if (!is.na(c2) && c2 > 0) c4 / c2^2 - 3 else NA_real_)
}

#' Gini index of non-negative responses
#'
#' @param x non-negative numeric vector
#'
#' @returns Gini index in [0, 1 - 1/N]

gini <- function(x) {
  N <- length(x)
  if (sum(x) <= 0) return(NA_real_)
  x <- sort(x)
  sum((2 * seq_len(N) - N - 1) * x) / (N * sum(x))
}

#' Plug-in and debiased sparsity of one site
#'
#' @param Y numeric matrix of stimuli x repetitions (one site)
#'
#' @returns named vector: each moment-based measure as `<measure>_plugin` and
#' `<measure>_debiased`, plus `gini_plugin`

sparsity_estimates <- function(Y) {
  N <- nrow(Y)
  plugin <- moment_sparsity(tuning_moments(Y, debias = F), N)
  debiased <- moment_sparsity(tuning_moments(Y, debias = T), N)
  c(setNames(plugin, paste0(names(plugin), "_plugin")),
    setNames(debiased, paste0(names(debiased), "_debiased")),
    gini_plugin = gini(rowMeans(Y)))
}

#' Sparsity of the true tuning curve
#'
#' @param mu numeric vector of true mean responses
#'
#' @returns named vector of the same measures as `sparsity_estimates()`
#' (suffix `_truth`)

sparsity_truth <- function(mu) {
  s <- c(moment_sparsity(tuning_moments(mu), length(mu)), gini = gini(mu))
  setNames(s, paste0(names(s), "_truth"))
}
