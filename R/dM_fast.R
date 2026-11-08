library(cli)
library(corpcor)

#' Mahalanobis distance between two matrices
#'
#' Calculates Mahalanobis distance between two groups in a data matrix.
#'
#' @param X1 num matrix for group 1 (p columns)
#' @param X2 num matrix for group 2 (p columns)
#' @param shrink logical apply the Ledoit-Wolf shrinkage covariance estimator,
#' defaults to `shrink = F`.
#' @param return_diff whether to return the components of the Mahalanobis
#' distance calculation, defaults to `return_diff = F`
#'
#' @returns a single distance value or a list of the calculated Mahalanobis
#' distance, the difference in means, and the inverse covariance matrix

dM <- function(X1, X2, shrink = F, return_diff = F) {
  n_1 <- nrow(X1)
  n_2 <- nrow(X2)
  p  <- ncol(X1)

  if (ncol(X2) != p)
    stop("X1 and X2 must have the same number of columns.")
  if (n_1 == 0 || n_2 == 0)
    stop("Input matrices have 0 rows.")

  mu_1 <- colMeans(X1)
  mu_2 <- colMeans(X2)
  mu_diff <- mu_2 - mu_1

  if (shrink) {
    # Center each matrix independently
    inv_cov <- rbind(
      sweep(X1, 2, mu_1, FUN = "-"),
      sweep(X2, 2, mu_2, FUN = "-")
    ) %>%
      corpcor::invcov.shrink(verbose = FALSE)
  } else {
    if (p > (n_1 + n_2))
      stop("p > n: The normal covariance estimator cannot be inverted. Apply Ledoit-Wolf correction with `shrink = T`.")
    inv_cov <- rbind(X1, X2) %>%
      Rfast::pooled.cov(ina = c(rep(1L, n_1), rep(2L, n_2))) %>%
      Rfast::spdinv() }

  dM <- as.numeric(sqrt(t(mu_diff) %*% inv_cov %*% mu_diff))

  if (return_diff) {
    return(list(
      dM = dM,
      mu_diff = mu_diff,
      inv_cov = inv_cov
    ))
  } else {
    return(dM)
  }
}

#' Pairwise Mahalanobis distances
#'
#' Given a data matrix and group labels, return the Mahalanobis distances for
#' all pairwise comparisons between the groups.
#'
#' @param dat numeric data matrix
#' @param grps vector of labels for the rows of `dat`
#' @param shrink whether to use the Ledoit-Wolf shrinkage covariance estimator,
#' defaults to `shrink = F`.
#' @param format output of "vector", "matrix", "dist", or "df"
#' @param progress logical, display `cli_progress_bar()` visualization,
#' defaults to `show_progress = F`.

dM_pairwise <- function(
  dat,
  grps,
  shrink = TRUE,
  format = c("vector", "matrix", "dist", "df"),
  progress = FALSE
) {
  format <- match.arg(format)
  if (!is.matrix(dat)) dat <- as.matrix(dat)

  # 1) split into list of group matrices
  if (is.factor(grps)) grps <- droplevels(grps)
  grps_ <- unique(grps)
  K <- length(grps_)
  grp_idx <- split(1:nrow(dat), grps)[as.character(grps_)]
  dat_grps <- lapply(grp_idx, function(idx) dat[idx, , drop = FALSE])

  # 2) enumerate all pairwise combinations
  pair_mat <- combn(K, 2)
  n_pairs <- ncol(pair_mat)

  dMs <- numeric(n_pairs)
  pair_labels <- character(n_pairs)

  # 3) iterative pairwise calculation
  # optional progress bar setup
  if (progress) {
    cli::cli_progress_bar(
      name = "Pairwise Mahalanobis distances",
      total = n_pairs,
      clear = FALSE )}

  for (i in 1:n_pairs) {
    g1 <- pair_mat[1, i]
    g2 <- pair_mat[2, i]
    dMs[i] <- tryCatch({
    dM(
      X1 = dat_grps[[g1]],
      X2 = dat_grps[[g2]],
      shrink = shrink,
      return_diff = FALSE)},
    error = function(e) {
      msg <- sprintf("Pair (%s, %s) failed: %s", g1, g2, e$message)
      if (progress) {
        cli::cli_warn(msg)
      } else {
        warning(msg, call. = FALSE)
      }
      return(NA_real_)
    })

    pair_labels[i] <- sprintf(
      "%03d_%03d", # e.g., 013_012
      as.numeric(grps_[g1]),
      as.numeric(grps_[g2]))

    if (progress) cli::cli_progress_update()
  }

  names(dMs) <- pair_labels

  # 4) format
  if (format == "vector") {
    return(dMs)
  } else if (format == "dist") {
    attr(dMs, "Size") <- K
    attr(dMs, "Labels") <- grps_
    attr(dMs, "Diag") <- FALSE
    attr(dMs, "Upper") <- FALSE
    class(dMs) <- "dist"
    return(dMs)
  } else if (format == "matrix") {
    mat <- matrix(0, nrow = K, ncol = K, dimnames = list(grps_, grps_))
    mat[lower.tri(mat)] <- dMs
    mat <- mat + t(mat)
    return(mat)
  } else if (format == "df") {
    return(data.frame(
      grp1 = grps_[pair_mat[1, ]],
      grp2 = grps_[pair_mat[2, ]],
      multi_d = dMs,
      stringsAsFactors = FALSE
    ))
  }
}
