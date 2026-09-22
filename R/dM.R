library(cli)
library(corpcor)

#' Mahalanobis distance between two groups
#'
#' Calculates Mahalanobis distance between two groups in a data matrix.
#'
#' @param dat numeric data matrix
#' @param grps a vector of binary group labels for the rows of `dat`
#' @param shrink whether to use the Ledoit-Wolf shrinkage covariance estimator,
#' defaults to `shrink = F`.
#' @param return_diff whether to return the components of the Mahalanobis
#' distance calculation, defaults to `return_diff = T`
#'
#' @returns a single distance value or a list of the calculated Mahalanobis
#' distance, the difference in means, and the inverse covariance matrix

dM <- function(dat, grps, shrink = F, return_diff = T) {
  grps <- as.numeric(as.factor(grps))
  mu_1 <- colMeans(dat[grps == 1, ])
  mu_2 <- colMeans(dat[grps == 2, ])
  mu_diff <- mu_2 - mu_1

  if (shrink) {
    inv_cov <- rbind(
      sweep(dat[grps == 1, ], 2, mu_1, FUN = "-"),
      sweep(dat[grps == 2, ], 2, mu_2, FUN = "-")) %>%
      corpcor::invcov.shrink(verbose = F)
  } else {
    inv_cov <- Rfast::spdinv(pooled.cov(dat, ina = grps)) }
  multi_d <- sqrt(t(mu_diff) %*% inv_cov %*% (mu_diff))

  if (return_diff) {
    return(list(
      multi_d = multi_d,
      mu_diff = mu_diff,
      inv_cov = inv_cov))
  } else {
    return(multi_d)
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
#' @param return_diff whether to return the components of the Mahalanobis
#' distance calculation, defaults to `return_diff = T`
#' @param show_progress whether to display `cli_progress_bar()` visualization,
#' defaults to `show_progress = F`.

dM_pairwise <- function(
  dat,
  grps,
  shrink = T,
  return_diff = T,
  show_progress = F
) {
  grp_pairs <- combn(unique(grps), 2, simplify = F)
  if (show_progress) {
    multi_d_list <- lapply(
      cli_progress_along(grp_pairs),
      function(i) {
        grp_ids <- grp_pairs[[i]]
        idx <- grps %in% grp_ids
        multi_d(
          dat = dat[idx, ],
          grps = grps[idx],
          shrink = shrink,
          return_diff = return_diff)})
    } else {
    multi_d_list <- lapply(
      grp_pairs,
      function(grp_ids) {
        idx <- grps %in% grp_ids
        dM(
          dat = dat[idx, ],
          grps = grps[idx],
          shrink = shrink,
          return_diff = return_diff)}) }

  if (return_diff) {
    multi_d_df <- data.frame(
      multi_d = sapply(multi_d_list, function(x) x[["multi_d"]]),
      grp1 = sapply(grp_pairs, "[[", 1),
      grp2 = sapply(grp_pairs, "[[", 2))
    mu_diff <- lapply(multi_d_list, function(x) x[["mu_diff"]]) %>%
      do.call(rbind, .)
    return(list(
      multi_d = multi_d_df,
      mu_diff = mu_diff
    ))
  } else {
    multi_d_df <- data.frame(
      multi_d = unlist(multi_d_list),
      grp1 = sapply(grp_pairs, "[[", 1),
      grp2 = sapply(grp_pairs, "[[", 2))
    return(multi_d_df)}
}
