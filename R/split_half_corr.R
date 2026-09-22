#' Diagonal of the correlation matrix
#'
#' A helper function that only provides the diagonal of the correlation matrix
#' between two given matrices.

diag_corr <- function(x1, x2) {
  colSums(scale(x1) * scale(x2)) / (nrow(x1) - 1)
}

#' Split-half correlation
#'
#' Given a list of two splits and labels for average groupings, calculates the
#' Pearson correlation between the two splits.
#'
#' @param splits list of two data matrices of matching dimension
#' @param avg_grp vector of group labels to average the observations
#' @param onset integer for the number of indices after stimulus presentation
#' for averaging to occur
#' @param l_window integer number of bins for averaging over
#' @param stim_idx integer bin number of stimulus presentation
#' @param incl_idx optional vector to indicate whether to subset the splits
#' before doing calculations

splits_corr <- function(splits, avg_grp, onset, l_window = 15, stim_idx = 50, incl_idx = NULL) {
  if (!is.null(incl_idx)) {
    splits <- lapply(splits, function(x) x[incl_idx, , ])
    avg_grp <- avg_grp[incl_idx]
  }

  res <- lapply(splits, function(split) {
    apply(split[, onset:(onset + l_window - 1) + stim_idx - 1, ], c(1, 3), mean) %>%
      data.frame() %>%
      mutate(avg_grp = avg_grp) %>%
      group_by(avg_grp) %>%
        summarize(across(where(is.numeric), mean)) %>%
        select(-avg_grp) %>%
        data.matrix() })
  diag_corr(res[[1]], res[[2]])
}
