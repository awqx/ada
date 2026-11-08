#' Spearman-Brown reliability correction
#'
#' Predicts the reliability of a measure after changing its length by a
#' factor of `scaling`.
#'
#' @param r initial correlation (e.g., split-half correlation)
#' @param scaling scaling factor for test length adjustment, defaults to
#' `scaling = 2` (doubling, i.e., halves to the full set of repetitions)
#'
#' @returns the corrected reliability, same shape as `r`

spearman_brown <- function(r, scaling = 2) {
  (scaling * r) / (1 + (scaling - 1) * r)
}

#' Balanced split pair indices
#'
#' Enumerates every unique way to split the repetitions into two groups. For
#' an even number of repetitions, the groups are of equal size and each
#' partition is listed once (i.e., {1, 2} vs. {3, 4} is not repeated as
#' {3, 4} vs. {1, 2}). For an odd number of repetitions, the first group has
#' `floor(n_reps / 2)` repetitions and the second has the remainder.
#'
#' @param n_reps number of repetitions (>= 2)
#'
#' @returns list of two matrices (repetitions x splits), where column k of
#' each matrix gives the repetition indices of the two groups in split k

balanced_splits <- function(n_reps) {
  stopifnot(n_reps >= 2)
  half <- floor(n_reps / 2)
  combos <- combn(n_reps, half)

  # combn() lists combinations lexicographically, so the complement of column
  # k is column (N + 1 - k); keeping the first half (all combinations
  # containing repetition 1) gives each partition exactly once
  if (n_reps %% 2 == 0)
    combos <- combos[, seq_len(ncol(combos) / 2), drop = F]

  grps_2 <- matrix(
    apply(combos, 2, function(x) seq_len(n_reps)[-x]),
    ncol = ncol(combos))

  list(combos, grps_2)
}

#' Diagonal of the correlation matrix
#'
#' Pearson correlations between matching columns of two matrices, without
#' computing the full correlation matrix.
#'
#' @param x1 numeric matrix
#' @param x2 numeric matrix of the same dimension as `x1`
#'
#' @returns vector of `ncol(x1)` correlations

diag_corr <- function(x1, x2) {
  colSums(scale(x1) * scale(x2)) / (nrow(x1) - 1)
}

#' Mean firing over a time window
#'
#' @param mua array of trials x time bins x sites
#' @param onset integer number of bins after stimulus presentation at which
#' the window starts
#' @param l_window integer number of bins to average over
#' @param stim_idx integer bin of stimulus presentation, defaults to
#' `stim_idx = 50`
#'
#' @returns matrix of trials x sites

window_mean <- function(mua, onset, l_window, stim_idx = 50) {
  bins <- onset:(onset + l_window - 1) + stim_idx - 1
  Reduce(`+`, lapply(bins, function(b) mua[, b, ])) / l_window
}

#' Repetition index from stimulus labels
#'
#' Assumes trials are presorted, so that the order in which a label appears
#' gives its repetition. Labels appearing more than `n_reps` times (e.g., a
#' single label shared by all blank images) are divided in order into `n_reps`
#' consecutive chunks of equal size, i.e., the first `count / n_reps`
#' appearances are repetition 1, and so on.
#'
#' @param labels vector of stimulus labels for each trial
#' @param n_reps number of repetitions, defaults to the smallest label count
#'
#' @returns integer vector of repetitions (1 to `n_reps`) for each trial

rep_index <- function(labels, n_reps = NULL) {
  counts <- table(labels)
  if (is.null(n_reps)) n_reps <- min(counts)
  if (any(counts %% n_reps != 0))
    stop("Every label count must be a multiple of `n_reps`.")

  occurrence <- ave(seq_along(labels), labels, FUN = seq_along)
  per_rep <- as.vector(counts[as.character(labels)]) / n_reps
  as.integer(ceiling(occurrence / per_rep))
}

#' Split-half correlation
#'
#' For each split of the repetitions, averages the trials in each half by
#' stimulus label and correlates the two halves across stimuli. Correlations
#' are averaged in Fisher z-space and back-transformed. Repetitions are
#' assigned by `rep_index()`, so trials must be presorted.
#'
#' @param dat numeric matrix of trials x sites (e.g., from `window_mean()`), or
#' a vector of trials for a single site
#' @param labels vector of stimulus labels for the rows of `dat`; trials
#' labeled `NA` are excluded
#' @param splits optional list of splits, as returned by `balanced_splits()`;
#' a single split may also be given as a list of two index vectors. Defaults
#' to all balanced splits.
#' @param n_reps number of repetitions, passed to `rep_index()`
#' @param reliability whether to apply the Spearman-Brown correction, defaults
#' to `reliability = F`
#' @param ... optional arguments passed to `spearman_brown()`
#'
#' @returns vector of per-site split-half correlations (or reliabilities)

split_half_corr <- function(
  dat,
  labels,
  splits = NULL,
  n_reps = NULL,
  reliability = F,
  ...
) {
  dat <- as.matrix(dat)
  keep <- !is.na(labels)
  dat <- dat[keep, , drop = F]
  stim <- as.integer(factor(labels[keep]))
  reps <- rep_index(stim, n_reps)

  if (is.null(splits)) splits <- balanced_splits(max(reps))
  if (is.vector(splits[[1]])) splits[[1]] <- matrix(splits[[1]], ncol = 1)
  if (is.vector(splits[[2]])) splits[[2]] <- matrix(splits[[2]], ncol = 1)
  stopifnot(ncol(splits[[1]]) == ncol(splits[[2]]))

  # stimuli x sites mean over the trials in a set of repetitions
  half_mean <- function(rep_set) {
    idx <- reps %in% rep_set
    rowsum(dat[idx, , drop = F], stim[idx], reorder = T) /
      tabulate(stim[idx])
  }

  n_splits <- ncol(splits[[1]])
  z <- vapply(
    seq_len(n_splits),
    function(k) {
      atanh(diag_corr(half_mean(splits[[1]][, k]), half_mean(splits[[2]][, k])))
    },
    numeric(ncol(dat)))
  res <- tanh(rowMeans(matrix(z, ncol = n_splits)))

  if (reliability) {
    return(spearman_brown(res, ...))
  } else {
    return(res)
  }
}
