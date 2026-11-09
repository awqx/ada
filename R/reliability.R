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

#' Nested random-effects ANOVA of stimulus responses
#'
#' Variance components for a balanced design with images nested in families
#' and repetitions crossed with both, fit for all sites at once from the
#' expected mean squares:
#'
#' y = mu + family + image(family) + rep + family:rep + error
#'
#' | source          | df              | E[MS]                               |
#' |-----------------|-----------------|-------------------------------------|
#' | family          | a - 1           | s2_e + b s2_fr + n s2_i + b n s2_f  |
#' | image(family)   | a (b - 1)       | s2_e + n s2_i                       |
#' | rep             | n - 1           | s2_e + b s2_fr + a b s2_r           |
#' | family:rep      | (a - 1)(n - 1)  | s2_e + b s2_fr                      |
#' | error           | a (b - 1)(n - 1)| s2_e                                |
#'
#' with a families, b images per family, and n repetitions (one trial per
#' image and repetition). Repetition terms absorb block-level drift, matching
#' the split-half setup in which halves are formed from repetition blocks.
#' Method-of-moments estimates are not truncated, so they can be negative.
#'
#' The simpler nested model without repetition terms,
#' y = mu + family + image(family) + error, is returned alongside (suffix
#' `_simple`). Its family and image mean squares are identical (balanced
#' design); its residual pools the rep, family:rep, and error sums of
#' squares, MS_E' = (SS_R + SS_FR + SS_E) / (a b (n - 1)). If repetition
#' effects exist, E[s2_f_simple] = s2_f + s2_fr / n and
#' E[s2_i_simple] = s2_i - (s2_r + s2_fr) / n.
#'
#' The intermediate model with a repetition main effect but no family:rep
#' term, y = mu + family + image(family) + rep + error, is also returned
#' (suffix `_rep`). Its residual pools the family:rep and error sums of
#' squares, MS_E'' = (SS_FR + SS_E) / ((a - 1)(n - 1) + a (b - 1)(n - 1)). Its
#' family component equals the simple model's; block offsets (s2_r) no longer
#' bias the image component, and only family:rep effects (s2_fr) would.
#'
#' @param dat numeric matrix of trials x sites (e.g., from `window_mean()`), or
#' a vector of trials for a single site
#' @param family vector of family labels for the rows of `dat`; trials
#' labeled `NA` are excluded
#' @param image vector of image labels, unique within each family
#' @param rep vector of repetition labels
#'
#' @returns data frame with one row per site: variance components of the full
#' model (`s2_f`, `s2_i`, `s2_fr`, `s2_e`, `s2_r`), the family share of
#' stimulus variance `rho = s2_f / (s2_f + s2_i)`, the mean squares (`ms_f`,
#' `ms_i`, `ms_r`, `ms_fr`, `ms_e`), and the simple model's components
#' (`s2_f_simple`, `s2_i_simple`, `s2_e_simple`, `rho_simple`) and the
#' repetition-main-effect model's (`s2_f_rep`, `s2_i_rep`, `s2_e_rep`,
#' `rho_rep`)

nested_anova <- function(dat, family, image, rep) {
  dat <- as.matrix(dat)
  keep <- !is.na(family)
  dat <- dat[keep, , drop = F]
  f <- as.integer(factor(family[keep]))
  fi <- as.integer(factor(paste(family[keep], image[keep], sep = "\r")))
  r <- as.integer(factor(rep[keep]))
  fr <- as.integer(factor(paste(f, r, sep = "\r")))

  a <- max(f)
  n <- max(r)
  b <- max(fi) / a
  if (b != round(b) || any(tabulate(fi) != n) || any(tabulate(fr) != b) ||
      nrow(dat) != a * b * n)
    stop("Design must be balanced: every image once per repetition, and the ",
         "same number of images in every family.")

  # cell means; rowsum() orders groups by integer code, so index by code
  grand <- colMeans(dat)
  m_f <- rowsum(dat, f, reorder = T) / (b * n)
  m_fi <- rowsum(dat, fi, reorder = T) / n
  m_r <- rowsum(dat, r, reorder = T) / (a * b)
  m_fr <- rowsum(dat, fr, reorder = T) / b
  f_of_fi <- f[match(seq_len(a * b), fi)]
  f_of_fr <- f[match(seq_len(a * n), fr)]
  r_of_fr <- r[match(seq_len(a * n), fr)]

  centered <- function(x, y) sweep(x, 2, y)
  ss_f <- b * n * colSums(centered(m_f, grand)^2)
  ss_i <- n * colSums((m_fi - m_f[f_of_fi, , drop = F])^2)
  ss_fr <- b * colSums(centered(
    m_fr - m_f[f_of_fr, , drop = F] - m_r[r_of_fr, , drop = F], -grand)^2)
  ss_e <- colSums((dat - m_fi[fi, , drop = F] - m_fr[fr, , drop = F] +
    m_f[f, , drop = F])^2)

  ms_f <- ss_f / (a - 1)
  ms_i <- ss_i / (a * (b - 1))
  ms_fr <- ss_fr / ((a - 1) * (n - 1))
  ms_e <- ss_e / (a * (b - 1) * (n - 1))

  s2_f <- (ms_f - ms_i - ms_fr + ms_e) / (b * n)
  s2_i <- (ms_i - ms_e) / n
  # repetition main effect (only needed for s2_r and the simple model)
  ss_r <- a * b * colSums(centered(m_r, grand)^2)
  ms_r <- ss_r / (n - 1)

  # simple model: rep, family:rep, and error pooled into the residual
  ms_e_simple <- (ss_r + ss_fr + ss_e) / (a * b * (n - 1))
  s2_f_simple <- (ms_f - ms_i) / (b * n)
  s2_i_simple <- (ms_i - ms_e_simple) / n

  # intermediate model: rep main effect kept, family:rep pooled with error
  ms_e_rep <- (ss_fr + ss_e) / ((a - 1) * (n - 1) + a * (b - 1) * (n - 1))
  s2_i_rep <- (ms_i - ms_e_rep) / n

  data.frame(
    s2_f = s2_f,
    s2_i = s2_i,
    s2_fr = (ms_fr - ms_e) / b,
    s2_e = ms_e,
    rho = s2_f / (s2_f + s2_i),
    s2_r = (ms_r - ms_fr) / (a * b),
    ms_f = ms_f,
    ms_i = ms_i,
    ms_r = ms_r,
    ms_fr = ms_fr,
    ms_e = ms_e,
    s2_f_simple = s2_f_simple,
    s2_i_simple = s2_i_simple,
    s2_e_simple = ms_e_simple,
    rho_simple = s2_f_simple / (s2_f_simple + s2_i_simple),
    s2_f_rep = s2_f_simple,
    s2_i_rep = s2_i_rep,
    s2_e_rep = ms_e_rep,
    rho_rep = s2_f_simple / (s2_f_simple + s2_i_rep))
}

#' Split-half correlation of one site's stimulus x repetition matrix
#'
#' Same estimator as `split_half_corr()` (Fisher z average over balanced
#' splits, optional Spearman-Brown correction), for data already arranged as
#' stimuli x repetitions. Faster when called many times, e.g., in simulations.
#'
#' @param Y numeric matrix of stimuli x repetitions
#' @param splits optional list of splits, as returned by `balanced_splits()`
#' @param reliability whether to apply the Spearman-Brown correction, defaults
#' to `reliability = F`
#' @param ... optional arguments passed to `spearman_brown()`
#'
#' @returns split-half correlation (or reliability)

split_half_reps <- function(Y, splits = NULL, reliability = F, ...) {
  if (is.null(splits)) splits <- balanced_splits(ncol(Y))
  z <- vapply(
    seq_len(ncol(splits[[1]])),
    function(k) atanh(cor(
      rowMeans(Y[, splits[[1]][, k], drop = F]),
      rowMeans(Y[, splits[[2]][, k], drop = F]))),
    numeric(1))
  res <- tanh(mean(z))
  if (reliability) spearman_brown(res, ...) else res
}
