#' Extract info from file names
#'
#' Gets the recorded region, age, and number of subjects from the
#' MATLAB file name.
#'
#' @param f str file name
#' @returns a list with values `region`, `age`, and `n_subjects`
file_info <- function(f) {
  region <- f %>%
    str_extract("^zz_[[:alnum:]]{2}") %>%
    str_remove("^zz_")
  age <- f %>%
    str_extract("Textures2_[[:digit:]]+") %>%
    str_remove("Textures2_") %>%
    as.integer()
  n_subjects <- f %>%
    str_extract("[[:digit:]]{1}subj\\.mat$") %>%
    str_remove("subj\\.mat") %>%
    as.integer()

  list(region = region, age = age, n_subjects = n_subjects)
}

#' Extract and clean the stimuli from MATLAB data
#'
#' Labels the repetitions. The `img_id`
#'
#' @param x readMat output of `textures2` data
#' @returns A data frame with columns coherence, family, id, id_blank, img_id
clean_stimuli <- function(x) {
  stims <- data.frame(x$save.label)
  colnames(stims) <- c(
    "coherence", "family", "id", "coherence_t3", "x5", "id_blank"
  )
  stims %>%
    select(-coherence_t3, -x5) %>%
    group_by(coherence, family, id, id_blank) %>%
    mutate(repetition = row_number(), .groups = "drop_last") %>%
    ungroup() %>%
    group_by(family, id) %>%
    mutate(famsamp = cur_group_id())
}

#' Calculate distance d' (d prime)
#'
#' Given two response distributions, calculates d' as described in
#' Physiological Data Analysis - Initial Analysis
#'
#' @param x1 a numeric vector
#' @param x2 a numeric vector
#' @param na.rm whether NAs should be removed, defaults to `TRUE`.
#'
#' @returns numeric d' value

dprime <- function(x1, x2, na.rm = TRUE) {
  (mean(x1, na.rm = na.rm) - mean(x2, na.rm = na.rm)) /
    sqrt(0.5 * (
      var(x1, na.rm = na.rm) + var(x2, na.rm = na.rm)
    ))
}

#' Apply d' calculation over matrices
#'
#' A modification to the `dprime` calculation that takes advantage of
#' `Rfast` for multiple calculations.
#'
#' @param x1 a numeric matrix
#' @param x2 a numeric matrix
#' @param nbins int for number of columns
#' @param na.rm defaults to `TRUE`
#'
#' @returns numeric vector of d' values

apply_dprime <- function(x1, x2, nbins, na.rm = TRUE) {
  if (is.null(dim(x1))) {
    x1 <- matrix(x1, ncol = nbins)
    x2 <- matrix(x2, ncol = nbins)
  }
  (colMeans(x1, na.rm = na.rm) - colMeans(x2, na.rm = na.rm)) /
    sqrt(0.5 * (
      Rfast::colVars(x1, na.rm = na.rm) + Rfast::colVars(x2, na.rm = na.rm)
    ))
}

#' Find onset time in single channel
#'
#' Given a vector of distances, find the first location where the
#' absolute value of d' exceeds 0.4 and the following 2 windows are also
#' suprathreshold (and of the same sign).
#'
#' @param dp numeric vector of distance values
#' @param threshold d prime threshold for onset, defaults to `0.4`.
#'
#' @returns integer index of onset time if conditions met, else `NA`

onset_channel <- function(dp, threshold) {
  n <- length(dp)
  s <- sign(dp)

  pass_thresh <- abs(dp) > threshold

  pass_thresh  <- pass_thresh[1:(n - 2)]
  pass_window1 <- pass_thresh[2:(n - 1)] &
    (s[1:(n - 2)] == s[2:(n - 1)])
  pass_window2 <- pass_thresh[3:n] &
    (s[1:(n - 2)] == s[3:n])

  which(pass_thresh & pass_window1 & pass_window2)[1]
}

#' Find session onset
#'
#' Given contrasting indices, return the median onset time in a session.
#'
#' @param firings numeric matrix of spiking values
#' @param indices two-element list containing contrasting indices for d' calculation
#' @param threshold d' threshold for onset, defaults to `0.4`.
#'
#' @return int of index of median onset

onset_session <- function(firings, indices, threshold = 0.4) {
  f1 <- firings[indices[[1]], , ]
  f2 <- firings[indices[[2]], , ]

  dprimes <- mapply(
    apply_dprime,
    split(f1, slice.index(f1, 3)),
    split(f2, slice.index(f2, 3)),
    MoreArgs = list(nbins = dim(firings)[2])
  )

  onsets <- apply(dprimes, 2, onset_channel, threshold = threshold)
  median(onsets, na.rm = T) - 2
}

#' Split repetitions
#'
#' Given the data frame containing stimuli details, return the indices
#' corresponding to the odd and even repetitions of the non-blank stimuli.
#'
#' @param stims data frame containing stimuli details
#'
#' @return list of the odd and even indices, respectively

odd_even_reps <- function(stims) {
  list(
    odd = (stims$repetition %% 2 == 0 & stims$coherence != 0),
    even = (stims$repetition %% 2 == 1 & stims$coherence != 0)
  )
}

#' Grouped row operations
#'
#' Apply an operation to rows given a vector of grouped indices.
#'
#' @param mat numeric matrix
#' @param grps integer vector of group labels
#' @param FUN operation to apply to rows
#' @param ... additional arguments passed to `FUN`
#'
#' @return matrix of dimension length(unique(grp)) x ncol(mat)

apply_rowgroups <- function(mat, grps, FUN, ...) {
  data.frame(mat, grps) %>%
    tidyr::gather("channel", "firing", -grps) %>%
    group_by(grps, channel) %>%
    summarize(firing = FUN(firing, ...)) %>%
    tidyr::pivot_wider(names_from = channel, values_from = firing) %>%
    ungroup() %>%
    select(-grps) %>%
    as.matrix()
}

#' Find responsive channels
#'
#' Average the response for even and odd repetitions of each stimulus,
#' then retain channels with a Pearson correlation exceeding some threshold.
#'
#' From what I can tell, repetitions within a sample and family are also averaged.
#'
#' @param firings numeric matrix of binned spiking data
#' @param indices two-element list containing indices for the groups to compare
#' @param subgroups indices of subgroups to average over
#' @param times indices of time points to compare
#' @param threshold numeric correlation threshold, default `0.4`
#'
#' @return int vector of responsive channels

find_channels <- function(firings, indices, subgroups, times, threshold = 0.4) {
  avg_firings <- lapply(
    indices,
    function(idx) {
      subgrps <- data.frame(idx = idx, grp = subgroups) %>%
        group_by(idx, grp) %>%
        mutate(subgrp = cur_group_id()) %>%
        .$subgrp
      subgrps <- subgrps[idx]
      apply(firings[idx, times, ], c(1, 3), mean, na.rm = T) %>%
        apply_rowgroups(subgrps, mean, na.rm = T)
    }
  )
  cors <- diag(cor(avg_firings[[1]], avg_firings[[2]]))
  which(cors > threshold)
}

#' Spike summary from binned data and stimulus data frame
#'
#' Given stimuli and the binned data, subset the spikes into
#' onset time and responsive channels.
#'
#' Currently assumes the standard groupings (TODO: elaborate)
#'
#' @param stims output of `clean_stimuli`
#' @param firings binned spiking data
#'
#' @return numeric matrix of summarized spike data

spike_summary <- function(
  stims,
  firings,
  keep_pre_onset = TRUE,
  keep_unresponsive = TRUE,
  keep_blanks = TRUE
) {
  # Get onset times
  onset_time <- onset_session(
    firings,
    list(stims$coherence > 0, stims$coherence == 0)
  )
  times <- onset_time:(onset_time + 14)
  message("Onset time: ", onset_time)

  # Get responsive channels
  rep_split <- odd_even_reps(stims)
  responsive <- find_channels(firings, rep_split, stims$famsamp, times)

  # Turn firings into tibble
  firings <- apply(firings[, times, responsive], c(1, 3), mean, na.rm = T) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("stimulus") %>%
    tidyr::gather("channel", "firing", -stimulus) %>%
    mutate(channel = as.factor(as.integer(gsub("V", "", channel))))

  # Average firings by family and sample
  avg_firings <- data.frame(stims, firings) %>%
    select(-id_blank, -repetition, -famsamp) %>%
    filter(coherence > 0) %>%
    group_by(coherence, family, id, channel) %>%
    summarize(firing = mean(firing)) %>%
    tidyr::pivot_wider(names_from = channel, values_from = firing) %>%
    ungroup() %>%
    select(-coherence:-id) %>%
    as.matrix()

  avg_firings
}
