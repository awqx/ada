

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
