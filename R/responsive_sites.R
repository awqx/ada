library(here)
source(file.path(here("R"), "onset_time.R"))

#' Split-half correlation helper
#'
#' Finds correlation of even and odd repetitions of the same stimulus
#'
#' @param channel the MUA readings for one channel
#' @param avg_grp the group labels for observations to average over
#' @param split_grp the group labels for splitting
#'
#' @return num Pearson correlation

split_half_corr <- function(channel, avg_grp, split_grp) {
  data.frame(channel, avg_grp, split_grp) %>%
    group_by(split_grp) %>%
    group_split(.keep = F) %>%
    sapply(., function(df)
      df %>%
        group_by(avg_grp) %>%
        summarize(channel = mean(channel)) %>%
        .$channel) %>%
    cor() %>%
    .[1, 2]
}

#' Find responsive channels
#'
#' Average the response for even and odd repetitions of each stimulus,
#' then retain channels with a Pearson correlation exceeding some threshold.
#'
#' @param mua numeric matrix of binned spiking data
#' @param onset_grp the group labels for contrasts for onset time
#' @param avg_grp the group labels for observations to average over
#' @param split_grp the group labels for splitting
#' @param stim_idx starting index for stimuli, default `50`
#' @param l_stim number of bins after `stim_idx` to look for `onset_time()`,
#'   default `20`
#' @param l_window number of bins after the calculated onset time to average for
#'   finding responsive sites
#'
#' @return int vector of responsive channels

responsive_sites <- function(
  mua,
  onset_grp,
  avg_grp,
  split_grp,
  onset = NULL,
  stim_idx = 50,
  l_stim = 20,
  l_window = 15
) {
  if (is.null(onset)) {
    onset <- onset_time(mua[, stim_idx:(stim_idx + l_stim), ], onset_grp) }

  mua_avg <- apply(
    mua[, onset:(onset + l_window - 1) + stim_idx - 1, ], c(1, 3), mean)
  sapply(
    1:dim(mua_avg)[2],
    function(x) split_half_corr(mua_avg[, x], avg_grp, split_grp))
}
