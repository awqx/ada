library(here)
source(file.path(here("R"), "dprime.R"))

#' Find onset time in single channel
#'
#' Given a vector of distances, find the first location where the
#' absolute value of d' exceeds 0.4 and the following 2 windows are also
#' suprathreshold (and of the same sign).
#'
#' @param x1 a numeric matrix
#' @param x2 a numeric matrix
#' @param threshold d prime threshold for onset, defaults to `0.4`.
#'
#' @returns integer index of onset time if conditions met, else `NA`

onset_channel <- function(x1, x2, threshold = 0.4) {
  dp <- dprime(x1, x2)
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

#' Find onset time in session
#'
#' Given contrasting groups, return the median onset time in a session.
#'
#' @param mua numeric matrix of MUA
#' @param onset_grp two-element list of contrasting groups for d' calculation
#' @param threshold d' threshold for onset, defaults to `0.4`.
#' @param all_times whether to return the entire vector, defaults to `F`
#'
#' @return int of index of median onset

onset_time <- function(
    firings, onset_grp, threshold = 0.4, default_onset = 10, all_times = F
) {
  f1 <- firings[onset_grp[[1]], , ]
  f2 <- firings[onset_grp[[2]], , ]

  onsets <- sapply(
    1:dim(firings)[3],
    function(x) onset_channel(f1[, , x], f2[, , x], threshold = threshold))

  if (all_times)
    return(onsets)

  onset <- median(onsets, na.rm = T) - 2
  ifelse(is.na(onset), default_onset, floor(onset))
}
