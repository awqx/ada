#' Average matrix rows by group
#'
#' Given an n x p matrix and an n-length vector, produce row averages based on
#' the categories of the vector. That is, the output should be a matrix with p
#' columns and with rows corresponding to the unique values of the vector.
#'
#' This functional is used elsewhere, but I made this a standalone function
#' later, so some functions may not implement it when appropriate.
#'
#' @param dat data frame
#' @param avg_grp vector of group labels to average the observations
#' @param incl_idx optional boolean vector to subset the data
#' @param return_grp logical (default `FALSE`) indicating whether to return the
#' column labeling the averaging group

avg_by_grp <- function(dat, avg_grp, incl_idx = NULL, return_grp = F) {
  if (!is.null(incl_idx)) {
    dat <- dat[incl_idx, ]
    if (length(avg_grp) > sum(incl_idx))
      avg_grp <- avg_grp[incl_idx]
  }

  res <- dat %>%
    mutate(avg_grp = avg_grp) %>%
    group_by(avg_grp) %>%
    summarize(across(where(is.numeric), mean))

  if (return_grp) return(res)

  select(res, -avg_grp)
}
