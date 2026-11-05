#' Label separate arrays within regions
#'
#' Given a data frame with the regions and the site numbers, labels the
#' individual arrays. The only special rule is the V1/V2 split; otherwise,
#' arrays are grouped by 96.
#'
#' @param df data frame with regions and site numbers

label_arrays <- function(df) {
  df %>%
    mutate(array = ceiling(site / 96)) %>%
    mutate(array = replace_when(array,
      region == "V1" & site > 61 ~ 2,
      region == "V2" & site > 35 ~ 2
    )) %>%
    mutate(array = factor(array))
}

#' Region relevel
#'
#' Helper to refactor the `region` column into the correct ordering.
#'
#' @param df data frame with region column.

order_regions <- function(df) {
  df %>%
    mutate(region = factor(region, levels = c("V1", "V2", "V4", "IT")))
}
