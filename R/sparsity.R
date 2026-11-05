#' Sparsity with ratio of moments
#'
#' Calculate the sparsity with the ratio of moments (similar to pq-mean)
#' (Rolls and Tovee, 1995; Vinje and Gallant, 2000).
#'
#' @param x numeric vector

rm_sparsity <- function(x) {
  if (any(x < 0))
    stop("Ratio of moments sparsity cannot be calculated on negative values.")

  scalar <- 1 / (1 - 1 / length(x))
  m1 <- mean(x) ^ 2
  m2 <- mean(x ^ 2)
  (1 - m1 / m2) * scalar
}

#' Lp norm
#'
#' Helper function to calculate the Lp norm for a vector.
#'
#' @param x numeric vector
#' @param p numerioc for Lp norm

Lp_norm <- function(x, p = 1) {
  if (is.infinite(p))
    return(max(abs(x)))
  (sum(abs(x) ^ p))^(1 / p)
}

#' Generalized p mean
#'
#' Helper function to calculate the generalized p mean for a vector.
#'
#' @param x numeric vector
#' @param p numeric for generalized p mean

p_mean <- function(x, p = 1) {
  (mean(x ^ p))^(1 / p)
}

#' pq-mean sparsity
#'
#' Calculate the normalized ratio of the generalized p and q mean, defaulting to p = 1, q = 2, similar to the L2/L1 ratio.
#'
#' @param x numeric vector
#' @param p numeric for generalized p-mean of numerator
#' @param q numeric for generalized q-mean of denominator

pq_sparsity <- function(x, p = 1, q = 2) {
  - p_mean(x, p) / p_mean(x, q)
}

#' Normalized Lp/Lq norm sparsity
#'
#' Calculate the normalized ratio of the Lp and Lq norm, defaulting to the L1 to
#' L2 ratio (Hoyer)
#'
#' @param x numeric vector
#' @param p int, where numerator takes Lp norm
#' @param q int, where denominator takes Lq norm
#' @param rescale logical, indicating whether to rescale to bound output between
#' 0 and 1. Defaults to `TRUE`.

Lp_sparsity <- function(x, p = 1, q = 2, rescale = T) {
  res <- Lp_norm(x, p) / Lp_norm(x, q)
  if (rescale) {
    n <- length(x)
    res <- (sqrt(n) - res) / (sqrt(n) - 1)
  }
  res
}

#' Calculate sparsity given a method
#'
#' Calculate the per-column sparsity.
#'
#' @param dat numeric matrix
#' @param sparsity_method function to calculate sparsity of a vector

# col_sparsities <- function(dat, sparsity_method) {
#   lapply(
#     regions,
#     function(rgn) {
#       dat[[rgn]] %>%
#         apply(2, sparsity_method) %>%
#         data.frame(sparsity = ., region = rgn) %>%
#         mutate(site = row_number()) }) %>%
#     do.call(rbind, .)
# }

col_sparsities <- function(dat, sparsity_method) {
  dat %>%
    apply(2, sparsity_method) %>%
    tibble::as_tibble_col(column_name = "sparsity") %>%
    mutate(site = row_number())
}
