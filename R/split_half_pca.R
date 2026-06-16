#' Split-half average helper
#'
#' Splits data into groups then averages within subgroups.
#'
#' @param channel the MUA readings for one channel
#' @param avg_grp the group labels for observations to average over
#' @param split_grp the group labels for splitting
#'
#' @return list of length `unique(split_grp)`

split_half_avg <- function(channel, avg_grp, split_grp) {
  data.frame(channel, avg_grp, split_grp) %>%
    group_by(split_grp) %>%
    group_split(.keep = F) %>%
    sapply(., function(df)
      df %>%
        group_by(avg_grp) %>%
        summarize(channel = mean(channel)) %>%
        .$channel)
}

#' Column centering helper
#'
#' Centers a matrix by subtracting column means.
#'
#' @param x numeric matrix
#'
#' @return matrix with same dimensions as `x`

center_cols <- function(x) {x - colMeans(x)[col(x)]}

#' Split half PCA
#'
#' Reports variance explained for training and test set PCA.
#'
#' @param x numeric matrix
#' @param avg_grp the group labels for observations to average over
#' @param split_grp the group labels for splitting
#'
#' @return list where `trn_var` reports the eigenvalues on training and `tst_rsq` reports R-squared of reconstruction on test.

split_half_pca <- function(x, avg_grp, split_grp) {
  # preprocessing
  x_splits <- data.frame(x) %>%
    mutate(avg_grp = avg_grp, split_grp = split_grp) %>%
    group_by(split_grp) %>%
    group_split(.keep = F) %>%
    lapply(., function(df) {
      df %>%
        group_by(avg_grp) %>%
        summarize(across(everything(), ~ mean(.x, na.rm = T))) %>%
        data.frame() %>%
        select(-avg_grp) %>%
        center_cols() %>%
        t() })

  # already centered in preprocessing above
  x_pca_trn <- prcomp(x_splits[[1]], center = F)
  trn_pcs <- x_pca_trn$rotation

  # test set scores
  tst_scores <- predict(x_pca_trn, x_splits[[2]])
  # test set reconstruction/loadings
  tst_loadings <- tst_scores %*% t(trn_pcs)

  trn_var <- x_pca_trn$sdev ^ 2
  tst_rsq <- sapply(
    1:ncol(trn_pcs),
    function(i) {
      tst_loadings <- tst_scores[, 1:i] %*% t(trn_pcs[, 1:i])
      data.frame(
        tst_true = matrix(x_splits[[2]], ncol = 1),
        tst_pred = matrix(tst_loadings, ncol = 1)) %>%
        lm(tst_true ~ tst_pred, .) %>%
        summary() %>%
        .$r.squared })

  list(trn_var = trn_var, tst_rsq = tst_rsq)
}

#' Split half PCA data frame formatting
#'
#' Formats output of `split_half_pca()` into useful data frame
#'
#' @param split_pca list output of `split_half_pca()`
#'
#' @return data frame with columns `n_pcs`, `trn_var`, `tst_rsq`, `trn_var_prop` (rescaled `trn_var`), `tst_rsq_prop` (rescaled `tst_rsq`), and `scaled_pcs` (rescaled `n_pcs`)

split_pca_df <- function(split_pca) {
  split_pca %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "n_pcs") %>%
    mutate(trn_var_prop = cumsum(trn_var) / sum(trn_var)) %>%
    mutate(tst_rsq_prop = tst_rsq / max(tst_rsq)) %>%
    mutate(n_pcs = as.numeric(n_pcs)) %>%
    mutate(scaled_pcs = (n_pcs - 1) / max(n_pcs - 1))
}
