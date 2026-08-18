nb_id <- "20260715"
here::i_am("scripts/pf_multi_d.R")

library(dplyr)
library(tidyr)
library(here)
library(qs2)
library(cli)
library(corpcor)
library(Rfast)

# 0 Setup #####################################################

set.seed(20260715)
args <- commandArgs(trailingOnly = TRUE)
n_sites <- ifelse(length(args) == 0, 30, args[1])
shrink <- T

# 0.1 Directories =============================================

pf_dir <- here("data", "processed", nb_id, "pairwise_family")
data_dir <- here("data", "processed", "20260625")

# 0.2 Functions ===============================================

pairwise_multi_d <- function(dat, shrink = T, return_diff = T) {
  multi_d_list <- lapply(
    cli_progress_along(family_pairs),
    function(i) {
      grp_ids <- family_pairs[[i]]
      idx <- stimuli$family %in% grp_ids

      multi_d(
        dat = dat[idx, ],
        grps = stimuli$family[idx],
        shrink = shrink,
        return_diff = return_diff)})
  multi_d_df <- data.frame(
    multi_d = sapply(multi_d_list, function(x) x[["multi_d"]]),
    grps = sapply(family_pairs, function(x) paste0(x, collapse = "vs")))

  if (return_diff) {
    mu_diff <- lapply(multi_d_list, function(x) x[["mu_diff"]]) %>%
      do.call(rbind, .)
    combvar <- lapply(multi_d_list, function(x) x[["combvar"]]) %>%
      do.call(rbind, .)

    return(list(
      multi_d = multi_d_df,
      mu_diff = mu_diff,
      combvar = combvar
    ))
  } else {
    return(multi_d_df)}

}

multi_d <- function(dat, grps, shrink = F, return_diff = T) {
  grps <- as.numeric(as.factor(grps))
  mu_1 <- colMeans(dat[grps == 1, ])
  mu_2 <- colMeans(dat[grps == 2, ])
  mu_diff <- mu_2 - mu_1

  if (shrink) {
    inv_cov <- rbind(
      sweep(dat[grps == 1, ], 2, mu_1, FUN = "-"),
      sweep(dat[grps == 2, ], 2, mu_2, FUN = "-")) %>%
      invcov.shrink(verbose = F)
  } else {
    inv_cov <- spdinv(pooled.cov(dat, ina = grps)) }
  multi_d <- sqrt(t(mu_diff) %*% inv_cov %*% (mu_diff))

  if (return_diff) {
    return(list(multi_d = multi_d, mu_diff = mu_diff, combvar = diag(inv_cov)))
  } else {
    return(multi_d)
  }
}

# 1 Data loading ##############################################

regions <- c("V1", "V2", "V4", "IT")

avgs <- qs_read(here(data_dir, "avgs.qs2"))
resp_sites <- readRDS(here(data_dir, "resp_sites.RDS"))
stimuli <- read.csv(here(data_dir, "stimuli.csv"))

# Limit to naturalistic stimuli and reliable sites
natur_idx <- stimuli$coherence == 2
avgs <- lapply(regions,
  function(x) avgs[[x]][natur_idx, resp_sites[[x]]])
names(avgs) <- regions

stimuli <- stimuli[natur_idx, ]

# 2 Pairwise family calculations ##############################

family_pairs <- combn(unique(stimuli$family), 2, simplify = F)
# family_pairs <- family_pairs[1:20]

eval_multi_d <- function(region, n_sites, shrink) {
  n_sites_f <- paste("n", n_sites, sep = "")
  if (shrink)
    multi_d_dir <- here(pf_dir, "multi_d_shrink", n_sites_f, region)
  else
    multi_d_dir <- here(pf_dir, "multi_d", n_sites_f, region)
  if (!dir.exists(multi_d_dir))
    dir.create(multi_d_dir, recursive = T)

  dat <- avgs[[region]]
  sites <- sample(1:ncol(dat), n_sites)
  dat <- dat[, sites]

  res <- pairwise_multi_d(dat, shrink = shrink)
  colnames(res$mu_diff) <- sites
  colnames(res$combvar) <- sites

  write.csv(res$multi_d, here(multi_d_dir, "dprime.csv"), row.names = F)
  write.csv(res$mu_diff, here(multi_d_dir, "mu_diff.csv"), row.names = F)
  write.csv(res$combvar, here(multi_d_dir, "combvar.csv"), row.names = F)
}


dummy <- lapply(regions, eval_multi_d, n_sites = n_sites, shrink = shrink)
