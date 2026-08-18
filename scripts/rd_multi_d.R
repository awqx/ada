nb_id <- "20260723"
here::i_am("scripts/rd_multi_d.R")

library(dplyr)
library(tidyr)
library(here)
library(qs2)
library(cli)
library(corpcor)
library(Rfast)

# 0 Setup #####################################################

set.seed(as.integer(nb_id))
shrink <- T
n_samp <- 15
n_dichot <- 2000

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  n_sites <- 30
  fam_id <- 13
} else {
  n_sites <- as.integer(args[1])
  fam_id <- as.integer(args[2])
}

# 0.1 Directories =============================================

rd_dir <- here(
  "data", "processed", nb_id,
  "rand_dichot", paste0("fam", sprintf("%03d", fam_id)))
if (!dir.exists(rd_dir)) dir.create(rd_dir, recursive = T)
data_dir <- here("data", "processed", "20260625")

# 0.2 Functions ===============================================

rdichot_multi_d <- function(dat, shrink = T, return_diff = T) {
  multi_d_list <- lapply(
    cli_progress_along(1:n_dichot),
    function(i) {
      dichot <- as.integer(dichotomies[i, ])
      vs_1 <- which(dichot == 1)
      vs_2 <- which(dichot == -1)
      grps <- ifelse(stimuli$id %in% vs_1, "1", "2")
      idx <- stimuli$id %in% c(vs_1, vs_2)

      multi_d(
        dat = dat[idx, ],
        grps = grps[idx],
        shrink = shrink,
        return_diff = return_diff)})
  multi_d_df <- data.frame(
    multi_d = sapply(multi_d_list, function(x) x[["multi_d"]]),
    grps = apply(dichotomies, MARGIN = 1, function(dichot) {
      paste(
        paste0(which(dichot == 1), collapse = "_"),
        which(dichot == 0),
        sep = "omit")}))

  if (return_diff) {
    mu_diff <- lapply(multi_d_list, function(x) x[["mu_diff"]]) %>%
      do.call(rbind, .)
    combvar <- lapply(multi_d_list, function(x) x[["combvar"]]) %>%
      do.call(rbind, .)

    return(list(
      multi_d = multi_d_df,
      mu_diff = mu_diff,
      combvar = combvar ))
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

# Limit to naturalistic, family of interest, and reliable sites
idx <- (stimuli$coherence == 2) & (stimuli$family == fam_id)
avgs <- lapply(regions,
  function(x) avgs[[x]][idx, resp_sites[[x]]])
names(avgs) <- regions
stimuli <- stimuli[idx, ]

dichotomies <- here("data", "processed", nb_id, "dichotomies.csv") %>%
  read.csv()

# 2 Randomized dichotomy calculations #########################

eval_multi_d <- function(region, n_sites, shrink) {
  n_sites_f <- paste("n", n_sites, sep = "")
  if (shrink)
    multi_d_dir <- here(rd_dir, "multi_d_shrink", n_sites_f, region)
  else
    multi_d_dir <- here(rd_dir, "multi_d", n_sites_f, region)
  if (!dir.exists(multi_d_dir))
    dir.create(multi_d_dir, recursive = T)

  dat <- avgs[[region]]
  sites <- sample(1:ncol(dat), n_sites)
  dat <- dat[, sites]

  res <- rdichot_multi_d(dat, shrink = shrink)
  colnames(res$mu_diff) <- sites
  colnames(res$combvar) <- sites

  write.csv(res$multi_d, here(multi_d_dir, "dprime.csv"), row.names = F)
  write.csv(res$mu_diff, here(multi_d_dir, "mu_diff.csv"), row.names = F)
  write.csv(res$combvar, here(multi_d_dir, "combvar.csv"), row.names = F)
}


dummy <- lapply(regions, eval_multi_d, n_sites = n_sites, shrink = shrink)
