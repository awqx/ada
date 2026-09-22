nb_id <- "20260826"
here::i_am("scripts/dM_rd_split_half_corr.R")

library(dplyr)
library(tidyr)
library(here)
library(qs2)
# library(corpcor)
# library(Rfast)

# 0 Setup #####################################################

set.seed(as.integer(nb_id))
shrink <- T
return_sites <- T
shrink_nm <- ifelse(shrink, "shrinkage", "ordinary")

l_corr <- 0.3
u_corr <- 0.6

dichotomies <- here("data", "processed", "20260723", "dichotomies.csv") %>%
  read.csv()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  n_sites <- 30
  fam_id <- 13
  coherence <- 2
} else {
  n_sites <- as.integer(args[1])
  fam_id <- as.integer(args[2])
  coherence <- as.integer(args[3])
}

coherence_nm <- ifelse(coherence == 1, "noisy", "natur")

# 0.1 Directories =============================================

rd_dir <- here(
  "data", "processed", nb_id,
  "dM_rd", coherence_nm, shrink_nm)
if (!dir.exists(rd_dir)) dir.create(rd_dir, recursive = T)
data_dir <- here("data", "processed", "20260625")

# 0.2 Functions ===============================================

rdichot_multi_d <- function(dat, shrink = T, return_diff = T) {
  multi_d_list <- lapply(
    1:nrow(dichotomies), # 1:10
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

  if (return_diff) {
    multi_d_df <- data.frame(
      multi_d = sapply(multi_d_list, function(x) x[["multi_d"]]))

    mu_diff <- lapply(multi_d_list, function(x) x[["mu_diff"]]) %>%
      do.call(rbind, .)

    return(list(
      multi_d = multi_d_df,
      mu_diff = mu_diff ))
  } else {
    multi_d_df <- data.frame(multi_d = multi_d_list)
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
      corpcor::invcov.shrink(verbose = F)
  } else {
    inv_cov <- Rfast::spdinv(pooled.cov(dat, ina = grps)) }
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
idx <- (stimuli$coherence == coherence) & (stimuli$family == fam_id)
avgs <- lapply(regions,
  function(x) avgs[[x]][idx, ])
names(avgs) <- regions
stimuli <- stimuli[idx, ]

# 2 Randomized dichotomy calculations #########################

eval_multi_d <- function(rgn) {
  if (!dir.exists(here(rd_dir, rgn))) dir.create(here(rd_dir, rgn))

  dat <- avgs[[rgn]]
  matched_sites <- read.csv(here("results", "20260826", "split_corrs.csv")) %>%
    filter(grp_type == "img", n_blanks == 1, region == rgn) %>%
    filter(corr > l_corr, corr < u_corr)
  sites <- sample(matched_sites$site, n_sites)
  dat <- dat[, sites]

  res <- rdichot_multi_d(dat, shrink)

  if (return_sites)
    data.frame(
      site = sites,
      corr = matched_sites %>% filter(site %in% sites) %>% .$corr) %>%
      write.csv(
        here(rd_dir, rgn, paste0("fam", sprintf("%03d", fam_id), "_site_corr.csv")),
        row.names = F)

  write.csv(
    res$multi_d,
    here(rd_dir, rgn, paste0("fam", sprintf("%03d", fam_id), "_dM.csv")),
    row.names = F)
}

dummy <- lapply(
  regions,
  eval_multi_d
)
