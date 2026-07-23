nb_id <- "20260715"
here::i_am("scripts/pf_lda.R")

library(dplyr)
library(tidyr)
library(here)
library(qs2)
library(cli)

# 0 Setup #####################################################

set.seed(20260715)
args <- commandArgs(trailingOnly = TRUE)
n_sites <- ifelse(length(args) == 0, 30, args[1])

# 0.1 Directories =============================================

out_dir <- here("data", "processed", nb_id)
t2_dir <- here("data", "raw", "6reps")
data_dir <- here("data", "processed", "20260625")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = T)

# 0.2 Functions ===============================================

dprime <- function(x, grps, return_components = T) {
  idx <- grps == unique(grps)[1]
  mu_diff <- mean(x[idx]) - mean(x[-idx])
  sd_comb <- sqrt(0.5 * (var(x[idx]) + var(x[-idx])))

  if (return_components)
    return(list(mu = mu_diff, sd = sd_comb))

  mu_diff / sd_comb }

dprime_tst <- function(dat, grps, trn_idx, return_trn = T) {
  fit <- MASS::lda(x = dat, grouping = grps, subset = trn_idx)

  tst_scores <- predict(fit, dat[-trn_idx, ])$x
  tst_grps <- grps[-trn_idx]

  if (return_trn) {
    trn_scores <- predict(fit, dat[trn_idx, ])$x
    trn_grps <- grps[trn_idx]

    return(list(
      trn = dprime(trn_scores, trn_grps),
      tst = dprime(tst_scores, tst_grps),
      scaling = as.numeric(fit$scaling)))
  } else {
    return(c(
      dprime(tst_scores, tst_grps),
      scaling = list(as.numeric(fit$scaling))))}
}

kfold_dprime <- function(dat, grps, cv_grps, return_trn = T) {
  k_tsts <- lapply(unique(cv_grps),
    function(grp) {
      trn_idx <- cv_grps != grp
      c(dprime_tst(dat, grps, trn_idx), k = grp) })

  # Reformat
  if (return_trn) {
    dprimes <- lapply(k_tsts,
      function(x) {
        rbind(
          x[["tst"]] %>% data.frame() %>% mutate(data_split = "tst"),
          x[["trn"]] %>% data.frame() %>% mutate(data_split = "trn")) %>%
          mutate(k = x[["k"]])})
  } else {
    dprimes <- lapply(k_tsts,
      function(x) data.frame(x[names(x) != "scaling"])) }

  dprimes <- do.call(rbind, dprimes) %>%
    mutate(dprime = mu / sd, .before = k)

  scalings <- lapply(k_tsts, function(x) x[["scaling"]]) %>%
    do.call(rbind, .)

  return(list(dprimes = dprimes, scalings = scalings))
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

pf_dir <- here(out_dir, "pairwise_family", "lda")
if (!dir.exists(pf_dir)) dir.create(pf_dir, recursive = T)

pairwise_dprimes <- function(region, n_sites) {

  dprime_f <- paste("dprime", n_sites, sep = "_")
  scalings_f <- paste("scalings", n_sites, sep = "_")

  if (!dir.exists(here(pf_dir, dprime_f, region)))
    dir.create(here(pf_dir, dprime_f, region), recursive = T)
  if (!dir.exists(here(pf_dir, scalings_f, region)))
    dir.create(here(pf_dir, scalings_f, region), recursive = T)

  dat <- avgs[[region]]
  sites <- sample(1:ncol(dat), n_sites)
  dat <- dat[, sites]

  res <- lapply(
    family_pairs,
    function(grp_ids) {
      idx <- stimuli$family %in% grp_ids

      # produce the 15-fold CV results for each pair
      res <- kfold_dprime(
        dat = dat[idx, ],
        grps = factor(stimuli$family[idx]),
        cv_grps = stimuli$id[idx])

      # save each result separately
      f_nm <- paste0(paste0(grp_ids, collapse = "vs"), ".csv")
      colnames(res$scalings) <- sites

      write.csv(res$dprimes, here(pf_dir, dprime_f, region, f_nm), row.names = F)
      write.csv(res$scalings, here(pf_dir, scalings_f, region, f_nm), row.names = F)
})}


dummy <- lapply(regions, pairwise_dprimes, n_sites = n_sites)
