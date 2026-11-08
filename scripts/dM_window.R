here::i_am("scripts/dM_window.R")
library(here)
shorten <- T
source(here("R", "dM.R"))

library(dplyr)
library(here)
library(cli)
library(qs2)

data_dir <- here("data", "processed", "20260625")
out_dir <- here("data", "processed", "20260920", "pf_dM")

args <- commandArgs(trailingOnly = TRUE)
STIM_IDX <- 50

if (length(args) > 0) {
  region <- args[1]
  onset <- as.numeric(args[2])
  l <- as.numeric(args[3])
} else {
  region <- "IT"
  onset <- 5
  l <- 10
}

if (!dir.exists(here(out_dir, "pf_dM")))
  dir.create(here(out_dir, "pf_dM"), recursive = T)

dM_window <- function(mua, onset, l, grps) {
  start_idx <- STIM_IDX + onset
  dat <- mua[, start_idx:(start_idx + l), ] %>%
    apply(c(1, 3), mean)
  dM_pairwise(dat, grps, shrink = T, return_diff = F)
}

mua <- here(data_dir, paste0(region, ".qs2")) %>%
  qs_read()
stimuli <- here("data", "processed", "20260826", "stimuli.csv") %>%
  read.csv()

incl_idx <- (stimuli$coherence == 2) & (stimuli$family != 0)
mua <- mua[incl_idx, , ]
stimuli <- stimuli[incl_idx, ]

# can change this later to sample
grps <- stimuli$family

dM_window(mua, onset, l, grps) %>%
  write.csv(here(out_dir, paste0(region, "_", onset, "_", l, ".csv")), row.names = F)

