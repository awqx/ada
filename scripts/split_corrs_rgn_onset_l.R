here::i_am("scripts/split_corrs_rgn_onset_l.R")
library(here)
shorten <- F
source(here("R", "split_half_corr.R"))

library(dplyr)
library(here)
library(cli)
library(qs2)

data_dir <- here("data", "processed", "20260826")
out_dir <- here("data", "processed", "20260920")

args <- commandArgs(trailingOnly = TRUE)
stim_idx <- 50

if (length(args) > 0) {
  region <- args[1]
  onset <- as.numeric(args[2])
  l_window <- as.numeric(args[3])
} else {
  region <- "V1"
  onset <- 3
  l_window <- 10
}

if (!dir.exists(here(out_dir, "split_corrs", region)))
  dir.create(here(out_dir, "split_corrs", region), recursive = T)

# keep only a single blank image
# for alternate schemes, see scripts/split_corrs.R
stimuli_grp <- read.csv(here(data_dir, "stimuli_grp.csv")) %>%
  mutate(grp = replace_when(family, coherence != 0 ~ grp)) %>%
  .$grp
nonblank <- read.csv(here(data_dir, "stimuli_grp.csv")) %>%
  .$coherence > 0

f_nms <- list.files(here(data_dir, "splits", region))
if (shorten) f_nms <- f_nms[1:2]

dummy <- lapply(
  f_nms,
  function(f) {
    qs_read(here(data_dir, "splits", region, f)) %>%
      splits_corr(
        splits = .,
        avg_grp = stimuli_grp,
        onset = onset,
        l_window = l_window,
        # incl_idx = nonblank,
        stim_idx = stim_idx) }) %>%
    do.call(cbind, .) %>%
    write.csv(
      here(out_dir, "split_corrs", paste0(region, "_", onset, "_", l_window, ".csv")), row.names = F)
