nb_id <- "20260826"
here::i_am("scripts/rep_splits.R")
local <- F
shorten <- T

library(dplyr)
library(here)
library(cli)
library(qs2)

out_dir <- here("data", "processed", nb_id)
scratch <- "/scratch/axin/rep_splits"
args <- commandArgs(trailingOnly = TRUE)

region <- ifelse(length(args) > 0, args[1], "V1")
scratch <- file.path(scratch, region)
if (local) scratch <- here(out_dir, "splits", region)

if (!dir.exists(scratch))
  dir.create(scratch, recursive = T)
if (!dir.exists(here(out_dir, "splits", region)))
  dir.create(here(out_dir, "splits"), recursive = T)

stimuli <- read.csv(here(out_dir, "stimuli.csv"))
img_grp <- stimuli %>%
  group_by(coherence, family, id) %>%
  mutate(grp = cur_group_id()) %>%
  .$grp

mua <- here("data", "processed", "20260625", paste0(region, ".qs2")) %>%
    qs_read()
if (shorten) mua <- mua[, 51:60, 1:4]
n <- dim(mua)[3]

rep_grps <- combn(6, 3)
if (shorten) rep_grps <- rep_grps[, 1:3]

split_half_avg <- function(channel, avg_grp, split_grp) {
  channel %>%
    data.frame() %>%
    mutate(avg_grp = avg_grp, split_grp = split_grp) %>%
    group_by(split_grp) %>%
    group_split(.keep = F) %>%
    lapply(., function(df)
      df %>%
        group_by(avg_grp) %>%
        summarize(across(where(is.numeric), mean)) %>%
        select(-avg_grp) %>%
        data.matrix())
}

avg_mua <- function() {
  X_1 <- X_2 <- array(0, dim = c(length(unique(img_grp)), dim(mua)[2:3]))
  cli_progress_bar(name = region, total = n * ncol(rep_grps))

  for (k in 1:ncol(rep_grps)) {
    rep_grp <- stimuli$repetition %in% rep_grps[, k]
    for (i in 1:n) {
      site_split <- split_half_avg(mua[, , i], img_grp, rep_grp)
      X_1[, , i] <- site_split[[1]]
      X_2[, , i] <- site_split[[2]]
      cli_progress_update()
    }
    qs_save(
      list(split_1 = X_1, split_2 = X_2),
      file.path(scratch, paste0(paste0(rep_grps[, k], collapse = ""), ".qs2"))) }

  cli_progress_done()
}

avg_mua()

