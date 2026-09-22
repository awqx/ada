nb_id <- "20260826"
here::i_am("scripts/split_corrs.R")
shorten <- F

library(dplyr)
library(here)
library(cli)
library(qs2)

out_dir <- here("data", "processed", nb_id)
args <- commandArgs(trailingOnly = TRUE)

region <- ifelse(length(args) > 0, args[1], "V1")

if (!dir.exists(here(out_dir, "split_corrs", region)))
  dir.create(here(out_dir, "split_corrs", region), recursive = T)

stimuli_grp <- read.csv(here(out_dir, "stimuli_grp.csv"))

# keeps 100 blank images
by_img_blank100 <- stimuli_grp$grp
# keeps only a single blank image
by_img_blank001 <- stimuli_grp %>%
  mutate(grp = replace_when(family, coherence != 0 ~ grp)) %>%
  .$grp
by_fam <- stimuli_grp$family

# helper for only calculating the columnwise correlations
diag_corr <- function(x1, x2) {
  colSums(scale(x1) * scale(x2)) / (nrow(x1) - 1)
}

onsets <- readRDS(file.path(out_dir, "onsets_all.RDS"))
onset <- median(onsets[[region]], na.rm = T)
l_window <- 15
stim_idx <- 50

splits_corr <- function(splits, avg_grp, onset, incl_idx = NULL) {
  if (!is.null(incl_idx)) {
    splits <- lapply(splits, function(x) x[incl_idx, , ])
    avg_grp <- avg_grp[incl_idx]
  }

  res <- lapply(splits, function(split) {
    apply(split[, onset:(onset + l_window - 1) + stim_idx - 1, ], c(1, 3), mean) %>%
      data.frame() %>%
      mutate(avg_grp = avg_grp) %>%
      group_by(avg_grp) %>%
        summarize(across(where(is.numeric), mean)) %>%
        select(-avg_grp) %>%
        data.matrix() })
  diag_corr(res[[1]], res[[2]])
}

save_corrs <- function(grp_nm, incl_idx = NULL) {
  lapply(
    f_nms,
    function(f) {
      qs_read(here(out_dir, "splits", region, f)) %>%
        splits_corr(grps[[grp_nm]], onset, incl_idx)
    }) %>%
    do.call(cbind, .) %>%
    write.csv(
      here(out_dir, "split_corrs", region, paste0(grp_nm, ".csv")), row.names = F)
}

f_nms <- list.files(here(out_dir, "splits", region))
if (shorten) f_nms <- f_nms[1:2]

grps <- list(by_img_blank100, by_img_blank001, by_fam)
names(grps) <- c("img_blank100", "img_blank001", "fam_blank001")

dummy <- lapply(names(grps), save_corrs)

grps <- list(by_img_blank100, by_fam)
names(grps) <- c("img_blank000", "fam_blank000")
nonblank <- stimuli_grp$coherence > 0

dummy <- lapply(names(grps), save_corrs, nonblank)
