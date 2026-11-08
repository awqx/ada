## Refined split-half correlations (see notebooks/20261104)
##
## For each region, onset scheme (earlier/median), and window length
## (10/15 bins = 100/150 ms), computes per-site split-half correlations over
## all 10 unique 3-vs-3 splits of the 6 repetitions. Correlations are averaged
## in Fisher z-space, back-transformed, and Spearman-Brown corrected.
##
## Data: raw MUA (trials x 10 ms bins x sites) from data/processed/20260625,
## averaged by image (coherence, family, id). Blanks are either excluded
## (none), averaged into a single image (aggregated), or treated as 100 exemplars
## with 6 repetitions each (separated).
##
## Usage: Rscript scripts/split_corrs_refined.R [region ...]

nb_id <- "20261104"
here::i_am("scripts/split_corrs_refined.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(qs2)
})
source(here("R", "reliability.R"))

out_dir <- here("data", "processed", nb_id)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
regions <- if (length(args) > 0) args else c("V1", "V2", "V4", "IT")

stim_idx <- 50 # bin of stimulus presentation
windows <- c(10, 15) # bins (100 ms, 150 ms)

# onset schemes (bins after stimulus presentation)
onsets_all <- readRDS(here("data", "processed", "20260826", "onsets_all.RDS"))
onset_schemes <- list(
  earlier = c(V1 = 3, V2 = 4, V4 = 5, IT = 6),
  median = sapply(onsets_all, function(x) floor(median(x, na.rm = TRUE))))

# ---- trial grouping ---------------------------------------------------------

# 20260826 and 20260625 stimuli share trial order; they differ only in blank
# labels. 20260826 labels blanks with id 1-100, each shown once per repetition.
stimuli <- read.csv(here("data", "processed", "20260826", "stimuli.csv"))
blank <- stimuli$coherence == 0
img <- paste(stimuli$coherence, stimuli$family, stimuli$id)

# stimulus labels per trial under each blank-handling scheme (NA = excluded).
# Trials are presorted into 6 blocks (one per repetition) of 1050 images and
# 100 blanks, so repetitions follow from label order (see rep_index()):
#   none:       blanks excluded
#   aggregated: all blanks share one label; its 600 trials split in order into
#               6 chunks of 100 (one per block), averaged as a single image
#   separated:  blank ids 1-100 (position within the blank block) are treated
#               as 100 images, each appearing once per block
blank_schemes <- list(
  none = ifelse(blank, NA, img),
  aggregated = ifelse(blank, "blank", img),
  separated = img)

splits <- balanced_splits(max(stimuli$repetition))

# ---- main loop --------------------------------------------------------------

res <- list()
for (region in regions) {
  t0 <- Sys.time()
  mua <- qs_read(here("data", "processed", "20260625", paste0(region, ".qs2")))
  stopifnot(dim(mua)[1] == nrow(stimuli))
  n_sites <- dim(mua)[3]

  for (scheme in names(onset_schemes)) {
    onset <- onset_schemes[[scheme]][[region]]
    for (l_window in windows) {
      win_mean <- window_mean(mua, onset, l_window, stim_idx)

      for (blanks in names(blank_schemes)) {
        corr <- split_half_corr(win_mean, blank_schemes[[blanks]], splits)
        res[[length(res) + 1]] <- data.frame(
          region = region,
          site = seq_len(n_sites),
          onset_type = scheme,
          onset = onset,
          ms = 10 * l_window,
          blanks = blanks,
          corr = corr,
          reliability = spearman_brown(corr))
      }
    }
  }
  rm(mua); gc(verbose = FALSE)
  message(region, ": ", n_sites, " sites, ",
          round(as.numeric(Sys.time() - t0, units = "secs"), 1), " s")
}

res <- do.call(rbind, res)
f_out <- if (length(args) > 0)
  paste0("split_corrs_", paste(regions, collapse = "_"), ".csv") else
  "split_corrs.csv"
write.csv(res, file.path(out_dir, f_out), row.names = FALSE)
message("Saved ", file.path(out_dir, f_out))
