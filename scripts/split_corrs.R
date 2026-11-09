## Per-site split-half correlations and Spearman-Brown reliabilities
##
## Computes every combination of region x onset x window x stimulus set
## (all/noisy/naturalistic) x averaging (image/family) x blank handling in a
## settings file (default config/split_corrs.yml) and writes one long table,
## plus a provenance record (<out>.json). Method: R/reliability.R. Background:
## notebooks/20261104.
##
## Usage:
##   Rscript scripts/split_corrs.R                         # full grid
##   Rscript scripts/split_corrs.R --regions V1,IT         # subset of regions
##   Rscript scripts/split_corrs.R --n_sites 5 --out x.csv # quick test
##   Rscript scripts/split_corrs.R --dry_run               # print settings only
##   Rscript scripts/split_corrs.R --help                  # all options

here::i_am("scripts/split_corrs.R")
suppressPackageStartupMessages({
  library(here)
  library(qs2)
})
source(here("R", "config.R"))
source(here("R", "reliability.R"))

# ---- options -----------------------------------------------------------------

opts <- optparse::parse_args(optparse::OptionParser(
  usage = "Rscript scripts/split_corrs.R [options]",
  option_list = list(
    optparse::make_option(
      "--config", default = here("config", "split_corrs.yml"),
      help = "settings file [default: %default]"),
    optparse::make_option(
      "--regions", default = NULL,
      help = "comma-separated regions, e.g., V1,IT [default: all in config]"),
    optparse::make_option(
      "--out", default = NULL,
      help = "output CSV [default: `out` in config]"),
    optparse::make_option(
      "--n_sites", type = "integer", default = 0L,
      help = "if > 0, only the first n sites per region, for testing [default: %default]"),
    optparse::make_option(
      "--dry_run", action = "store_true", default = F,
      help = "print the settings and exit"))))
opts$help <- NULL

cfg <- config::get(file = opts$config)
opts$regions <- if (is.null(opts$regions)) unlist(cfg$regions) else
  strsplit(opts$regions, ",", fixed = T)[[1]]
if (is.null(opts$out)) opts$out <- resolve_path(cfg$out)
regions <- opts$regions
out <- opts$out

# ---- settings ----------------------------------------------------------------

onsets_all <- readRDS(data_path("onsets"))
onset_bins <- function(scheme, region) {
  spec <- cfg$onsets[[scheme]]
  if (identical(spec, "median")) floor(median(onsets_all[[region]], na.rm = T))
  else spec[[region]]
}

grid <- expand.grid(
  window = unlist(cfg$windows),
  onset_type = names(cfg$onsets),
  stringsAsFactors = F)

# Stimulus labels per trial (NA = excluded) for each combination of stimulus
# set, averaging, and blank handling. Trials are presorted into one block per
# repetition, so labels and n_reps determine the repetitions (see
# rep_index()); n_reps is given explicitly because family labels appear
# several times per repetition.
stimuli <- read.csv(data_path("stimuli"))
n_reps <- max(stimuli$repetition)
blank <- stimuli$coherence == 0

stim_label <- function(stim_set, average, blanks) {
  keep <- switch(
    stim_set,
    all = stimuli$coherence > 0,
    noisy = stimuli$coherence == 1,
    naturalistic = stimuli$coherence == 2)
  lab <- switch(
    average,
    image = paste(stimuli$coherence, stimuli$family, stimuli$id),
    family = paste(stimuli$coherence, "family", stimuli$family))
  lab[!keep] <- NA
  lab[blank] <- switch(
    blanks,
    none = NA,
    aggregated = "blank",
    separated = paste("blank", stimuli$id[blank]))
  lab
}

schemes <- do.call(rbind, lapply(names(cfg$averages), function(average)
  expand.grid(
    stimuli = unlist(cfg$stimuli),
    average = average,
    blanks = unlist(cfg$averages[[average]]),
    stringsAsFactors = F)))
labels <- lapply(seq_len(nrow(schemes)), function(j)
  stim_label(schemes$stimuli[j], schemes$average[j], schemes$blanks[j]))

# check: label order reproduces the recorded repetitions in every scheme
for (lab in labels) {
  keep <- !is.na(lab)
  stopifnot(identical(rep_index(lab[keep], n_reps), stimuli$repetition[keep]))
}

splits <- balanced_splits(n_reps)

message(
  "Regions: ", paste(regions, collapse = ", "), "\n",
  "Onsets: ", paste(names(cfg$onsets), collapse = ", "),
  " | windows: ", paste(cfg$windows, collapse = ", "), " bins\n",
  "Schemes (stimuli / average / blanks): ",
  paste(apply(schemes, 1, paste, collapse = "/"), collapse = ", "), "\n",
  "Output: ", out)
if (isTRUE(opts$dry_run)) quit(save = "no")

# ---- compute -----------------------------------------------------------------

res <- list()
for (region in regions) {
  t0 <- Sys.time()
  mua <- qs_read(data_path("mua", region = region))
  stopifnot(dim(mua)[1] == nrow(stimuli))
  if (opts$n_sites > 0) mua <- mua[, , seq_len(min(opts$n_sites, dim(mua)[3])), drop = F]

  for (i in seq_len(nrow(grid))) {
    onset <- onset_bins(grid$onset_type[i], region)
    win_mean <- window_mean(mua, onset, grid$window[i], cfg$stim_idx)

    for (j in seq_len(nrow(schemes))) {
      corr <- split_half_corr(win_mean, labels[[j]], splits, n_reps = n_reps)
      res[[length(res) + 1]] <- data.frame(
        region = region,
        site = seq_along(corr),
        onset_type = grid$onset_type[i],
        onset = onset,
        ms = 10 * grid$window[i],
        stimuli = schemes$stimuli[j],
        average = schemes$average[j],
        blanks = schemes$blanks[j],
        corr = corr,
        reliability = spearman_brown(corr))
    }
  }
  rm(mua); invisible(gc())
  message(region, ": ", dim(win_mean)[2], " sites, ",
          round(as.numeric(Sys.time() - t0, units = "secs"), 1), " s")
}

# ---- save --------------------------------------------------------------------

dir.create(dirname(out), recursive = T, showWarnings = F)
write.csv(do.call(rbind, res), out, row.names = F)
write_provenance(
  out,
  script = here("scripts", "split_corrs.R"),
  opts = opts,
  inputs = c(
    config = opts$config,
    stimuli = data_path("stimuli"),
    onsets = data_path("onsets"),
    sapply(regions, function(r) data_path("mua", region = r))))
message("Saved ", out)
