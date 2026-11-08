## Per-site split-half correlations and Spearman-Brown reliabilities
##
## Computes every combination of region x onset x window x blank handling in
## a settings file (default config/split_corrs.yml) and writes one long table,
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

# stimulus labels per trial for each blank scheme (NA = excluded). Trials are
# presorted into one block per repetition, so labels alone determine the
# repetitions (see rep_index()).
stimuli <- read.csv(data_path("stimuli"))
blank <- stimuli$coherence == 0
img <- paste(stimuli$coherence, stimuli$family, stimuli$id)
labels <- list(
  none = ifelse(blank, NA, img),
  aggregated = ifelse(blank, "blank", img),
  separated = img)[unlist(cfg$blanks)]

splits <- balanced_splits(max(stimuli$repetition))

message(
  "Regions: ", paste(regions, collapse = ", "), "\n",
  "Onsets: ", paste(names(cfg$onsets), collapse = ", "),
  " | windows: ", paste(cfg$windows, collapse = ", "), " bins",
  " | blanks: ", paste(names(labels), collapse = ", "), "\n",
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

    for (blanks in names(labels)) {
      corr <- split_half_corr(win_mean, labels[[blanks]], splits)
      res[[length(res) + 1]] <- data.frame(
        region = region,
        site = seq_along(corr),
        onset_type = grid$onset_type[i],
        onset = onset,
        ms = 10 * grid$window[i],
        blanks = blanks,
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
