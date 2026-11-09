## Simulation: sparsity inflation at low reliability, plug-in vs. debiased
##
## For a single site, true mean counts mu_i (stimuli i = 1..N) are drawn from
## a tuning shape with mean count m and squared coefficient of variation CV^2
## across stimuli (true activity ratio a = m1^2 / m2 = 1 / (1 + CV^2)), and K
## repetitions of Poisson counts are drawn per stimulus. The reliability of the
## K-repetition mean is
##
##   R = CV^2 m / (CV^2 m + F / K),
##
## with Fano factor F (1 for Poisson), so varying m at fixed CV^2 traces bias
## vs. reliability at fixed true sparsity.
##
## Each replication records, per sparsity measure, the truth (computed on the
## realized mu), the plug-in estimate (on trial averages), and the debiased
## estimate (R/sparsity_debiased.R), plus the estimated split-half
## reliability. Settings: config/sim_sparsity.yml. Background:
## notebooks/20261104.
##
## Usage:
##   Rscript scripts/sim_sparsity.R --dry_run                 # design only
##   Rscript scripts/sim_sparsity.R --replications 5 --n_conditions 10 --out x
##   Rscript scripts/sim_sparsity.R                           # full run

here::i_am("scripts/sim_sparsity.R")
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(qs2)
  library(SimDesign)
})
source(here("R", "config.R"))

# ---- options -----------------------------------------------------------------

opts <- optparse::parse_args(optparse::OptionParser(
  usage = "Rscript scripts/sim_sparsity.R [options]",
  option_list = list(
    optparse::make_option(
      "--config", default = here("config", "sim_sparsity.yml"),
      help = "settings file [default: %default]"),
    optparse::make_option(
      "--replications", type = "integer", default = NULL,
      help = "replications per condition [default: from config]"),
    optparse::make_option(
      "--n_conditions", type = "integer", default = 0L,
      help = "if > 0, run only n conditions spread evenly across the design, for testing [default: %default]"),
    optparse::make_option(
      "--out", default = NULL,
      help = "output path without extension [default: `out` in config]"),
    optparse::make_option(
      "--dry_run", action = "store_true", default = F,
      help = "print the design and exit"))))
opts$help <- NULL

cfg <- config::get(file = opts$config)
if (is.null(opts$replications)) opts$replications <- cfg$replications
if (is.null(opts$out)) opts$out <- resolve_path(cfg$out)

# ---- tuning and noise --------------------------------------------------------

#' Spike-and-slab elevation: with a fraction p of stimuli at rate c * b and
#' the rest at b, solve c > 1 so that the squared CV across stimuli is cv2
spike_slab_elevation <- function(p, cv2) {
  q <- 1 - p
  A <- (1 + cv2) * p^2 - p
  B <- 2 * (1 + cv2) * q * p
  C <- q * ((1 + cv2) * q - 1)
  if (A >= 0) stop("CV^2 = ", cv2, " is not reachable with p = ", p, ".")
  max((-B + c(-1, 1) * sqrt(B^2 - 4 * A * C)) / (2 * A))
}

#' True mean counts per stimulus with mean m and squared CV cv2
draw_tuning <- function(shape, N, m, cv2, p_active = 0.1) {
  switch(
    shape,
    constant = rep(m, N),
    lognormal = {
      s2 <- log1p(cv2)
      m * exp(sqrt(s2) * rnorm(N) - s2 / 2)
    },
    gamma = rgamma(N, shape = 1 / cv2, scale = m * cv2),
    spike_slab = {
      n_active <- max(1, round(p_active * N))
      p <- n_active / N
      elev <- spike_slab_elevation(p, cv2)
      b <- m / (1 - p + p * elev)
      sample(c(rep(b * elev, n_active), rep(b, N - n_active)))
    },
    stop("Unknown tuning shape '", shape, "'."))
}

#' Counts for stimuli x repetitions given true means; `block_sd` adds a
#' multiplicative gain per repetition (lognormal, mean 1)
draw_counts <- function(mu, K, noise = "poisson", block_sd = 0) {
  gain <- if (block_sd > 0) exp(block_sd * rnorm(K) - block_sd^2 / 2) else rep(1, K)
  rate <- outer(mu, gain)
  switch(
    noise,
    poisson = matrix(rpois(length(rate), rate), nrow = length(mu)),
    stop("Noise model '", noise, "' is not implemented."))
}

# ---- design ------------------------------------------------------------------

fano <- switch(cfg$noise, poisson = 1)

# note: dplyr verbs (mutate, filter, bind_rows) silently drop rows of SimDesign
# `Design` objects, so only base-R assignment, which(), and rbind() are used
signal <- createDesign(
  shape = unlist(cfg$shapes),
  cv2 = unlist(cfg$cv2),
  m = unlist(cfg$rates),
  N = unlist(cfg$n_stimuli),
  K = unlist(cfg$n_reps))
null <- createDesign(
  shape = "constant",
  cv2 = 0,
  m = unlist(cfg$rates),
  N = unlist(cfg$n_stimuli),
  K = unlist(cfg$n_reps))
Design <- rbind(signal, null)
Design$a <- 1 / (1 + Design$cv2)
Design$R_target <- Design$cv2 * Design$m / (Design$cv2 * Design$m + fano / Design$K)
if (opts$n_conditions > 0)  # evenly spaced rows, so a smoke test spans the design
  Design <- Design[unique(round(seq(1, nrow(Design), length.out = opts$n_conditions))), ]

message(
  "Conditions: ", nrow(Design), " (", sum(Design$shape == "constant"), " null)\n",
  "Replications: ", opts$replications, " per condition\n",
  "Output: ", opts$out, ".{qs2,csv}")
if (isTRUE(opts$dry_run)) {
  print(as.data.frame(Design) %>% count(shape, N, K) %>%
          tidyr::pivot_wider(names_from = N, values_from = n, names_prefix = "N = "), n = Inf)
  message("Target reliability by CV^2 (rows) and mean count (columns), K = 6:")
  print(round(xtabs(R_target ~ cv2 + m, data = as.data.frame(Design)[
    Design$K == 6 & Design$N == min(Design$N) & Design$shape == Design$shape[1], ]), 3))
  quit(save = "no")
}

# ---- simulation --------------------------------------------------------------

# estimator functions in their own environment so parallel workers get them
# through fixed_objects
fns <- new.env()
sys.source(here("R", "reliability.R"), envir = fns)
sys.source(here("R", "sparsity_debiased.R"), envir = fns)
fns$draw_tuning <- draw_tuning
fns$draw_counts <- draw_counts
fns$spike_slab_elevation <- spike_slab_elevation
environment(fns$draw_tuning) <- fns
environment(fns$draw_counts) <- fns
environment(fns$spike_slab_elevation) <- fns

# balanced splits per K, capped at max_splits (fixed random subset)
split_subset <- function(K, max_splits, seed) {
  sp <- fns$balanced_splits(K)
  if (ncol(sp[[1]]) <= max_splits) return(sp)
  set.seed(seed + K)
  keep <- sort(sample(ncol(sp[[1]]), max_splits))
  lapply(sp, function(x) x[, keep, drop = F])
}

fixed <- list(
  fns = fns,
  splits = lapply(setNames(nm = unique(Design$K)),
                  split_subset, max_splits = cfg$max_splits, seed = cfg$seed),
  noise = cfg$noise,
  block_sd = cfg$block_sd,
  p_active = cfg$spike_slab_p)

Generate <- function(condition, fixed_objects) {
  f <- fixed_objects$fns
  mu <- f$draw_tuning(condition$shape, condition$N, condition$m, condition$cv2,
                      fixed_objects$p_active)
  list(mu = mu, Y = f$draw_counts(mu, condition$K, fixed_objects$noise,
                                  fixed_objects$block_sd))
}

Analyse <- function(condition, dat, fixed_objects) {
  f <- fixed_objects$fns
  c(f$sparsity_truth(dat$mu),
    f$sparsity_estimates(dat$Y),
    reliability = f$split_half_reps(
      dat$Y, fixed_objects$splits[[as.character(condition$K)]], reliability = T))
}

Summarise <- function(condition, results, fixed_objects) {
  measures <- c("activity", "tr", "hoyer", "kurtosis", "gini")
  out <- c(
    reliability_mean = mean(results[, "reliability"], na.rm = T),
    reliability_sd = sd(results[, "reliability"], na.rm = T))
  for (meas in measures) {
    truth <- results[, paste0(meas, "_truth")]
    out[paste0(meas, "_truth")] <- mean(truth, na.rm = T)
    for (version in c("plugin", "debiased")) {
      col <- paste0(meas, "_", version)
      if (!col %in% colnames(results)) next
      err <- results[, col] - truth
      out[paste0(col, "_bias")] <- mean(err, na.rm = T)
      out[paste0(col, "_sd")] <- sd(results[, col], na.rm = T)
      out[paste0(col, "_rmse")] <- sqrt(mean(err^2, na.rm = T))
      out[paste0(col, "_na")] <- mean(is.na(results[, col]))
    }
  }
  out
}

res <- runSimulation(
  design = Design,
  replications = opts$replications,
  generate = Generate,
  analyse = Analyse,
  summarise = Summarise,
  fixed_objects = fixed,
  seed = genSeeds(Design, iseed = cfg$seed),
  parallel = isTRUE(cfg$parallel),
  ncores = cfg$ncores,
  store_results = isTRUE(cfg$store_results),
  save = F,
  control = list(allow_na = T, allow_nan = T))

# ---- save --------------------------------------------------------------------

dir.create(dirname(opts$out), recursive = T, showWarnings = F)
qs_save(res, paste0(opts$out, ".qs2"))
write.csv(as.data.frame(res), paste0(opts$out, ".csv"), row.names = F)
write_provenance(
  paste0(opts$out, ".qs2"),
  script = here("scripts", "sim_sparsity.R"),
  opts = c(opts, list(settings = unclass(cfg))),
  inputs = c(config = opts$config))
message("Saved ", opts$out, ".{qs2,csv}")
