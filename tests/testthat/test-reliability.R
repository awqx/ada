# Simulations with known answers for R/reliability.R

test_that("balanced_splits gives each 3-vs-3 partition of 6 reps once", {
  s <- balanced_splits(6)
  expect_equal(ncol(s[[1]]), choose(6, 3) / 2)
  expect_true(all(apply(rbind(s[[1]], s[[2]]), 2, setequal, 1:6)))
  expect_false(anyDuplicated(apply(s[[1]], 2, paste, collapse = "")) > 0)
  # no split is the mirror of another
  mirrored <- apply(s[[2]], 2, function(x)
    any(apply(s[[1]], 2, setequal, x)))
  expect_false(any(mirrored))
})

test_that("balanced_splits handles odd repetitions", {
  s <- balanced_splits(5)
  expect_equal(ncol(s[[1]]), choose(5, 2))
  expect_equal(nrow(s[[1]]), 2)
  expect_equal(nrow(s[[2]]), 3)
  expect_error(balanced_splits(1))
})

test_that("rep_index uses order of appearance", {
  expect_identical(rep_index(rep(c("a", "b", "c"), times = 4)), rep(1:4, each = 3))
})

test_that("rep_index splits a shared label into consecutive chunks", {
  # one block per repetition: two images, then two blanks sharing a label
  lab <- rep(c("a", "b", "blank", "blank"), times = 3)
  expect_identical(rep_index(lab), rep(1:3, each = 4))
})

test_that("rep_index with explicit n_reps handles labels repeated within a block", {
  # family-style labels: each family appears twice per repetition, 3 reps
  lab <- rep(c("f1", "f1", "f2", "f2", "blank"), times = 3)
  expect_identical(rep_index(lab, n_reps = 3), rep(1:3, each = 5))
  # the default (smallest count) would wrongly infer 3 reps from "blank" here,
  # but 6 from a set without it
  expect_identical(max(rep_index(lab[lab != "blank"])), 6L)
})

test_that("split_half_corr with family-style labels matches pre-averaged images", {
  set.seed(3)
  n_fam <- 300
  per_fam <- 2
  fam_signal <- rnorm(n_fam)
  # one block per rep; each family shown as 2 images per block
  labels <- rep(rep(seq_len(n_fam), each = per_fam), 6)
  trials <- fam_signal[labels] + rnorm(length(labels))
  by_family <- split_half_corr(trials, labels, n_reps = 6)
  # equivalent: average the 2 images within each block first
  block_mean <- colMeans(matrix(trials, nrow = per_fam))
  expect_equal(by_family, split_half_corr(block_mean, rep(seq_len(n_fam), 6)))
})

test_that("rep_index errors when counts are not multiples of n_reps", {
  expect_error(rep_index(c("a", "a", "b", "b", "b")), "multiple")
})

# response = signal (var 1) + noise (var 1) per trial, so a mean over k reps
# has reliability k / (k + 1): halves of 6 reps (k = 3) correlate at 3/4, and
# Spearman-Brown to 6 reps gives 6/7
simulate_trials <- function(n_stim = 2000, n_reps = 6, n_sites = 3, signal_sd = 1) {
  signal <- matrix(rnorm(n_stim * n_sites, sd = signal_sd), n_stim)
  list(
    trials = signal[rep(seq_len(n_stim), n_reps), , drop = F] +
      matrix(rnorm(n_stim * n_reps * n_sites), ncol = n_sites),
    labels = rep(seq_len(n_stim), n_reps)) # presorted, one block per rep
}

test_that("split_half_corr recovers known reliability", {
  set.seed(20261104)
  sim <- simulate_trials()
  expect_equal(split_half_corr(sim$trials, sim$labels), rep(3 / 4, 3), tolerance = 0.02)
  expect_equal(
    split_half_corr(sim$trials, sim$labels, reliability = T),
    rep(6 / 7, 3), tolerance = 0.02)
})

test_that("split_half_corr is near 0 for pure noise", {
  set.seed(1)
  sim <- simulate_trials(signal_sd = 0)
  expect_true(all(abs(split_half_corr(sim$trials, sim$labels)) < 0.05))
})

test_that("split_half_corr input handling", {
  set.seed(2)
  sim <- simulate_trials(n_stim = 200)
  r <- split_half_corr(sim$trials, sim$labels)

  # single site vector matches the matrix column
  expect_equal(split_half_corr(sim$trials[, 2], sim$labels), r[2])

  # NA labels are excluded
  keep <- sim$labels <= 100
  expect_equal(
    split_half_corr(sim$trials, replace(sim$labels, !keep, NA)),
    split_half_corr(sim$trials[keep, ], sim$labels[keep]))

  # order within a repetition block does not matter
  perm <- unlist(lapply(0:5, function(b) b * 200 + sample(200)))
  expect_equal(split_half_corr(sim$trials[perm, ], sim$labels[perm]), r)

  # a single split given as two index vectors
  one <- split_half_corr(sim$trials, sim$labels, splits = list(c(1, 2, 3), c(4, 5, 6)))
  expect_length(one, 3)
})

test_that("spearman_brown", {
  expect_equal(spearman_brown(0.5), 2 / 3)
  expect_equal(spearman_brown(0.3, scaling = 1), 0.3)
})

# nested design: a families x b images x n reps, one trial per image and rep,
# presorted with one block per repetition
simulate_nested <- function(a = 35, b = 15, n = 6, n_sites = 4,
                            s2_f = 1, s2_i = 1, s2_r = 0.5, s2_fr = 0.2, s2_e = 4) {
  grid <- expand.grid(image = seq_len(b), family = seq_len(a), rep = seq_len(n))
  y <- sapply(seq_len(n_sites), function(s) {
    fam <- rnorm(a, sd = sqrt(s2_f))
    img <- matrix(rnorm(a * b, sd = sqrt(s2_i)), b, a)
    rep_eff <- rnorm(n, sd = sqrt(s2_r))
    fr_eff <- matrix(rnorm(a * n, sd = sqrt(s2_fr)), a, n)
    fam[grid$family] + img[cbind(grid$image, grid$family)] + rep_eff[grid$rep] +
      fr_eff[cbind(grid$family, grid$rep)] + rnorm(nrow(grid), sd = sqrt(s2_e))
  })
  list(y = y, grid = grid)
}

test_that("nested_anova recovers variance components and rho", {
  set.seed(4)
  sim <- simulate_nested(n_sites = 300, s2_f = 1, s2_i = 3)
  est <- nested_anova(sim$y, sim$grid$family, sim$grid$image, sim$grid$rep)
  expect_equal(nrow(est), 300)
  # means over sites are unbiased for the true components
  expect_equal(mean(est$s2_f), 1, tolerance = 0.1)
  expect_equal(mean(est$s2_i), 3, tolerance = 0.05)
  expect_equal(mean(est$s2_fr), 0.2, tolerance = 0.25)
  expect_equal(mean(est$s2_e), 4, tolerance = 0.02)
  expect_equal(median(est$rho), 1 / 4, tolerance = 0.15)
})

test_that("nested_anova matches lm() sums of squares on one site", {
  set.seed(5)
  sim <- simulate_nested(a = 6, b = 4, n = 3, n_sites = 1)
  d <- data.frame(y = sim$y[, 1], f = factor(sim$grid$family),
                  i = factor(paste(sim$grid$family, sim$grid$image)), r = factor(sim$grid$rep))
  ms <- anova(lm(y ~ f + r + f:r + f/i, data = d))[["Mean Sq"]]
  names(ms) <- rownames(anova(lm(y ~ f + r + f:r + f/i, data = d)))
  est <- nested_anova(d$y, d$f, d$i, d$r)
  expect_equal(est$s2_e, ms[["Residuals"]])
  expect_equal(est$s2_i, (ms[["f:i"]] - ms[["Residuals"]]) / 3)
  expect_equal(est$s2_f, (ms[["f"]] - ms[["f:i"]] - ms[["f:r"]] + ms[["Residuals"]]) / (4 * 3))
})

test_that("nested_anova rejects unbalanced designs", {
  sim <- simulate_nested(a = 4, b = 3, n = 2, n_sites = 1)
  expect_error(nested_anova(sim$y[-1], sim$grid$family[-1], sim$grid$image[-1], sim$grid$rep[-1]), "balanced")
})

test_that("nested_anova simple model matches lm() with family/image only", {
  set.seed(6)
  sim <- simulate_nested(a = 6, b = 4, n = 3, n_sites = 1)
  d <- data.frame(y = sim$y[, 1], f = factor(sim$grid$family),
                  i = factor(paste(sim$grid$family, sim$grid$image)), r = factor(sim$grid$rep))
  tab <- anova(lm(y ~ f/i, data = d))
  ms <- setNames(tab[["Mean Sq"]], rownames(tab))
  est <- nested_anova(d$y, d$f, d$i, d$r)
  expect_equal(est$s2_e_simple, ms[["Residuals"]])
  expect_equal(est$s2_i_simple, (ms[["f:i"]] - ms[["Residuals"]]) / 3)
  expect_equal(est$s2_f_simple, (ms[["f"]] - ms[["f:i"]]) / (4 * 3))
  # rep main effect mean square from the full model
  full <- anova(lm(y ~ f + r + f:r + f/i, data = d))
  expect_equal(est$ms_r, full["r", "Mean Sq"])
})

test_that("simple and full models agree without repetition effects", {
  set.seed(7)
  sim <- simulate_nested(n_sites = 300, s2_f = 1, s2_i = 1, s2_r = 0, s2_fr = 0)
  est <- nested_anova(sim$y, sim$grid$family, sim$grid$image, sim$grid$rep)
  expect_equal(mean(est$s2_f_simple), mean(est$s2_f), tolerance = 0.05)
  expect_equal(mean(est$s2_i_simple), mean(est$s2_i), tolerance = 0.02)
  expect_equal(mean(est$s2_r), 0, tolerance = 0.02)
})

test_that("simple model is biased in the predicted direction with repetition effects", {
  set.seed(8)
  n <- 6
  sim <- simulate_nested(n_sites = 400, s2_f = 1, s2_i = 1, s2_r = 2, s2_fr = 1.5, s2_e = 2)
  est <- nested_anova(sim$y, sim$grid$family, sim$grid$image, sim$grid$rep)
  # E[s2_f_simple] = s2_f + s2_fr / n; E[s2_i_simple] = s2_i - (s2_r + s2_fr) / n
  expect_equal(mean(est$s2_f_simple), 1 + 1.5 / n, tolerance = 0.1)
  expect_equal(mean(est$s2_i_simple), 1 - (2 + 1.5) / n, tolerance = 0.1)
  # full model stays unbiased
  expect_equal(mean(est$s2_f), 1, tolerance = 0.1)
  expect_equal(mean(est$s2_i), 1, tolerance = 0.05)
  expect_gt(median(est$rho_simple), median(est$rho))
})

test_that("nested_anova repetition-main-effect model matches lm()", {
  set.seed(9)
  sim <- simulate_nested(a = 6, b = 4, n = 3, n_sites = 1)
  d <- data.frame(y = sim$y[, 1], f = factor(sim$grid$family),
                  i = factor(paste(sim$grid$family, sim$grid$image)), r = factor(sim$grid$rep))
  tab <- anova(lm(y ~ r + f/i, data = d))
  ms <- setNames(tab[["Mean Sq"]], rownames(tab))
  est <- nested_anova(d$y, d$f, d$i, d$r)
  expect_equal(est$s2_e_rep, ms[["Residuals"]])
  expect_equal(est$s2_i_rep, (ms[["f:i"]] - ms[["Residuals"]]) / 3)
  expect_equal(est$s2_f_rep, (ms[["f"]] - ms[["f:i"]]) / (4 * 3))
})

test_that("repetition-main-effect model is unbiased with block offsets only", {
  set.seed(10)
  sim <- simulate_nested(n_sites = 400, s2_f = 1, s2_i = 1, s2_r = 2, s2_fr = 0, s2_e = 2)
  est <- nested_anova(sim$y, sim$grid$family, sim$grid$image, sim$grid$rep)
  expect_equal(mean(est$s2_f_rep), 1, tolerance = 0.1)
  expect_equal(mean(est$s2_i_rep), 1, tolerance = 0.05)
  # the simple model is not: its image component shrinks by s2_r / n
  expect_equal(mean(est$s2_i_simple), 1 - 2 / 6, tolerance = 0.1)
})
