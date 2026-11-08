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
