# R/sparsity_debiased.R

test_that("rep_power_mean is unbiased for mu^p under Poisson noise", {
  set.seed(11)
  mu <- 3
  for (K in c(4, 6, 10)) {
    Y <- matrix(rpois(2e5 * K, mu), ncol = K)
    for (p in 1:4) {
      expect_equal(mean(rep_power_mean(Y, p)), mu^p, tolerance = 0.02,
                   label = paste0("K = ", K, ", p = ", p))
    }
  }
  # plug-in is biased: E[xbar^2] = mu^2 + mu / K
  Y <- matrix(rpois(2e5 * 4, mu), ncol = 4)
  expect_equal(mean(rowMeans(Y)^2), mu^2 + mu / 4, tolerance = 0.01)
  expect_error(rep_power_mean(Y[, 1:3], 4), "at least 4")
})

test_that("rep_power_mean matches brute-force products of distinct reps", {
  Y <- matrix(c(1, 4, 2, 0, 3, 5, 2, 2, 1, 7), nrow = 2)  # 2 stimuli x 5 reps
  brute <- function(y, p) mean(apply(combn(length(y), p), 2, function(i) prod(y[i])))
  for (p in 2:4) expect_equal(rep_power_mean(Y, p), apply(Y, 1, brute, p = p))
})

test_that("moment measures agree with existing sparsity functions on true vectors", {
  set.seed(12)
  mu <- rgamma(50, shape = 0.7, rate = 0.2)
  s <- moment_sparsity(tuning_moments(mu), length(mu))
  expect_equal(s[["tr"]], rm_sparsity(mu))
  expect_equal(s[["hoyer"]], Lp_sparsity(mu))
  expect_equal(s[["activity"]], mean(mu)^2 / mean(mu^2))
  expect_equal(s[["kurtosis"]], mean((mu - mean(mu))^4) / mean((mu - mean(mu))^2)^2 - 3)
})

test_that("gini bounds", {
  expect_equal(gini(rep(2, 10)), 0)
  expect_equal(gini(c(5, rep(0, 9))), 1 - 1 / 10)
})

test_that("debiased activity ratio removes trial-noise bias; plug-in does not", {
  set.seed(13)
  mu <- 2 * rgamma(200, shape = 1)  # fixed true tuning, a ~ 0.5
  truth <- sparsity_truth(mu)
  est <- t(replicate(400, sparsity_estimates(
    matrix(rpois(200 * 6, mu), ncol = 6))))
  bias_plugin <- mean(est[, "activity_plugin"]) - truth[["activity_truth"]]
  bias_debiased <- mean(est[, "activity_debiased"]) - truth[["activity_truth"]]
  # Poisson noise adds m1 / K to the second moment: E[plug-in] ~ m1^2 / (m2 + m1 / K)
  m <- tuning_moments(mu)
  predicted <- m[["m1sq"]] / (m[["m2"]] + m[["m1"]] / 6) - m[["m1sq"]] / m[["m2"]]
  expect_equal(bias_plugin, predicted, tolerance = 0.1)
  expect_lt(bias_plugin, -0.015)
  expect_lt(abs(bias_debiased), 0.005)
  expect_lt(abs(mean(est[, "kurtosis_debiased"]) - truth[["kurtosis_truth"]]),
            abs(mean(est[, "kurtosis_plugin"]) - truth[["kurtosis_truth"]]))
})

test_that("split_half_reps matches split_half_corr on the same data", {
  set.seed(14)
  Y <- matrix(rnorm(300 * 6, mean = rnorm(300)), ncol = 6)
  labels <- rep(seq_len(300), 6)  # column-major: block of reps
  expect_equal(split_half_reps(Y), split_half_corr(as.vector(Y), labels))
  expect_equal(split_half_reps(Y, reliability = T),
               split_half_corr(as.vector(Y), labels, reliability = T))
})
