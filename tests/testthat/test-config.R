test_that("data_path fills placeholders and checks names", {
  expect_match(data_path("mua", region = "IT"), "IT\\.qs2$")
  expect_equal(data_path("stimuli"), here::here(config::get("stimuli", file = here::here("config", "paths.yml"))))
  expect_error(data_path("mua"), "Unfilled placeholder")
  expect_error(data_path("not_a_source"), "Unknown data source")
})

test_that("data_path respects config sections and absolute paths", {
  f <- withr::local_tempfile(fileext = ".yml")
  writeLines(c(
    "default:",
    "  mua: data/{region}.qs2",
    "  stimuli: data/stimuli.csv",
    "server:",
    "  mua: /scratch/ada/{region}.qs2"), f)

  # default section, relative to the project root
  expect_equal(data_path("mua", region = "V1", file = f), here::here("data", "V1.qs2"))

  # section chosen explicitly or via R_CONFIG_ACTIVE; absolute paths kept
  expect_equal(data_path("mua", region = "V1", file = f, config = "server"), "/scratch/ada/V1.qs2")
  withr::local_envvar(R_CONFIG_ACTIVE = "server")
  expect_equal(data_path("mua", region = "V1", file = f), "/scratch/ada/V1.qs2")

  # keys missing from a section are inherited from default
  expect_equal(data_path("stimuli", file = f), here::here("data", "stimuli.csv"))
})
