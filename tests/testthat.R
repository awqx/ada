## Run all tests: Rscript tests/testthat.R
## (or, in an R session, testthat::test_dir("tests/testthat"))
here::i_am("tests/testthat.R")
testthat::test_dir(here::here("tests", "testthat"), stop_on_failure = TRUE)
