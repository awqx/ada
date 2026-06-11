library(stringr)

#' Extract info from file names
#'
#' Gets the recorded region, age, and number of subjects from the
#' MATLAB file name.
#'
#' @param f str file name
#' @returns a list with values `region`, `age`, and `n_subjects`

file_info <- function(f) {
  region <- f %>%
    str_extract("^zz_[[:alnum:]]{2}") %>%
    str_remove("^zz_")
  age <- f %>%
    str_extract("Textures2_[[:digit:]]+") %>%
    str_remove("Textures2_") %>%
    as.integer()
  n_subjects <- f %>%
    str_extract("[[:digit:]]{1}subj\\.mat$") %>%
    str_remove("subj\\.mat") %>%
    as.integer()

  list(region = region, age = age, n_subjects = n_subjects)
}

#' Extract info from ID names
#'
#' Gets the recorded region, age, and number of subjects from the
#' shortened ID label used in processed data
#'
#' @param id str ID
#' @returns a list with values `region`, `age`, and `subj`

id_info <- function(id) {
  region <- str_extract(id, "^[[:alnum:]]{2}")
  age <- id %>%
    str_extract("[[:digit:]]{3}") %>%
    str_remove("_M") %>%
    as.integer()
  subj <- id %>%
    str_extract("_M[[:digit:]]{1}[[:lower:]]") %>%
    str_remove("_")

  list(region = region, age = age, subj = subj)
}
