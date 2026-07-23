nb_id <- "20260715"
here::i_am("scripts/pf_lda_avgs.R")

library(dplyr)
library(tidyr)
library(here)
library(qs2)
library(cli)

# 0 Setup #####################################################

set.seed(20260715)
args <- commandArgs(trailingOnly = TRUE)
n_sites <- ifelse(length(args) == 0, 30, args[1])

# 0.1 Directories =============================================

pf_dir <- here("data", "processed", nb_id, "pairwise_family", "lda")
res_dir <- here("results", nb_id)

if (!dir.exists(res_dir)) dir.create(res_dir, recursive = T)

# 0.2 Functions ===============================================

avg_folds <- function(dprimes) {
  dprimes %>%
    group_by(data_split) %>%
    summarize(across(c("mu", "sd", "dprime"), mean))
}

avg_file <- function(f_nm, dir_nm = ".") {
  group_ids <- f_nm %>%
    basename() %>%
    tools::file_path_sans_ext() %>%
    strsplit("vs") %>%
    unlist()

  file.path(dir_nm, f_nm) %>%
    read.csv() %>%
    avg_folds() %>%
    mutate(grp1 = group_ids[1], grp2 = group_ids[2])
}

avg_region_dir <- function(region, dir_nm = "data", n_sites) {
  dprime_f <- paste("dprime", n_sites, sep = "_")
  dprime_dir <- here(dir_nm, dprime_f, region)
  f_nms <- list.files(dprime_dir)

  lapply(f_nms, avg_file, dir_nm = dprime_dir) %>%
    do.call(rbind, .) %>%
    mutate(region = region)
}

# 1 Data loading ##############################################

regions <- c("V1", "V2", "V4", "IT")

dprime_avgs <- lapply(
  regions,
  avg_region_dir,
  dir_nm = pf_dir,
  n_sites = n_sites) %>%
  do.call(rbind, .) %>%
  mutate(n_sites = n_sites)

res_f <- paste0(paste("pf_lda_avgs", n_sites, sep = "_"), ".csv")

write.csv(dprime_avgs, here(res_dir, res_f), row.names = F)
