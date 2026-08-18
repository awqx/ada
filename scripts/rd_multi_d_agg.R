nb_id <- "20260723"
here::i_am("scripts/rd_multi_d_agg.R")

library(dplyr)
library(here)
library(cli)

# 0 Setup #####################################################

# 0.1 Directories =============================================

rd_dir <- here("data", "processed", nb_id, "rand_dichot")
if (!dir.exists(rd_dir)) dir.create(rd_dir, recursive = T)
data_dir <- here("data", "processed", "20260625")
res_dir <- here("results", nb_id)

# 0.2 Functions ===============================================

stimuli <- read.csv(here(data_dir, "stimuli.csv")) %>%
  filter(coherence > 0)
families <- unique(stimuli$family)
agg_rd <- function(rd_dir, fam_ids, n_sites, regions, shrink) {
  dir_nms <- list.files(rd_dir)
  f_nms <- expand.grid(fam_ids, shrink, n_sites, regions)
  lapply(
    1:nrow(f_nms),
    function(i) {
      x <- f_nms[i, ]
      md_dir <- ifelse(x[[2]], "multi_d_shrink", "multi_d")
      f_nm <- here(
        rd_dir, paste0("fam", sprintf("%03d", x[[1]])),
        md_dir, paste0("n", x[[3]]), x[[4]], "dprime.csv")
      if (!file.exists(f_nm)) return()
      f_nm %>%
        read.csv() %>%
        dplyr::select(-grps) %>%
        mutate(
          family = x[[1]],
          shrink = x[[2]],
          n_sites = x[[3]],
          region = x[[4]])}) %>%
    do.call(rbind, .)
}

# 1 Data loading ##############################################

regions <- c("V1", "V2", "V4", "IT")
n_sites <- 10 * 1:4
shrink <- c(T, F)
fam_ids <- read.csv(here(data_dir, "stimuli.csv")) %>%
  filter(coherence > 0) %>%
  .$family %>%
  unique()

# 2 Aggregation over loop #####################################

res <- agg_rd(rd_dir, fam_ids, n_sites, regions, shrink)
write.csv(res, here(res_dir, "rd_infam.csv"), row.names = F)
