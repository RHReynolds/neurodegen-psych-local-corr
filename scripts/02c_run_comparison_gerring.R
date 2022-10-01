# Description: compare to results of Gerring 2022 (https://pubmed.ncbi.nlm.nih.gov/35525699/)

# Load packages -----------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

args <-
  list(
    url = "https://ars.els-cdn.com/content/image/1-s2.0-S000632232201071X-mmc2.zip",
    download_file_name = here::here("raw_data", "02_univar_bivar", "gerring2022_supp_materials.zip"),
    results = here::here("results", "02_univar_bivar_test", "results_summary.rds")
  )

# Download supplementary --------------------------------------------------
dir.create(here::here("raw_data", "02_univar_bivar"), showWarnings = T)

# check if file exists
if (!file.exists(args$download_file_name)) {

  download.file(
    url = args$url,
    destfile = args$download_file_name
  )

  unzip(
    args$download_file_name,
    files = "supplementary_table_3.xlsx",
    exdir = here::here("raw_data", "02_univar_bivar")
  )

}

# Load data ---------------------------------------------------------------

gerring <-
  readxl::read_xlsx(
    here::here("raw_data", "02_univar_bivar", "supplementary_table_3.xlsx"),
    skip = 1
    )

lava <-
  readRDS(args$results)

# Main --------------------------------------------------------------------

overlap <-
  lava$bivar %>%
  dplyr::filter(
    phen1 %in% c("BIP", "MDD", "SCZ"),
    phen2 %in% c("BIP", "MDD", "SCZ")
  ) %>%
  dplyr::inner_join(
    gerring %>%
      dplyr::rename(
        n_snps = n.snps,
        n_pcs = n.pcs
      ) %>%
      dplyr::mutate(
        phen1 =
          stringr::str_replace(phen1, "DEP", "MDD"),
        phen2 =
          stringr::str_replace(phen2, "DEP", "MDD")
      ),
    by = c("chr", "start", "stop"),
    suffix = c(".rhr", ".gerring")
  ) %>%
  dplyr::mutate(
    overlapping_phen =
      case_when(
        phen1.rhr == phen1.gerring & phen2.rhr == phen2.gerring ~ 2,
        phen1.rhr == phen1.gerring & phen2.rhr != phen2.gerring |
          phen1.rhr != phen1.gerring & phen2.rhr == phen2.gerring ~ 1,
        TRUE ~ 0
      )
  )

# sort column names, so that they alternate by suffix
cols_sorted <-
  overlap %>%
  dplyr::select(-chr, -start, -stop, -fdr, -overlapping_phen) %>%
  colnames() %>%
  sort()

overlap <-
  overlap %>%
  dplyr::select(chr, start, stop, overlapping_phen, cols_sorted, -fdr)

# Save data ---------------------------------------------------------------
out_dir <- here::here("results", "02_univar_bivar_test")
dir.create(out_dir, showWarnings = T)
saveRDS(
  overlap,
  file = file.path(out_dir, "results_overlap_with_gerring.rds")
)
