# Description: tidy local genetic correlation results

# Load packages -----------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

args <-
  list(
    results_dir = here::here("results", "02_univar_bivar_test")
  )

# Load data ---------------------------------------------------------------

files <-
  list.files(
    path =
      args$results_dir,
    pattern = ".lava.rds",
    full.names = T
  )

results <-
  setNames(
    object = files %>%
      lapply(., function(file){

        list <-
          file %>%
          readRDS()

        list %>%
          purrr::discard(is.null) %>%
          qdapTools::list_df2df() %>%
          dplyr::select(-X1)

      }),
    nm = files %>%
      basename() %>%
      str_remove(".lava.rds") %>%
      str_remove(".*:") %>%
      str_remove(".*\\.")
  )

# Main --------------------------------------------------------------------

# Tidy phenotypes
results$univ <-
  results$univ %>%
  dplyr::mutate(
    phen = phen %>%
      str_replace_all("[:digit:]", "") %>%
      str_remove("\\..*") %>%
      fct_relevel(
        fct_disease
      )
  )

results$bivar <-
  results$bivar %>%
  dplyr::mutate(
    phen1 = phen1 %>%
      str_replace_all("[:digit:]", "") %>%
      str_remove("\\..*"),
    phen2 = phen2 %>%
      str_replace_all("[:digit:]", "") %>%
      str_remove("\\..*"),
    fdr = p.adjust(p, method = "fdr")
  )


# Save data ---------------------------------------------------------------
out_dir <- here::here("results", "02_univar_bivar_test")
dir.create(out_dir, showWarnings = T)
saveRDS(
  results,
  file = file.path(out_dir, "results_summary.rds")
)
