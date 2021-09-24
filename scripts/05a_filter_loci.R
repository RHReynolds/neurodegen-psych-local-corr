# Description: get LD blocks with significant bivariate corr between AD/PD

# Load packages -----------------------------------------------------------

library(GenomicRanges)
library(stringr)
library(tidyverse)

# Set arguments -----------------------------------------------------------

args <-
  list(
    phenotypes = c(c("AD2019", "PD2019.meta5.ex23andMe", "LBD2020")),
    loc_file = here::here("results", "01_input_prep", "gwas_filtered.loci"),
    bivar_results = here::here("results", "02_univar_bivar_test"),
    # Bivariate threshold from previous analyses (02_run_univar_bivar_test.rmd)
    # Used to filter bivariate results to include only phenos with significant local rg
    # These phenos will then be used together with the eQTL
    bivar_threshold = 0.05/1603, # n bivariate tests,
    out_dir = here::here("results", "05_noproxy_univar_bivar")
  )

# Load data ---------------------------------------------------------------

loci <- LAVA::read.loci(args$loc_file)

# Load bivariate results
results <-
  list.files(
    path = args$bivar_results,
    pattern = "bivar.lava.rds",
    full.names = T
  ) %>%
  lapply(., function(file){

    list <-
      file %>%
      readRDS()

    list %>%
      purrr::discard(is.null) %>%
      qdapTools::list_df2df() %>%
      dplyr::select(-X1)

  }) %>%
  qdapTools::list_df2df(col1 = "list_name")

# Main --------------------------------------------------------------------

# Filter bivariate results to include only loci with significant local rg of phenotypes
locus_of_interest <-
  results %>%
  dplyr::filter(
    p < args$bivar_threshold
  ) %>%
  dplyr::select(locus, contains("phen")) %>%
  tidyr::pivot_longer(
    cols = !locus
  ) %>%
  dplyr::select(-name, phenotype = value) %>%
  dplyr::distinct(locus, phenotype) %>%
  dplyr::mutate(
    phenotype_of_interest =
      case_when(phenotype %in% args$phenotypes ~ TRUE,
                TRUE ~ FALSE)
  ) %>%
  dplyr::group_by(locus) %>%
  dplyr::summarise(
    true_sum = sum(phenotype_of_interest)
  ) %>%
  dplyr::filter(
    true_sum == length(args$phenotypes)
  ) %>%
  .[["locus"]]

filtered_loci <-
  loci %>%
  dplyr::filter(
    LOC %in% locus_of_interest
  )

# Save data ---------------------------------------------------------------

dir.create(args$out_dir, showWarnings = T)
write_delim(
  filtered_loci,
  delim = "\t",
  file = file.path(args$out_dir, "filtered.loci")
)

