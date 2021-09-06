# Description: run multiple regressions for GWAS traits

# Load packages -----------------------------------------------------------

library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

phenotypes <- c("AD2019", "LBD2020", "PD2019.meta5.ex23andMe", "BIP2021", "MDD2019", "SCZ2018") %>% sort()

args <-
  list(
    ref_prefix = "/tools/MAGMA/reference_files/g1000_eur/g1000_eur",
    loc_file = here::here("results", "01_input_prep", "gwas_filtered.loci"),
    info_file = here::here("results", "01_input_prep", "input.info.txt"),
    sample_overlap_file = here::here("results", "01_input_prep", "sample_overlap.txt"),
    phenotypes = phenotypes,
    output_filename = str_c(phenotypes, collapse = ":")
  )

print(args)

# Load data ---------------------------------------------------------------

loci <- LAVA::read.loci(args$loc_file)
input <-
  LAVA::process.input(
    input.info.file = args$info_file,
    sample.overlap.file = args$sample_overlap_file,
    ref.prefix = args$ref_prefix,
    phenos = args$phenotypes
  )

# Identify loci of interest -----------------------------------------------

bivar_results <-
  readRDS(
    list.files(
      here::here("results", "02_univar_bivar_test"),
      pattern = "bivar.lava.rds",
      full.names = T
      )
    ) %>%
  purrr::discard(is.null) %>%
  qdapTools::list_df2df() %>%
  dplyr::select(-X1)

# Filter for bivariate correlations that are significant
bivar_results_filt <-
  bivar_results %>%
  dplyr::filter(
    p < 0.05/nrow(bivar_results)
    )

# For each locus, count:
# (i) number of times a phenotype occurs in a bivariate correlation
# (ii) number of bivariate correlations overall
# (iii) number of distinct phenotypes (and list these phenotypes)
# Then filter for loci:
# (i) > 1 bivariate correlation
# (ii) where a phenotype occurs > 1 times (i.e. involved in multiple bivariate correlations and therefore should be modelled as outcome)
# And select phenotype that occurs the most times within a locus (can result in ties)
loci_of_interest <-
  bivar_results_filt %>%
  dplyr::select(locus, contains("phen")) %>%
  tidyr::pivot_longer(
    cols = contains("phen"),
    names_to = "phen_number",
    values_to = "phen"
  ) %>%
  dplyr::count(locus, phen, name = "n_phen") %>%
  dplyr::inner_join(
    bivar_results_filt %>%
      dplyr::count(locus, name = "n_bivar")
  ) %>%
  dplyr::filter(
    n_bivar > 1,
    n_phen > 1
  ) %>%
  dplyr::group_by(locus) %>%
  dplyr::slice_max(order_by = n_phen, n = 1)

# From these loci of interest, want to get our predictor phenotypes
# This should only be phenotypes correlated with the outcome phenotype
models <-
  bivar_results_filt %>%
  dplyr::filter(locus %in% loci_of_interest$locus) %>%
  dplyr::mutate(index = row_number()) %>%
  dplyr::select(locus, index, contains("phen")) %>%
  tidyr::pivot_longer(
    cols = contains("phen"),
    names_to = "phen_number",
    values_to = "phen"
  ) %>%
  # For each bivariate correlation, list phenotypes
  dplyr::group_by(
    locus, index
  ) %>%
  dplyr::summarise(
    phen_list = list(phen)
  ) %>%
  dplyr::inner_join(
    loci_of_interest %>%
      dplyr::rename(outcome = phen)
  ) %>%
  # Once loci of interest joined, filter for bivariate correlations that contain outcome of interest
  dplyr::group_by(locus, index) %>%
  dplyr::mutate(
    outcome_present =
      any(outcome %in% unlist(phen_list))
  ) %>%
  dplyr::filter(
    outcome_present == TRUE
  ) %>%
  # List phenotypes by locus and outcome (outcome necessary as some loci may have multiple outcome phenotypes where ties have occurred in previous df)
  dplyr::group_by(locus, outcome) %>%
  dplyr::summarise(
    phen_list = unlist(phen_list) %>%
      list() %>%
      lapply(., unique)
  ) %>%
  # Remove outcome phenotype from phenotype list, leaving predictor phenotypes
  dplyr::group_by(locus, outcome) %>%
  dplyr::mutate(
    predictors =
      phen_list %>%
      lapply(., function(vec){

        vec[vec != outcome]

      })
  )

# Main --------------------------------------------------------------------

results <-
  setNames(
    object = vector(mode = "list", length = length(unique(loci_of_interest$locus))),
    nm = unique(loci_of_interest$locus)
  ) %>%
  lapply(., function(x) x <- list())

for(i in 1:length(results)){

  # Process locus
  locus <-
    LAVA::process.locus(
      locus =
        loci %>%
        dplyr::filter(LOC == names(results)[i]),
      input = input
    )

  filtered_models <-
    models %>%
    dplyr::filter(
      locus == names(results)[i]
    )

  for(j in 1:nrow(filtered_models)){

    # Multi-reg
    results[[i]][[j]] <-
      run.multireg(
        locus,
        phenos =
          c(filtered_models$predictors[j] %>% unlist(),
            filtered_models$outcome[j])
        )

  }

  names(results[[i]]) <- filtered_models$outcome

}

# Save data ---------------------------------------------------------------
out_dir <- here::here("results", "03_partial_corr_multi_reg")
dir.create(out_dir, showWarnings = T)
saveRDS(
  results,
  file = file.path(out_dir, str_c(args$output_filename, ".multireg.lava.rds"))
)

print(str_c("Done! Analysis output written to ", out_dir, "/", args$output_filename, ".multireg.lava.rds"))



