# Description: run multiple regressions for GWAS traits with no by-proxy cases

# Load packages -----------------------------------------------------------

library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

phenotypes <-  c("AD2019.Kunkle", "PD2019.ex23andMe.exUKBB", "LBD2020") %>% sort()

args <-
  list(
    ref_prefix = "/tools/MAGMA/reference_files/g1000_eur/g1000_eur",
    loc_file = here::here("results", "05_noproxy_univar_bivar", "filtered.loci"),
    info_file = here::here("results", "05_noproxy_univar_bivar", "input.info.txt"),
    sample_overlap_file = here::here("results", "05_noproxy_univar_bivar", "sample_overlap.txt"),
    phenotypes = phenotypes,
    loci_of_interest = c(2351),
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

# Locus of interest--------------------------------------------------------

models <-
  tibble(
    locus = c(2351),
    outcome = c("LBD2020"),
    phen_list = list(c("AD2019.Kunkle", "PD2019.ex23andMe.exUKBB", "LBD2020")),
    predictors = list(c("AD2019.Kunkle", "PD2019.ex23andMe.exUKBB")),
  )

# Main --------------------------------------------------------------------

results <-
  setNames(
    object = vector(mode = "list", length = length(args$loci_of_interest)),
    nm = unique(args$loci_of_interest)
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
out_dir <- here::here("results", "05_noproxy_univar_bivar")
dir.create(out_dir, showWarnings = T)
saveRDS(
  results,
  file = file.path(out_dir, str_c(args$output_filename, ".multireg.lava.rds"))
)

print(str_c("Done! Analysis output written to ", out_dir, "/", args$output_filename, ".multireg.lava.rds"))



