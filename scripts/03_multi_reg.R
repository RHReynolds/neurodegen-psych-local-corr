# Load packages -----------------------------------------------------------

library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

phenotypes <- c("AD2019", "LBD2020", "PD2019.meta5.ex23andMe", "PD2018.AOO", "RBD2020")
outcome <- "PD2019.meta5.ex23andMe"

args <- 
  list(
    ref_prefix = "/tools/MAGMA/reference_files/g1000_eur/g1000_eur",
    loc_file = here::here("results", "01_input_prep", "ad_pd.loci"),
    info_file = here::here("results", "01_input_prep", "input.info.txt"),
    sample_overlap_file = here::here("results", "01_input_prep", "sample_overlap.txt"),
    phenotypes = phenotypes,
    output_filename = str_c(phenotypes, collapse = ":")
  )

print(args)

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

# WARNING: Still need to add filter to ensure outcome is present for locus
loci_of_interest <-
  bivar_results %>% 
  dplyr::filter(p < 0.05/nrow(bivar_results)) %>% 
  dplyr::group_by(locus) %>% 
  dplyr::summarise(
    distinct_phen = 
      list(c(phen1, phen2)) %>% 
      lapply(., unique),
    n_distinct_phen =
      distinct_phen %>% 
      unlist() %>% 
      length()
    ) %>% 
    dplyr::filter(n_distinct_phen > 2)

# Re-order phenotypes to ensure outcome is at the end
loci_of_interest$distinct_phen <- 
  loci_of_interest$distinct_phen %>% 
  lapply(., function(vec){
    
    c((vec[vec != outcome]), outcome)
    
  })
    
# Load data ---------------------------------------------------------------

loci <- LAVA::read.loci(args$loc_file)
input <- 
  LAVA::process.input(
    input.info.file = args$info_file,
    sample.overlap.file = args$sample_overlap_file,
    ref.prefix = args$ref_prefix,
    phenos = args$phenotypes
  )

# Main --------------------------------------------------------------------

# for(i in 1:nrow(loci_of_interest)){}

locus_of_interest <- 
  loci_of_interest %>% 
  dplyr::slice(i)

# Process locus
locus <- 
  LAVA::process.locus(
    locus = 
      loci %>% 
      dplyr::filter(LOC == locus_of_interest$locus), 
    input = input
  ) 

# Multi-reg
results <- 
  run.multireg(locus, phenos = c(phenotypes[phenotypes != outcome], outcome))



