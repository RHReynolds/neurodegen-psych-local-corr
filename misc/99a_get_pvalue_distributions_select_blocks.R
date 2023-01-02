# Description: extract p-values across LD blocks where no
# GWAS significant hits, but rg detected

# Load packages -----------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(vroom)

# Set arguments -----------------------------------------------------------

args <-
  list(
    ld_blocks = here::here("results", "01_input_prep", "gwas_filtered_loci_w_genes.rds"),
    results = here::here("results", "02_univar_bivar_test", "results_summary.rds"),
    out_dir = here::here("results", "02_univar_bivar_test")
  )

gwas_list <-
  setNames(
    object = list(
      vroom::vroom("/data/LDScore/GWAS/AD2019/AD_sumstats_Jansenetal_2019sept.lava.gz"),
      vroom::vroom("/data/LDScore/GWAS/LBD2020/LBD2020_hg19_rsids.lava.gz"),
      vroom::vroom("/data/LDScore/GWAS/PD2019_meta5_ex23andMe/PD2019_ex23andMe_hg19.lava.gz"),
      vroom::vroom("/data/LDScore/GWAS/BIP2021/BIP2021.lava.gz"),
      vroom::vroom("/data/LDScore/GWAS/MDD2019_ex23andMe/MDD2019.lava.gz"),
      vroom::vroom("/data/LDScore/GWAS/SCZ2018/SCZ2018.lava.gz")
    ),
    nm = c(
      "ad", "lbd", "pd", "bip", "mdd", "scz"
      )
  )


# Load data ---------------------------------------------------------------

results <-
  readRDS(args$results)

# LD blocks
ld_blocks <-
  readRDS(args$ld_blocks) %>%
  dplyr::ungroup()


# Main --------------------------------------------------------------------

# Determine number of traits in pair with genome-wide significant SNPs overlapping LD block
p_threshold <- 0.05/nrow(results$bivar)

bivar_pair_type <-
  results$bivar %>%
  dplyr::inner_join(
    ld_blocks %>%
      dplyr::select(locus, assoc_gwas)
  ) %>%
  dplyr::filter(
    p < p_threshold
  ) %>%
  dplyr::group_by(locus, phen1, phen2) %>%
  dplyr::mutate(
    phen1_overlapping_assoc_gwas =
      any(phen1 %in% unlist(assoc_gwas)),
    phen2_overlapping_assoc_gwas =
      any(phen2 %in% unlist(assoc_gwas)),
    sum_phen_overlapping_assoc_gwas =
      sum(phen1_overlapping_assoc_gwas, phen2_overlapping_assoc_gwas)
  )

# For those traits where there were no overlapping genome-wide significant SNPs
# What was the p-value of the lead SNP?
# What was the p-value distribution across the blocks?
# What was p-value of lead SNP across other blocks that were non-significant (but not tested)?
ld_blocks_phen_of_interest <-
  bivar_pair_type %>%
  dplyr::select(-c("phen1_overlapping_assoc_gwas", "phen2_overlapping_assoc_gwas")) %>%
  dplyr::filter(sum_phen_overlapping_assoc_gwas <= 1) %>%
  tidyr::pivot_longer(
    cols = c("phen1", "phen2"),
    names_to = "phen_col",
    values_to = "phen"
  ) %>%
  dplyr::group_by(locus, phen) %>%
  dplyr::mutate(
    phen_overlapping_assoc_gwas =
      any(phen %in% unlist(assoc_gwas))
  ) %>%
  dplyr::filter(phen_overlapping_assoc_gwas == FALSE) %>%
  dplyr::distinct(
    locus, chr, start, stop, phen
  )

# Filter LD blocks
ld_blocks_filtered <-
  ld_blocks %>%
  dplyr::filter(
    locus %in% ld_blocks_phen_of_interest$locus
  ) %>%
  dplyr::select(contains("locus"), chr)

ld_blocks_filtered <-
  setNames(
  ld_blocks_filtered %>%
    dplyr::group_by(locus) %>%
    dplyr::group_split(),
  nm =
    sort(ld_blocks_filtered$locus)
)

p_values <-
  setNames(
    vector(mode = "list", length = length(ld_blocks_filtered)),
    nm = names(ld_blocks_filtered)
  )

for(locus in names(ld_blocks_filtered)){

  print(stringr::str_c("Processing LD block: ", locus))

  ld_block <-
    ld_blocks_filtered[[locus]]

  phen <-
    ld_blocks_phen_of_interest %>%
    dplyr::filter(locus == ld_block$locus) %>%
    dplyr::pull(phen) %>%
    stringr::str_to_lower()

  p_values[[locus]] <-
    gwas_list[phen] %>%
    lapply(., function(df){
      df %>%
        dplyr::select(CHR, BP, SNP, P) %>%
        dplyr::filter(
          CHR == ld_block$chr,
          BP >= ld_block$locus_start,
          BP <= ld_block$locus_end
        )
    })

}

# Save data ---------------------------------------------------------------

saveRDS(
  p_values,
  file = file.path(args$out_dir, "ld_blocks_p.rds")
)

print(str_c(Sys.time(), " - Analysis done!"))

