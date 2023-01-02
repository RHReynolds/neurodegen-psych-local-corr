# Description: extract top 50 p-values across each untested LD block

# Load packages -----------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(vroom)

# Set arguments -----------------------------------------------------------

args <-
  list(
    ld_blocks_all = "https://github.com/cadeleeuw/lava-partitioning/raw/main/LAVA_s2500_m25_f1_w200.blocks",
    ld_blocks_tested = here::here("results", "01_input_prep", "gwas_filtered_loci_w_genes.rds"),
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
ld_blocks_all <-
  read_delim(args$ld_blocks_all, delim = "\t")

ld_blocks_tested <-
  readRDS(args$ld_blocks_tested) %>%
  dplyr::ungroup()


# Main --------------------------------------------------------------------

# Remove tested blocks
ld_blocks <-
  ld_blocks_all %>%
  dplyr::mutate(
    locus = dplyr::row_number()
  ) %>%
  dplyr::anti_join(ld_blocks_tested, by = "locus")

ld_blocks <-
  setNames(
    ld_blocks %>%
      dplyr::group_by(locus) %>%
      dplyr::group_split(),
    nm =
      sort(ld_blocks$locus)
  )

p_values <-
  setNames(
    vector(mode = "list", length = length(ld_blocks)),
    nm = names(ld_blocks)
  )

for(locus in names(ld_blocks)){

  print(stringr::str_c("Processing LD block: ", locus))

  ld_block <-
    ld_blocks[[locus]]

  p_values[[locus]] <-
    gwas_list %>%
    lapply(., function(df){
      df %>%
        dplyr::select(CHR, BP, SNP, P) %>%
        dplyr::filter(
          CHR == ld_block$chr,
          BP >= ld_block$start,
          BP <= ld_block$stop
        ) %>%
        dplyr::slice_min(n = 50, order_by = P)
    })

}

# Save data ---------------------------------------------------------------

saveRDS(
  p_values,
  file = file.path(args$out_dir, "ld_blocks_p_min_50.rds")
)

print(str_c(Sys.time(), " - Analysis done!"))

