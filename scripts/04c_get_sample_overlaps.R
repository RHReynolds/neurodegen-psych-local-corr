# Load packages -----------------------------------------------------------

library(tidyverse)
library(stringr)
library(qdapTools)

# Set arguments -----------------------------------------------------------

out_dir <- here::here("results", "04_eqtl_univar_bivar")

gwas <- c("AD2019", "LBD2020", "PD2019.meta5.ex23andMe", "BIP2021", "MDD2019", "SCZ2018") %>% sort()
eqtl_gene <-
  list.files(
    path = here::here("results", "04_eqtl_univar_bivar", "qtl_files"),
    pattern = "lava.gz"
  ) %>%
  basename() %>%
  str_remove(".lava.gz")

phenotypes <- c(gwas, eqtl_gene)

# Load ldsc results -------------------------------------------------------

gwas_rg <-
  read_delim(
    file = here::here("results", "01_input_prep", "ldsc_corr", "ldsc_correlations.txt"),
    delim = "\t"
  )

# Main --------------------------------------------------------------------

# Generate all possible permutations with repetitions
permutations <-
  tibble(
    p1 = c(gwas, eqtl_gene),
    p2 = c(gwas, eqtl_gene)
  ) %>%
  tidyr::expand(p1, p2)

# LDSC already run for gwas phenotypes
# Just need to add 0 in where eqtl phenotypes are used
# We are assuming that no sample overlaps exist between gwas and eqtl cohorts
all_rg <-
  permutations %>%
  dplyr::left_join(
    gwas_rg
  ) %>%
  dplyr::select(p1, p2, gcov_int) %>%
  dplyr::mutate(
    gcov_int =
      case_when(
        p1 %in% eqtl_gene & !p2 %in% eqtl_gene ~ 0,
        !p1 %in% eqtl_gene & p2 %in% eqtl_gene ~ 0,
        p1 %in% eqtl_gene & p2 %in% eqtl_gene & p1 != p2 ~ 0,
        p1 %in% eqtl_gene & p2 %in% eqtl_gene & p1 == p2 ~ 1,
        TRUE ~ gcov_int
      )
  )

###### Creating sample overlap matrix ######

n <- length(phenotypes)
covar_matrix <- matrix(NA,n,n)
rownames(covar_matrix) <- colnames(covar_matrix) <- phenotypes

for(i in phenotypes) {
  for(j in phenotypes) {

    covar_matrix[i,j] <-
      all_rg %>%
      dplyr::filter(p1 == i, p2 == j) %>%
      .[["gcov_int"]]

  }
}

# Sometimes there might be small differences in gcov_int depending on which phenotype was analysed as the outcome / predictor
if (!all(t(covar_matrix)==covar_matrix)) {
  covar_matrix[lower.tri(covar_matrix)] <- t(covar_matrix)[lower.tri(covar_matrix)]
}

# Standardise the matrix
covar_matrix <-
  covar_matrix %>%
  cov2cor() %>%
  round(digits = 5)

# Save data ---------------------------------------------------------------

write.table(
  covar_matrix,
  file = file.path(out_dir,
                   "sample_overlap.txt"),
  quote = F
)
