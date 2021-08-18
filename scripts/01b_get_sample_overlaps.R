# Load packages -----------------------------------------------------------

library(tidyverse)
library(stringr)
library(qdapTools)

# Set arguments -----------------------------------------------------------

gwas_dir <- "/data/LDScore/GWAS/"

phenotypes <- c("AD2019", "LBD2020", "PD2019.meta5.ex23andMe", "BIP2021", "MDD2019", "SCZ2018") %>% sort()

out_dir <- here::here("results", "01_input_prep", "ldsc_corr")
dir.create(out_dir, showWarnings = T)

ldsc_args <-
  list(
    ldsc = "/tools/ldsc/ldsc.py",
    ref_ld = "/data/LDScore/Reference_Files/eur_w_ld_chr/"
    )

# Load files paths --------------------------------------------------------

gwas_details <-
  tibble(
    full_path = list.files(path = gwas_dir,
                           recursive = T,
                           full.names = T,
                           pattern = ".sumstats.gz")) %>%
  dplyr::mutate(
    file_name = basename(full_path) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_replace_all("_", "\\.")) %>%
  dplyr::filter(file_name %in% phenotypes)

# Main --------------------------------------------------------------------

###### Running bivariate LDSC ######

# Generate all possible permutations with repetitions of file paths to run against each other
permutations <-
  gwas_details %>%
  dplyr::mutate(
    full_path_1 = full_path %>%
      as.character(),
    full_path_2 = full_path %>%
      as.character()
  ) %>%
  tidyr::expand(full_path_1, full_path_2) %>%
  dplyr::mutate(
    phen_1 =
      basename(full_path_1) %>%
      str_remove("\\..*"),
    phen_2 =
      basename(full_path_2) %>%
      str_remove("\\..*")
  ) %>%
  dplyr::arrange(full_path_1)

# Run loop across permutations, with external call to ldsc
# If many permutations, might consider parallelising process
for(i in 1:nrow(permutations)){

  args <-
    list(
      phen_1 = permutations$full_path_1[i],
      phen_2 = permutations$full_path_2[i],
      out_prefix =
        stringr::str_c(
          permutations$phen_1[i],
          "_",
          permutations$phen_2[i]
        )
    )

  rg_arg <-
    stringr::str_c(
      " ", ldsc_args$ldsc,
      " --rg ", args$phen_1, ",", args$phen_2,
      " --ref-ld-chr ", ldsc_args$ref_ld,
      " --w-ld-chr ", ldsc_args$ref_ld,
      " --out ", out_dir, "/", args$out_prefix, "_rg"
    )

  # Run command
  system2(command = "python", args = rg_arg)

}

###### Extracting LDSC results ######

# Extract all outputs into single file
file_paths <-
  list.files(
    out_dir,
    pattern = "_rg.log",
    full.names = T
  )

files <-
  vector(mode = "list",
         length = length(file_paths))

for(i in 1:length(file_paths)){

  files[[i]] <-
    read_table(
      file = file_paths[i],
      skip = 60, # Number of lines to skip in log file
      n_max = 1 # Number of lines to read
    )

}

all_rg <-
  files %>%
  qdapTools::list_df2df(col1 = "list_name") %>%
  dplyr::select(-list_name) %>%
  dplyr::mutate(
    p1 = basename(p1) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_replace_all("_", "\\."),
    p2 = basename(p2) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_replace_all("_", "\\.")
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

write_delim(
  all_rg,
  delim = "\t",
  file = file.path(out_dir, "ldsc_correlations.txt")
)

write.table(
  covar_matrix,
  file = file.path(here::here("results", "01_input_prep"),
                   "sample_overlap.txt"),
  quote = F
)
