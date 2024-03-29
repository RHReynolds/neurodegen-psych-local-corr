# Description: get sample overlaps for AD/PD GWASs without by-proxy cases

# Load packages -----------------------------------------------------------

library(tidyverse)
library(stringr)
library(qdapTools)

# Set arguments -----------------------------------------------------------

out_dir <- here::here("results", "05_noproxy_univar_bivar", "ldsc_corr")
dir.create(out_dir, showWarnings = T, recursive = T)

gwas_dir <- "/data/LDScore/GWAS/"
phenotypes <- c("AD2019", "AD2019.Kunkle", "PD2019.ex23andMe.exUKBB", "PD2019.meta5.ex23andMe", "LBD2020") %>% sort()

ldsc_args <-
  list(
    ldsc = "/tools/ldsc/ldsc.py",
    ref_ld = "/data/LDScore/Reference_Files/eur_w_ld_chr/"
  )

# Load files paths --------------------------------------------------------

gwas_details <-
  tibble(
    full_path = list.files(
      path = gwas_dir,
      recursive = T,
      full.names = T,
      pattern = ".sumstats.gz"
    )
  ) %>%
  dplyr::mutate(
    file_name = basename(full_path) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_replace_all("_", "\\.")
  ) %>%
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
      str_remove(".sumstats.gz"),
    phen_2 =
      basename(full_path_2) %>%
      str_remove(".sumstats.gz")
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

  rg_warning <-
    any(
      read_lines(file_paths[i]) %>%
        stringr::str_detect("rg out of bounds")
    )

  h2_warning <-
    any(
      read_lines(file_paths[i]) %>%
        stringr::str_detect("h2 out of bounds")
    )

  if(rg_warning){

    print(stringr::str_c("File ", i, "WARNING: rg out of bounds"))

    files[[i]] <-
      read_table(
        file = file_paths[i],
        skip = 62, # Number of lines to skip in log file
        n_max = 1 # Number of lines to read
      )

  }

  if(h2_warning){

    print("ERROR: h2 out of bounds for one phenotype; no rg calculated")

  }

  if(rg_warning == FALSE & h2_warning == FALSE){

    files[[i]] <-
      read_table(
        file = file_paths[i],
        skip = 60, # Number of lines to skip in log file
        n_max = 1 # Number of lines to read
      )

  }

}

# Rename AD to AD2013 to fit with LAVA file
all_rg <-
  files %>%
  qdapTools::list_df2df(col1 = "list_name") %>%
  dplyr::select(-list_name) %>%
  dplyr::mutate(
    p1 = basename(p1) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_replace_all("_", "\\."),
      # stringr::str_replace("AD$", "AD2013"),
    p2 = basename(p2) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_replace_all("_", "\\.")
      # stringr::str_replace("AD$", "AD2013")
  )

###### Creating sample overlap matrix ######

# Rename AD to AD2013 to fit with LAVA file
# phenotypes <- str_replace(phenotypes, "AD$", "AD2013")

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
  file = file.path(here::here("results", "05_noproxy_univar_bivar"),
                   "sample_overlap.txt"),
  quote = F
)
