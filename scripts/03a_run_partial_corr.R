# Description: run partial correlations for GWAS traits

# Load packages -----------------------------------------------------------

library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

neuropsych_phen <- c("BIP2021", "MDD2019", "SCZ2018")

args <-
  list(
    ref_prefix = "/tools/MAGMA/reference_files/g1000_eur/g1000_eur",
    loc_file = here::here("results", "01_input_prep", "gwas_filtered.loci"),
    info_file = here::here("results", "01_input_prep", "input.info.txt"),
    sample_overlap_file = here::here("results", "01_input_prep", "sample_overlap.txt"),
    phenotypes = neuropsych_phen,
    output_filename = str_c(neuropsych_phen, collapse = ":")
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

# Filter for significant loci where:
# (i) > 2 distinct phenotypes
# (ii) all neuropsychiatric phenotypes are present
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
      length(),
    all_neuropsych_present =
      all(neuropsych_phen %in% unlist(distinct_phen))
  ) %>%
  dplyr::filter(
    all_neuropsych_present == TRUE
  )

# Combinations
combn <-
  neuropsych_phen %>%
  combn(m = 2) %>%
  t() %>%
  as_tibble() %>%
  dplyr::rename(
    phen1 = V1,
    phen2 = V2
  )

# Main --------------------------------------------------------------------

# Create empty lists for results
results <-
  vector(mode = "list", length = nrow(loci_of_interest)) %>%
  lapply(., function(x) x <- list())

for(i in 1:nrow(loci_of_interest)){

  locus_of_interest <-
    loci_of_interest %>%
    dplyr::slice(i)

  print(Sys.time())
  print(str_c("Locus: ", locus_of_interest$locus))

  # Process locus
  locus <-
    LAVA::process.locus(
      locus =
        loci %>%
        dplyr::filter(LOC == locus_of_interest$locus),
      input = input
    )

  # Extract some general locus info for the output
  loc_info <-
    data.frame(
      locus = locus$id,
      chr = locus$chr,
      start = locus$start,
      stop = locus$stop,
      n_snps = locus$n.snps,
      n_pcs = locus$K
    )

  for(j in 1:nrow(combn)){

    phen_to_corr <- c(combn$phen1[j], combn$phen2[j])
    phen_to_adj <- neuropsych_phen[!neuropsych_phen %in% phen_to_corr]

    print(Sys.time())
    print(str_c("Correlating: ", phen_to_corr[1], " ", phen_to_corr[2]))
    print(str_c("Adjusting for: ", phen_to_adj))

    # Partial corr
    loc_out <-
      run.pcor(locus, phenos = c(phen_to_corr, phen_to_adj))

    results[[i]][[j]] <-
      loc_info %>%
      dplyr::bind_cols(loc_out)

  }

}

# Save data ---------------------------------------------------------------
out_dir <- here::here("results", "03_partial_corr_multi_reg")
dir.create(out_dir, showWarnings = T)
saveRDS(
  results,
  file = file.path(out_dir, str_c(args$output_filename, ".partcorr.lava.rds"))
)

print(str_c("Done! Analysis output written to ", out_dir, "/", args$output_filename, ".partcorr.lava.rds"))
