# Description: run partial correlations for gwas-eqtl traits

# # Run with:
# cd /home/rreynolds/misc_projects/neurodegen-psych-local-corr
#
# nohup Rscript \
# /home/rreynolds/misc_projects/neurodegen-psych-local-corr/scripts/04f_run_partial_corr.R \
# &>/home/rreynolds/misc_projects/neurodegen-psych-local-corr/logs/04f_run_partial_corr.log&

# Load packages -----------------------------------------------------------

library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

args <-
  list(
    ref_prefix = "/tools/MAGMA/reference_files/g1000_eur/g1000_eur",
    loc_file = here::here("results", "04_eqtl_univar_bivar", "window_100000", "gene_filtered.loci"),
    info_file = here::here("results", "04_eqtl_univar_bivar", "input.info.txt"),
    sample_overlap_file = here::here("results", "04_eqtl_univar_bivar", "sample_overlap.txt"),
    loci_of_interest = c(
      "ENSG00000168078", # PBK
      "ENSG00000168079", # SCARA5
      "ENSG00000170209", # ANKK1
      "ENSG00000069399", # BCL3
      "ENSG00000130202", # PVRL2
      "ENSG00000130204" # TOMM40
    ),
    tidy_bivar = readRDS(here::here("results", "04_eqtl_univar_bivar", "results_summary.rds")),
    output_filename = c("eqtl_pcor")
  )

print(args)

# Load data ---------------------------------------------------------------

loci <- LAVA::read.loci(args$loc_file)

# Identify loci of interest -----------------------------------------------

bivar_results <-
  setNames(
    object =
      list.files(
        path =
          here::here(
            "results",
            "04_eqtl_univar_bivar",
            "window_100000",
            "bivar"
          ),
        pattern = ".lava.rds",
        full.names = T
      ) %>%
      lapply(., function(file) readRDS(file)),
    nm =
      list.files(
        path =
          here::here(
            "results",
            "04_eqtl_univar_bivar",
            "window_100000",
            "bivar"
          ),
        pattern = ".lava.rds"
      ) %>%
      str_remove(., "\\..*.lava.rds")
  ) %>%
  purrr::discard(is.null) %>%
  qdapTools::list_df2df(col1 = "eqtl_gene") %>%
  dplyr::rename(
    gene_locus = locus
  ) %>%
  tidyr::separate(
    col = eqtl_gene,
    into = c("eqtl_dataset", NA),
    sep = ":"
  ) %>%
  dplyr::mutate(
    fdr = p.adjust(p, method = "fdr")
  )

# Filter significant loci with triangles of correlation with a gene eQTL:
loci_of_interest <-
  bivar_results %>%
  dplyr::filter(gene_locus %in% args$loci_of_interest, fdr < 0.05) %>%
  dplyr::filter(
    !(eqtl_dataset == "PSYCHENCODE" & phen1 == "BIP2021") &
      !(eqtl_dataset == "PSYCHENCODE" & gene_locus %in% c("ENSG00000130202", "ENSG00000130204")) &
      !(eqtl_dataset == "EQTLGEN" & gene_locus %in% c("ENSG00000130202", "ENSG00000130204") & phen2 == "PD2019.meta5.ex23andMe") &
      !(eqtl_dataset == "EQTLGEN" & gene_locus %in% c("ENSG00000069399") & phen2 == "LBD2020")
  ) %>%
  dplyr::group_by(eqtl_dataset, gene_locus) %>%
  dplyr::summarise(
    distinct_phen =
      list(c(phen1, phen2)) %>%
      lapply(., unique)
  ) %>%
  dplyr::ungroup()

loci_of_interest <-
  loci_of_interest %>%
  dplyr::inner_join(
    args$tidy_bivar$bivar$window_100000 %>%
      dplyr::distinct(gene_locus, gene_name)
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

  phenotypes <-
    locus_of_interest$distinct_phen %>%
    unlist()

  # Combinations
  combn <-
    phenotypes %>%
    combn(m = 2) %>%
    t() %>%
    as_tibble() %>%
    dplyr::rename(
      phen1 = V1,
      phen2 = V2
    )

  print(Sys.time())
  print(str_c("Locus: ", locus_of_interest$gene_locus))

  input <-
    LAVA::process.input(
      input.info.file = args$info_file,
      sample.overlap.file = args$sample_overlap_file,
      ref.prefix = args$ref_prefix,
      phenos = phenotypes
    )

  # Process locus
  locus <-
    LAVA::process.locus(
      locus =
        loci %>%
        dplyr::filter(LOC == locus_of_interest$gene_locus),
      input = input
    )

  # Extract some general locus info for the output
  loc_info <-
    data.frame(
      gene_locus = locus$id,
      gene_name = locus_of_interest$gene_name,
      chr = locus$chr,
      start = locus$start,
      stop = locus$stop,
      n_snps = locus$n.snps,
      n_pcs = locus$K
    )

  for(j in 1:nrow(combn)){

    phen_to_corr <- c(combn$phen1[j], combn$phen2[j])
    phen_to_adj <- phenotypes[!phenotypes %in% phen_to_corr]

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
out_dir <- here::here("results", "04_eqtl_univar_bivar", "window_100000", "pcor")
dir.create(out_dir, showWarnings = T)
saveRDS(
  results,
  file = file.path(out_dir, str_c(args$output_filename, ".partcorr.lava.rds"))
)

print(str_c("Done! Analysis output written to ", out_dir, "/", args$output_filename, ".partcorr.lava.rds"))
