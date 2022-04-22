# Description: tidy local genetic correlation results for gwas-eqtl traits

# Load packages -----------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

args <-
  list(
    results_dir = here::here("results", "04_eqtl_univar_bivar"),
    gene_filtered_loci = readRDS(
      here::here("results", "04_eqtl_univar_bivar", "window_100000", "gene_filtered_loci.rds")
    )
  )

# Load data ---------------------------------------------------------------

# LAVA results
results <-
  setNames(
    vector(mode = "list", length = 2),
    nm = c("univ", "bivar")
  ) %>%
  lapply(., function(x){

    setNames(
      vector(mode = "list", length = 2),
      nm = c("window_100000", "window_50000")
    )

  })

for(test in names(results)){

  for(window in names(results[[test]]))

    results[[test]][[window]] <-
      setNames(
        object =
          list.files(
            path =
              file.path(
                args$results_dir,
                window,
                test
              ),
            pattern = ".lava.rds",
            full.names = T
          ) %>%
          lapply(., function(file) readRDS(file)),
        nm =
          list.files(
            path =
              file.path(
                args$results_dir,
                window,
                test
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
      dplyr::inner_join(
        args$gene_filtered_loci %>%
          dplyr::select(
            ld_block = locus, gene_id, gene_name
          ),
        by = c("gene_locus" = "gene_id")
      ) %>%
      dplyr::select(
        ld_block, eqtl_dataset, gene_locus, gene_name, everything()
      ) %>%
      dplyr::arrange(ld_block, gene_locus, eqtl_dataset) %>%
      as_tibble()

}

# Main --------------------------------------------------------------------

# Tidy phenotypes and apply multiple test correction
results$univ <-
  results$univ %>%
  lapply(., function(df){

    df %>%
      dplyr::mutate(
        phen =
          case_when(
            !str_detect(phen, "ENSG") ~ phen %>%
              str_replace_all("[:digit:]", "") %>%
              str_remove("\\..*"),
            str_detect(phen, "ENSG") ~ str_remove(phen, ".*_")
          )
      )

  })

results$bivar <-
  results$bivar %>%
  lapply(., function(df){

    df %>%
      dplyr::mutate(
        across(
          .cols = contains("phen"),
          ~ case_when(
            !str_detect(.x, "ENSG") ~ .x %>%
              str_replace_all("[:digit:]", "") %>%
              str_remove("\\..*"),
            str_detect(.x, "ENSG") ~ str_remove(.x, ".*_")
          )
        )
      ) %>%
      dplyr::mutate(
        fdr = p.adjust(p, method = "fdr"),
        bonferroni = p.adjust(p, method = "bonferroni")
      )

  })

# Save data ---------------------------------------------------------------
out_dir <- here::here("results", "04_eqtl_univar_bivar")
dir.create(out_dir, showWarnings = T)
saveRDS(
  results,
  file = file.path(out_dir, "results_summary.rds")
)
