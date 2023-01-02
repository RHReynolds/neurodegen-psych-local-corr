# Description: tidy partial correlations for gwas-eqtl traits

# Load packages -----------------------------------------------------------

library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

args <-
  list(
    loci_of_interest = c(
      "ENSG00000168078", # PBK
      "ENSG00000168079", # SCARA5
      "ENSG00000170209", # ANKK1
      "ENSG00000069399", # BCL3
      "ENSG00000130202", # PVRL2
      "ENSG00000130204" # TOMM40
    )
  )

print(args)

# Load data ---------------------------------------------------------------
source(here::here("R", "theme_rhr.R"))

pcorr <-
  readRDS(here::here("results", "04_eqtl_univar_bivar", "window_100000", "pcor", "eqtl_pcor.partcorr.lava.rds"))

bivar_results <-
  readRDS(here::here("results", "04_eqtl_univar_bivar", "results_summary.rds"))$bivar$window_100000 %>%
  dplyr::filter(gene_locus %in% args$loci_of_interest)

# Main --------------------------------------------------------------------

pcorr <-
  pcorr %>%
  lapply(., function(list){

    list %>%
      qdapTools::list_df2df(col1 = "list_name") %>%
      dplyr::mutate(
        eqtl_dataset =
          case_when(
            stringr::str_detect(phen1, "PSYCH") |
              stringr::str_detect(phen2, "PSYCH") |
              stringr::str_detect(z, "PSYCH") ~ "PSYCHENCODE",
            stringr::str_detect(phen1, "EQTLGEN") |
              stringr::str_detect(phen2, "EQTLGEN") |
              stringr::str_detect(z, "EQTLGEN") ~ "EQTLGEN",
          )
      ) %>%
      dplyr::mutate(
        across(
          .cols = phen1:z,
          ~ case_when(
            !str_detect(.x, "ENSG") ~ .x %>%
              str_replace_all("[:digit:]", "") %>%
              str_remove("\\..*"),
            str_detect(.x, "ENSG") ~ str_remove(.x, ".*_")
          )
        )
      )

  }) %>%
  qdapTools::list_df2df(col1 = "list_name_2") %>%
  dplyr::select(
    eqtl_dataset,
    gene_locus = locus,
    everything(),
    -contains("list_name")
    ) %>%
  as_tibble()

# Tile --------------------------------------------------------------------

data_to_plot <-
  pcorr %>%
  dplyr::filter(!is.na(pcor)) %>%
  dplyr::select(
    eqtl_dataset, gene_locus, phen1, phen2, pcor, p_cond = p
  ) %>%
  dplyr::inner_join(
    bivar_results %>%
      dplyr::select(
        ld_block,
        eqtl_dataset,
        gene_locus,
        gene_name,
        contains("phen"),
        rho,
        p_uncond = p
      ),
    by = c("eqtl_dataset", "gene_locus", "phen1", "phen2")
  ) %>%
  tidyr::pivot_longer(
    cols = c("pcor", "p_cond", "rho", "p_uncond")
  ) %>%
  dplyr::mutate(
    model =
      case_when(
        name %in% c("pcor", "p_cond") ~ "Conditioned",
        TRUE ~ "Unconditioned"
      ) %>%
      fct_rev(),
    coef =
      case_when(
        name %in% c("pcor", "rho") ~ "corr",
        TRUE ~ "p"
      )
  ) %>%
  dplyr::select(-name) %>%
  tidyr::pivot_wider(
    names_from = c("coef"),
    values_from = c("value")
  ) %>%
  dplyr::inner_join(
    bivar_results %>%
      dplyr::select(
        eqtl_dataset, gene_locus, phen1, phen2, fdr
      )
  )

data_to_plot %>%
  dplyr::mutate(
    rg_fill =
      dplyr::case_when(
        # Different p-value cut-offs depending on conditioning
        p < 0.05 & model == "Conditioned" ~ round(corr, 2),
        fdr < 0.05 & model == "Unconditioned" ~ round(corr, 2)
      ),
    locus =
      stringr::str_c(
        eqtl_dataset, ": ", gene_name
      )
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = phen1,
      y = phen2,
      fill = rg_fill,
      label = round(corr, 2)
    )
  ) +
  ggplot2::geom_tile(colour = "black") +
  ggplot2::geom_text(
    size = 3
  ) +
  ggplot2::facet_wrap(vars(ld_block, locus, model), scales = "free") +
  ggplot2::labs(x = "", y = "", fill = "Local genetic correlation (rg)") +
  ggplot2::scale_fill_distiller(
    palette = "RdYlBu",
    limits = c(-1, 1)
  ) +
  ggplot2::scale_colour_manual(values = c("black", "white")) +
  ggplot2::guides(colour = F) +
  # theme_rhr +
  ggplot2::theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  )

data_to_plot %>%
  dplyr::mutate(
    rg_fill =
      dplyr::case_when(
        # Different p-value cut-offs depending on conditioning
        p < 0.05 & model == "Conditioned" ~ round(corr, 2),
        fdr < 0.05 & model == "Unconditioned" ~ round(corr, 2)
      ),
    locus =
      stringr::str_c(
        eqtl_dataset, ": ", gene_name
      )
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = phen1,
      y = phen2,
      fill = rg_fill,
      label = round(p, 6)
    )
  ) +
  ggplot2::geom_tile(colour = "black") +
  ggplot2::geom_text(
    size = 3
  ) +
  ggplot2::facet_wrap(vars(ld_block, locus, model), scales = "free") +
  ggplot2::labs(x = "", y = "", fill = "Local genetic correlation (rg)") +
  ggplot2::scale_fill_distiller(
    palette = "RdYlBu",
    limits = c(-1, 1)
  ) +
  ggplot2::scale_colour_manual(values = c("black", "white")) +
  ggplot2::guides(colour = F) +
  # theme_rhr +
  ggplot2::theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  )

# Pval ---------------------------------------------------------

data_to_plot <-
  pcorr %>%
  dplyr::filter(!is.na(pcor)) %>%
  dplyr::select(
    eqtl_dataset, gene_locus, phen1, phen2, pcor, p_cond = p
  ) %>%
  dplyr::inner_join(
    bivar_results %>%
      dplyr::select(
        ld_block,
        eqtl_dataset,
        gene_locus,
        gene_name,
        contains("phen"),
        rho,
        p_uncond = p
      ),
    by = c("eqtl_dataset", "gene_locus", "phen1", "phen2")
  ) %>%
  tidyr::pivot_longer(
    cols = c("pcor", "p_cond", "rho", "p_uncond")
  ) %>%
  dplyr::mutate(
    model =
      case_when(
        name %in% c("pcor", "p_cond") ~ "Conditioned",
        TRUE ~ "Unconditioned"
      ) %>%
      fct_rev(),
    coef =
      case_when(
        name %in% c("pcor", "rho") ~ "corr",
        TRUE ~ "p"
      )
  )

data_to_plot %>%
  dplyr::filter(
    !str_detect(phen2, "ENSG"),
    coef == "p"
  ) %>%
  dplyr::mutate(
    log10p = -log10(value),
    label = str_c(gene_name, "\n", eqtl_dataset),
    trait = str_c(phen1, "-", phen2)
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = label,
      y = log10p,
      fill = model
    )
  ) +
  ggplot2::geom_col(
    position = position_dodge2(preserve = "single")
  ) +
  ggplot2::facet_wrap(
    vars(ld_block, trait),
    scales = "free",
    nrow = 1
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  ggplot2::labs(x = "", y = "-log10(p)", fill = "Correlation") +
  theme_rhr
