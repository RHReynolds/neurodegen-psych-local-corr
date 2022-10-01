# Description: compare to results of Smeland 2021 (https://pubmed.ncbi.nlm.nih.gov/32201043/)

# Load packages -----------------------------------------------------------

library(colochelpR)
library(here)
library(LAVA)
library(janitor)
library(tidyverse)
library(stringr)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# Set arguments -----------------------------------------------------------

args <-
  list(
    # table 3 downloaded from
    # https://www.sciencedirect.com/science/article/pii/S0006322320300615?via%3Dihub#tbl1
    smeland = here::here("raw_data", "02_univar_bivar", "smeland2021_table3.xlsx"),
    results = here::here("results", "02_univar_bivar_test", "results_summary.rds"),
    loci = here::here("results", "01_input_prep", "gwas_filtered.loci")
  )

# Load data ---------------------------------------------------------------

input_data <-
  list(
    smeland =
      readxl::read_xlsx(
        args$smeland
      ) %>%
      janitor::clean_names(),
    lava = readRDS(args$results),
    loci = LAVA::read.loci(args$loci)
  )

dbSNP <- SNPlocs.Hsapiens.dbSNP144.GRCh37

# Main --------------------------------------------------------------------

# turn rs id into GRCh37 location and create granges
input_data$smeland <-
  input_data$smeland %>%
  colochelpR::convert_rs_to_loc(SNP_column = "lead_snp", dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37) %>%
  dplyr::select(-chr) %>%
  tidyr::separate(
    col = loc,
    into = c("chr", "bp"),
    sep = ":"
  )

smeland_gr <-
  GenomicRanges::makeGRangesFromDataFrame(
    .,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = NULL,
    seqnames.field = "chr",
    start.field = "bp",
    end.field = "bp"
  )

# construct granges object with LD blocks for overlap with reference gtf
loci_gr <-
  input_data$loci %>%
  GenomicRanges::makeGRangesFromDataFrame(
    .,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = NULL,
    seqnames.field = "CHR",
    start.field = "START",
    end.field = "STOP"
  )

overlap <-
  GenomicRanges::findOverlaps(loci_gr, smeland_gr) %>%
  tibble::as_tibble()

loci_of_interest <-
  input_data$loci[overlap$queryHits, ] %>%
  bind_cols(
    input_data$smeland[overlap$subjectHits, ] %>%
      dplyr::mutate(
        snp_loc = str_c(chr, ":", bp),
        same_direction_of_effect =
          case_when(
            z_score_pd < 0 & odds_ratio_scz < 1 |
              z_score_pd > 0 & odds_ratio_scz > 1~ TRUE,
            TRUE ~ FALSE
          )
      ) %>%
      dplyr::select(lead_snp, snp_loc, a1_a2, nearest_gene, p_value_pd, z_score_pd, p_value_scz, odds_ratio_scz, same_direction_of_effect, conj_fdr) %>%
      select_all(., list(~ paste0("smeland_", .)))
    )

# overlap bivariate results with smeland results
# compare rho with consistency across effect directions in smeland
# i.e if both SNPs have same effect directions, might expect positive rho in local rg
# but if both SNPs have different effect directions, might expect negative rho in local rg
# note that SCZ is OR -- OR > 1 is equivalent to positive z-score
# while OR < 1 is equivalent to negative z-score
overlap_smeland <-
  input_data$lava$bivar %>%
  dplyr::filter(
    phen1 == "PD",
    phen2 == "SCZ"
  ) %>%
  dplyr::inner_join(
    loci_of_interest,
    by =
      c(
        "locus" = "LOC",
        "chr" = "CHR",
        "start" = "START",
        "stop" = "STOP"
      )
  ) %>%
  dplyr::mutate(
    rho_consistent_w_effect_consistency =
      case_when(
        rho < 0 & smeland_same_direction_of_effect == FALSE |
          rho > 0 & smeland_same_direction_of_effect == TRUE ~ TRUE,
        TRUE ~ FALSE
      )
  )

# Save data ---------------------------------------------------------------
out_dir <- here::here("results", "02_univar_bivar_test")
dir.create(out_dir, showWarnings = T)
saveRDS(
  overlap_smeland,
  file = file.path(out_dir, "results_overlap_with_smeland.rds")
)
