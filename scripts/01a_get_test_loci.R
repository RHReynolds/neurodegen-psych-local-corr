# Description: get test loci from GWAS traits

# Load packages -----------------------------------------------------------

library(data.table)
library(GenomicRanges)
library(qdapTools)
library(stringr)
library(tidyverse)

# Load data ---------------------------------------------------------------

ld_blocks <-
  read_delim("https://github.com/cadeleeuw/lava-partitioning/raw/main/LAVA_s2500_m25_f1_w200.blocks",
             delim = "\t")

gwas_list <-
  setNames(
    object = list(
      fread("/data/LDScore/GWAS/AD2019/AD_sumstats_Jansenetal_2019sept.lava.gz"),
      fread("/data/LDScore/GWAS/LBD2020/LBD2020_hg19_rsids.lava.gz"),
      fread("/data/LDScore/GWAS/PD2019_meta5_ex23andMe/PD2019_ex23andMe_hg19.lava.gz"),
      fread("/data/LDScore/GWAS/BIP2021/BIP2021.lava.gz"),
      fread("/data/LDScore/GWAS/MDD2019_ex23andMe/MDD2019.lava.gz"),
      fread("/data/LDScore/GWAS/SCZ2018/SCZ2018.lava.gz")
    ),
    nm = c("ad", "lbd", "pd", "bip", "mdd", "scz")
  )

# Main --------------------------------------------------------------------

# Filter for only genome-wide significant loci and convert to granges
gr_list <-
  gwas_list %>%
  lapply(., function(gwas){
    gwas %>%
      dplyr::filter(P < 5e-8) %>%
      GenomicRanges::makeGRangesFromDataFrame(
        .,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqinfo = NULL,
        seqnames.field = "CHR",
        start.field = "BP",
        end.field = "BP"
      )
  })

# Add locus id and rename remaining columns to fit LAVA requirements
ld_blocks <-
  ld_blocks %>%
  dplyr::rename_with(
    .fn = stringr::str_to_upper,
    .cols = everything()
  ) %>%
  dplyr::mutate(
    LOC = dplyr::row_number()
  ) %>%
  dplyr::select(LOC, everything())

# Convert to granges
ld_blocks_gr <-
  ld_blocks %>%
  GenomicRanges::makeGRangesFromDataFrame(
    .,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = NULL,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "stop"
  )

# Overlap granges objects
overlap_list <-
  gr_list %>%
  lapply(., function(gr){
    GenomicRanges::findOverlaps(gr, ld_blocks_gr, type = "within") %>%
      as_tibble()
  })


# Extract relevant rows from ld_blocks using overlap indices
# Rename remaining columns to fit with LAVA requirements
loci <-
  ld_blocks %>%
  dplyr::slice(
    overlap_list %>%
      qdapTools::list_df2df(col1 = "gwas") %>%
      .[["subjectHits"]] %>%
      unique()
    ) %>%
  dplyr::arrange(LOC)

# Generate df of GWAS loci that overlap which LD blocks
overlap_df_list <- vector(mode = "list", length = length(gr_list))

for(i in 1:length(gr_list)){

  gr <- gr_list[[i]]

  names(overlap_df_list)[i] <- names(gr_list)[i]

  if(names(gr_list)[i] == c("mdd")){

    overlap_df_list[[i]] <-
      gr %>%
      as_tibble() %>%
      dplyr::mutate(
        BETA = NA,
        OR = NA
      ) %>%
      dplyr::select(seqnames, start, end, SNP, A1, A2, MAF, BETA, OR, logOdds, SE, P, N)


  } else if(names(gr_list)[i] == c("scz")){

    overlap_df_list[[i]] <-
      gr %>%
      as_tibble() %>%
      dplyr::mutate(
        BETA = NA,
        logOdds = NA
      ) %>%
      dplyr::select(seqnames, start, end, SNP, A1, A2, MAF, BETA, OR, logOdds, SE, P, N)

  } else{

    overlap_df_list[[i]] <-
      gr %>%
      as_tibble() %>%
      dplyr::mutate(
        OR = NA,
        logOdds = NA
      ) %>%
      dplyr::select(seqnames, start, end, SNP, A1, A2, MAF, BETA, OR, logOdds, SE, P, N)

  }

  overlap_df_list[[i]] <-
    overlap_df_list[[i]] %>%
    dplyr::rename_with(
      ~ stringr::str_c("GWAS", .x, sep = "_")
    ) %>%
    dplyr::slice(overlap_list[[i]]$queryHits) %>%
    dplyr::bind_cols(
      ld_blocks_gr %>%
        as_tibble() %>%
        dplyr::rename_with(
          ~ stringr::str_c("LD", .x, sep = "_")
        ) %>%
        dplyr::slice(overlap_list[[i]]$subjectHits)
    ) %>%
    dplyr::select(-contains("strand")) %>%
    dplyr::rename_with(
      ~ stringr::str_replace(.x,
                             pattern = "seqnames",
                             replacement = "CHR"),
      .col = dplyr::ends_with("seqnames"))

}

overlap_df <-
  overlap_df_list %>%
  qdapTools::list_df2df(col1= "GWAS")

# Save data ---------------------------------------------------------------

out_dir <- here::here("results", "01_input_prep")
dir.create(out_dir, showWarnings = T)
write_delim(
  loci,
  delim = "\t",
  file = file.path(out_dir, "gwas_filtered.loci")
  )
saveRDS(
  overlap_df,
  file = file.path(out_dir, "gwas_filtered_loci.rds")
  )
