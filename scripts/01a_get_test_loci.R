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

ref <- import("/data/references/ensembl/gtf_gff3/v87/Homo_sapiens.GRCh37.87.gtf")
ref <- ref %>% keepSeqlevels(c(1:22), pruning.mode = "coarse")
ref <- ref[ref$type == "gene"]

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

# Generate df of associated GWAS and genes that overlap investigated LD blocks
loci_gr <-
  loci %>%
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
  GenomicRanges::findOverlaps(loci_gr, ref) %>%
  tibble::as_tibble()

loci_genes <-
  tibble::tibble(
    locus = loci_gr[overlap$queryHits]$LOC,
    chr = loci_gr[overlap$queryHits] %>% GenomeInfoDb::seqnames() %>% as.character(),
    locus_start = loci_gr[overlap$queryHits] %>% BiocGenerics::start(),
    locus_end = loci_gr[overlap$queryHits] %>% BiocGenerics::end(),
    gene_id = ref[overlap$subjectHits]$gene_id,
    gene_name =
      ref[overlap$subjectHits]$gene_name %>%
      stringr::str_replace_all("-", "_") %>%
      stringr::str_replace_all("\\.", "_")
  ) %>%
  dplyr::group_by(locus, chr, locus_start, locus_end) %>%
  dplyr::summarise(
    n_overlapping_genes = n(),
    overlapping_gene_id = list(gene_id),
    overlapping_gene_name = list(gene_name)
      )

loci_genes_gwas_df <-
  loci_genes %>%
  dplyr::inner_join(
    overlap_df %>%
      dplyr::distinct(GWAS, LD_LOC) %>%
      dplyr::group_by(LD_LOC) %>%
      dplyr::summarise(
        n_assoc_gwas = n(),
        assoc_gwas = list(str_to_upper(GWAS))
        ),
    by = c("locus" = "LD_LOC")
  )

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
saveRDS(
  loci_genes_gwas_df,
  file = file.path(out_dir, "gwas_filtered_loci_w_genes.rds")
)
