# Description: get gene loci from LD blocks of interest

# Load packages -----------------------------------------------------------

library(GenomicRanges)
library(stringr)
library(tidyverse)

# Load data ---------------------------------------------------------------

ld_blocks <-
  read_delim("https://github.com/cadeleeuw/lava-partitioning/raw/main/LAVA_s2500_m25_f1_w200.blocks",
             delim = "\t")

gwas_overlap_df <-
  readRDS(
    file = here::here("results", "01_input_prep", "gwas_filtered_loci.rds")
  )

loci_of_interest <- c(681, 1273, 1719, 2281, 2351)

window_size <- 100000

ref <- rtracklayer::import("/data/references/ensembl/gtf_gff3/v87/Homo_sapiens.GRCh37.87.gtf")
ref <- ref %>% GenomeInfoDb::keepSeqlevels(c(1:22), pruning.mode = "coarse")
ref <- ref[ref$type == "gene"]

# Main --------------------------------------------------------------------

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
  dplyr::select(LOC, everything()) %>%
  dplyr::filter(LOC %in% loci_of_interest)

# Convert to granges
ld_blocks_gr <-
  ld_blocks %>%
  GenomicRanges::makeGRangesFromDataFrame(
    .,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = NULL,
    seqnames.field = "CHR",
    start.field = "START",
    end.field = "STOP"
  )

# Overlap granges objects
overlap <-
  GenomicRanges::findOverlaps(ld_blocks_gr, ref) %>%
  tibble::as_tibble()

# Generate df of genes that overlap which LD blocks
overlap_df <-
  tibble::tibble(
    locus = ld_blocks_gr[overlap$queryHits]$LOC,
    gene_id = ref[overlap$subjectHits]$gene_id,
    gene_name =
      ref[overlap$subjectHits]$gene_name %>%
      stringr::str_replace_all("-", "_") %>%
      stringr::str_replace_all("\\.", "_"),
    gene_biotype = ref[overlap$subjectHits]$gene_biotype %>% as.factor(),
    chr = ref[overlap$subjectHits] %>% GenomeInfoDb::seqnames() %>% as.character(),
    strand = ref[overlap$subjectHits] %>% BiocGenerics::strand() %>% as.character(),
    start = c(ref[overlap$subjectHits] %>% BiocGenerics::start()),
    end = ref[overlap$subjectHits] %>% BiocGenerics::end(),
    width = ref[overlap$subjectHits] %>% BiocGenerics::width()
  ) %>%
  # If there are two versions of a gene, choose the longest
  dplyr::group_by(gene_name) %>%
  dplyr::top_n(1, wt = width) %>%
  dplyr::filter(
    gene_biotype %in% c("protein_coding", "antisense", "lincRNA")
  ) %>%
  # Add 100 kb window
  dplyr::mutate(
    start_100kb =
      case_when(
        start - window_size <= 0 ~ 0,
        TRUE ~ start - window_size
      ),
    end_100kb = end + window_size
  )

# Create locus file
loci <-
  overlap_df %>%
  dplyr::ungroup() %>%
  dplyr::distinct(
    gene_id, chr, start_100kb, end_100kb
  ) %>%
  dplyr::select(
    LOC = gene_id,
    CHR = chr,
    START = start_100kb,
    STOP = end_100kb
  )

# Generate df of gwas loci that overlap genes/ld blocks
gwas_gr <-
  gwas_overlap_df %>%
  GenomicRanges::makeGRangesFromDataFrame(
    .,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = NULL,
    seqnames.field = "GWAS_CHR",
    start.field = "GWAS_start",
    end.field = "GWAS_end"
  )

locus_gr <-
  overlap_df %>%
  dplyr::ungroup() %>%
  dplyr::distinct(
    locus, gene_id, chr, start_100kb, end_100kb
  ) %>%
  GenomicRanges::makeGRangesFromDataFrame(
    .,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = NULL,
    seqnames.field = "chr",
    start.field = "start_100kb",
    end.field = "end_100kb"
  )

overlap <-
  GenomicRanges::findOverlaps(gwas_gr, locus_gr) %>%
  tibble::as_tibble()

gwas_gene_overlap_df <-
  gwas_overlap_df %>%
  dplyr::slice(overlap$queryHits) %>%
  dplyr::select(contains("GWAS")) %>%
    dplyr::bind_cols(
      locus_gr %>%
        as_tibble() %>%
        dplyr::rename_with(
          ~ stringr::str_c("GENE", .x, sep = "_")
        ) %>%
        dplyr::slice(overlap$subjectHits)
    ) %>%
  dplyr::select(contains("GENE"), contains("GWAS"), -contains("strand")) %>%
  dplyr::rename_with(
    ~ stringr::str_replace(.x,
                           pattern = "seqnames",
                           replacement = "CHR"),
    .col = dplyr::ends_with("seqnames")
    )


# Save data ---------------------------------------------------------------

out_dir <- here::here("results", "04_eqtl_univar_bivar")
dir.create(out_dir, showWarnings = T)
write_delim(
  loci,
  delim = "\t",
  file = file.path(out_dir, "gene_filtered.loci")
  )
saveRDS(
  overlap_df,
  file = file.path(out_dir, "gene_filtered_loci.rds")
  )
saveRDS(
  gwas_gene_overlap_df,
  file = file.path(out_dir, "gene_gwas_filtered_loci.rds")
)
