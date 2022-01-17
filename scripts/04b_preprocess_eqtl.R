# Description: pre-process eQTLs in gene loci of interest

# Load packages -----------------------------------------------------------

library(data.table)
library(tidyverse)
library(stringr)
library(qdapTools)

# Load data ---------------------------------------------------------------

# Genes run will be the same, just with different window sizes, so only need to load one locus file
# Used genes from 100 kb window
loci <-
  read_delim(
    file = here::here("results", "04_eqtl_univar_bivar", "window_100000", "gene_filtered.loci"),
    delim = "\t"
  )


setNames(
  object =
    list.files(
      path =
        here::here("results",
                   "04_eqtl_univar_bivar"),
      pattern = ".loci$",
      full.names = T,
      recursive = T
    ) %>%
    lapply(., function(file) read_delim(file, delim = "\t")),
  nm =
    list.files(
      path =
        here::here("results",
                   "04_eqtl_univar_bivar"),
      pattern = ".loci$",
      recursive = T
    ) %>%
    str_remove(., "/.*.loci")
)

eqtlgen <-
  fread("/data/eQTLGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")

psychencode <-
  fread("/data/psychencode/Full_hg19_cis-eQTL.txt.gz")
colnames(psychencode) <- colnames(fread("/data/psychencode/DER-08a_hg19_eQTL.significant.txt"))[!str_detect(colnames(fread("/data/psychencode/DER-08a_hg19_eQTL.significant.txt")), "FDR")]
snp_info <- fread("/data/psychencode/SNP_Information_Table_with_Alleles.txt")

# Main --------------------------------------------------------------------

# Filter for genes from loci
# Format for LAVA
eqtlgen_filtered <-
  eqtlgen %>%
  dplyr::filter(
    Gene %in% loci$LOC
  ) %>%
  dplyr::select(
    GENE = Gene,
    CHR = SNPChr,
    BP = SNPPos,
    SNP,
    A1 = AssessedAllele,
    A2 = OtherAllele,
    Z = Zscore,
    P = Pvalue,
    N = NrSamples
  )

psychencode_filtered <-
  psychencode %>%
  dplyr::filter(str_detect(gene_id, str_c(loci$LOC, collapse = "|"))) %>%
  dplyr::select(
    GENE = gene_id,
    chr = SNP_chr,
    position = SNP_start,
    PEC_id = SNP_id,
    BETA = regression_slope,
    P = nominal_pval
  ) %>%
  dplyr::inner_join(
    snp_info
  ) %>%
  dplyr::mutate(
    GENE =
      str_remove(GENE, "\\..*"),
    chr =
      str_remove(chr, "chr"),
    N = 1387
  ) %>%
  dplyr::select(
    GENE,
    CHR = chr,
    BP = position,
    SNP = Rsid,
    # As Psychencode is based on GTEx pipeline, and GTEx betas refer to ALT allele, A1 = ALT here
    A1 = ALT,
    A2 = REF,
    BETA,
    P,
    N
  )

eqtlgen_list <-
  setNames(
    eqtlgen_filtered %>%
      dplyr::group_split(
        GENE
      ),
    nm =
      str_c(
        "EQTLGEN_",
        eqtlgen_filtered %>%
          .[["GENE"]] %>%
          unique() %>%
          sort()
      )
  )

psychencode_list <-
  setNames(
    psychencode_filtered %>%
      dplyr::group_split(
        GENE
      ),
    nm =
      str_c(
        "PSYCHENCODE_",
        psychencode_filtered %>%
          .[["GENE"]] %>%
          unique() %>%
          sort()
      )
  )

# Save data ---------------------------------------------------------------

out_dir <- here::here("results", "04_eqtl_univar_bivar", "qtl_files")
dir.create(out_dir, showWarnings = T)

lava_ext <-  ".lava.gz"

for(i in 1:length(eqtlgen_list)){

  fwrite(
    eqtlgen_list[[i]],
    file =
      stringr::str_c(
        out_dir, "/",
        names(eqtlgen_list)[i],
        lava_ext
      ),
    sep = "\t"
  )

}

for(i in 1:length(psychencode_list)){

  fwrite(
    psychencode_list[[i]],
    file =
      stringr::str_c(
        out_dir, "/",
        names(psychencode_list)[i],
        lava_ext
      ),
    sep = "\t"
  )

}
