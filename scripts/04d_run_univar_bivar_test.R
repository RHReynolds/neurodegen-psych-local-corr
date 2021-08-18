# Load packages -----------------------------------------------------------

library(doParallel)
library(foreach)
library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

eqtls <- c("EQTLGEN", "PSYCHENCODE") %>% sort()

out_dir <-
  setNames(
    c(
      here::here("results", "04_eqtl_univar_bivar", "univ"),
      here::here("results", "04_eqtl_univar_bivar", "bivar")
    ),
    nm = c("univ", "bivar")
  )

args <-
  list(
    cores = 10,
    ref_prefix = "/tools/MAGMA/reference_files/g1000_eur/g1000_eur",
    loc_file = here::here("results", "04_eqtl_univar_bivar", "gene_filtered.loci"),
    info_file = here::here("results", "04_eqtl_univar_bivar", "input.info.txt"),
    sample_overlap_file = here::here("results", "04_eqtl_univar_bivar", "sample_overlap.txt"),
    phenotypes = c(""),
    output_filename = c("")
  )

print(args)

# Bivariate threshold from previous analyses (02_run_univar_bivar_test.rmd)
# Used to filter bivariate results to include only phenos with significant local rg
# These phenos will then be used together with the eQTL
bivar_threshold <- 0.05/1603 # n bivariate tests

# Set univariate threshold to liberal 0.05
univar_threshold <-
  0.05

# Set up for parallel run -------------------------------------------------

# Create out directories
out_dir %>%
  lapply(., function(x) dir.create(x, showWarnings = T))

# Run in parallel
cl <- parallel::makeCluster(args$cores)

# Register clusters
doParallel::registerDoParallel(cl)
foreach::getDoParWorkers()

# Load data ---------------------------------------------------------------

loci <-
  LAVA::read.loci(args$loc_file) %>%
  dplyr::arrange(LOC)

gene_filtered_loci <-
  readRDS(
    here::here("results", "04_eqtl_univar_bivar", "gene_filtered_loci.rds")
  )

input_info <-
  read_delim(
    args$info_file,
    delim = "\t"
  )

# Load and filter bivariate results to include only phenos with significant local rg
results <-
  list.files(
    path =
      here::here("results",
                 "02_univar_bivar_test"),
    pattern = "bivar.lava.rds",
    full.names = T
  ) %>%
  lapply(., function(file){

    list <-
      file %>%
      readRDS()

    list %>%
      purrr::discard(is.null) %>%
      qdapTools::list_df2df() %>%
      dplyr::select(-X1)

  }) %>%
  qdapTools::list_df2df(col1 = "list_name") %>%
  dplyr::filter(
    locus %in% gene_filtered_loci$locus,
    p < bivar_threshold
  )

# Main --------------------------------------------------------------------

print(str_c(Sys.time(), " - Starting LAVA analysis for ", nrow(loci), " genes"))

# LAVA looks for common SNPs across all phenotypes when processing input
# Thus, must load separately for each gene/eqtl dataset otherwise no common SNPs will be found
foreach::foreach(
  i = 1:nrow(loci),
  .verbose = TRUE,
  .packages = c("LAVA", "tidyverse", "stringr")
) %dopar% {

  gene <- loci$LOC[i]

  # Filter results by gene/locus to narrow phenotypes used
  results_filtered <-
    results %>%
    dplyr::filter(
      locus == c(gene_filtered_loci %>%
                   dplyr::filter(
                     gene_id == gene
                   ) %>%
                   .[["locus"]])
    )

  gwas <-
    c(results_filtered$phen1, results_filtered$phen2) %>% unique()

  for(j in 1:length(eqtls)){

    eqtl_gene <- str_c(eqtls[j], "_", gene)

    if(!eqtl_gene %in% input_info$phenotype) next

    # Update args
    args$phenotypes <-
      c(
        gwas,
        eqtl_gene
      )

    args$output_filename <-
      str_c(eqtls[j], ":", gene)

    # Load input
    input <-
      LAVA::process.input(
        input.info.file = args$info_file,
        sample.overlap.file = args$sample_overlap_file,
        ref.prefix = args$ref_prefix,
        phenos = args$phenotypes
      )

    # Process locus
    locus <-
      LAVA::process.locus(
        loci[i,],
        input
      )

    # It is possible that the locus cannot be defined for various reasons (e.g. too few SNPs),
    # The !is.null(locus) check is necessary before calling the analysis functions.
    if (!is.null(locus)) {

      # extract some general locus info for the output
      loc_info <-
        data.frame(
          locus = locus$id,
          chr = locus$chr,
          start = locus$start,
          stop = locus$stop,
          n_snps = locus$n.snps,
          n_pcs = locus$K
        )

      # Create empty list
      loc_out <-
        setNames(
          vector(mode = "list", length = 2),
          nm = c("univ", "bivar")
        )

      # Run the univariate and bivariate tests seperately
      # Only run bivariate if p-value for eqtl_gene is < univar_threshold
      loc_out[["univ"]] <-
        LAVA::run.univ(locus)

      if(loc_out[["univ"]] %>%
         dplyr::filter(phen == eqtl_gene) %>%
         .[["p"]] < univar_threshold){

        loc_out[["bivar"]] <-
          LAVA::run.bivar(locus)

      }

      # Bind locus info
      loc_out <-
        loc_out %>%
        lapply(., function(df){

          if(!is.null(df)){
            loc_info %>% dplyr::bind_cols(df)
          }

        })

      # Save files
      for(k in 1:length(loc_out)){

        file_type <- names(loc_out)[k]

        if(!is.null(loc_out[[k]])){

          saveRDS(
            loc_out[[k]],
            file = file.path(out_dir[file_type], str_c(args$output_filename, ".", file_type, ".lava.rds"))
          )

        }

      }

    }

  }

}

# Stop cluster
parallel::stopCluster(cl)

print(str_c(Sys.time(), " - Analysis done!"))

