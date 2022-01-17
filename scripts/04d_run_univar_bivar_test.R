# Description: run univariate and bivariate tests for eQTL/GWAS traits

# Load packages -----------------------------------------------------------

library(doSNOW)
library(foreach)
library(gtools)
library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

eqtls <- c("EQTLGEN", "PSYCHENCODE") %>% sort()

args <-
  list(
    cores = 10,
    ref_prefix = "/tools/MAGMA/reference_files/g1000_eur/g1000_eur",
    loc_file =
      setNames(
        list.files(
          here::here("results", "04_eqtl_univar_bivar"),
          recursive = T,
          full.names = T,
          pattern = ".loci$"
        ),
        nm =
          list.files(
            here::here("results", "04_eqtl_univar_bivar"),
            recursive = T,
            pattern = ".loci$"
          ) %>%
          str_remove("/.*")
      ),
    info_file = here::here("results", "04_eqtl_univar_bivar", "input.info.txt"),
    sample_overlap_file = here::here("results", "04_eqtl_univar_bivar", "sample_overlap.txt"),
    phenotypes = c(""),
    output_filename = c(""),
    path_to_log = here::here("logs", "04d_run_univar_bivar_test_multwindows_foreach.log"),
    # Bivariate threshold from previous analyses (02_run_univar_bivar_test.rmd)
    # Used to filter bivariate results to include only phenos with significant local rg
    # These phenos will then be used together with the eQTL
    bivar_threshold = 0.05/1603, # n bivariate tests
    # Will be set to 0.05/n_loci
    univar_threshold = c()
  )

out_dir <-
  setNames(
    names(args$loc_file),
    names(args$loc_file)
  ) %>%
  lapply(., function(x){

    setNames(
      c(
        here::here("results", "04_eqtl_univar_bivar", x, "univ"),
        here::here("results", "04_eqtl_univar_bivar", x, "bivar")
      ),
      nm = c("univ", "bivar")
    )

  })

# Set up for parallel run -------------------------------------------------

# Create out directories
out_dir %>%
  lapply(., function(list){

    list %>%
      lapply(., function(dir_path) dir.create(dir_path, showWarnings = T))

  })

# Run in parallel
cl <- parallel::makeCluster(args$cores, outfile = args$path_to_log)

# Register clusters
doSNOW::registerDoSNOW(cl)

# Load common files -------------------------------------------------------

input_info <-
  read_delim(
    args$info_file,
    delim = "\t"
  )

# Main --------------------------------------------------------------------

for(window in names(args$loc_file)){

  loci <-
    read.loci(args$loc_file[[window]]) %>%
    dplyr::arrange(LOC)

  gene_filtered_loci <-
    readRDS(
      here::here("results", "04_eqtl_univar_bivar", window, "gene_filtered_loci.rds")
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
      p < args$bivar_threshold
    )

  # Update univariate threshold
  args$univar_threshold <- 0.05/nrow(loci)

  print(
    str_c(Sys.time(),
          " - Starting LAVA analysis for ",
          nrow(loci),
          " genes using a window size of ",
          str_remove(window, ".*_"),
          " bp")
  )

  print(str_c(Sys.time(), " - Log file of foreach parallelisation outputted here: ", args$path_to_log))

  foreach::foreach(
    i = 1:nrow(loci),
    .verbose = TRUE,
    .packages = c("LAVA", "tidyverse", "stringr"),
    .errorhandling = "pass",
    .combine = "list"
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

      if(!eqtl_gene %in% input_info$phenotype){

        cat("Chunk", i," -- no gene-eqtl for", eqtl_gene, "\n")
        next

      }

      cat("Chunk", i, "-- processing locus for", eqtl_gene,"\n")

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

        cat("Chunk", i, "-- running univariate for", eqtl_gene,"\n")

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
        # Only run bivariate if eQTL p-value is not NULL/NA and is < univar_threshold
        loc_out[["univ"]] <-
          LAVA::run.univ(locus)

        eqtl_univ_p <-
          loc_out[["univ"]] %>%
          dplyr::filter(phen == eqtl_gene) %>%
          .[["p"]]

        if(!gtools::invalid(eqtl_univ_p) &&
           is.numeric(eqtl_univ_p) &&
           eqtl_univ_p < args$univar_threshold){

          cat("Chunk", i, " -- running bivariate for", eqtl_gene,"\n")

          loc_out[["bivar"]] <-
            LAVA::run.bivar(
              locus,
              # Subset for other phenotypes with univ p < univar_threshold
              phenos =
                loc_out[["univ"]] %>%
                dplyr::filter(
                  p < args$univar_threshold
                ) %>%
                .[["phen"]] %>%
                as.character()
            )

        } else{

          cat("Chunk ", i, " -- no bivariate test run for ", eqtl_gene, ". Univariate p-value =", eqtl_univ_p, "\n", sep = "")

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
              file = file.path(out_dir[[window]][file_type], str_c(args$output_filename, ".", file_type, ".lava.rds"))
            )

          }

        }

      } else {

        cat("Chunk", i, "-- after processing locus is null for", eqtl_gene, "\n")

      }

    }

  }


}

# Stop cluster
parallel::stopCluster(cl)

print(str_c(Sys.time(), " - Analysis done!"))





