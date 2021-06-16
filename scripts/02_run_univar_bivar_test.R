# Load packages -----------------------------------------------------------

library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

phenotypes <- c("AD2019", "LBD2020", "PD2019.meta5.ex23andMe", "BIP2021", "MDD2019", "SCZ2018") %>% sort()

args <- 
  list(
    ref_prefix = "/tools/MAGMA/reference_files/g1000_eur/g1000_eur",
    loc_file = here::here("results", "01_input_prep", "gwas_filtered.loci"),
    info_file = here::here("results", "01_input_prep", "input.info.txt"),
    sample_overlap_file = here::here("results", "01_input_prep", "sample_overlap.txt"),
    phenotypes = phenotypes,
    output_filename = str_c(phenotypes, collapse = ":")
  )

print(args)

# Load data ---------------------------------------------------------------

loci <- LAVA::read.loci(args$loc_file)
n_loci <- nrow(loci)
input <- 
  LAVA::process.input(
    input.info.file = args$info_file,
    sample.overlap.file = args$sample_overlap_file,
    ref.prefix = args$ref_prefix,
    phenos = args$phenotypes
  )

# Main --------------------------------------------------------------------

# Print progress
print(str_c("Starting LAVA analysis for ", n_loci, " loci"))
progress <- 
  quantile(
    x = 1:n_loci, 
    probs = seq(.05,1,.05)
    ) %>% 
  ceiling()

# Set univariate threshold to 0.05/n_loci
univar_threshold <- 
  0.05/n_loci

univar = bivar = list()

for (i in 1:n_loci) {
  
  if (i %in% progress) print(str_c("..", names(progress[which(progress==i)])))     # (printing progress)
  
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
    
    # Run the univariate and bivariate tests
    loc_out <- 
      LAVA::run.univ.bivar(
        locus, 
        univ.thresh = univar_threshold
      )
    
    # Bind 
    univar[[i]] <- 
      loc_info %>% 
      dplyr::bind_cols(loc_out$univ)
    
    if(!is.null(loc_out$bivar)){
      
      bivar[[i]] <- 
        loc_info %>% 
        dplyr::bind_cols(loc_out$bivar)

      } 
  
    }

  }

# Save data ---------------------------------------------------------------

out_dir <- here::here("results", "02_univar_bivar_test")
dir.create(out_dir, showWarnings = T)
saveRDS(
  univar, 
  file = file.path(out_dir, str_c(args$output_filename, ".univ.lava.rds"))
)
saveRDS(
  bivar, 
  file = file.path(out_dir, str_c(args$output_filename, ".bivar.lava.rds"))
)

print(str_c("Done! Analysis output written to ", out_dir, "/", args$output_filename,".*.lava"))
