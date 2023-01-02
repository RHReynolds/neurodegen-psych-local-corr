# Description: plot p-value distributions

# Load packages -----------------------------------------------------------

library(here)
library(tidyverse)

# Set arguments -----------------------------------------------------------

args <-
  list(
    tested = here::here("results", "02_univar_bivar_test", "ld_blocks_p.rds"),
    untested = here::here("results", "02_univar_bivar_test", "ld_blocks_p_min_50.rds"),
    n_snps = 50
  )

# Load data ---------------------------------------------------------------

tested <-
  readRDS(args$tested)

untested <-
  readRDS(args$untested)

# Main --------------------------------------------------------------------

tested <-
  tested %>%
  lapply(., function(list){
    list %>%
      lapply(., function(df){
        df %>%
          dplyr::slice_min(n = args$n_snps, order_by = P)
      }) %>%
      qdapTools::list_df2df(col1 = "phen")
  }) %>%
  qdapTools::list_df2df(col1 = "locus")  %>%
  dplyr::mutate(
    logP = -log10(P),
    type = "tested"
  )

untested <-
  untested %>%
  lapply(., function(list){
    list %>%
      lapply(., function(df){
        df %>%
          dplyr::slice_min(n = args$n_snps, order_by = P)
      }) %>%
      qdapTools::list_df2df(col1 = "phen")
  }) %>%
  qdapTools::list_df2df(col1 = "locus")  %>%
  dplyr::mutate(
    logP = -log10(P),
    type = "untested"
  )

combined <-
  tested %>%
  bind_rows(untested)

med_df <-
  combined %>%
  dplyr::group_by(phen, type) %>%
  dplyr::summarise(median = median(logP))

combined %>%
  ggplot(aes(
    x = logP,
    fill = type
  )) +
  geom_density(alpha = 0.5)  +
  geom_vline(
    data = med_df,
    aes(
      xintercept = median,
      colour = type
    )) +
  facet_wrap(vars(phen))
