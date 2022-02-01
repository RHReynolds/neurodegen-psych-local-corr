# Description: generate N column for PD2019 summary statistics without 23andMe and UKBB

# Load packages -----------------------------------------------------------

library(data.table)
library(janitor)
library(stringr)
library(tidyverse)
library(readxl)

# Arguments ---------------------------------------------------------------

args <-
  list(
    file_path_pd = "/data/LDScore/GWAS/PD2019_meta5_ex23andMe_exUKBB/META_no_UKBB_no_231.tbl.gz",
    file_path_cohort = here::here("raw_data", "05_noproxy_check", "Table_S1_Study_descriptive_statistics.xlsx"),
    out_path = "/data/LDScore/GWAS/PD2019_meta5_ex23andMe_exUKBB/META_no_UKBB_no_231_N.tbl.gz"
  )

print(args)

# Load data ---------------------------------------------------------------

# Load cohort sample sizes
cohort_n <-
  readxl::read_excel(
    path = args$file_path_cohort
  )

# Main --------------------------------------------------------------------

# Tidy data and add ordering
# --> Input File 1 : ../individual_cohort_sumstats_meta5/toMeta.dbgapNeuroX.tab
# --> Input File 2 : ../individual_cohort_sumstats_meta5/toMeta.SHUL.tab
# --> Input File 3 : ../individual_cohort_sumstats_meta5/toMeta.FINLAND_no_age.tab
# --> Input File 4 : ../individual_cohort_sumstats_meta5/toMeta.HBS.tab
# --> Input File 5 : ../individual_cohort_sumstats_meta5/toMeta.GILL_PD_C.tab
# --> Input File 6 : ../individual_cohort_sumstats_meta5/toMeta.OSLO.tab
# --> Input File 7 : ../individual_cohort_sumstats_meta5/toMeta.PDBP.tab
# --> Input File 8 : ../individual_cohort_sumstats_meta5/toMeta.PPMI.tab
# --> Input File 9 : ../individual_cohort_sumstats_meta5/toMeta.queensland.tab
# --> Input File 10 : ../individual_cohort_sumstats_meta5/toMeta.SPAIN3.tab
# --> Input File 11 : ../individual_cohort_sumstats_meta5/toMeta.TUBI_no_overlap.tab
# --> Input File 12 : ../individual_cohort_sumstats_meta5/toMeta.VANCE.tab
# --> Input File 13 : ../individual_cohort_sumstats_meta5/toMeta.COURAGE_UK.tab
cohort_n <-
  cohort_n %>%
  janitor::clean_names() %>%
  dplyr::filter(!is.na(study) & study != "Total" & included == "X") %>%
  dplyr::mutate(
    order =
      case_when(
        str_detect(study, "NeuroX") ~ 1,
        str_detect(study, "Baylor") ~ 2,
        str_detect(study, "Finnish") ~ 3,
        str_detect(study, "HBS") ~ 4,
        str_detect(study, "McGill") ~ 5,
        str_detect(study, "Oslo") ~ 6,
        str_detect(study, "PDBP") ~ 7,
        str_detect(study, "PPMI") ~ 8,
        str_detect(study, "SGPD") ~ 9,
        str_detect(study, "Spanish") ~ 10,
        str_detect(study, "Tubingen") ~ 11,
        str_detect(study, "Vance") ~ 12,
        str_detect(study, "UK PDMED") ~ 13
      )
  ) %>%
  dplyr::arrange(order)

print(stringr::str_c(Sys.time(), " - generating N column"))

# Generate N column
pd_n <-
  fread(args$file_path_pd) %>%
  dplyr::select(MarkerName, Direction) %>%
  dplyr::mutate(
    Direction =
      stringr::str_split(Direction, pattern = "")
  ) %>%
  tidyr::unnest_wider(
    col = Direction,
    names_sep = "_"
  ) %>%
  tidyr::pivot_longer(
    cols = Direction_1:Direction_13,
    names_to = "order"
  ) %>%
  dplyr::filter(
    value %in% c("+", "-")
  ) %>%
  dplyr::inner_join(
    cohort_n %>%
      dplyr::mutate(
        order =
          stringr::str_c("Direction_", order)
      ),
    by = c("order")
  ) %>%
  dplyr::group_by(MarkerName) %>%
  dplyr::summarise(
    N_cases = sum(cases_n),
    N_controls = sum(controls_n),
    N = sum(total_n)
  )

print(stringr::str_c(Sys.time(), " - N column generated"))

# Save data ---------------------------------------------------------------

data.table::fwrite(
  pd_n,
  args$out_path,
  sep = "\t"
)

print(stringr::str_c(Sys.time(), " - Done!"))
