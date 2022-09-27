# Description: get h2 for GWAS traits

# Load packages -----------------------------------------------------------

library(janitor)
library(tidyverse)
library(stringr)
library(qdapTools)

# Set arguments -----------------------------------------------------------

gwas_dir <- "/data/LDScore/GWAS/"

phenotypes <- c("AD2019", "AD2019.Kunkle", "LBD2020", "PD2019.meta5.ex23andMe",  "PD2019.ex23andMe.exUKBB", "BIP2021", "MDD2019", "SCZ2018") %>% sort()

out_dir <- here::here("results", "00_trait_h2")
dir.create(out_dir, showWarnings = T)

ldsc_args <-
  list(
    ldsc = "/tools/ldsc/ldsc.py",
    ref_ld = "/data/LDScore/Reference_Files/eur_w_ld_chr/"
  )

# Load files paths --------------------------------------------------------

gwas_details <-
  tibble(
    full_path = list.files(path = gwas_dir,
                           recursive = T,
                           full.names = T,
                           pattern = ".sumstats.gz")) %>%
  dplyr::mutate(
    file_name = basename(full_path) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_replace_all("_", "\\.")) %>%
  dplyr::filter(file_name %in% phenotypes)

# Main --------------------------------------------------------------------

# Run loop across rows, with external call to ldsc
for(i in 1:nrow(gwas_details)){

  args <-
    list(
      phen = gwas_details$full_path[i],
      out_prefix =
          gwas_details$file_name[i]
    )

  h2_arg <-
    stringr::str_c(
      " ", ldsc_args$ldsc,
      " --h2 ", args$phen,
      " --ref-ld-chr ", ldsc_args$ref_ld,
      " --w-ld-chr ", ldsc_args$ref_ld,
      " --out ", out_dir, "/", args$out_prefix, "_h2"
    )

  # Run command
  system2(command = "python", args = h2_arg)

}

###### Extracting LDSC results ######

# Extract all outputs into single file
file_paths <-
  list.files(
    out_dir,
    pattern = "_h2.log",
    full.names = T
  )

files <-
  vector(mode = "list",
         length = length(file_paths))

for(i in 1:length(file_paths)){

  files[[i]] <-
    read_lines(
      file = file_paths[i],
      skip = 25, # Number of lines to skip in log file
      n_max = 4 # Number of lines to read
    ) %>%
    stringr::str_split(": ") %>%
    lapply(., function(vector){

      if(str_detect(vector[2], "\\(")){

        tibble(
          name =
            c(
              vector[1],
              stringr::str_c(vector[1], " se")
            ),
          value = vector[2] %>%
          stringr::str_remove(".*: ") %>%
          stringr::str_split(" ") %>%
          unlist() %>%
          readr::parse_number()
          )

      } else{

        tibble(
          name = vector[1],
          value = vector[2] %>%
            readr::parse_number()
        )

      }


    }) %>%
    qdapTools::list_df2df(col1 = "list_name") %>%
    as_tibble() %>%
    bind_cols(
      phen = basename(file_paths[i]) %>% stringr::str_remove("_h2.log")
    ) %>%
    dplyr::select(-list_name)

}

h2 <-
  files %>%
  qdapTools::list_df2df(col1 = "list_name") %>%
  dplyr::select(-list_name) %>%
  tidyr::pivot_wider(
    names_from = name,
    values_from = value
  ) %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    z = total_observed_scale_h2/total_observed_scale_h2_se
  )

# Save data ---------------------------------------------------------------

write_delim(
  h2,
  delim = "\t",
  file = file.path(out_dir, "h2.txt")
)


