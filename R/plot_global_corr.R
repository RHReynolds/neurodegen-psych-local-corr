source(here::here("R", "theme_rhr.R"))

#' Plot global correlations
#'
#' @description This function will plot global genetic correlations between
#'   phenotypes (as determined using LDSC) as a lower triangle heatmap.
#'   Significant correlations will be determined with multiple test corrections
#'   applied (Bonferroni). Number of tests = to number of unique combinations
#'   between phenotypes (excluding comparisons between the same phenotype).
#'
#' @param global_corr a `data.frame` or [tibble][tibble::tbl_df-class] object,
#'   with the following columns (all of which are output by LDSC genetic
#'   correlation analyses):
#'  \itemize{
#'  \item `p1`: name of phenotype 1
#'  \item `p2`: name of phenotype 2
#'  \item `rg`: the estimated genetic correlation
#'  \item `p`: p-value of the genetic correlation
#'  }
#' @param n_phenotypes `integer` vector indicating number of phenotypes run
#' @param tidy_label `logical` vector indicating whether phenotype labels should
#'   be "tidied" i.e. remove digits and anything after `.`. Default is TRUE.
#'
#' @return `ggplot` displaying the genetic correlations between phenotypes.
#' \itemize{
#' \item x and y-axis display the phenotypes.
#' \item Significant negative and positive correlations are indicated by blue
#' and red fill, respectively. Non-significant correlations (p >= 0.05/n_tests)
#' have a grey fill.
#' }
#'
#' @export
#'
#' @references \itemize{ \item Bulik-Sullivan et al. (2015) An atlas of genetic
#'   correlations across human diseases and traits \emph{Nature Genetics}, 2015
#'   Nov;47(11):1236-41. \url{https://www.nature.com/articles/ng.3406} PMID:
#'   26414676 }

plot_global_corr <-
  function(
    global_corr,
    n_phenotypes,
    tidy_label = TRUE
  ) {

    # Only need first half of matrix, thus must extract appropriate rows from dataframe
    indices <-
      .generate_global_corr_indices(n_phenotypes = n_phenotypes)

    # Determine number of combinations for multiple test correction
    n_combn <-
      length(indices) - n_phenotypes

    if(tidy_label == TRUE){

      global_corr <-
        global_corr %>%
        dplyr::mutate(
          p1 = p1 %>%
            stringr::str_replace_all("[:digit:]", "") %>%
            stringr::str_remove("\\..*"),
          p2 = p2 %>%
            stringr::str_replace_all("[:digit:]", "") %>%
            stringr::str_remove("\\..*")
        )

    }

    global_corr %>%
      dplyr::mutate(
        rg_fill =
          dplyr::case_when(
            p < 0.05 / n_combn ~ round(rg, 2)
          )
      ) %>%
      dplyr::slice(indices) %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = p1,
          y = p2 %>% fct_rev(),
          fill = rg_fill,
          label = round(rg, 2)
        )
      ) +
      ggplot2::geom_tile(colour = "black") +
      ggplot2::geom_text(
        size = 3
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(x = "", y = "", fill = "Global correlation (rg)") +
      ggplot2::scale_fill_distiller(
        palette = "RdYlBu",
        limits = c(-1, 1)
        ) +
      ggplot2::scale_colour_manual(values = c("black", "white")) +
      ggplot2::guides(colour = "none") +
      ggplot2::theme_bw(
        base_family = "Helvetica",
        base_size = 10
      ) +
      ggplot2::theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(vjust = 0.6),
        panel.spacing = unit(0.1, "lines")
      )

  }

#' Plot global and local correlations as heatmap
#'
#' @description This function will plot global genetic correlations between
#'   phenotypes (as determined using LDSC) as a lower triangle heatmap, and
#'   local genetic correlations as an upper triangle heatmap. Significant global
#'   correlations will be determined with multiple test corrections applied
#'   (Bonferroni), and indicated with an asterisk. Number of tests = to number
#'   of unique combinations between phenotypes (excluding comparisons between
#'   the same phenotype). Significant local correlations will be determined with
#'   multiple test corrections applied (Bonferroni), and the number of
#'   significant local genetic correlations will be indicated within each tile.
#'
#' @param global_corr a `data.frame` or [tibble][tibble::tbl_df-class] object,
#'   with the following columns (all of which are output by LDSC genetic
#'   correlation analyses):
#'  \itemize{
#'  \item `p1`: name of phenotype 1
#'  \item `p2`: name of phenotype 2
#'  \item `rg`: the estimated genetic correlation
#'  \item `p`: p-value of the genetic correlation
#'  }
#' @param bivar_corr a `data.frame` or [tibble][tibble::tbl_df-class] object,
#'   with the following columns (all of which are output by LAVA's bivariate test):
#'  \itemize{
#'  \item `phen1`: name of phenotype 1
#'  \item `phen2`: name of phenotype 2
#'  \item `rho`: the estimated genetic correlation
#'  \item `p`: p-value of the genetic correlation
#'  }
#' @param n_phenotypes `integer` vector indicating number of phenotypes run
#' @param tidy_label `logical` vector indicating whether phenotype labels should
#'   be "tidied" i.e. remove digits and anything after `.`. Default is TRUE.
#'
#' @return `ggplot` displaying the genetic correlations between phenotypes.
#' \itemize{
#' \item x and y-axis display the phenotypes.
#' \item The fill of each tile in the lower triangle indicates the global
#' genetic correlation, while in the upper triangle it represents the mean rho
#' across all tested LD blocks.
#' \item Significant global correlations are indicated with an asterisk.
#' \item The number of significant local genetic correlations will be indicated
#' within each tile.
#' }
#'
#' @export
#'
#' @references \itemize{ \item Bulik-Sullivan et al. (2015) An atlas of genetic
#'   correlations across human diseases and traits \emph{Nature Genetics}, 2015
#'   Nov;47(11):1236-41. \url{https://www.nature.com/articles/ng.3406} PMID:
#'   26414676 }

plot_global_bivar <-
  function(
    global_corr,
    n_phenotypes,
    bivar_corr
  ) {

    # Only need first half of matrix, thus must extract appropriate rows from dataframe
    indices <-
      .generate_global_corr_indices(n_phenotypes = n_phenotypes)

    # Determine number of combinations for global multiple test correction
    n_combn <-
      length(indices) - n_phenotypes

    # Prepare global correlation df
    global_corr <-
      global_corr %>%
      dplyr::select(
        p1, p2, rg, p
      ) %>%
      dplyr::mutate(
        label =
          dplyr::case_when(
            p < 0.05 / n_combn & p1 != p2 ~ "*"
          )
      ) %>%
      dplyr::slice(indices)

    # Prepare local correlation df
    # Average correlation across tested LD blocks
    bivar_corr <-
      bivar_corr %>%
      dplyr::group_by(phen1, phen2) %>%
      dplyr::summarise(
        rg = mean(rho)
      ) %>%
      dplyr::left_join(
        bivar_corr %>%
          dplyr::filter(
            p < 0.05/nrow(bivar_corr)
          ) %>%
          dplyr::group_by(phen1, phen2) %>%
          dplyr::summarise(
            label = n() %>% as.character()
          )
      ) %>%
      # Reverse phenotypes to ensure bivar appear on top half of triangle
      dplyr::rename(
        p2 = phen1,
        p1 = phen2
      )

    # Plot
    global_corr %>%
      dplyr::bind_rows(
        bivar_corr
      ) %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = p1,
          y = p2 %>% fct_rev(),
          fill = rg,
          label = label
        )
      ) +
      ggplot2::geom_tile(colour = "black") +
      ggplot2::geom_text(
        size = 4
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(
        x = "", y = "",
        fill = "Correlation (rg)",
        title = "LDSC (bottom) vs LAVA (top)"
        ) +
      ggplot2::scale_fill_distiller(
        palette = "RdBu",
        limits = c(-1, 1)
      ) +
      ggplot2::theme_bw(
        base_family = "Helvetica",
        base_size = 10
      ) +
      ggplot2::theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(vjust = 0.6),
        panel.spacing = unit(0.1, "lines")
      )

  }

#' Extract appropriate rows from global correlations dataframe
#'
#'@description This function will extract the appropriate rows from the global
#'  genetic correlations dataframe, ensuring that it will be plotted as a lower
#'  triangle heatmap.
#'
#'@param n_phenotypes `integer` vector indicating number of phenotypes run
#'
#'@return vector of indices to extract
#'

.generate_global_corr_indices <-
  function(
    n_phenotypes
  ){

    for (i in 1:n_phenotypes) {

      # General formula for extracting appropriate indices
      index <- (i * n_phenotypes - (n_phenotypes - i)):(i * n_phenotypes)

      if (i == 1) {
        indices <- index
      } else {
        indices <- c(indices, index)
      }
    }

    return(indices)

  }
