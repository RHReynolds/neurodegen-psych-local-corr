#' Plot chord diagram with bivariate correlations.
#'
#' @description This function will plot significant bivariate local genetic
#'   correlations as a chord diagram. Connections between phenotypes will be
#'   coloured by the direction of the correlation i.e. positive correlations =
#'   red, negative correlations = blue.
#'
#' @param bivar_corr a `data.frame` or [tibble][tibble::tbl_df-class] object,
#'   with the following columns (all of which are output by LAVA's bivariate test):
#'  \itemize{
#'  \item `phen1`: name of phenotype 1
#'  \item `phen2`: name of phenotype 2
#'  \item `rho`: the estimated genetic correlation
#'  \item `p`: p-value of the genetic correlation
#'  }
#' @param p_threshold `numeric` vector indicating p-value threshold to filter
#'   results by. Default is NULL.
#' @param fct_phen `character` vector indicating order in which phenotypes
#'   should be displayed in plot.
#' @param palette `character` vector indicating palette to be used. Default is
#'   NULL.
#'
#' @return chord diagram displaying the bivariate local genetic correlations
#'   between phenotypes.
#' \itemize{
#' \item Outer circle displays phenotypes, with correlations between phenotypes
#' indicated by a line connecting the two
#' \item Negative and positive
#' correlations are indicated by blue and red, respectively.
#'   }
#' @export
#'
#' @references \itemize{
#' \item Gu et al. (2014) circlize Implements and enhances circular
#' visualization in R \emph{Bioinformatics}, 2014 Oct;30(19):2811-2.
#' \url{https://pubmed.ncbi.nlm.nih.gov/24930139/} PMID: 24930139 }

plot_bivar_chord_diagram <-
  function(
    bivar_corr,
    p_threshold = NULL,
    fct_phen,
    palette = NULL
  ) {

    # If p-value threshold set, filter by threshold
    if (!is.null(p_threshold)) {
      bivar_corr <-
        bivar_corr %>%
        dplyr::filter(
          p < p_threshold
        )
    }

    # Number of unique phenotypes
    n_phenotype <- c(bivar_corr$phen1, bivar_corr$phen2) %>% unique() %>% length()

    # Set palette to number of unique phenotypes
    if (is.null(palette)) {
      palette <-
        scales::hue_pal()(n_phenotype)
    }

    data_to_plot <-
      bivar_corr %>%
      dplyr::mutate(
        phen1 = phen1 %>%
          forcats::fct_relevel(
            fct_phen
          ),
        phen2 = phen2 %>%
          forcats::fct_relevel(
            fct_phen
          )
      ) %>%
      dplyr::arrange(phen1, rho) %>%
      dplyr::select(phen1, phen2, rho)

    circlize::circos.par(
      gap.after =
        stats::setNames(
          object = rep(x = 3, n_phenotype),
          nm = fct_phen
        )
    )

    data_to_plot %>%
      dplyr::select(phen1, phen2) %>%
      dplyr::mutate(n = 1) %>%
      dplyr::rename(
        from = phen1,
        to = phen2,
        value = n
      ) %>%
      circlize::chordDiagram(
        list(
          track.height =
            c(
              rownames(data_to_plot),
              colnames(data_to_plot)
            ) %>%
            length()
        ),
        annotationTrack = "grid",
        preAllocateTracks = 1,
        grid.col = palette,
        col = ifelse(data_to_plot$rho > 0, "#de2d26", "#3182bd")
      )

    circlize::circos.track(
      track.index = 1, panel.fun = function(x, y) { # Add text labels
        circlize::circos.axis(
          h = "top",
          major.at = seq(from = 0, to = nrow(data_to_plot), by = 1),
          labels = FALSE,
          minor.ticks = 0,
          major.tick.percentage = 1,
          track.index = 2
        )
        circlize::circos.text(CELL_META$xcenter,
                              CELL_META$ylim[1] + mm_y(5),
                              CELL_META$sector.index,
                              facing = "clockwise", # Direction of text
                              niceFacing = TRUE,
                              adj = c(0, 0.5),
                              cex = 0.7
        )

      },
      bg.border = NA
    )
  }
