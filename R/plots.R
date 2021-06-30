# Set theme
theme_rhr <-
  ggplot2::theme_bw(
    base_family = "Helvetica",
    base_size = 10
  ) +
  ggplot2::theme(
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(vjust = 0.6),
    panel.spacing = unit(0.1, "lines")
  )

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
           n_phenotypes) {

    # Only need first half of matrix, thus must extract appropriate rows from dataframe
    n <- n_phenotypes

    for (i in 1:n) {

      # General formula for extracting appropriate indices
      index <- (i * n - (n - i)):(i * n)

      if (i == 1) {
        indices <- index
      } else {
        indices <- c(indices, index)
      }
    }

    # Determine number of combinations for multiple test correction
    n_combn <-
      length(indices) - n_phenotypes


    global_corr %>%
      dplyr::mutate(
        fdr = p.adjust(p, method = "fdr"),
        rg_fill =
          dplyr::case_when(
            p < 0.05 / n_combn ~ round(rg, 2)
          ),
        p1 = p1 %>%
          stringr::str_replace_all("[:digit:]", "") %>%
          stringr::str_remove("\\..*"),
        p2 = p2 %>%
          stringr::str_replace_all("[:digit:]", "") %>%
          stringr::str_remove("\\..*")
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
      ggplot2::guides(colour = F) +
      theme_rhr +
      ggplot2::theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank()
      )
  }

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
        sector.name <- circlize::get.cell.meta.data("sector.index")
        circlize::circos.text(CELL_META$xcenter,
                              CELL_META$ylim[1],
                              CELL_META$sector.index,
                              facing = "clockwise", # Direction of text
                              niceFacing = TRUE,
                              adj = c(-0.1, 0.5),
                              cex = 1
        )
        circlize::circos.axis(
          h = "top",
          major.at = seq(from = 0, to = nrow(data_to_plot), by = 1),
          labels = FALSE,
          labels.cex = 0.5,
          minor.ticks = 0,
          sector.index = sector.name,
          major.tick.percentage = 1,
          labels.niceFacing = TRUE,
          track.index = 2
        )
      },
      bg.border = NA
    )
  }

#' Plot edge diagrams for bivariate correlations
#'
#' @description This function will plot significant bivariate local genetic
#'   correlations as an edge diagram for each locus. Connections between
#'   phenotypes will be coloured by the direction of the correlation i.e.
#'   positive correlations = red, negative correlations = blue.
#'
#' @param bivar_corr a `data.frame` or [tibble][tibble::tbl_df-class] object,
#'   with the following columns (all of which are output by LAVA's bivariate test):
#'  \itemize{
#'  \item `phen1`: name of phenotype 1
#'  \item `phen2`: name of phenotype 2
#'  \item `rho`: the estimated genetic correlation
#'  \item `rho.lower`: lower 95% confidence estimate for rho
#'  \item `rho.upper`: upper 95% confidence estimate for rho
#'  \item `p`: p-value of the genetic correlation
#'  }
#' @param p_threshold `numeric` vector indicating p-value threshold to filter
#'   results by. Default is NULL.
#' @param phen `character` vector indicating phenotypes present.
#' @param ncol `integer` vector indicating number of columns in facet.
#'
#' @return edge diagram displaying the bivariate local genetic correlations
#'   between phenotypes at each locus with more than one bivariate correlation.
#'   Negative and positive´correlations are indicated by blue and red,
#'   respectively. }
#' @export
#'

plot_edge_diagram <-
  function(
    bivar_corr,
    p_threshold = NULL,
    phen,
    ncol = 3
  ) {

    # If p-value threshold set, filter by threshold
    if (!is.null(p_threshold)) {
      bivar_corr <-
        bivar_corr %>%
        dplyr::filter(
          p < p_threshold
        )
    }

    # Filter for loci with more than one bivariate correlation
    bivar_corr <-
      bivar_corr %>%
      dplyr::group_by(locus) %>%
      dplyr::filter(n() > 1)

    edges <-
      bivar_corr %>%
      # Add confidence intervals in case want to use these in the figure
      dplyr::mutate(
        display_rho =
          sprintf(
            "%.2f [%.2f, %.2f]" ,
            rho,
            rho.lower,
            rho.upper
          )
      ) %>%
      dplyr::select(
        contains("phen"), rho, display_rho, locus
      )

    nodes <-
      tibble::tibble(
        name =
          factor(
            phen,
            levels = phen),
        id = c(1:length(phen))
      )

    edge_tbl_graph <-
      tidygraph::as_tbl_graph(
        x = edges,
        nodes = nodes,
        directed = T
      )

    # Set seed so that same graph drawn everytime
    set.seed(89)

    ggraph::ggraph(
      edge_tbl_graph,
      layout = "igraph",
      algorithm = "graphopt"
    ) +
      ggraph::geom_edge_link(
        ggplot2::aes(
          colour = rho,
          label = round(rho, 2)
        ),
        angle_calc = 'along',
        label_dodge = unit(1.5, 'mm'),
        label_size = 2,
        width = 1
      ) +
      ggraph::geom_node_label(
        ggplot2::aes(
          label = name
        ),
        force = 0,
        repel = T,
        color = 'black',
        fill = 'white',
        size = 3
      ) +
      ggplot2::scale_y_reverse() +
      ggraph::scale_edge_colour_distiller(
        palette = "RdYlBu",
        limits = c(-1,1),
        direction = -1
      ) +
      ggraph::facet_edges(
        vars(locus),
        ncol = ncol
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top"
      )

  }

#' Plot genes within a LAVA LD block.
#'
#' This function will plot all genes of the gene biotypes "protein_coding",
#' "antisense", and "lincRNA" within the start and stop of an LD block.
#'
#' @param locus_gr gr. Granges object with coordinates of each locus, together
#'   with at least 1 additional metadata column:
#'   \itemize{
#'   \item locus - locus identifier.
#'   }
#' @param ref gr. Granges object containing reference ensembl gtf.
#'
#' @return ggplot2 plot of genes (of biotype "protein_coding", "antisense", and
#'   "lincRNA") within locus.
#' @export
#'
#'

plot_locus <-
  function(
           locus_gr,
           ref) {

    overlap <-
      GenomicRanges::findOverlaps(locus_gr, ref) %>%
      tibble::as_tibble()

    coords_to_plot <-
      tibble::tibble(
        group_type = as.factor("genes"),
        gene_id = ref[overlap$subjectHits]$gene_id,
        gene_name =
          ref[overlap$subjectHits]$gene_name %>%
            stringr::str_replace_all("-", "_") %>%
            stringr::str_replace_all("\\.", "_"),
        gene_biotype = ref[overlap$subjectHits]$gene_biotype %>% as.factor(),
        chr = ref[overlap$subjectHits] %>% GenomeInfoDb::seqnames() %>% as.character(),
        strand = ref[overlap$subjectHits] %>% BiocGenerics::strand() %>% as.character(),
        start = ref[overlap$subjectHits] %>% BiocGenerics::start(),
        end = ref[overlap$subjectHits] %>% BiocGenerics::end(),
        width = ref[overlap$subjectHits] %>% BiocGenerics::width()
      ) %>%
      # If there are two versions of a gene, choose the longest
      dplyr::group_by(gene_name) %>%
      dplyr::top_n(1, wt = width) %>%
      dplyr::filter(
        gene_biotype %in% c("protein_coding", "antisense", "lincRNA")
      )

    # Add locus coordinates
    coords_to_plot <-
      coords_to_plot %>%
      dplyr::bind_rows(
        tibble(
          group_type = as.factor("locus"),
          gene_name = stringr::str_c("Locus_", locus_gr$locus),
          chr = locus_gr %>% GenomeInfoDb::seqnames() %>% as.character(),
          start = locus_gr %>% BiocGenerics::start(),
          end = locus_gr %>% BiocGenerics::end(),
          gene_biotype = as.factor("locus")
        )
      )

    # Re-factor groups
    coords_to_plot <-
      coords_to_plot %>%
      dplyr::mutate(
        gene_biotype =
          gene_biotype %>%
            fct_relevel(
              .,
              c("locus", "protein_coding", "antisense", "lincRNA")
            )
      )

    ggplot2::ggplot(
      data = coords_to_plot
    ) +
      ggplot2::geom_linerange(
        ggplot2::aes(
          x = start,
          y = forcats::fct_reorder(
            gene_name,
            -start
          ),
          xmin = start,
          xmax = end
          # colour = gene_biotype
        ),
        size = 1
      ) +
      ggplot2::labs(title = stringr::str_c("Locus: ", locus_gr$locus)) +
      ggplot2::scale_fill_manual(
        values = c(
          "locus" = "#868686FF",
          "protein_coding" = "#0073C2FF",
          "antisense" = "#EFC000FF",
          "lincRNA" = "#CD534CFF"
        )
      ) +
      ggplot2::scale_x_continuous(labels = scales::comma) +
      ggplot2::labs(
        x = stringr::str_c("Chromosome ", c(coords_to_plot$chr %>% unique()))
      ) +
      ggplot2::facet_grid(
        rows = vars(gene_biotype),
        scales = "free_y",
        space = "free_y"
      ) +
      ggplot2::theme_bw(
        base_size = 8
      ) +
      ggplot2::theme(
        strip.text.y = element_text(angle = 0)
      )
  }
