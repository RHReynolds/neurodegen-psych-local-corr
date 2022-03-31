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
#'  \item `locus`: locus id
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
#' @param locus_labels named vector indicating labels that should be used for
#'   facets. Default is to use the locus names in the supplied `data.frame`.
#' @param multiple_corr logical vector indicating whether to filter loci such
#'   that only those with more than one bivariate correlation are displayed.
#'   Default is TRUE.
#' @param geom_node_size `integer` vector indicating node size. Default is 3.
#' @param geom_label_size `integer` vector indicating edge link label size.
#'   Default is 3.
#' @param base_size `integer` vector indicating base font size in pts. Default
#'   is 11, i.e. default of \code{ggplot2::\link[ggplot2:theme_bw]{theme_bw}}.
#' @param ncol `integer` vector indicating number of columns in facet.
#' @param seed `integer` vector indicating seed to be used for generating the
#'   layout of the graph (i.e. how nodes are placed on the plot), which is
#'   random by default. Setting a seed ensures the same output each time.
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
    locus_labels = NULL,
    multiple_corr = TRUE,
    geom_node_size = 3,
    geom_label_size = 2,
    base_size = 11,
    ncol = 3,
    seed = 89
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
    if(multiple_corr == TRUE){

      bivar_corr <-
        bivar_corr %>%
        dplyr::group_by(locus) %>%
        dplyr::filter(n() > 1)

    }

    # Create locus labels vector with original names if argument set NULL
    if(is.null(locus_labels)){

      locus_labels <-
        setNames(
          unique(bivar_corr$locus),
          nm = unique(bivar_corr$locus)
        )

    }

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
    set.seed(seed)

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
        label_size = geom_label_size,
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
        size = geom_node_size
      ) +
      ggplot2::scale_y_reverse() +
      ggraph::scale_edge_colour_distiller(
        palette = "RdYlBu",
        limits = c(-1,1),
        direction = -1
      ) +
      ggraph::facet_edges(
        vars(locus),
        ncol = ncol,
        labeller = as_labeller(locus_labels)
      ) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top"
      )

  }

#' Plot edge diagrams for bivariate GWAS-eQTL correlations
#'
#' @description This function will plot significant bivariate local genetic
#'   correlations between GWAS and eQTLs as an edge diagram for each gene
#'   locus/eQTL dataset. Connections between phenotypes will be coloured by the
#'   direction of the correlation i.e. positive correlations = red, negative
#'   correlations = blue. Nodes will be coloured by the phenotype type i.e.
#'   gwas/disease trait = grey, eQTL = white.
#'
#' @param bivar_corr a `data.frame` or [tibble][tibble::tbl_df-class] object,
#'   with the following columns (most of which are output by LAVA's bivariate
#'   test):
#'   \itemize{
#'   \item `ld_block`: name of LD block in which gene locus is located
#'   \item `eqtl_dataset`: name of eQTL dataset from which gene locus is derived
#'   \item `gene_locus`: ensembl gene id (used as gene locus identifier)
#'   \item `gene_name`: HGNC gene symbol
#'   \item `phen1`: name of phenotype 1
#'   \item `phen2`: name of phenotype 2
#'   \item `rho`: the estimated genetic correlation
#'   \item `rho.lower`: lower 95% confidence estimate for rho
#'   \item `rho.upper`: upper 95% confidence estimate for rho
#'   \item `p`: p-value of the genetic correlation
#'   }
#' @param p_threshold `numeric` vector indicating p-value threshold to filter
#'   results by. Default is NULL.
#' @param phen `character` vector indicating GWAS phenotypes present.
#' @param seed `integer` vector indicating seed to be used for generating the
#'   layout of the graph (i.e. how nodes are placed on the plot), which is
#'   random by default. Setting a seed ensures the same output each time.
#'   Default is 89.
#'
#' @return list of plots with edge diagrams displaying the bivariate local
#'   genetic correlations between gwas and eQTL phenotypes at each locus with
#'   (i) more than one bivariate correlation and (ii) with a bivariate
#'   correlation that includes an eQTL. Negative and positive correlations are
#'   indicated by blue and red, respectively. GWAS and eQTL nodes are indicated
#'   by grey and white fill, respectively.
#' @export
#'

plot_qtl_edge_diagram <-
  function(
    bivar_corr_qtl,
    p_threshold = NULL,
    phen,
    seed = 89
  ) {

    # If p-value threshold set, filter by threshold
    if (!is.null(p_threshold)) {
      bivar_corr_qtl <-
        bivar_corr_qtl %>%
        dplyr::filter(
          p < p_threshold
        )
    }

    # Filter for loci where there is a bivariate correlation with an eQTL
    bivar_corr_qtl <-
      bivar_corr_qtl %>%
      dplyr::inner_join(
        bivar_corr_qtl %>%
          dplyr::filter(
            str_detect(phen2, "ENSG")
          ) %>%
          dplyr::distinct(eqtl_dataset, gene_locus),
        by = c("eqtl_dataset", "gene_locus")
      ) %>%
      dplyr::mutate(
        phen_type =
          case_when(
            str_detect(phen2, "ENSG") ~ "eQTL",
            TRUE ~ "GWAS"
          ),
        phen2 =
          case_when(
            str_detect(phen2, "ENSG") ~ gene_name,
            TRUE ~ phen2
          )
      ) %>%
      dplyr::ungroup()

    bivar_corr_list <-
      setNames(
        object =
          bivar_corr_qtl %>%
          dplyr::group_split(gene_name, eqtl_dataset),
        nm =
          bivar_corr_qtl %>%
          dplyr::arrange(gene_name, eqtl_dataset) %>%
          dplyr::mutate(
            list_name = str_c(gene_name, ":", eqtl_dataset)
          ) %>%
          .[["list_name"]] %>%
          unique()
      )

    plots <-
      bivar_corr_list %>%
      lapply(., function(x){

        edges <-
          x %>%
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
          # phen1 and phen2 (i.e. to and from) need to be first two columns
          dplyr::select(
            contains("phen"), ld_block, rho, display_rho, eqtl_dataset, gene_name
          )

        nodes <-
          tibble::tibble(
            name =
              unique(c(edges$phen1, edges$phen2))
          ) %>%
          dplyr::mutate(
            id = row_number()
          )

        edge_tbl_graph <-
          tidygraph::as_tbl_graph(
            x = edges,
            nodes = nodes,
            directed = T
          )

        fill <-
          as_tibble(
            edge_tbl_graph
          ) %>%
          dplyr::mutate(
            fill =
              case_when(
                name %in% phen ~ "grey",
                TRUE ~ "white"
              )
          )

        # Set seed so that same graph drawn everytime
        set.seed(seed)

        ggraph::ggraph(
          edge_tbl_graph,
          layout = "igraph",
          algorithm = "fr"
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
            fill = fill$fill,
            size = 3
          ) +
          ggplot2::scale_y_reverse() +
          ggraph::scale_edge_colour_distiller(
            palette = "RdYlBu",
            limits = c(-1,1),
            direction = -1
          ) +
          ggraph::facet_edges(
            vars(eqtl_dataset, gene_name)
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            legend.position = "top"
          )

      })

    return(plots)

  }
