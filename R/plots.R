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
      ggplot2::scale_fill_distiller(palette = "RdYlBu", limits = c(-1, 1)) +
      ggplot2::scale_colour_manual(values = c("black", "white")) +
      ggplot2::guides(colour = F) +
      theme_rhr +
      ggplot2::theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank()
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
