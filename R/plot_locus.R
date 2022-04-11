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
#' @param window `integer` vector indicating base pairs to be added around gene
#'   bodies.
#' @param highlight_gene `character` vector with genes (as ensembl IDs) to
#'   highlight in locus. Highlighted genes will be coloured blue and labelled
#'   "TRUE".
#' @param highlight_gene_label `character` vector with label to use for
#'   `highlight_gene` legend.
#'
#' @return ggplot2 plot of genes (of biotype "protein_coding", "antisense", and
#'   "lincRNA") within locus.
#' @export
#'
#'

plot_locus <-
  function(
    locus_gr,
    ref,
    window = NULL,
    highlight_gene = NULL,
    highlight_gene_label = NULL
  ) {

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

    if(!is.null(window)){

      coords_to_plot <-
        coords_to_plot %>%
        dplyr::mutate(
          start =
            case_when(
              start - window < 0 ~ 0,
              TRUE ~ start - window
            ),
          end = end + window
        )

    }

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

    if(!is.null(highlight_gene)){

      coords_to_plot <-
        coords_to_plot %>%
        dplyr::mutate(
          gene_highlight =
            case_when(
              gene_id %in% highlight_gene ~ TRUE,
              TRUE ~ FALSE
            )
        )

      plot <-
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
            xmax = end,
            colour = gene_highlight
          ),
          size = 1
        ) +
        ggplot2::scale_colour_manual(
          values =
            c(
              "TRUE" = "#00BFC4",
              "FALSE" = "black"
            )
        ) +
        ggplot2::labs(
          colour = highlight_gene_label
        ) +
        theme(
          legend.position = "top"
        )


    } else{

      plot <-
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
          ),
          size = 1
        )

    }

    plot +
      ggplot2::labs(title = stringr::str_c("Locus: ", locus_gr$locus)) +
      ggplot2::scale_x_continuous(labels = scales::comma) +
      ggplot2::labs(
        x = stringr::str_c("Chromosome ", c(coords_to_plot$chr %>% unique())),
        y = ""
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
