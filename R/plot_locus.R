#' Plot genes within a LAVA LD block.
#'
#' This function will plot all genes of the gene biotypes "protein_coding",
#' "antisense", and "lincRNA" within the start and stop of an LD block.
#'
#' @param locus_gr gr. Granges object with coordinates of each locus, together
#'   with at least 1 additional metadata column: \itemize{ \item locus - locus
#'   identifier. }
#' @param ref gr. Granges object containing reference ensembl gtf.
#'
#' @return ggplot2 plot of genes (of biotype "protein_coding", "antisense", and
#'   "lincRNA") within locus.
#' @export
#' 


plot_locus <- function(locus_gr, ref){
  
  overlap <- 
    findOverlaps(locus_gr, ref) %>% 
    as_tibble()
  
  coords_to_plot <-
    tibble(
      group_type = as.factor("genes"),
      gene_id = ref[overlap$subjectHits]$gene_id,
      gene_name = 
        ref[overlap$subjectHits]$gene_name %>% 
        str_replace_all("-", "_") %>% 
        str_replace_all("\\.", "_"),
      gene_biotype = ref[overlap$subjectHits]$gene_biotype %>% as.factor(),
      chr = ref[overlap$subjectHits] %>% seqnames() %>% as.character(),
      strand = ref[overlap$subjectHits] %>% strand() %>% as.character(),
      start = ref[overlap$subjectHits] %>% start(),
      end = ref[overlap$subjectHits] %>% end(),
      width = ref[overlap$subjectHits] %>% width()
    )  %>%
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
        gene_name = str_c("Locus_", locus_gr$locus),
        chr = locus_gr %>% seqnames() %>% as.character(),
        start = locus_gr %>% start(),
        end = locus_gr %>% end(),
        gene_biotype = as.factor("locus")
      )
    )
  
  # Re-factor groups
  coords_to_plot <- 
    coords_to_plot %>% 
    dplyr::mutate(
      gene_biotype = 
        gene_biotype %>%
        forcats::fct_relevel(.,
                             c("locus", "protein_coding", "antisense", "lincRNA")
        )
    ) 
  
  ggplot2::ggplot(
    data = coords_to_plot
  ) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        x = start,
        y = fct_reorder(
          gene_name,
          -start
        ),
        xmin = start,
        xmax = end
        # colour = gene_biotype
      ),
      size = 1
    )+
    ggplot2::labs(title = str_c("Locus: ", locus_gr$locus)) +
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
      x = str_c("Chromosome ", c(coords_to_plot$chr %>% unique()))
    ) +
    ggplot2::facet_grid(
      rows = vars(gene_biotype), 
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::theme_bw(
      base_size = 8
    ) +
    theme(
      strip.text.y = element_text(angle = 0)
    )
  
}