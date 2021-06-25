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
#'   \item x and y-axis display the phenotypes. 
#' \item Significant negative and positive correlations are indicated by blue
#' and red fill, respectively. Non-significant correlations (p >= 0.05/n_tests)
#' have a grey fill.
#'   }
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile geom_text coord_fixed labs
#'   scale_fill_distiller scale_colour_manual guides theme_bw theme unit
#' @importFrom forcats fct_relevel
#' @importFrom stats setNames
#'
#' @references \itemize{ \item Bulik-Sullivan et al. (2015) An atlas of genetic
#'   correlations across human diseases and traits \emph{Nature Genetics}, 2015
#'   Nov;47(11):1236-41. \url{https://www.nature.com/articles/ng.3406} PMID:
#'   26414676 }

plot_global_corr <- 
  function(
    global_corr,
    n_phenotypes
  ){
    
    # Only need first half of matrix, thus must extract appropriate rows from dataframe
    n <- n_phenotypes
    
    for(i in 1:n){
      
      # General formula for extracting appropriate indices
      index <- (i*n- (n-i)):(i*n)
      
      if(i == 1){
        
        indices <- index
        
      } else{
        
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
          case_when(
            p < 0.05/n_combn ~ round(rg, 2)
          ),
        p1 = p1 %>% 
          str_replace_all("[:digit:]", "") %>% 
          str_remove("\\..*"),
        p2 = p2 %>% 
          str_replace_all("[:digit:]", "") %>% 
          str_remove("\\..*")
      ) %>% 
      dplyr::slice(indices) %>% 
      ggplot(
        aes(
          x = p1,
          y = p2 %>% fct_rev(),
          fill = rg_fill,
          label = round(rg, 2)
        )
      ) +
      geom_tile(colour = "black") + 
      geom_text(
        size = 3
      ) +
      coord_fixed() +
      labs(x = "", y = "", fill = "Global correlation (rg)") +
      scale_fill_distiller(palette = "RdYlBu", limits = c(-1,1)) +
      scale_colour_manual(values = c("black", "white")) +
      guides(colour = F) +
      theme_rhr +
      theme( 
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

plot_locus <- 
  function(
    locus_gr, 
    ref
    ){
  
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