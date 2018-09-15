#' helper function: Plots taxonomy pie chart
#'
#'
#' @param x data frame or tibble containing at least the columns:
#' + Taxonomy as character vector (1st clomn)
#' + number of taxa ('n')
#' + average abundance acoss data set ('avg')
#' + color values ('col')
#' @param count select the reference: otu counts (n) or otu average abundances (avg)
#' @export
#' @examples
#' @return ggplot pie chart
#' @keywords internal


plot_pie = function(x, # long format otu table with col, tax ('level')
                    count = c("n","avg"),
                    level = "Phylum",
                    taxa_groups = NULL # a list of char vectors containing grou taxa names
                    ) {

   # if no groups are provided/wanted
   # set group colum to ""
   if (is.null(taxa_groups) == F || is.null(x$group) == TRUE) {
      x$group = ""
      # message("proceeding without groups")
   }

   # prune the taxa provided
   if(is.null(taxa_groups) == F) {
      x =
         x %>%
         dplyr::filter(otu %in% taxa_groups)
   }

   # factorize groups to preserve order
   x$group = factor (x$group, levels = unique(x$group))

   # summarize
   x %>%
      dplyr::arrange_(paste(level)) %>%
         dplyr::group_by_(level, "group") %>%
         dplyr::summarise(n=n(),
                          avg=sum(avg),
                          col= unique(col)) -> x1

   # get the color vector for the plot
   tax_col =
      x1 %>%
      dplyr::group_by_(level) %>%
      dplyr::summarise(col= unique(col))

   # plot according to selection in 'count' (n or avg)
   x1 %>%
      ggplot2::ggplot(ggplot2::aes_string( y = count[1] , x = "0", fill= paste(level)) ) +
      ggplot2::geom_bar(stat="sum",
                        color="black",
                        size=0.5,
                        position = "fill") +
      ggplot2::coord_polar("y", start = 0) +
      # use the color vector previously generated for color mapping
      ggplot2::scale_fill_manual(values = tax_col$col) +
      ggplot2::facet_wrap("group") +
      ggplot2::guides(fill = ggplot2::guide_legend(title = paste0(level," (", count[1], ")"))) -> p
   return(p)
}

