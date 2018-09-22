#' helper function: Plots taxonomy pie chart
#'
#' @param x data frame or tibble containing at least the columns:
#' + Taxonomy as character vector (1st cloumn)
#' + number of taxa ('n')
#' + summed abundance acoss data set ('sum')
#' + color values ('col')
#' @param count select the summary method: otu counts (n) or otu summed abundances (sum)
#' @param ncol Number of columns in facets
#' @param taxa_groups A optional list of character vector containing grou taxa names for facetting
#' @param level Taxonomic level
#' @param key_col Number of columns of the color key
#' @export
#' @examples
#' @return ggplot pie chart, and a summary table
#' @keywords visualization


plot_pie = function(x, # long format otu table with col, tax ('level')
                    count = c("n","sum"),
                    level = "Phylum",
                    taxa_groups = NULL, # a list of char vectors containing grou taxa names
                    ncol = 5,
                    key_col = 2
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


   # summarize
   x %>%
      dplyr::arrange_(paste(level)) %>%
         dplyr::group_by_(level, "group") %>%
         dplyr::summarise(n=n(),
                          sum = sum(sum),
                          col= unique(col)) -> x1

   # get the color vector for the plot
   tax_col =
      x1 %>%
      dplyr::group_by_(level) %>%
      dplyr::summarise(col= unique(col))

   # plot according to selection in 'count' (n or sum)
   x1 %>%
      ggplot2::ggplot(ggplot2::aes_string( y = count[1] , x = "0", fill= paste(level)) ) +
      ggplot2::geom_bar(stat="sum",
                        color="black",
                        size=0.5,
                        position = "fill") +
      ggplot2::coord_polar("y", start = 0) +
      # use the color vector previously generated for color mapping
      ggplot2::scale_fill_manual(values = tax_col$col) +
      ggplot2::facet_wrap("group",
                          ncol = ncol) +
      ggplot2::guides(fill = ggplot2::guide_legend(
         title = paste0(level," (", count[1], ")"),
         ncol = key_col)
         ) -> p
   return(
      list(
         "pie" = p,
         "tab" = x1
         )
   )

}


