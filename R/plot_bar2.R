#' helper function: Plots taxonomy bar chart
#'
#'
#' @param x data frame or tibble containing the columns:
#' 1) Taxonomy (character)
#' 2) number of taxa (integers)
#' 3) color values ('col')
#' @export
#' @examples
#' @return ggplot bar chart
#' @keywords internal


plot_bar2 = function(x, count = c("n","avg")) {
   x %>%
      ggplot2::ggplot(ggplot2::aes_string( y = count[1] , x = "0", fill=x[[1]]) ) +
      ggplot2::geom_bar(stat="identity") +
      #ggplot2::coord_polar("y", start = 0) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = colnames(x[1]) ) ) +
      ggplot2::scale_fill_manual(values = x$col) +
      ggplot2::scale_y_continuous(expand = c(0,0)) +
      ggplot2::theme(aspect.ratio = 4) -> p
   return(p)
}



