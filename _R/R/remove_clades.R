#' Removing taxonomic groups from Phyloseq objects
#'
#' This function removes taxonomic groups from Phyloseq objects. Clades can be of different ranks.
#' @param phy Phyloseq object
#' @param clades Character vector with taxonomic groups to be removed
#' @export
#' @examples
#' @keywords remove groups taxonomic


remove_clades = function(phy, clades) {
 all.taxa = phyloseq::taxa_names (phy)
  x = as.data.frame (phy@tax_table, stringsAsFactors = F)
  ranks = colnames(x)
  x$otu = all.taxa
  a=character()
  for (i in ranks) {
    x %>%
    dplyr::filter(get(i) %in% clades) -> y
   head(y)
   a = c (a, y$otu)
  }
  b = all.taxa[which(all.taxa %!in% a)]
  return (phyloseq::prune_taxa (b, phy))
}

