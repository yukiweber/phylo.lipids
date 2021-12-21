#' trim Phyloseq onject by threshold sum across all samples
#'
#'
#' @param phy Phyloseq object containing an OTU and taxonomy table.
#' @param min Threshold value (minimum) for taxa to be returned.
#' Context dependent: reads, concentrations, relative abundanece, etc
#' @export
#' @examples
#'

##
trim_taxa_sum = function( phy, min) {
  phy1 = phyloseq::filter_taxa (phy, function(x) sum(x) >= min, prune = TRUE)
  return (phy1)
}
