#' trim by occurrence in min number/fraction of samples
#'
#'
#' @param phy Phyloseq object containing an OTU and taxonomy table.
#' @param min threshold for taxa to return
#' @param fraction if TRUE, a relative threshold is used
#' i.e., 'min' das to be a frction between 0 and 1
#' @export
#' @examples
#'

trim_taxa_samp = function(phy, min, fraction = F) {
  x = as.data.frame(phy@otu_table)
  nsamples = ncol(x)
  rownames(x) -> otu.names
  # make an 'occurrence matrix'
  x %>%
    dplyr::mutate_all (dplyr::funs (ifelse (. > 0, 1, .))) -> y
  rownames(y) = otu.names

  if (fraction == T) {
   if (min <= 0 | min >= 1) {
     stop("'min' is not a fraction!")
   }
   taxa.to.keep = names(
     rowSums (y))[which(
       rowSums(y)/nsamples >= min)]
  }
  if (fraction == F) {
    taxa.to.keep = names ( rowSums (y))[which (rowSums (y) >= min)]
  }
  return( phyloseq::prune_taxa( taxa.to.keep, phy))
}
