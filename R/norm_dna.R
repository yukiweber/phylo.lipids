#' normalizes read numbers
#'
#' Returns either relative abundances (relative == TRUE), or OTU-specific DNA abundances in [ng] (!)
#' @param phy Phyloseq object containing an OTU and taxonomy table.
#' @param dna numeric vector containing absolute DNA abundance in µg,
#' which will be multipled with the relative OTU abudances
#' @param relative If TRUE, the relative abundances will be returned, without multiplication with 'dna.
#' Set to FASLE for absolute read counts. (default = FASLE)
#' @keywords normalize DNA concentrations OTU-specific
#' @export
#' @examples
#' @family Main functions

norm_dna = function(phy, dna = 1, relative = T) {
  # make relative abundances
  x = phyloseq::transform_sample_counts(phy, function (x) x / sum (x))
  # get otu table and modify
  otu = as.data.frame (x@otu_table, stringsAsFactors = F)
  otu.names = rownames (otu)
  # multiply with total dna abundance ()
  otu.norm = data.frame (mapply ('*', otu, dna * 1000)) ### from micro gram to ng /L !!!
  rownames (otu.norm) = otu.names
  phy1 =
   phyloseq::phyloseq (
     phyloseq::otu_table (as.matrix(otu.norm), taxa_are_rows = T),
     phy@tax_table,
     phy@sam_data,
     phy@refseq )
  if ( relative == F) {
    return (phy1)
    warning("output: abundances in ng (if DNA input in micro gram) !!!")
  }
  if (relative == T) {
    return (x)
  }
}
