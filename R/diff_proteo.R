#' Replaces the phylum Proteobacteria with the classes (i.e., Aplpha-, Beta-, Deltaproteobacteria, etc.)
#'
#' Aussumes that second taxonomic rank is 'Phylum'
#' @param phy Phyloseq object containing an OTU and taxonomy table.
#' @export
#' @examples
#'

diff_proteo = function(phy) {
  phy@tax_table %>%
    as.data.frame (stringsAsFactors=FALSE) -> x
  x %>%
    dplyr::mutate_at (dplyr::vars (2), dplyr::funs (ifelse(. == "Proteobacteria", Class, .))) -> y
  rownames(y) = rownames(x)
  phyloseq::tax_table (as.matrix (y)) -> tax
  phy2 = phyloseq::phyloseq( try(phy@refseq) ,   try(phy@sam_data) , try(phy@phy_tree) , phy@otu_table, tax )
  return(phy2)
}
