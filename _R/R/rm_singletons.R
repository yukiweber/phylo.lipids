#' Removes singletons from phyloseq object
#'
#' i.e., OTUs with only one read across all samples
#' @param phy Phyloseq object containing an OTU and taxonomy table.
#' @export
#' @examples
#'

### remove singletons (== row sums are exactly 1 == occurr)
rm_singletons = function(phy) {
  all_taxa = phyloseq::taxa_names(phy)
  test = as.data.frame (phy@otu_table)
  test_sums = as.data.frame (rowSums(test))
  single_tab = subset (test_sums, test_sums[,1] == 1)
  singles = rownames(single_tab)
  not_singles = all_taxa[which(all_taxa %!in% singles)]
  return (phyloseq::prune_taxa( not_singles, phy))
}
