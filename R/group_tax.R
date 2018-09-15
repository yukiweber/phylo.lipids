#' Finds groups of OTUs/taxa with similar spatial distribution
#'
#' Aussumes that second taxonomic rank is 'Phylum'
#' @param phy Phyloseq object containing at least an OTU and a taxonomy table.
#' @param h Height passed to 'cuttree'. Will be overridden by 'k' if k is numeric.
#' @param k Number of branched for 'cuttree'.
#' Overrides 'h' setting if provided (i.e., numeric).
#' @param dist Distance matrix type passed to 'dist()' to be used in clustering
#' @param clust Hclust method passed to 'hclust()'
#' @export
#' @examples
#' @return A list with i elements, containing the OTU/taxa names in the i-th group

group_tax = function(phy,
                     h = FALSE,
                     k = 20,
                     dist = c( "euclidean","bray","..."),
                     clust = c("ward.D2", "ward.D", "single", "complete", "average", "mcquitty", "median", "centroid")
                     ) {

   if (is.numeric(h)==T & is.numeric(k)==T){stop("Provide either h or k to cut the cluster tree!")}

   # CLUSTER helper function----
   clust_tax = function (phy, h = F, k = F) {

      otu = tibble::as.tibble(as.data.frame(phy@otu_table, strngsAsFactors = F))

      # make PROFILES
      otu = sapply(otu,`/`,rowSums(otu)) ### profiled PTU abundances

      # do the clustering
      hc = hclust( dist(otu, method = dist[1] ), method = clust[1])


      # cut tree gives a table with otu name vs group affiliation
      if (is.numeric(h) == T) {hc1 = cutree(hc, h = h)}
      if (is.numeric(k) == T) {hc1 = cutree(hc, k = k)}
      if ( sum( is.numeric(h) , is.numeric(k)) == 0)  {
         stop("Provide either h or k to cut the cluster tree!")}

      # plot the tree
      # as to be written... or not

      # format that table
      hc2 =
         hc1 %>%
      as.data.frame() %>%
      dplyr::mutate(otu = names(hc1)) %>%
      dplyr::rename(group = ".") %>%
      tibble::as.tibble()

      return(hc2)
   }

   clust_tax(phy, h=h, k=k) -> C

   # number of taxa in each group
   C %>%
      dplyr::group_by(group) %>%
      dplyr::tally() -> Cn

   # prepare a lis
   L = list()
   # for each group
   for (i in 1: max(C$group)){
      L[paste("G",i,sep = "")] = C["otu"][which(C["group"] == i), ]
   }
return(L)
}
