#' Finds groups of OTUs/taxa with similar spatial distribution
#'
#' Aussumes that second taxonomic rank is 'Phylum'
#' @param phy Phyloseq object containing at least an OTU and a taxonomy table.
#' @param h Height passed to 'cuttree'. Will be overridden by 'k' if k is numeric.
#' @param k Number of branched for 'cuttree'.
#' Overrides 'h' setting if provided (i.e., numeric).
#' @param dist Distance matrix type passed to 'dist()' to be used in clustering
#' @param clust Hclust method passed to 'hclust()'
#' @param standardize Method of standardization of the otu abundances,
#' either to the sum or the maximum across all samples
#' @param count Determines whether the statistics are based on the NUMBER of otus or their ABUNDANCE SUMS
#' @export
#' @examples
#' @return A list with i elements, containing the OTU/taxa names in the i-th group

group_tax = function(phy,
                     h = FALSE,
                     k = FALSE,
                     dist = c( "euclidean","bray","..."),
                     clust = c("ward.D2", "ward.D", "single", "complete", "average", "mcquitty", "median", "centroid"),
                     standardize = c("sum","max"),
                     count = c("n","sum")
                     ) {

   if (is.numeric(h)==T & is.numeric(k)==T){stop("Provide either h or k to cut the cluster tree!")}

   # CLUSTER helper function----
   clust_tax = function (phy, h = F, k = F) {

      otu = tibble::as.tibble(as.data.frame(phy@otu_table, strngsAsFactors = F))

      # compute abundance SUMs across all samples ----
      sum = Reduce('+', otu) # / ncol(x) # prev: averages

      # make PROFILES
      # otu = sapply(otu,`/`,rowSums(otu)) ### deprecated as it is done in the step before

      # standardize abundances of each otu relative to  ----
      # the max value across all samples
      otu0 = standardize_tax(phy, method = standardize[1])
      otu = otu0 %>%
         dplyr::select(-sum)

      # do the clustering
      hc = hclust( dist(otu, method = dist[1] ), method = clust[1])

      # cut tree gives a table with otu name vs group affiliation
      if (is.numeric(h) == T) {hc1 = cutree(hc, h = h)}
      if (is.numeric(k) == T) {hc1 = cutree(hc, k = k)}
      if ( sum( is.numeric(h) , is.numeric(k)) == 0)  {
         stop("Provide either h or k to cut the cluster tree!")}

      # plot the tree
      # as to be written... or not

      # format that group table
      hc2 =
         hc1 %>%
      as.data.frame() %>%
      dplyr::mutate(otu = names(hc1),
                    sum = otu0$sum) %>%
      dplyr::rename(group = ".") %>%
      tibble::as.tibble()

      return(hc2)
   }

   clust_tax(phy, h=h, k=k) -> C

   # compute group count/abundance stats of otus
   C %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(n = n(),
                       sum=sum(sum)) %>%
      # depending on selectioj in 'count'
      dplyr::mutate(rel = get(count)/sum( get(count) )
                    ) -> gStat

   # prepare a lis
   L = list()
   # for each group
   for (i in 1: max(C$group)){
      L[paste("G",i,sep = "")] = C["otu"][which(C["group"] == i), ]
   }



      # add statistics to group names
   names(L) =
      paste0(
         names(L),
         " (",
         sapply(gStat$rel, function(x) signif(x,2)) *100,
         "%)")

   # order groups by abundance
   ord =  order(gStat[[count[1]]])
   L = L[rev(ord)]

   return(L)
}


