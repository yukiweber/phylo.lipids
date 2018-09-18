#' Standardized otu abundances
#'
#' @return Returns a data standardized otu table with a 'sum' column containg the sums of the original data
#' @param phy Phyloseq object containing at least an OTU and a taxonomy table.
#' @param method Method of standardization of the otu abundances,
#' either to the sum or the maximum across all samples
#' @export


standardize_tax = function (phy, method) {

   otu = as.data.frame(unclass(phy@otu_table))

   # compute abundance SUMs across all samples ----
   sum = Reduce('+', otu) # / ncol(x) # prev: averages


   if (method[1] == "max") {
            otu =
               as.data.frame(
                  lapply ( otu , `/` , # apply to each row
                          apply(otu, 1, max) # max values in reach row (otu)
                        )
               )
         }
         # OR the abundance sum across all samples
         if (method[1] == "sum") {
            otu = as.data.frame(sapply(otu,`/`,rowSums(otu)))
         }


   # add the sums to the table
   otu$sum = sum
   return(otu)
}
