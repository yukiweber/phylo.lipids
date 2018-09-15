#' Summarizes taxonomy at a given taxonomic level
#'
#' a cutoff can be fifined for subsummation as 'other, as well as expemtions to this rule
#' @param phy Phyloseq object containing at least an OTU and a taxonomy table.
#' @export
#' @examples
#' @return A list with i elements, containing the OTU/taxa names in the i-th group


sum_tax = function (phy,
                    cutoff = 0.03,
                    level= "Phylum",
                    exempt = FALSE
                    ) {

     ## extract data from phyloseq object ----
      # OTU table
      x = as.data.frame (phy@otu_table, stringsAsFactors = F) %>% # get otu table
        dplyr::mutate(avg = Reduce('+', .)/ ncol(.), # compute average abundance across samples
                      otu = rownames(.)) # add otu column
      # TAX table
      y = as.data.frame (phy@tax_table[,level], stringsAsFactors = F) # get tax table, no factors!
      y1 = tibble::as.tibble (y) %>%
       dplyr::select (level) %>% # get taxonomy assignments at the specified level
       dplyr::mutate (otu = rownames (.)) %>% # add otu column
       dplyr::mutate (!!level := as.character (.[[level]])) %>% # make sure taxonomy column it's not a factor
       dplyr::left_join (x, by ="otu") # add abundance data


  ## check if 'exempt' taxon is in the tax table ----
   if (is.character (exempt)) {
     if (sum(exempt %in% unique(y[[level]]) ) < length(exempt)) {stop ("You have provided a wrong 'exempt' taxon! Check spelling!")
     }
   }

  #### summarize

   # helper function to summarize
   # summaraizes abundace and counts,
   # replaces NA with 'ukn'
   sum_tax_help = function (y) { # y: taxonomy table with abundance data
   z =
      y %>%
      dplyr::group_by_(paste(level)) %>%
      dplyr::summarise(n = n(),
                       avg = sum(avg)) %>%
      # relative counts
      dplyr::mutate(rel.cnt = n/ sum(n)) %>%
      # replace unassigned taxonomy 'NA' with 'ukn.'
      dplyr::mutate_at(dplyr::vars(level),
                       dplyr::funs(ifelse( is.na(.) | .=="NA",
                                         "ukn."
                                         , .)
                      ))
   return (z) }

     # summarize first
     y2 = sum_tax_help(y1)

     ## determine the cutoff taxa ----
     # list of low-abundance taxa
     y2[[level]][which(y2$rel.cnt < cutoff)] -> cu
     # remove the 'ukn' from low abundance list, because we always want to see those
     if ("ukn." %in% cu ) { # if ukn.s are present, remove them (without 'if', returns c)
       cu = cu[-which (cu == "ukn.")] }
     cu = cu[order (cu)] # order that list alphabetically

     ## apply exemption to cutoff ----
     # only keep taxa in the cutoff list that are not in exempt
     cu = cu[which(cu %!in% exempt)]

     # replace below cutoff taxa values in with 'other'----
     y1 %>%
        dplyr::mutate_at(dplyr::vars(level),
                         dplyr::funs(ifelse( . %in% cu, "other", .)
                         )
        )-> y3
     # summarize again
     y4 = sum_tax_help(y3)

   # make taxonomy a factor and re-level
   tax = as.factor(y4[[level]])
   if ("other" %in% levels(tax) ) {
      tax = relevel(tax ,"other") }
   if ("ukn." %in% levels(tax) ) {
      tax = relevel(tax ,"ukn.") }
   y4[level] = tax

   return(y4)
}

