#' Correlate taxa abundances with lipid abundances
#'
#' This function allows you to extract the phylogeneic affiliations of OTUs, using a theshold correlation coefficient
#' Returns a table with each row representing one OTU, and each column one lipid compound. Values are affiliations of above-threshold OTUs.
#' 'lower' = below threshold, 'other' = sumsumed according to cutoff fraction
#'
#' @param phy Phyloseq object containing an OTU and taxonomy table.
#' @param lipids Data frame with concentration of lipids (or any other variables/parameters) to correlate with phylogenetic data. Rows are observations/samples and must match the samples in 'phy'. Used by 'phlo_lipids'
#' @param method  correaltion method for 'cor' {stats}, i.e., "pearson", "kendall", or "spearman"
#' @param thresh threshold correlation coefficient, above which OTUs will be considered (default = 0.75)
#' @param level character string with taxonomic level of analysis, e.g., "Phylum", "Class", etc
#' @param cutoff cutoff below which OTUs will be subsumed as "other" (applies for each compound/lipid separately)
#' @param exempt taxonomic unit that is expemt from the cutoff, i.e., is always displayed regardless of count numbers
#' @examples
#' @return A data frame with each row representing one OTU, and each column one lipid compound. Values are character strings (phylogenetic affiliations) of above-threshold OTUs. 'lower' = below threshold, 'other' = sumsumed according to cutoff fraction
#' @export
#' @family Main functions


cor_tax = function (phy,
                    lipids,
                    thresh = 0.75,
                    cutoff = 0.03,
                    method= "pearson",
                    level= "Phylum",
                    exempt = FALSE
                    ) {

  ## extract data from phyloseq object ----
  x = as.data.frame (t (phy@otu_table), stringsAsFactors = F) # get otu table
  y = as.data.frame (phy@tax_table[,level], stringsAsFactors = F) # get tax table, no factors!
  y1 = tibble::as.tibble (y) %>%
    dplyr::select (level) %>% # get taxonomy assignments at the specified level
    dplyr::mutate (otu= rownames (.)) %>% # add otu column
    dplyr::mutate (!!level := as.character (.[[level]]))


  ## check if row numbers match ----
  if ( nrow(lipids) != nrow(x) ) {
    stop("Number of rows do not match between DNA and lipid data")}

  ## check if 'exempt' taxon is in the tax table ----
   if (is.character (exempt)) {
     if (try (exempt %in% unique (y [[level]])) == FALSE) {stop ("You have provided a wrong 'exempt' taxon! Check spelling!")
     }
   }

  #### correlate each lipid/compound in 'lipids' with all OTUs in 'phy' ----
  stats::cor(x, lipids, method = method[1]) %>%
    as.data.frame() -> cormat
  cormat %>%
     tibble::as.tibble() %>%
     # add otu column
     dplyr::mutate(otu = rownames(cormat),
     # add taxonomy of sepcified level
     !!level := y[[level]]) -> a

  b = a ## working copy within the the loop
  b1 = a ## for substitution

  #### FOR EACH LIPID ----
  ## loop through all lipids/compounds and gather taxonomy
  for (j in colnames(lipids)) { # for each lipid...

    b %>%
      dplyr::select(j, level) %>%
      # replace below threshold taxa with "low.cor"
      dplyr::mutate_at(dplyr::vars(level),
                       dplyr::funs(ifelse(b[[j]] < thresh, "low.cor", .)
                         )) %>%
      # replace unassigned taxonomy 'NA' with 'ukn.'
      dplyr::mutate_at(dplyr::vars(level),
                       dplyr::funs(ifelse( is.na(.),
                                         "ukn."
                                         , .)
                      )) -> c

     # determine the cutoff taxa for this lipid ----
     c %>%
        dplyr::select(j, level) %>%
        dplyr::filter( .[[j]] > thresh) %>%
        dplyr::group_by_(paste(level)) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::mutate(rel.cnt = n/ sum(n)) -> d
unique(d$Phylum)
     if (nrow(d) > 0) {
       if (sum(d$rel.cnt) != 1) {
         stop("false rel.cnt")
       }
     }

     # list of low-abundance taxa
     d[[level]][which(d$rel.cnt < cutoff)] -> cu
     # remove the 'ukn' from low abundance list, because we always want to see those
     if ("ukn." %in% cu ) { # if ukn.s are present, remove them (without 'if', returns c)
       cu = cu[-which (cu == "ukn.")] }
     cu = cu[order (cu)] # order that list alphabetically

     ## apply exemption to cutoff ----
     # check if exempt taxa/taxon are/is among the low abundance taxa
     # if so, remove this from the list of low abundance taxa
     if (try (exempt %in% cu) == TRUE ) { cu = cu[-which(cu == exempt)]}

     # replace below cutoff taxa values in 'c' with 'other'----
     c %>%
        dplyr::mutate_at(dplyr::vars(level),
                         dplyr::funs(ifelse( . %in% cu, "other", .)
                         )
        ) -> c1

     # finally replace data in cor.matrix with filtered taxonomy
     b1[[j]] = c1[[level]]
  }
## END OF LOOP

  # remove the rank (level) column
  b1 %>% dplyr::select( colnames(lipids), otu) -> b2

  return(b2)
}

