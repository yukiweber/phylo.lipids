#' Phylogenetic diversity of correlated taxa
#'
#' This function correlates absolute abundances of OTUs vs lipids, and summarizes the phylogenetic diversity of correlated taxa for each lipid.
#' @param phy Phyloseq object containing an OTU and taxonomy table.
#' @param lipids Data frame with concentration of lipids (or any other variables/parameters) to correlate with phylogenetic data. Rows are observations/samples and must match the samples in 'phy'.
#' @param method correaltion coefficient used in 'cor' {stats}, i.e., "pearson", "kendall", or "spearman"
#' @param thresh threshold correlation coefficient, above which OTUs will be considered (default = 0.75)
#' @param level character string with taxonomic level of analysis, e.g., "Phylum", "Class", etc
#' @param cutoff cutoff below which OTUs will be subsumed as "other" (applies for each compound/lipid separately)
#' @param exempt taxonomic unit(s) that is/are expemt from the cutoff, i.e., is always displayed on top regardless of cutoff value and sorting
#' @param n.col number of columns for pie charts (default = 3)
#' @param sort.tax sort taxa aphabetically or by counts. Either "alpha" or "counts" (default = "alpha")
#' @param plot.theme a ggplot theme object which overrides defualt plot theme
#' @param text.size global text size, titles and labels are relative (default = 12)
#' @param test randomply subsample otus for testing large data sets (default = FALSE)
#' @param test.size number of random otus to draw (default = 100)
#' @param guide.col number of legend columns
#' @keywords correlation lipids DNA pie charts
#' @examples
#' @importFrom magrittr %>%
#' @return A list, containing a pie chart (.$pies), a table in long format (.$table), and a summary text (.$summary).
#' @export
#' @family Main functions
#' @seealso vignette("phylo.lipids")
#' @seealso ?phylo.lipids

phylo_lipids = function (
  phy,
  lipids,
  method = c("pearson", "spearman", "kendall"),
  thresh = 0.75,
  level = "Phylum",
  cutoff = 0.0 ,
  exempt = FALSE,
  n.col = 3,
  sort.tax = c("alpha","counts"),
  plot.theme = FALSE,
  text.size = 12 ,
  test = FALSE ,
  test.size = 102,
  guide.col = 1
  ) {

  ## disable warnings temporarily in this function (turn back on after function execution)
  oldw <- getOption("warn") ## safe old warning setting
  options(warn = -1)

  #### check inputs ----
  if (phyloseq::ntaxa (phy) < 3) {
    stop ("Less than 3 taxa provided !!!")
  }

  #### check tax level ----
  colnames(phy@tax_table) -> ranks
  if ( level %in% ranks == F) {
    stop ("You provided the wrong taxonomic level! ",
       paste ("Rank '", level, "' not found in taxonomy table.", sep=""))
  }

  ## subsample for testing? ----
  if (test == TRUE) {
    if (try (is.numeric(test.size)) == FALSE ) {
      stop("Provide sample size for subsampling ('test.size')")
    }
    sample (split (phyloseq::taxa_names (phy), phyloseq::taxa_names (phy)) , size = test.size) -> p
    phy = phyloseq::prune_taxa (as.character(p) , phy)
  }


  ## Correlate all lipids with all OTUs and summarize the OTUs' taxonomic  affiliations at given taxonomic level ----
  COR = # matrix with cor coeffs for each lipid x each otu; + taxonomic afilliation
    cor_tax (phy = phy, lipids = L, method = method[1], thresh = thresh, level = level, cutoff = cutoff, exempt = exempt) %>%
    tibble::as.tibble()

    ## get total number of OTUs that correlated above threshold with at least one of the lipids ---
  cor = COR %>%
    tidyr::gather (key = variable, value = value, -otu) %>%
    tibble::as.tibble() %>%
    dplyr::filter (value != "low.cor") # remove all 'low' correlations with any of the lipids
  ## if correlation was avove thresh for any of the lipids, the otu name will be retained
  n.cor =  cor$otu %>% unique() %>% length()

  ## summarize number of OTUs for each lipid by taxonomy ---
  tax =
    COR %>%
    dplyr::mutate(otu = NULL) %>% # delete otu name column before long format
    tidyr::gather (key= variable , value= !!level, colnames(.)[colnames(.)%in%colnames(L)] ) %>% # same as melt: make long format
    dplyr::select (!!level , variable) %>%
    dplyr::filter (paste(get(level)) != "low.cor") %>% # omit 'low.cor' below threshold correlation
    dplyr::group_by_ (paste(level), "variable") %>% # standard evaluation
    dplyr::summarise (n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange (variable) %>% # sort by compound names
    dplyr::mutate (variable = as.factor (variable)) # factorize lipid names
     # dplyr::mutate(!!level := as.factor(.[[!!level]])) # but not taxa names yet

  ## list of well-correlating taxonomic groups across all lipids ----
  unique (tax[[level]]) %>%
     tibble::tibble() -> tax.levels # character vector
  t = list() # put it in a list
  # fill the list
  for (i in tax.levels) {
   t[i]=i
   }

  ## get number of highly correalted OTUs for each compound
  ntax=
    tax %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(n = sum(n, na.rm = TRUE)) %>%
    ## fill randomply with taxa names for ggplot aesthetics
    dplyr::mutate(value=1, # for ggplot aesthetics
           !!level := as.character ( sample(t, size = nrow(.) , replace = TRUE)
                                     )
           ) %>%
    dplyr::mutate(variable = as.factor(variable)) %>% # make compound column a factor
    dplyr::arrange(variable) %>% # sort by compound names
    dplyr::mutate(rel.tot = .$n/phyloseq::ntaxa(phy)) %>% ## compute fractions relative to total OTUS
    dplyr::arrange(variable)

  ## get order of taxa count numbers across all lipids
  tax.ord =
    tax %>%
      dplyr::group_by_(paste(level) ) %>%
      dplyr::summarise(n = sum(n) )  %>% # sum of counts for each taxon
      dplyr::arrange(-n) %>% # most common taxa at the top
      dplyr::ungroup()

  ## factorize taxa in 'tax'
  tax[[level]] = as.factor(tax[[level]])

  ## sort taxa levels: alphabetically or by counts ----
  if (sort.tax[1] == "alpha") {
    tax=
      tax %>%
        dplyr::mutate(!!level := factor (as.character (.[[!!level]]),
                                         levels = as.character (levels(.[[!!level]])[order(as.character(levels(.[[!!level]])))])
                                         )
                      )
  }

  if( sort.tax[1] =="count") {
    tax=
      tax %>%
      dplyr::mutate (!!level := factor (as.character (.[[!!level]] ),
                                        levels = as.character (tax.ord[[level]])
                                        )
                     )
  }


  ### prepare data for pie charts

  ## normalize OTU counts for each compound (count/concentration-sums) for pie charts ----
#  tax = # to be used by ggplot
#    tax %>%
#      dplyr::group_by (variable) %>%
#      dplyr::mutate_at (dplyr::vars (n), dplyr::funs (./sum (., na.rm = TRUE))) %>%
#      dplyr::ungroup () %>%
      # factorize lipid variable
#      dplyr::mutate (variable = as.factor (variable))

  ## check if normalization was correct ----
#  if (sum (x$n, na.rm = T) !=  length (unique (tax$variable))) {
#    stop ("Normalization per compound failed!")
#  }

  ## relevel so that expemt taxa are always on top ----
  if (is.character(exempt)) {
     y = tax[[level]] # factor with taxa
     for (i in length(exempt):1) { # reverse the order
        if (exempt[i] %in% levels(y)) {
           y = relevel(y, exempt[i])
        }
     }
  tax[[level]] = y
  }


  ## relevel so that 'other' and 'NA' (ukn) are first ----
  if ("other" %in%  levels (tax[[paste (level)]]) == TRUE) {
    tax[[level]] %>%
      relevel(. , "other" ) -> y
    tax[level] = y
  }
  if ("ukn." %in% levels (tax[[paste (level)]]) == TRUE) {
    tax[[level]] %>%
      relevel (., "ukn.") -> y
    tax[level] = y
  }

  ####  make pie charts ####

  ## check for varaible mismatch ----
  if (sum(levels(tax$variable) == levels(ntax$variable)) <
      length(levels(tax$variable))) {
    stop ("Variables jumbled!!")}

  ## make color palett ----
  colourCount = length( unique(tax[[level]]))
  ptb =
    tibble::tibble (randomcoloR::distinctColorPalette (colourCount), levels(tax[[level]]))
  ptb[[1]][which(ptb[[2]] == "other")] = "gray95" ## replace color assignment for other
  ptb[[1]][which(ptb[[2]] == "ukn.")] = "black" ## replace color assignment for ukn.

  ## define plot tile + info ----
  string1 =
    paste (
      n.cor, " out of ", phyloseq::ntaxa(phy), " OTUs analysed yield ", method[1], " correlation coefficents larger than ", thresh,
      "\n" , "with at least one of the compounds/lipids ", "(i.e., ",
      round(n.cor/phyloseq::ntaxa(phy)*100, digits = 1) , "% of all OTUs).",
      "\n" ,
      "Individual compounds correlated well (> ", thresh, ") with ",
      round(min(ntax$rel.tot), digits = 2)*100, "-", round(max(ntax$rel.tot),digits = 2)*100 ,
      "% of the OTUs." ,
      "\n" ,  sep ="")

  if (try (cutoff > 0) == T) {
    string2 =
      paste("Taxa < ", cutoff*100, "% are included in 'other'" ,sep = "")
  } else {
    string2 =""
  }


  if (try (is.character (exempt)) == T) {
    string3 = paste(", except ", paste(exempt, collapse=" and "), ".", sep = "")
  } else {
    string3 = ifelse (cutoff > 0, ".", "")
  }


  ## plot theme ---
  mytheme =
    ggplot2::theme(
      legend.position = "right",
      axis.title = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = text.size*1.25),
      axis.line = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank() ,
      panel.border = ggplot2::element_rect(fill = NA, size = 1),
      axis.ticks.length = ggplot2::unit(.15, "cm") ,
      axis.ticks.y = ggplot2::element_blank() ,
      axis.text = ggplot2::element_blank() ,
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      aspect.ratio = 1,
      plot.title = ggplot2::element_text(lineheight=1, face="plain", size = text.size*1),
      legend.text = ggplot2::element_text(size = text.size)
      )


  ## bar graph
  P2 =
    tax %>%
      ggplot2::ggplot () +
      ggplot2::geom_bar (width = 0.85, stat = "identity",
                         mapping = ggplot2::aes_string( "variable" , "n", fill=paste(level)),
                         position = "fill"
                         ) +
      #ggplot2::facet_wrap ("variable", ncol = n.col) +
      mytheme +
      ggplot2::theme(axis.text = ggplot2::element_text(),
                     panel.border = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_line(),
                     axis.title = ggplot2::element_text(size = text.size)) +
      ggplot2::scale_y_continuous(expand = c(0.0,0)) +
      ggplot2::scale_fill_manual (values = c (ptb[[1]])) +
      ggplot2::geom_text (ggplot2::aes ( variable, value),
                          data=ntax,
                          label = paste ("n=", ntax$n, sep = ""),
                          vjust=-0.7,
                          hjust=0.5,
                          size = text.size/2.75 ) +
      ggplot2::labs (title = paste (string1, string2, string3, "\n\n", sep = "")) +
      ggplot2::coord_cartesian(clip = 'off') +
      ggplot2::xlab("compounds") +
      ggplot2::ylab("fraction")
      #ggplot2::labs(caption = paste (string1, string2, string3, "\n\n", sep = ""))



  ## plot pie chart ----
  P1=
    tax %>%
      ggplot2::ggplot () +
      ggplot2::geom_bar (width = 1, stat = "identity",
                         mapping = ggplot2::aes_string( '""' , "n", fill=paste(level)),
                         position = "fill" )  +
      ggplot2::coord_polar ("y", start=0) +
      ggplot2::facet_wrap ("variable", ncol = n.col) +
      ggplot2::guides (fill=FALSE) +
      ggplot2::scale_fill_manual (values = c (ptb[[1]])) +
      ggplot2::geom_text (ggplot2::aes ( "", value), data=ntax, label = paste ("n=", ntax$n, sep = ""), vjust=0, hjust=0, size = text.size/3 ) +
      ggplot2::labs (title = paste (string1, string2, string3, sep = "")) +
      mytheme # apply the default theme

  ## if specified, override with custom theme
  if (ggplot2::is.theme (plot.theme) == TRUE) {
    P1 = P1 + plot.theme
  }

  ## Return Plots and raw data in a list
  return (
   list(
     bars = P2 + ggplot2::guides(fill = ggplot2::guide_legend (ncol = guide.col)),
     pies = P1 + ggplot2::guides(fill = ggplot2::guide_legend (ncol = guide.col)), # plot with legend in 2 columns
     table = tibble::as.tibble(tax),
     summary = gsub("\n", " ",   paste (string1, string2, string3, sep = ""))
     )
   )

  #### revocer old warning setting ----
  options(warn = oldw)
  print("DONE")

  }

