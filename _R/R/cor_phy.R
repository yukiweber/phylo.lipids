#' Correlates each lipid with all OTUs
#'
#' This function allows you to extract the phylogeneic affiliations of OTUs, using a theshold correlation coefficient
#' Returns a table with each row representing one OTU, and each column one lipid compound. Values are affiliations of above-threshold OTUs.
#' 'lower' = below threshold, 'other' = sumsumed according to cutoff fraction
#' @name cor_phy
#' @param phy Phyloseq object containing an OTU and taxonomy table.
#' @param lipids Data frame with concentration of lipids (or any other variables/parameters) to correlate with phylogenetic data. Rows are observations/samples and must match the samples in 'phy'. Used by 'phlo_lipids'
#' @param method  correaltion coefficient used in 'cor' {stats}, i.e., "pearson", "kendall", or "spearman"
#' @param level taxonomic level to evaluate, if desplay is set to "all"
#' @param display taxa to be displayed/evaluated. The remining taxa will be subsumed as '_other'. If 'display = "all"',
#' all taxa of the specified 'level' will be displayed (can be slow for lower taxonomic levels).
#' @param shape plot shape of 'display' taxa
#' @param shape.other plot shape of '_other' taxa
#' @param alpha opacitiy of 'display' taxa
#' @param alpha.other opacity of '_other' taxa
#' @param col.other color of '_other' taxa
#' @param colors type of color palette from 'randomcolorR::randomColor',
#' @param luminosity luminosity setting of solor palette, if 'colors = "random"'
# // deactivated// @param runTsne if TRUE, improves color distinction if 'colors = "distinct". may be slower
#' @param flip flips axis; if TRUE, horizontal plot
#' @param box.plot if TRUE, a boxplot will be drawn instead of a jitter plot
#' @param width width of jitter (or boxplot)
#' @param aspect aspect ratio of plot (direction depends of 'flip')
#' @return A list with two data frames containing summaries (max and average coefficients) by taxa (.$taxa) and lipids (.$lipids), respectively, as well as a jitter plot (.$plot).
#' @examples
#' @family Main functions

#'
#' @export

cor_phy = function (phy,
                    lipids,
                    method = c("pearson", "spearman", "kendall"),
                    level = "Phylum",
                    display = "all",
                    shape = 19,
                    shape.other = 1,
                    alpha = 0.75,
                    alpha.other = 0.2,
                    col.other = "gray75",
                    colors = c("random", "distinc"),
                    luminosity = c("random", "light", "bright", "dark"),
                    # runTsne = FALSE, // deactivated//
                    flip = TRUE,
                    box.plot = FALSE,
                    width = 0.25,
                    aspect = 1
                    ) {

  # get taxonomy table ----
  y = tibble::as.tibble (as.data.frame (phy@tax_table, stringsAsFactors = F))
  rownames(phy@tax_table) -> all.taxa
  # add otu names
  y$otu = all.taxa

  # get transposed otu table ----
  x = tibble::as.tibble (t (as.data.frame (phy@otu_table, stringsAsFactors = F)))
  if (sum (rownames (phy@tax_table) == colnames(x)) < length (all.taxa)) {stop ("otu names mismatch")}

  ## input checks ----
  # check if input tax level is valid
  colnames (phy@tax_table) -> ranks
  if (display == "all" & level %in% ranks == F) {
    stop ("You provided the wrong taxonomic level! ",
       paste ("Rank '", level, "' not found in taxonomy table.", sep = ""))
  }

  # check if dimensions match
  if (nrow (lipids) < nrow (x) ) {
    stop ("different number of rows!")
    }

  # Correlate ----
  # get correlation coefficients
  cormat = tibble::as.tibble (cor (x, lipids, method = method[1]))

    # join with tax table
  cormat %>%
    dplyr::mutate (otu = all.taxa) %>%
    dplyr::left_join (y , by = "otu") -> joined

  # long format cormat
  long =
    cormat %>%
      dplyr::mutate (otu = all.taxa) %>%
      tidyr::gather (key = lipids, value = coef , -otu)

  # long format taxonomy
  tax =
    y %>%
      tidyr::gather (key = rank, value = tax , -otu )

  ## join taxonomy and correlations
  DF0 = dplyr::full_join( tax, long, by = "otu")


  # check if taxa name input is correct ----
  display[display %!in% DF0[["tax"]]] -> not.found

  if (display[1] != "all" & length (not.found > 0)) {
  # make string of not found taxa
  a=character()
  count = 1
    while (count <= length(not.found)) {
     a = c(a,
           paste("'",
                 not.found[count],
                 "'",
                 switch(count < length(not.found)-1, ", "),
                 switch(count == length(not.found)-1, ", and "),
                 sep = "")
     )
    count=count+1
    }

  # stop with error message
  stop (
        c (
          round(length (not.found)/length(display)*100, digits = 1),
          "% of the taxa you provided were not found!",
          expression("\n"),
          a,
          " do not exist!"
        )
  )
  }


   ## plot theme ----
  mytheme = ggplot2::theme(
    legend.position = "right", axis.title.x = ggplot2::element_text(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank() ,
    panel.border = ggplot2::element_rect(fill = NA, size = 1),
    axis.ticks.length = ggplot2::unit(.15, "cm") ,
    axis.ticks.y = ggplot2::element_line() ,
    axis.text = ggplot2::element_text(size = 10) ,
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    aspect.ratio = aspect,
    axis.text.x = ggplot2::element_text(angle = 90),
    legend.title = ggplot2::element_blank()
    )

  ## filter taxonomy depending on input ----
  X =
   DF0 %>%
   dplyr::mutate (show = tax)

  # replace with label '_other' if display not set to "all"
  if (display[1] != "all") {
    X =
     X %>%
      dplyr::mutate (show = ifelse(tax %in% display, tax, "_other" ))
  }

  # filter out all other levels/ranks, if level is given display is set to "all"
  if (unique (is.character (level) & display[1] == "all")) {
   X =
     X %>%
      dplyr::filter(., rank == level)
  }

  # order alphabetically to enforce plot order
  X =
    X %>%
     dplyr::arrange(show) %>%
     tibble::as.tibble() -> DF1

  ## assgin plot shapes conditionally ----
  DF2 =
   DF1 %>%
   dplyr::mutate ( plot.cat = ifelse(.["show"] == "_other", "_other", "displayed")) %>%
   tibble::as.tibble()

  # subset to control layer order in ggplot ----
  # factor levels are ignored with multiple layers. therefore use alphabetic order
    l1 = ggplot2::geom_jitter(data = DF2[which(DF2[ ,"plot.cat"] == "_other"), ], width = width, stat = "unique")
    l2 = ggplot2::geom_jitter(data = DF2[which(DF2[ ,"plot.cat"] == "displayed"), ], width = width, stat = "unique")
  # for optional boxplot
  b1 = ggplot2::geom_boxplot(data = DF2 , position= "dodge", width = width)

  # make ggplot ----
  DF2 %>%
   ggplot2::ggplot( ggplot2::aes( lipids,
                                  coef,
                                  color = show,
                                  fill = show,
                                  shape = show,
                                  alpha = show) ) -> P1

  # boxplots instead of jitter ?
  if (box.plot == TRUE) {
    P1 = P1 + b1
  } else {
  # add the jitter layers in order
  P1 = P1 + l1 + l2
  }

  # add other components
  P1 = P1 +
   mytheme  +
   ggplot2::scale_y_continuous (expand = c(0,0),
                                breaks = seq (from = -1, to = 1, by = 0.25),
                                limits = c(-1,1)) +
   ggplot2::xlab("") +
   ggplot2::ylab("corr coef") +
   ggplot2::geom_hline(yintercept = 0, color="gray20") +
   ggplot2::geom_hline(yintercept = 0.5, color="gray20", linetype="dashed") -> P1


  # conditional formatting of plots ----
  # flip axis
  if (flip == TRUE) {
    P1 =
     P1 + ggplot2::coord_flip() -> P1
  }

  ## assign color and shapes and alpha values---
  rand.cols =
    randomcoloR::randomColor (
          ifelse(display[1] == "all",
                 length(unique(DF2[["tax"]])),
                 # else
                 length(display)
          ),
          luminosity = luminosity[1]
        )

  dist.cols =
    randomcoloR::distinctColorPalette (
          ifelse(display[1] == "all",
                 length(unique(DF2[["tax"]])),
                 # else
                 length(display)
          )
        )

  alpha.vals =
   rep (alpha,
        ifelse(display[1] == "all",
               length(unique(DF2[["tax"]])),
               length(display)))

  shape.vals =
   rep (shape,
        ifelse(display[1] == "all",
               length(unique(DF2[["tax"]])),
               length(display)))

    # add '*.other' value to beginning of color/alpha/shape list if display is not "all" ----
  if (display[1] != "all" & "_other" %in% DF2[["show"]]) {
   rand.cols = c(col.other, rand.cols)
   dist.cols = c(col.other, dist.cols)
   alpha.vals = c(alpha.other, alpha.vals)
   shape.vals = c(shape.other, shape.vals)
  }

  # add manual scales ----
  P1 =
    P1 +
      ggplot2::scale_color_manual (
        values =
          dplyr::if_else(
            colors[1] == "random",
            list(rand.cols), # using lists to match length of consition
            # else
            list(dist.cols)
          )[[1]] # get vector from element 1 of list
      ) +
   ggplot2::scale_shape_manual (
          values = shape.vals
   ) +
   ggplot2::scale_alpha_manual (
          values = alpha.vals
   )

  # make summary tables ----
  summary_lipids =
   DF2 %>%
   dplyr::group_by(lipids) %>%
   dplyr::summarise( # n = n(),
                     # min = min(coef),
                     max = max(coef),
                     avg = mean(coef)
   )

  coef_taxa =
   DF2 %>%
   dplyr::group_by(show) %>%
   dplyr::summarise( # n = n(),
                     # min = min(coef),
                     max = round(max(coef), 3),
                     avg = round(mean(coef), 3)) %>%
   dplyr::rename(tax = show)

  # get number of OTUs in for each taxon
  ntaxa =
     tax %>%
        dplyr::mutate_at(dplyr::vars(tax),
                         dplyr::funs(ifelse(
                            . %in% display,
                            .,
                            "_other"))) %>%
        dplyr::filter(tax != "_other") %>%
        dplyr::group_by(tax) %>%
        dplyr::summarise(n = n()) %>%
        tibble::add_row(tax = "_other", n = length(all.taxa) - sum (.$n))

  # merge with coeficcient stats
  summary_taxa =
     dplyr::full_join(ntaxa, coef_taxa, by = "tax")


  # return a list ----
  return(
    list(
     plot = P1,
     lipids = summary_lipids,
     taxa = summary_taxa
     )
  )
}
