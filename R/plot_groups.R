#' Plots vertical profiles of OTU groups with taxonomy summary
#'
#'
#' @param phy Phyloseq object containing at least an OTU and a taxonomy table.
#' @param label A column of 'phy'at'samdata' that is used to label observations (character string),
#' or a vector of same lengh as observations
#' @param cutoff cutoff abundance for assignment of 'other'
#' @param exempt character vector of taxa exempt from cutoff
#' @param level taxonomic level
#' @param ncol Number of columns in the arranged plots/facets
#' @param key_col Number of columns of the color key
#' @param clust Hclust method passed to hclust() in group_tax()
#' @param use_groups Use a previously generated list of otu groups
#' @param prune_groups A character vector with the grop names to be displayed. If NULL, all are shown.
#' @param highlight Taxa to be highlighted
#' @param standardize Method of standardization of the otu abundances,
#' either to the sum or the maximum across all samples
#' @param count Determines whether the statistics are based on the NUMBER of otus or their ABUNDANCE SUMS
#' @export
#' @return A list of summary plots and data


plot_groups = function (phy,
                        level="Phylum",
                        label ="NULL", # alternative labe to sample names
                        use_groups = F,
                        k = 10,
                        h = F,
                        cutoff = 0.03,
                        exempt = FALSE,
                        otus = T, # plot individual otus
                        ribbon =T, # plot ribbons
                        prune_groups = NULL, # plot only selevted groups
                        count = c("n","sum"), # select count type
                        fancy = T, # plots overlay plots
                        highlight = "NULL",
                        standardize = c("sum","max"),
                        ncol = NULL,
                        key_col = 2,
                        clust = c("ward.D2", "ward.D", "single", "complete", "average", "mcquitty", "median", "centroid")
                        ) {


   # input checks ----

   # check if label value is valid
   # and set a flag 'LAB'
   LAB = F
   if( stringr::str_detect(paste(colnames(phy@sam_data), collapse = " "),
                           stringr::regex(label[1], ignore_case = TRUE )
                           ) == T) {
      LAB = T
   } else {
      stop ("Label not found in phy@sam_data!")
      }



   # cluster groups ----
   # if no groups were provided
   if (is.list(use_groups) == F) {
   phy %>%
      group_tax(k = k,
                h = h,
                clust = clust,
                standardize = standardize,
                count = count) -> g
   }

   # otherwise use provided group list
   if (is.list(use_groups) == T) {
      g = use_groups
   }

   # if groups are provided, prune these groups from the list
   if (is.null(prune_groups) == F ) {
      g = g[prune_groups]
      }

   # set ncol to default
   if ( is.null(ncol) == T) {
      ncol = length(g)
   }

   ## extract data from phyloseq object
   # OTU table (rows = taxa)
   x = as.data.frame ( unclass (phy@otu_table) , stringsAsFactors = F) #
      #tibble::as.tibble()
   colnames(x) -> sam_names
   # class 'phyloseq caused 'subset out of bound' error >> unclass
   s = as.data.frame ( unclass( phy@sam_data), stringsAsFactors = F )

   # if alternative label was defined
   # replace sample names with sample data (e.g. "depth)
   # use stringr to ignore case
   sam_var="sample"
   if( LAB == T) {
      sam_var = label[1]
      colnames(x) =
         s[[grep( label[1], colnames(s), ignore.case = T, value = T)]]
   }
head(x)

   class(x)
class(   s[[label]] )
as.data.frame( unclass(s)) ->ss
class(ss)

class(phy@sam_data) ->> cp
colnames(s)
is.data.frame(phy@sam_data)

   # compute abundance SUMs across all samples ----
   sum = Reduce('+', x) # / ncol(x) # prev: averages

   # add the SUMs and otu names to the DF
   x = x %>%
      dplyr::mutate(sum = sum, # previously computed abundance SUMs
                    otu = rownames(.)) %>% # add otu column
      tibble::as.tibble()  # standardize abundances

   # TAX table ----
   # no factors!
   # use 'select' when phyloseq package is not loaded!
   y = as.data.frame (phy@tax_table , stringsAsFactors = F) %>%
      dplyr::select(level)

   # add everything together ----
   x1 = tibble::as.tibble (y) %>%
    dplyr::select(level) %>% # get taxonomy assignments at the specified level
    dplyr::mutate (otu = rownames (.)) %>% # add otu column
    dplyr::mutate (!!level := as.character (.[[level]])) %>% # make sure taxonomy column it's not a factor
    dplyr::left_join (x, by ="otu") # add abundance data

   # gather with 'sam_var'
   x1 %>%
      tidyr::gather(key = !!sam_var, value = "abundance", -otu, -sum, -!!level ) %>%
      tibble::as.tibble() -> x2

   # assign groups to otus
   x2$group=""
   for (i in names (g) ) {
      x2$group[which(x2$otu %in% g[[i]])] = i }

   # drop all otus without group assignement
   x2 = x2 %>%
      dplyr::filter(group != "")

   # determine cutoff taxa for each group ----
   for (i in names(g)) {
      a = x2 %>%
         dplyr::filter(group==i) %>%
         dplyr::group_by_(level) %>%
         dplyr::summarise(n = n(),
                          sum = sum(sum)) %>%
         # relative counts depending on the count type
         dplyr::mutate(rel.cnt = get(count[1]) / sum(get(count[1])))

      # get the below-cutoff taxa as character vector
      a[[level]][which(a$rel.cnt < cutoff) ] -> cut
      # only keep taxa in the cutoff list that are not in exempt
      cut = cut[which(cut %!in% exempt)]
      # replace tax assignements of below cutoff taxa with 'other'
      x2[level][which(x2$group == i & x2[[level]] %in% cut ), ] = "other"
   }

   # replace NA with 'ukn.'
   x2 %>%
      dplyr::mutate_at(dplyr::vars(level),
                       dplyr::funs(ifelse( is.na(.) | .=="NA",
                                         "ukn."
                                         , .)
                                   )
                       ) -> x3

   # make taxonomy a factor and re-level ----
   tax = as.factor(x3[[level]])
   # put expemt taxa on the top
   # if they are present after cutoff
   if (is.character(exempt)) {
      for(i in rev(exempt)) {
         if ( i %in% levels(tax) ){
            tax = relevel(tax , i) } } }
   # show 'other' second
   if ("other" %in% levels(tax) ) {
      tax = relevel(tax ,"other") }
   # show ukn first
   if ("ukn." %in% levels(tax) ) {
      tax = relevel(tax ,"ukn.") }
   # replace factor in x3
   x3[level] = tax


   # if the label data is numeric, convert to numeric (otherwise character!)
   if (is.numeric(phy@sam_data[[label[1]]]) == T) {
      x3[label] = as.numeric(x3[[label]])
   }

   # color assignments ----
   colourCount = length( unique(x3[[level]]))
   ptb =
      tibble::tibble (randomcoloR::distinctColorPalette (colourCount), levels(x3[[level]]))
   ptb[[1]][which(ptb[[2]] == "other")] = "gray90" ## replace color assignment for other
   ptb[[1]][which(ptb[[2]] == "ukn.")] = "black" ## replace color assignment for ukn.
   colnames(ptb) = c("col",paste(level))
   # make sure the factor levels are the same order
   ptb[level] = factor(ptb[[level]], levels = levels(x3[[level]]))
   # add color assignements to otu data
   x3 %>%
      dplyr::left_join(ptb, by= paste(level)) -> x3

   # if highlight was defined
   # and all provided taxa are valid at the given level
   # replace all other tax colors with gray
   if (highlight[1] != "NULL" &
       sum ( highlight %in% y[[level]]) == length(highlight)) {
      x3["col"][which(
         ! x3[[level]] %in% c( highlight, "other" , "ukn.") ), ] = "white"
   }


   ## THEMES ----
   # balnk theme for pies
   blank_theme = ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
   )

   # balnk theme for line plots
   blank_theme_x = ggplot2::theme(
      axis.line.x =  ggplot2::element_blank(),
      axis.line.y =  ggplot2::element_line(),
      axis.text.x =  ggplot2::element_blank(),
      axis.title.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "gray95"),
      plot.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", size=1, fill = NA),
      strip.text = ggplot2::element_text(size = 9),
      aspect.ratio = 2.5 )

   ### fancy overlay plots ----

   # prepare a list to hold the plots
   L = list()

   if (fancy == TRUE) {
   # for each group:
   for (i in names(g)) {

      # make a distribution line plot
      plot_a_group(taxnames = g[[i]],
                   phy = phy,
                   label = label,
                   otus = otus,
                   ribbon = ribbon,
                   is.numeric.label = num) +
      blank_theme_x +
      ggplot2::labs(title = i)   -> a

      # and a pie chart
      x3 %>%
         plot_pie(x = .,
                  taxa_groups = g[[i]],
                  level = level,
                  count = count[1]) -> bb
      bb$pie +
         blank_theme +
         ggplot2::guides(fill=F) -> b

      # put the two together
      cowplot::ggdraw() +
      cowplot::draw_plot(a) +
         cowplot::draw_plot(b,
                            x = 0.5, #The x location of the lower left corner of the plot.
                            y = 0.3, #The y location of the lower left corner of the plot.
                            width = 0.5,
                            height = 0.5,
                            scale = 1
                            ) -> c

      # add the grobs to a list
      L[[i]] =  c
   }
   }  # END of fancy branch

   ## color legend ----
   # get color values of ALL taxa
   x3 %>% # all data
   dplyr::group_by_(level) %>%
   dplyr::summarise(n=n(),
                    sum = sum(sum),
                    col = unique(col)) -> x5

   # dummy plot for color legend, from whole data set
   x5 %>%
   ggplot2::ggplot(ggplot2::aes_string(x= paste(level), y = count[1] , fill = level)) +
      ggplot2::scale_fill_manual(values = x5$col ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::guides(fill = ggplot2::guide_legend (ncol = key_col,
                                                    title = paste0(level," (", count[1], ")"))
                      ) ->  d

   ## extract legend from dummy plot
   # helper function to extract only the legend from a ggplot
   g_legend <- function(a.gplot){
     tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
     leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
     legend <- tmp$grobs[[leg]]
   return(legend)
   }

   # extract legend ----
   leg = g_legend(d) %>%
      # convert to ggplot
      ggplotify::as.ggplot()


   ### simple facet wraps ----
   Fg =
      plot_a_group(taxnames = g,
                   phy = phy,
                   label = label,
                   otus = otus,
                   ribbon = ribbon,
                   ncol = ncol,
                   is.numeric.label = num) +
         blank_theme_x

   Fp =
      x3 %>%
         plot_pie(x = .,
                  count = count,
                  level = level,
                  taxa_groups = NULL,
                  ncol = ncol,
                  key_col = key_col
                  )
   Fp = Fp$pie + blank_theme


   ### outputs ----

   # arrage fancy plots ----
   A = ""
   B = ""
   if (fancy==T) {
      # plots only
      gridExtra::grid.arrange(grobs = c(L),
                              ncol = ncol,
                              newpage = F
                              ) %>%
         ggplotify::as.ggplot() -> A

      ## as a single figure with legend
      gridExtra::grid.arrange(A,leg,
                              newpage = F) %>%
         ggplotify::as.ggplot() -> B
            }

   # make a summary stat list ----
   summary = list()
   for (i in names(g)) {
     summary[[i]] =
       x3 %>%
         dplyr::filter(group == i) %>%
         dplyr::group_by_(level) %>%
         dplyr::summarise(n = n(),
                          sum = sum(sum)) %>%
        dplyr::mutate(rel.n = n/sum(n),
                      rel.sum = sum/sum(sum),
                      group = i) %>%
        dplyr::select(-sum)
   }

   ## return a list ----
   return(
      list(
         "groups_pies_only" = A,
         "legend_only" = leg,
         "groups_pies_full" = B,
         "group_dat" = g,
         "facet_pies" = Fp,
         "facet_groups" = Fg,
         "summary" = summary
         )
   )
}

