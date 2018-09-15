#' Plots vertical profiles of OTU groups with taxonomy summary
#'
#'
#' @param phy Phyloseq object containing at least an OTU and a taxonomy table.
#' @param x list of n groups containing otu names of each group
#' @param label A column of 'phy'at'samdata' that is used to label observations (character string), or a vector of same lengh as observations
#' @export
#' @examples
#' @return a summary plot


plot_groups = function (phy,
                        level="Phylum",
                        label ="NULL", # alternative labe to sample names
                        k = 10,
                        cutoff = 0.03,
                        exempt = FALSE,
                        otus = T, # plot individual otus
                        ribbon =T, # plot ribbons
                        prune_groups = NULL, # plot only selevted groups
                        count = c("n","avg"), # select count type
                        fancy = T, # plots overlay plots
                        highlight = "NULL"
                        ) {


   # cluster groups
   phy %>%
      group_tax(k = k) -> g0
   g = g0

   # if groups are pecified, prune these groups from the list
   if (is.null(prune_groups) == F ) {
      g = g0[prune_groups]
      }

   ## extract data from phyloseq object
   # OTU table (rows = taxa)
   x = as.data.frame (phy@otu_table, stringsAsFactors = F)
   colnames(x) -> sam_names

   # replace sample names sample data (e.g. "depth)
   sam_var="sample"
   if(label[1] %in% colnames(phy@sam_data)){
         colnames(x) = phy@sam_data[[label]]
         sam_names = phy@sam_data[[label]]
         sam_var = label[1]
   }

   # compute average abundance across all samples
   avg = Reduce('+', x) / ncol(x)
   # standardize abundances
   x = as.data.frame(sapply(x,`/`,rowSums(x)))
   x = x %>%
      dplyr::mutate(avg = avg, # compute average abundance across samples
                    otu = rownames(.)) %>% # add otu column
      tibble::as.tibble()  # standardize abundances

   # TAX table
   # no factors!
   # use 'select' when phyloseq package is not loaded!
   y = as.data.frame (phy@tax_table , stringsAsFactors = F) %>%
      dplyr::select(level)

   # add everything together
   x1 = tibble::as.tibble (y) %>%
    dplyr::select(level) %>% # get taxonomy assignments at the specified level
    dplyr::mutate (otu = rownames (.)) %>% # add otu column
    dplyr::mutate (!!level := as.character (.[[level]])) %>% # make sure taxonomy column it's not a factor
    dplyr::left_join (x, by ="otu") # add abundance data

   # gather
   x1 %>%
      tidyr::gather(key = !!sam_var, value = "abundance", -otu, -avg, -!!level ) %>%
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
                          avg = sum(avg)) %>%
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
   ptb[[1]][which(ptb[[2]] == "other")] = "gray95" ## replace color assignment for other
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
         ! x3[[level]] %in% c( highlight, "other" , "ukn.") ), ] = "gray95"
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
                   ribbon = ribbon) +
      blank_theme_x +
      ggplot2::labs(title = i)   -> a

      # and a pie chart
      x3 %>%
         plot_pie(x = .,
                  taxa_groups = g[[i]],
                  level = level,
                  count = count[1]) +
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
                    avg = sum(avg),
                    col = unique(col)) -> x5

   # dummy plot for color legend, from whole data set
   x5 %>%
   ggplot2::ggplot(ggplot2::aes_string(x= paste(level), y = "avg" , fill = paste(level))) +
      ggplot2::scale_fill_manual(values = x5$col ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::guides(fill = ggplot2::guide_legend (ncol = 3,
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

   # extract legend
   dd = g_legend(d)
   # convert to ggplot
   ggplotify::as.ggplot(dd) -> D
   ggpubr::as_ggplot(dd) -> D1


   ### simple facet wraps
   Fg =
      plot_a_group(taxnames = g,
                   phy = phy,
                   label = label,
                   otus = otus,
                   ribbon = ribbon) +
         blank_theme_x

   Fp =
      x3 %>%
         plot_pie(x = .,
                  count = count,
                  level = level,
                  taxa_groups = NULL
                  ) +
         blank_theme


   ### outputs ----

   # arrage plots
   A = ""
   B = ""
   if (fancy==T) {
      # plots only
      gridExtra::grid.arrange(grobs = c(L), ncol = round(length(L)/1.8, 0)) %>%
         ggplotify::as.ggplot() -> A

      ## as a single figure with legend
      # get the table grob from legend
      dd$grobs[1] -> ddd
      # place it at the end of grob list
      L1 = c(L,ddd)
      # plot together
      gridExtra::grid.arrange(
         ggplotify::as.grob(gridExtra::grid.arrange(grobs = L1)),
         ncol = round(length(L1)/1.8, 0)
         ) %>%
         ggplotify::as.ggplot() -> B
      }


   return(
      list(
         "groups_pies_only" = A,
         "legend_only" = D,
         "groups_pies_full" = B,
         "group_dat" = g,
         "facet_pies" = Fp,
         "facet_groups" = Fg
         )
   )
}

