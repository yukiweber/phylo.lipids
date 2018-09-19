#' helper function: plots standardized otu abundances
#'
#' @param x character vector with otu names
#' @param phy Phyloseq object
#' @param taxnames A (list of) char vector(s) containing otu names to be plotted
#' If list, output will be facetted by the elements of the list.
#' @param label a element of phy'at'sam_data to replace the original sample names
#' @param otus display individual otus?
#' @param ribbon Plot error bounds (+/- 2*SD)
#' @param ncol Number of columns in facets
#' @param wt.means Plot weighted means by abundance of the otus (!!! under development !!!)
#' @param standardize Method of standardization of the otu abundances,
#' either to the sum or the maximum across all samples
#' @param is.numeric.label If TRUE, force coertion of labels to numeric
#' @export
#' @return ggplot line plot

# x character vector with otu names
plot_a_group = function (phy,
                         taxnames,
                         label = "NULL",
                         otus = TRUE,
                         ribbon = TRUE,
                         ncol = 5,
                         is.numeric.label = F,
                         standardize = c("sum","max"),
                         wt.means = F # under development
                         ) {

   # check if taxa names are valid
   if (mean(Reduce(c, taxnames) %in% phyloseq::taxa_names(phy)) != 1) {
      stop("wrong taxa names provided!") }

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

   # prune all taxa in taxnames
   phyloseq::prune_taxa(Reduce(c, taxnames), phy) -> P

   # get otu table
   as.data.frame(unclass(P@otu_table), strngsAsFactors = F) -> otu
   # standardize abundances and sums
   OT = standardize_tax(P, method = standardize) #%>%
      #dplyr::select(-sum)

   # get sam_data
   # class 'phyloseq may cause 'subset out of bound' error >> unclass
   s = as.data.frame ( unclass( P@sam_data), stringsAsFactors = F )

   # replace sample names sample data (e.g. "depth)
   sam_var="sample"
   if(LAB == T){
      sam_var = label
      # disregard the sum column
      colnames(OT)[1:ncol(OT)-1] =
         s[[grep(label[1], colnames(s), ignore.case = T, value = T)]]
      sam_labs = s[[grep(label[1], colnames(s), ignore.case = T, value = T)]]
      }

   # add otu column to otu table
   OT$otu=rownames(OT)

   ## assign groups to otus
   OT$group=""
   # if a list of groups was provided
   if (is.list(taxnames)==T) {
   for (i in names (taxnames) ) {
      OT$group[which(OT$otu %in% taxnames[[i]])] = i
   }
   # drop all otus without group assignement
   OT = OT %>%
      dplyr::filter(group != "")
   }

   # gather data ----
   OT %>%
      tidyr::gather(key = !!sam_var, value = "abundance", -otu, -group, -sum) %>%
      tibble::as.tibble() -> OT1

   # if alternative labels are numeric, coerce to muneric
   if (is.numeric(sam_labs) || is.numeric.label == T ) {
      OT1[[sam_var]] = as.numeric(OT1[[sam_var]])
   }

   # factorize groups to preserve order
   OT1$group = factor (OT1$group, levels = names(taxnames))

   OT0=OT1

   ### add weighted abundances ----
   # if taxnames was not a list, set dummy group to 'x'

#   if(wt.means == T) {

   if (length(levels(OT1$group)) == 0) {
      OT1$group = factor ("x", levels=c("x"))
   }
   # get group sums
   mm = as.data.frame(matrix(nrow = 0,ncol = 8))
   colnames(mm) = c("otu"   ,    "sum"   ,    "group"  ,   "depth"  ,   "abundance", "g.sum"   ,  "g.max.abu", "g.rel.sum")
   for (i in levels(OT1$group)) {
      OTa = OT1 %>%
            dplyr::filter(group==i)
      OTx =
         OT1 %>%
            dplyr::filter(group==i) %>%
            dplyr::group_by_("otu") %>%
            dplyr::summarise(g.sum = sum(sum),  # absolute sum of sums within group
                             g.max.abu = max(abundance)
            ) %>%
            dplyr::mutate(g.rel.sum = g.sum / sum(g.sum))
      mmx= merge(OTa,OTx, by = "otu")
      #dim(mmx)
      #dim(OTa)
      mm = rbind(mm,mmx)
   }

   mm = tibble::as.tibble(mm)
   ## levels(mm$group)

   OT1 = mm
   #OT1 ->> OT1

   # compute weighted abundances within groups
   OT1 =
      OT1 %>%
      dplyr::mutate(wt.abu = abundance * sum / g.sum / g.max.abu)
 #  }


   #########################
   ## initialize a plot ----
   OT1 %>%
   ggplot2::ggplot( ggplot2::aes_string(y= "abundance",
                                        x = sam_var,
                                        group = "otu"
                                        )
                    ) -> p

   # plus/minus SD Ribbon
   if (ribbon == TRUE) {
      p = p +
         ggplot2::stat_summary(fun.ymax =  function(y) mean(y) + mean(2*sd(y)) ,
                               fun.ymin =  function(y) mean(y) - mean(2*sd(y)) ,
                               geom = "ribbon",
                               ggplot2::aes(group = 1), # character vector for samples!
                               fill= "red",
                               alpha = 0.5)
   }

   # individual otus ----
   if (otus == TRUE) {
      p = p +
         ggplot2::geom_line(color="gray30",
                            size = 0.5,
                            ggplot2::aes(alpha = g.rel.sum)
                            #ggplot2::aes(alpha=sum)
                            )
   }

   # mean line ----
   p = p +
      ggplot2::stat_summary(fun.y="mean",
                            geom = "line",
                            ggplot2::aes(group = 1),
                            size=1,
                            alpha = 0.5,
                            color="red") +
      #ggplot2::coord_cartesian(ylim=c(0, 1), expand = c(0,0)) +
      ggplot2::scale_y_continuous(expand = c(0,0))+
      ggplot2::coord_flip(ylim = c(0, max(OT1$abundance)), clip = "on")

   # weighted mean line ----
   if (wt.means == T) {
   p = p +
      ggplot2::stat_summary(ggplot2::aes_string(y= "wt.abu", x = sam_var, group =1),
                            fun.y= function(y) mean(y),
                            geom = "line",
                            size=1.2,
                            color="darkblue",
                            data=OT1)
   }


   # if another label was specified
   if( stringr::str_detect(paste(colnames(phy@sam_data), collapse = " "),
                              stringr::regex(label[1], ignore_case = TRUE )
                              ) == T) {
         p = p + ggplot2::xlab(label=label)
         }

   # if observations are depths and numeric, reverse scale
   if(length(grep("depth", sam_var, ignore.case = T)) >= 1
      & is.numeric(OT1[[sam_var]]) == TRUE) {
      p = p + ggplot2::scale_x_reverse(expand = c(0,0),
                                       limits = c(max(OT1[[sam_var]]), 0))
      }

   p + ggplot2::ylab(label="standardized abundacne") -> p

   # if a list was provided facet by groups
   if (is.list(taxnames) == T) {
      p = p +
         ggplot2::facet_wrap("group",
                             scales = "fixed",
                             ncol = ncol )
   }
   return(p)
}

