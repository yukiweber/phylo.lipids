# phylo.lipids

This R package was developed to investigate quantitative links between lipid biomarker concentration and meta taxonomic (16S rRNA gene) data.

### Version 0.0.0.900 

This is a development version and was used to generate some of the results published in:  

Yuki Weber, Jaap S. Sinninghe Damsté, Jakob Zopfi, Cindy De Jonge, Adrian Gili, Carsten J. Schubert, Fabio Lepori, Moritz F. Lehmann, Helge Niemann (2018): *"Redox-dependent niche differentiation of tetraether lipid producing bacteria: Evidence for multiple branched GDGT sources in lakes"*.  

More features will be added soon. Please report any bugs to: yuki.t.weber@gmail.com


#### Installation of latest source package
```{r}
install.packages("devtools")
devtools::install_github("yukiweber/phylo.lipids", build_vignettes = F)
```

#### Reproducing published analysis
To reproduce the results published in the above paper:

1.) Install 'R', 'RStudio', and 'phylo.lipids' (see above) on your machine. 

2.) Download Supplementary Information Datasets 1 and 2 from the publisher's website and place files into the root of the RStudio project folder.

3.) Run the code in the R markdown file *'Lake_Lugano.Rmd'* in Rstudio.
