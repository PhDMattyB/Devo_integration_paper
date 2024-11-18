##############################
## wild vrel interlandmark dist
##
## Matt Brachmann (PhDMattyB)
##
## 18.11.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(tidyverse)
library(geomorph)


# wild traits mag integration ---------------------------------------------

wild_lmk_dist = read_csv('Wild_Univariate_traits.csv')

wild_dist = wild_lmk_dist %>% 
  select(2:29)

# lmk_dist = geomorph.data.frame(lmk_dist)
wild_lmk_matrix = as.matrix(wild_dist)
wild_lmk_array = arrayspecs(wild_lmk_matrix, 14, 2)

wild_lmk_sub = coords.subset(wild_lmk_array, 
              wild_lmk_dist$Lake_morph)


vrel_wild_lmkdist = Map(function(x) integration.Vrel(x), 
                       wild_lmk_sub)

ASHN_compare = compare.ZVrel(vrel_wild_lmkdist$ASHNC, 
                             vrel_wild_lmkdist$ASHNW)
MYV_compare = compare.ZVrel(vrel_wild_lmkdist$MYVC, 
                            vrel_wild_lmkdist$MYVW)

SKR_compare = compare.ZVrel(vrel_wild_lmkdist$SKRC, 
                            vrel_wild_lmkdist$SKRW)

GTS_CSWY_compare = compare.ZVrel(vrel_wild_lmkdist$CSWY, 
                                 vrel_wild_lmkdist$GTS)


# F1 effect mag integration -----------------------------------------------



# F2 effect mag integration -----------------------------------------------


