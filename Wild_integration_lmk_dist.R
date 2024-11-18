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

wild_lmk_dist = read_csv('Wild_Univariate_traits.csv')

lmk_dist = wild_lmk_dist %>% 
  select(2:28)

test = as.matrix(lmk_dist)

arrayspecs(test, 27, 2)

plethodon$land    #original data in the form of 3D array

two.d.array(plethodon$land)   

x2 = matrix(rnorm(18), ncol = 2) 
arrayspecs(x2, 3, 2)


coord_sub = coords.subset(lmk_dist, 
                          wild_lmk_dist$Lake_morph)

subset_wild_coords = coords.subset(wild_gpa$coords, 
                                   wild_identifiers$Lake_morph)

vrel_wild_coords = Map(function(x) integration.Vrel(x), 
                       subset_wild_coords)

ASHN_compare = compare.ZVrel(vrel_wild_coords$ASHNC, 
                             vrel_wild_coords$ASHNW)

MYV_compare = compare.ZVrel(vrel_wild_coords$MYVC, 
                            vrel_wild_coords$MYVW)

SKR_compare = compare.ZVrel(vrel_wild_coords$SKRC, 
                            vrel_wild_coords$SKRW)

# STN_compare = compare.ZVrel(vrel_wild_coords$STNC, 
#                             vrel_wild_coords$STNW)
# 
# RKLT_compare = compare.ZVrel(vrel_wild_coords$RKLTC, 
#                             vrel_wild_coords$RKLTW)

GTS_CSWY_compare = compare.ZVrel(vrel_wild_coords$CSWY, 
                                 vrel_wild_coords$GTS)


