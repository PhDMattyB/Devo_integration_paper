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

# lmk_dist = geomorph.data.frame(lmk_dist)

## data needs to be a 3d array to subset the data
arrayspecs(wild_univariate_traits, 13, 2)

## the data going in needs to be a 3d array
coords.subset(wild_univariate_traits, 
              wild_lmk_dist$Lake_morph)

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


