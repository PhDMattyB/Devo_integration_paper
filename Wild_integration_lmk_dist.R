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

ASHN_compare_wild = compare.ZVrel(vrel_wild_lmkdist$ASHNC, 
                             vrel_wild_lmkdist$ASHNW)
MYV_compare_wild = compare.ZVrel(vrel_wild_lmkdist$MYVC, 
                            vrel_wild_lmkdist$MYVW)

SKR_compare_wild = compare.ZVrel(vrel_wild_lmkdist$SKRC, 
                            vrel_wild_lmkdist$SKRW)

GTS_CSWY_compare_wild = compare.ZVrel(vrel_wild_lmkdist$CSWY, 
                                 vrel_wild_lmkdist$GTS)


# F2 Uncorrected integration ----------------------------------------------


F2_raw_lmk_dist = read_csv('F2_Original_univariate_traits.csv')
F2_raw_dist = F2_raw_lmk_dist %>% 
  select(2:29)

# lmk_dist = geomorph.data.frame(lmk_dist)
F2_raw_lmk_matrix = as.matrix(F2_raw_dist)
F2_raw_lmk_array = arrayspecs(F2_raw_lmk_matrix, 14, 2)

F2_raw_lmk_sub = coords.subset(F2_raw_lmk_array, 
                           F2_raw_lmk_dist$Lake_morph)


vrel_F2_raw_lmkdist = Map(function(x) integration.Vrel(x), 
                      F2_raw_lmk_sub)

ASHN_compare_raw = compare.ZVrel(vrel_F2_raw_lmkdist$ASHNC, 
                             vrel_F2_raw_lmkdist$ASHNW)
MYV_compare_raw = compare.ZVrel(vrel_F2_raw_lmkdist$MYVC, 
                            vrel_F2_raw_lmkdist$MYVW)

SKR_compare_raw = compare.ZVrel(vrel_F2_raw_lmkdist$SKRC, 
                            vrel_F2_raw_lmkdist$SKRW)

GTS_CSWY_compare_raw = compare.ZVrel(vrel_F2_raw_lmkdist$CSWY, 
                                 vrel_F2_raw_lmkdist$GTS)


##
# F1 effect mag integration -----------------------------------------------
F1_lmk_dist = read_csv('F1_Plasticity_Corrected.csv')
F1_dist = F1_lmk_dist %>% 
  select(2:29)

# lmk_dist = geomorph.data.frame(lmk_dist)
F1_lmk_matrix = as.matrix(F1_dist)
F1_lmk_array = arrayspecs(F1_lmk_matrix, 14, 2)

F1_lmk_sub = coords.subset(F1_lmk_array, 
                             F1_lmk_dist$Lake_morph)


vrel_F1_lmkdist = Map(function(x) integration.Vrel(x), 
                        F1_lmk_sub)

ASHN_compare_f1 = compare.ZVrel(vrel_F1_lmkdist$ASHNC, 
                             vrel_F1_lmkdist$ASHNW)
MYV_compare_f1 = compare.ZVrel(vrel_F1_lmkdist$MYVC, 
                            vrel_F1_lmkdist$MYVW)

SKR_compare_f1 = compare.ZVrel(vrel_F1_lmkdist$SKRC, 
                            vrel_F1_lmkdist$SKRW)

GTS_CSWY_compare_f1 = compare.ZVrel(vrel_F1_lmkdist$CSWY, 
                                 vrel_F1_lmkdist$GTS)


# F2 effect mag integration -----------------------------------------------

F2_lmk_dist = read_csv('F2_Corrected_F2_temp_only.csv')
F2_dist = F2_lmk_dist %>% 
  select(2:29)

# lmk_dist = geomorph.data.frame(lmk_dist)
F2_lmk_matrix = as.matrix(F2_dist)
F2_lmk_array = arrayspecs(F2_lmk_matrix, 14, 2)

F2_lmk_sub = coords.subset(F2_lmk_array, 
                           F2_lmk_dist$Lake_morph)


vrel_F2_lmkdist = Map(function(x) integration.Vrel(x), 
                      F2_lmk_sub)

ASHN_compare_f2 = compare.ZVrel(vrel_F2_lmkdist$ASHNC, 
                             vrel_F2_lmkdist$ASHNW)
MYV_compare_f2 = compare.ZVrel(vrel_F2_lmkdist$MYVC, 
                            vrel_F2_lmkdist$MYVW)

SKR_compare_f2 = compare.ZVrel(vrel_F2_lmkdist$SKRC, 
                            vrel_F2_lmkdist$SKRW)

GTS_CSWY_compare_f2 = compare.ZVrel(vrel_F2_lmkdist$CSWY, 
                                 vrel_F2_lmkdist$GTS)
