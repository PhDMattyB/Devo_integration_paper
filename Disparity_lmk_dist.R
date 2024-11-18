##############################
## Disparity analysis linear distances
##
## Matt Brachmann (PhDMattyB)
##
## 18.11.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(tidyverse)
library(geomorph)


# wild disparity ----------------------------------------------------------

wild_lmk_dist = read_csv('Wild_Univariate_traits.csv') %>% 
  filter(Lake %in% c('ASHN', 
                     'MYV', 
                     'SKR', 
                     'GTS', 
                     'CSWY'))

wild_dist = wild_lmk_dist %>% 
  select(2:29)

wild_mat = as.matrix(wild_dist)

Wild_disparity = morphol.disparity(wild_mat ~ 1,
                                       groups = ~Lake_morph,
                                       data = wild_lmk_dist,
                                       iter = 999)

wild_pval = Wild_disparity$PV.dist.Pval
wild_disparity_dist = Wild_disparity$PV.dist
wild_proc_var = Wild_disparity$Procrustes.var

