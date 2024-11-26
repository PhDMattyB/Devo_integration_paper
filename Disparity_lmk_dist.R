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
  dplyr::select(2:29)

wild_mat = as.matrix(wild_dist)

Wild_disparity = morphol.disparity(wild_mat ~ 1,
                                       groups = ~Lake_morph,
                                       data = wild_lmk_dist,
                                       iter = 999)

wild_pval = Wild_disparity$PV.dist.Pval
wild_disparity_dist = Wild_disparity$PV.dist
wild_proc_var = Wild_disparity$Procrustes.var

wild_disp_pval = wild_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()


wild_disp_dist = wild_disparity_dist %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  # as_tibble() %>% 
  melt(id.vars = c('rowname')) %>% 
  as_tibble()

wild_disp_data = bind_cols(wild_disp_dist,
                           wild_disp_pval) %>% 
  dplyr::select(1:3, 
         6) %>% 
  rename(Ecotype1 = 1, 
         Ecotype2 = 2, 
         zscore = 3, 
         pvalue = 4) %>% 
  mutate(across(where(is.numeric),
                ~ round(., 3))) %>% 
  mutate(label = 'Grandparental (wild) generation')



# F2 uncorrected data -----------------------------------------------------

F2_raw_lmk_dist = read_csv('F2_Original_univariate_traits.csv')
F2_raw_dist = F2_raw_lmk_dist %>% 
  dplyr::select(2:29)

F2_raw_mat = as.matrix(F2_raw_dist)

F2_raw_disparity = morphol.disparity(F2_raw_mat ~ 1,
                                   groups = ~Lake_morph,
                                   data = F2_raw_lmk_dist,
                                   iter = 999)

F2_raw_pval = F2_raw_disparity$PV.dist.Pval
F2_raw_disparity_dist = F2_raw_disparity$PV.dist
F2_raw_proc_var = F2_raw_disparity$Procrustes.var

# TGP effect on disparity --------------------------------------------------
F1_lmk_dist = read_csv('F1_Plasticity_Corrected.csv')
F1_dist = F1_lmk_dist %>% 
  dplyr::select(2:29)

TGP_mat = as.matrix(F1_dist)

TGP_disparity = morphol.disparity(TGP_mat ~ 1,
                                     groups = ~Lake_morph,
                                     data = F2_raw_lmk_dist,
                                     iter = 999)

TGP_pval = TGP_disparity$PV.dist.Pval
TGP_disparity_dist = TGP_disparity$PV.dist
TGP_proc_var = TGP_disparity$Procrustes.var


# WGP effect on disparity -------------------------------------------------
F2_lmk_dist = read_csv('F2_Corrected_F2_temp_only.csv')
F2_dist = F2_lmk_dist %>% 
  dplyr::select(2:29)

WGP_mat = as.matrix(F2_dist)

WGP_disparity = morphol.disparity(WGP_mat ~ 1,
                                  groups = ~Lake_morph,
                                  data = F2_raw_lmk_dist,
                                  iter = 999)

WGP_pval = WGP_disparity$PV.dist.Pval
WGP_disparity_dist = WGP_disparity$PV.dist
WGP_proc_var = WGP_disparity$Procrustes.var




