setwd("~/Parsons_postdoc/Stickleback_Morphometric_data/Updated Landmarks/")

library(tidyverse)
library(lme4)
library(lmerTest)
library(brms)
library(paran)
library(broom.mixed)
library(MatrixCorrelation)
library(vegan)
library(evolqg)


### Need to compare the F2, WGP, and TGP covariance matrices 

ashnc_f2_tgp_12_12_proc_test = protest(ashnc_f2_12_12_cov_mat, 
                                       ashnc_TGP_12_12_cov_mat, 
                                       permutations = 999)

ashnc_f2_wgp_12_12_proc_test = protest(ashnc_f2_12_12_cov_mat, 
                                       ashnc_WGP_12_12_cov_mat, 
                                       permutations = 999)

ashnc_tgp_wgp_12_12_proc_test = protest(ashnc_TGP_12_12_cov_mat, 
                                       ashnc_WGP_12_12_cov_mat, 
                                       permutations = 999)


ashnw_f2_tgp_12_12_proc_test = protest(ashnw_f2_12_12_cov_mat, 
                                       ashnw_TGP_12_12_cov_mat, 
                                       permutations = 999)

ashnw_f2_wgp_12_12_proc_test = protest(ashnw_f2_12_12_cov_mat, 
                                       ashnw_WGP_12_12_cov_mat, 
                                       permutations = 999)

ashnw_tgp_wgp_12_12_proc_test = protest(ashnw_TGP_12_12_cov_mat, 
                                        ashnw_WGP_12_12_cov_mat, 
                                        permutations = 999)



# kendall test three matrices ---------------------------------------------



ashn_f2_12_12_kendal = ashn_f2_12_12_residuals %>% 
  dplyr::rename(ASHN_F2 = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

# ashn_WGP_12_12_kendal = ashn_WGP_12_12_residuals %>% 
#   dplyr::rename(ASHN = residual) %>% 
#   # rownames_to_column() %>% 
#   # mutate(Pop_pair = 'ASHN') %>% 
#   dplyr::select(-resid_10) %>% 
#   column_to_rownames('trait')

ashn_TGP_12_12_kendal = ashn_TGP_12_12_residuals %>% 
  dplyr::rename(ASHN_TGP = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'ASHN') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')



f2_12_12_kendal_df = bind_cols(ashn_f2_12_12_kendal, 
                               ashn_TGP_12_12_kendal)

group = c('ASHN_F2','ASHN_TGP')
group = c(1,2)


f2_12_12_kendal = kendall.global(f2_12_12_kendal_df,
                                 # group = group,
                                 nperm = 999, 
                                 mult = 'hochberg')
f2_12_12_kendal_post = kendall.post(f2_12_12_kendal_df, 
                                    mult = 'hochberg', 
                                    nperm = 999)




myv_f2_12_12_kendal = myv_f2_12_12_residuals %>% 
  dplyr::rename(myv_F2 = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'myv') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')

myv_WGP_12_12_kendal = myv_WGP_12_12_residuals %>%
  dplyr::rename(myv_WGP = residual) %>%
  # rownames_to_column() %>%
  # mutate(Pop_pair = 'myv') %>%
  dplyr::select(-resid_10) %>%
  column_to_rownames('trait')

myv_TGP_12_12_kendal = myv_TGP_12_12_residuals %>% 
  dplyr::rename(myv_TGP = residual) %>% 
  # rownames_to_column() %>% 
  # mutate(Pop_pair = 'myv') %>% 
  dplyr::select(-resid_10) %>% 
  column_to_rownames('trait')



myv_f2_12_12_kendal_df = bind_cols(myv_f2_12_12_kendal,
                               myv_WGP_12_12_kendal,
                               myv_TGP_12_12_kendal)

myv_f2_12_12_kendal = kendall.global(myv_f2_12_12_kendal_df,
                                 # group = group,
                                 nperm = 999, 
                                 mult = 'hochberg')
myv_f2_12_12_kendal_post = kendall.post(myv_f2_12_12_kendal_df, 
                                    mult = 'hochberg', 
                                    nperm = 999)

myv_f2_12_12_nemeny_df = myv_f2_12_12_residuals %>% 
  mutate(pop_pair = 'F2') %>% 
  dplyr::select(-resid_10)

myv_TGP_12_12_nemeny_df = myv_TGP_12_12_residuals %>% 
  mutate(pop_pair = 'TGP') %>% 
  dplyr::select(-resid_10)

myv_WGP_12_12_nemeny_df = myv_WGP_12_12_residuals %>% 
  mutate(pop_pair = 'WGP') %>% 
  dplyr::select(-resid_10)


myv_comp_12_12_nemeny_df = bind_rows(myv_f2_12_12_nemeny_df, 
                                     myv_TGP_12_12_nemeny_df, 
                                     myv_WGP_12_12_nemeny_df)


myv_comp_12_12_nemeny_test = frdAllPairsNemenyiTest(y = myv_comp_12_12_nemeny_df$residual, 
                                              groups = myv_comp_12_12_nemeny_df$trait, 
                                              blocks = myv_comp_12_12_nemeny_df$pop_pair, 
                                              p.adjust.method = 'hochberg')

myv_comp_12_12_nemeny_test$p.value %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(trait = rowname) %>% 
  filter(if_any(c(2:33), function(x) x<0.05)) %>% 
  View()
  # write_csv('F2_1212_residual_NemenyTest_results.csv')



# MYV 12_12 COV matrix comparisons ----------------------------------------

krz_F2_TGP = KrzCor()


# Covariance matrices
F2_cov  <- cov(F2_traits)
TGP_cov <- cov(TGP_traits)
WGP_cov <- cov(WGP_traits)

# Compare to the overall F2 phenotype
krz_F2_TGP <- KrzCor(F2_cov, TGP_cov)
krz_F2_WGP <- KrzCor(F2_cov, WGP_cov)

# Difference in similarity
delta_krz <- krz_F2_TGP - krz_F2_WGP

delta_krz