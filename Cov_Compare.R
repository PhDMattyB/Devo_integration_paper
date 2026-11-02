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

library(geomorph)



# ASHN --------------------------------------------------------------------



ashn_f2_cov_comp = ashn_f2_eigen_residuals %>% 
  dplyr::select(eigen_residual)

ashn_wgp_cov_comp = ashn_wgp_eigen_residuals %>% 
  dplyr::select(eigen_residual)

ashn_tgp_cov_comp = ashn_tgp_eigen_residuals %>% 
  dplyr::select(eigen_residual)

ashn_f2_wgp_cor = cor(ashn_f2_cov_comp, 
                     ashn_wgp_cov_comp, 
                     method = 'pearson')

ashn_f2_tgp_cor = cor(ashn_f2_cov_comp, 
                     ashn_tgp_cov_comp, 
                     method = 'pearson')

ashn_tgp_wgp_cor = cor(ashn_tgp_cov_comp, 
                     ashn_wgp_cov_comp, 
                     method = 'pearson')

ashn_delta_cor = ashn_f2_tgp_cor - ashn_f2_wgp_cor


cor.test(ashn_f2_cov_comp$eigen_residual, 
         ashn_wgp_cov_comp$eigen_residual)

cor.test(ashn_f2_cov_comp$eigen_residual, 
         ashn_tgp_cov_comp$eigen_residual)

cor.test(ashn_tgp_cov_comp$eigen_residual, 
         ashn_wgp_cov_comp$eigen_residual)


# MYV ---------------------------------------------------------------------



myv_f2_cov_comp = myv_f2_eigen_residuals %>% 
  dplyr::select(eigen_residual)

myv_wgp_cov_comp = myv_wgp_eigen_residuals %>% 
  dplyr::select(eigen_residual)

myv_tgp_cov_comp = myv_tgp_eigen_residuals %>% 
  dplyr::select(eigen_residual)

myv_f2_wgp_cor = cor(myv_f2_cov_comp, 
                      myv_wgp_cov_comp, 
                      method = 'pearson')

myv_f2_tgp_cor = cor(myv_f2_cov_comp, 
                      myv_tgp_cov_comp, 
                      method = 'pearson')

myv_tgp_wgp_cor = cor(myv_tgp_cov_comp, 
                       myv_wgp_cov_comp, 
                       method = 'pearson')

myv_delta_cor = myv_f2_tgp_cor - myv_f2_wgp_cor


cor.test(myv_f2_cov_comp$eigen_residual, 
         myv_wgp_cov_comp$eigen_residual)

cor.test(myv_f2_cov_comp$eigen_residual, 
         myv_tgp_cov_comp$eigen_residual)

cor.test(myv_tgp_cov_comp$eigen_residual, 
         myv_wgp_cov_comp$eigen_residual)



# SKR ---------------------------------------------------------------------




skr_f2_cov_comp = skr_f2_eigen_residuals %>% 
  dplyr::select(eigen_residual)

skr_wgp_cov_comp = skr_wgp_eigen_residuals %>% 
  dplyr::select(eigen_residual)

skr_tgp_cov_comp = skr_tgp_eigen_residuals %>% 
  dplyr::select(eigen_residual)

skr_f2_wgp_cor = cor(skr_f2_cov_comp, 
                      skr_wgp_cov_comp, 
                      method = 'pearson')

skr_f2_tgp_cor = cor(skr_f2_cov_comp, 
                      skr_tgp_cov_comp, 
                      method = 'pearson')

skr_tgp_wgp_cor = cor(skr_tgp_cov_comp, 
                       skr_wgp_cov_comp, 
                       method = 'pearson')

skr_delta_cor = skr_f2_tgp_cor - skr_f2_wgp_cor


cor.test(skr_f2_cov_comp$eigen_residual, 
         skr_wgp_cov_comp$eigen_residual)

cor.test(skr_f2_cov_comp$eigen_residual, 
         skr_tgp_cov_comp$eigen_residual)

cor.test(skr_tgp_cov_comp$eigen_residual, 
         skr_wgp_cov_comp$eigen_residual)


# GTS/GAR -----------------------------------------------------------------



gts_gar_f2_cov_comp = gts_gar_f2_eigen_residuals %>% 
  dplyr::select(eigen_residual)

gts_gar_wgp_cov_comp = gts_gar_wgp_eigen_residuals %>% 
  dplyr::select(eigen_residual)

gts_gar_tgp_cov_comp = gts_gar_tgp_eigen_residuals %>% 
  dplyr::select(eigen_residual)

gts_gar_f2_wgp_cor = cor(gts_gar_f2_cov_comp, 
                      gts_gar_wgp_cov_comp, 
                      method = 'pearson')

gts_gar_f2_tgp_cor = cor(gts_gar_f2_cov_comp, 
                      gts_gar_tgp_cov_comp, 
                      method = 'pearson')

gts_gar_tgp_wgp_cor = cor(gts_gar_tgp_cov_comp, 
                       gts_gar_wgp_cov_comp, 
                       method = 'pearson')

gts_gar_delta_cor = gts_gar_f2_tgp_cor - gts_gar_f2_wgp_cor


cor.test(gts_gar_f2_cov_comp$eigen_residual, 
         gts_gar_wgp_cov_comp$eigen_residual)

cor.test(gts_gar_f2_cov_comp$eigen_residual, 
         gts_gar_tgp_cov_comp$eigen_residual)

cor.test(gts_gar_tgp_cov_comp$eigen_residual, 
         gts_gar_wgp_cov_comp$eigen_residual)

