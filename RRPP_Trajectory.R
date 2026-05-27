##############################
## RRPP trajectory analysis 
##
## Matt Brachmann (PhDMattyB)
##
## 27.05.2026
##
##############################


setwd("~/Parsons_postdoc/Stickleback_Morphometric_data/Updated Landmarks/")

library(tidyverse)
library(paran)
library(RRPP)


F2_pca_data = read_csv('F2_orig_PCA_sig_axes.csv')

WGP_pca_data = read_csv('WGP_PCA_sig_axes.csv')

TGP_pca_data = read_csv('TGP_PCA_sig_axes.csv')

wild_pca_data = read_csv('wild_PCA_sig_axes.csv')

F2_PC_axes = F2_pca_data %>% dplyr::select(PC1:PC7) %>% 
  as.matrix()

F2_RRPP_mod <- lm.rrpp(F2_PC_axes ~ Morph*Offspring_temp*Parent_temp , 
                       SS.type = "I", 
               data = F2_pca_data, 
               print.progress = FALSE, 
               iter = 999,
               verbose = TRUE) 
summary(F2_RRPP_mod, 
        formula = FALSE)
anova(F2_RRPP_mod) 

coef(F2_RRPP_mod, test = TRUE)
