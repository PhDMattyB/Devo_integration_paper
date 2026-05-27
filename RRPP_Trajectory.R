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


