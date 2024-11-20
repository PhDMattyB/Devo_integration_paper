##############################
## Family effects - integration paper
##
## Matt Brachmann (PhDMattyB)
##
## 20.11.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(tidyverse)

F1_lmk_dist = read_csv('F1_Plasticity_Corrected.csv')

F1_lmk_dist %>% 
  select(individualID)

