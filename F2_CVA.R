##############################
##  F2 integration CVA analysis
##
## Matt Brachmann (PhDMattyB)
##
## 13.05.2024
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')


library(geomorph)
library(RRPP)
library(MASS)
library(tidyverse)

F2_craniofacial = readland.tps('F2_Craniofacial_LM.TPS',
                               specID = 'imageID')

identifiers = read_csv('F2_metadata.csv') %>% 
  unite('Ecotype_Pair_Full_Temp', 
        Ecotype_pair, 
        Full_temp, 
        sep = '_', 
        remove = F)

identifiers = read_csv('F2_metadata.csv') %>% 
  unite('lake_morph_Pair_Full_Temp', 
        Lake_morph, 
        Full_temp, 
        sep = '_', 
        remove = F)

## perform gpa on craniofacial data
F2_craniofacial_gpa = gpagen(F2_craniofacial,
                             print.progress = F)


