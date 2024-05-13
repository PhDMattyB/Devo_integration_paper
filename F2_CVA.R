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

# identifiers = read_csv('F2_metadata.csv') %>% 
#   unite('Ecotype_Pair_Full_Temp', 
#         Ecotype_pair, 
#         Full_temp, 
#         sep = '_', 
#         remove = F)

identifiers = read_csv('F2_metadata.csv') %>% 
  unite('lake_morph_Pair_Full_Temp', 
        Lake_morph, 
        Full_temp, 
        sep = '_', 
        remove = F) %>% 
  unite('Ecotype_Pair_Full_Temp', 
        Ecotype_pair, 
        Full_temp, 
        sep = '_', 
        remove = F)

## perform gpa on craniofacial data
F2_craniofacial_gpa = gpagen(F2_craniofacial,
                             print.progress = F)



F2_geo_df = geomorph.data.frame(coords = two.d.array(F2_craniofacial_gpa$coords), 
                                Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                parent_temp = identifiers$Parent_temp, 
                                offspring_temp = identifiers$Offspring_temp,
                                grand_temp = identifiers$Grand_temp,
                                morph = identifiers$Morph, 
                                population = identifiers$Lake, 
                                lake_morph = identifiers$lake_morph_Pair_Full_Temp)

F2_off_temp = lm.rrpp(coords ~ offspring_temp, 
                     data = F2_geo_df, 
                     turbo = F, 
                     verbose = T)
