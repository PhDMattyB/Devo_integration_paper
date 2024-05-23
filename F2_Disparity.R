#################################################
## Disparity analyses - ecotypes and environments
##
## Matt Brachmann (PhDMattyB)
##
##  23.05.2024
## 
################################################


setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')


library(geomorph)
library(RRPP)
# library(MASS)
# library(ppcor)
# library(igraph)
library(tidyverse)

identifiers = read_csv('F2_metadata.csv') %>% 
  rename(individualID = Names) %>% 
  unite('lake_morph_Pair_Full_Temp', 
        Lake_morph, 
        Full_temp, 
        sep = '_', 
        remove = F) %>% 
  unite('Ecotype_Pair_Full_Temp', 
        Ecotype_pair, 
        Full_temp, 
        sep = '_', 
        remove = F) %>% 
  mutate(across(c('Lake_morph',
                  'Offspring_temp',
                  'Parent_temp',
                  'Grand_temp'),
                factor)) %>% 
  arrange(individualID)

## Craniofacial traits
F2_craniofacial = readland.tps('F2_Craniofacial_LM.TPS',
                               specID = 'imageID')
F2_craniofacial_gpa = gpagen(F2_craniofacial,
                             print.progress = F)
F2_cranio_geo_df = geomorph.data.frame(coords = two.d.array(F2_craniofacial_gpa$coords), 
                                       Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                       parent_temp = identifiers$Parent_temp, 
                                       offspring_temp = identifiers$Offspring_temp,
                                       grand_temp = identifiers$Grand_temp,
                                       morph = identifiers$Morph, 
                                       population = identifiers$Lake,
                                       lake_morph = identifiers$Lake_morph,
                                       lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)

## body shape data
F2_body = readland.tps('F2_Body_LM.TPS', 
                       specID = 'imageID')
F2_body_gpa = gpagen(F2_body,
                     print.progress = F)
F2_body_geo_df = geomorph.data.frame(coords = two.d.array(F2_body_gpa$coords), 
                                     Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                     parent_temp = identifiers$Parent_temp, 
                                     offspring_temp = identifiers$Offspring_temp,
                                     grand_temp = identifiers$Grand_temp,
                                     morph = identifiers$Morph, 
                                     population = identifiers$Lake,
                                     lake_morph = identifiers$Lake_morph,
                                     lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)

## 4 bar linkage
F2_4bar = readland.tps('F2_4bar_linkage.TPS', 
                       specID = 'imageID')
F2_4bar_gpa = gpagen(F2_4bar,
                     print.progress = F)
F2_4bar_geo_df = geomorph.data.frame(coords = two.d.array(F2_4bar_gpa$coords), 
                                     Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                     parent_temp = identifiers$Parent_temp, 
                                     offspring_temp = identifiers$Offspring_temp,
                                     grand_temp = identifiers$Grand_temp,
                                     morph = identifiers$Morph, 
                                     population = identifiers$Lake,
                                     lake_morph = identifiers$Lake_morph,
                                     lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)

## eye shape 
multi_shape_345 = readland.tps('multi_shape_345.TPS', 
                               specID = 'imageID')
multi_shape_345 = gpagen(multi_shape_345)
F2_multi_345 = geomorph.data.frame(coords = two.d.array(multi_shape_345$coords), 
                                   Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                   parent_temp = identifiers$Parent_temp, 
                                   offspring_temp = identifiers$Offspring_temp,
                                   grand_temp = identifiers$Grand_temp,
                                   morph = identifiers$Morph, 
                                   population = identifiers$Lake,
                                   lake_morph = identifiers$Lake_morph,
                                   lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)


## operculum shape
multi_shape_789 = readland.tps('multi_shape_789.TPS', 
                               specID = 'imageID')
multi_shape_789 = gpagen(multi_shape_789)
F2_multi_789 = geomorph.data.frame(coords = two.d.array(multi_shape_789$coords), 
                                   Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                   parent_temp = identifiers$Parent_temp, 
                                   offspring_temp = identifiers$Offspring_temp,
                                   grand_temp = identifiers$Grand_temp,
                                   morph = identifiers$Morph, 
                                   population = identifiers$Lake,
                                   lake_morph = identifiers$Lake_morph,
                                   lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)


# ashn = landmarks %>%
#   filter(POP_only == 'ASHN')
# 
# ashn_pheno = as.matrix(ashn[which(names(ashn) == 'LM1X'):
#                               which(names(ashn) == 'LM22Y')])   
# 
# rownames(ashn_pheno) = ashn$ID
# 
# # dim(phenotypes)
# ashn_array = arrayspecs(ashn_pheno, 
#                         p = 22, 
#                         k = 2)
# 
# ashn_gpa = gpagen(ashn_array)
# 
# ashn_df = geomorph.data.frame(ashn_gpa, 
#                               # species = plethodon$species, 
#                               species = ashn$POP)
# 
# morphol.disparity(coords ~ 1,
#                   groups = ~species,
#                   data = ashn_df, 
#                   iter = 999)
