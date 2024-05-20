#########################################
## Pheno Integration F2 CVA + univariate
## 
## Matt Brachmann (PhDMattyB)
##
## 
## 20.05.2023
########################################


setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')


library(geomorph)
library(RRPP)
library(MASS)
library(ppcor)
library(igraph)
library(tidyverse)


# F2 Data -----------------------------------------------------------------

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


F2_whole_body = readland.tps('F2_No_GT.TPS', 
                      specID = 'imageID')
F2_whole_body_gpa = gpagen(F2_whole_body,
                             print.progress = F)
F2_WB_geo_df = geomorph.data.frame(coords = two.d.array(F2_whole_body_gpa$coords), 
                                       Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                       parent_temp = identifiers$Parent_temp, 
                                       offspring_temp = identifiers$Offspring_temp,
                                       grand_temp = identifiers$Grand_temp,
                                       morph = identifiers$Morph, 
                                       population = identifiers$Lake,
                                       lake_morph = identifiers$Lake_morph,
                                       lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)


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

multi_shape_1011 = readland.tps('multi_shape_1011.TPS', 
                               specID = 'imageID')
multi_shape_1011 = gpagen(multi_shape_1011)
F2_multi_1011 = geomorph.data.frame(coords = two.d.array(multi_shape_1011$coords), 
                                   Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                   parent_temp = identifiers$Parent_temp, 
                                   offspring_temp = identifiers$Offspring_temp,
                                   grand_temp = identifiers$Grand_temp,
                                   morph = identifiers$Morph, 
                                   population = identifiers$Lake,
                                   lake_morph = identifiers$Lake_morph,
                                   lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)

multi_shape_12 = readland.tps('multi_shape_12.TPS', 
                               specID = 'imageID')
multi_shape_12 = gpagen(multi_shape_12)
F2_multi_12 = geomorph.data.frame(coords = two.d.array(multi_shape_12$coords), 
                                   Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                   parent_temp = identifiers$Parent_temp, 
                                   offspring_temp = identifiers$Offspring_temp,
                                   grand_temp = identifiers$Grand_temp,
                                   morph = identifiers$Morph, 
                                   population = identifiers$Lake,
                                   lake_morph = identifiers$Lake_morph,
                                   lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)


# CVA multivariate traits -------------------------------------------------


# eye landmarks 345 -------------------------------------------------------

multi_shape345_mod = procD.lm(coords ~ offspring_temp * lake_morph, 
                       data = F2_multi_345)
summary(multi_shape345_mod)
prep.lda(multi_shape345_mod, 
         inherent.groups = TRUE) # see groups available

lda.args = prep.lda(multi_shape345_mod) 
multi_shape345_CVA = do.call(lda, lda.args)

## CVS scores for each individual
multi_shape345_cva_scores = predict(multi_shape345_CVA)

multi_shape345_cva_scores = multi_shape345_cva_scores$x %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(individualID = rowname) %>% 
  arrange(individualID)

multi_shape345_cva = bind_cols(multi_shape345_cva_scores, 
                            identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)

multi_shape345_cva %>% 
  write_csv('F2_eye_shape_cva_per_pop.csv')


# Inter-LM distances - get univariate traits ------------------------------------------------------


lmks = data.frame(jaw_length = c(1, 2), 
                  fbar_23_24 = c(23, 24), 
                  fbar_8_24 = c(8, 24), 
                  fbar_8_27 = c(8, 27), 
                  fbar_23_27 = c(23, 27), 
                  fbar_25_26 = c(25, 26), 
                  body_width = c(12, 21), 
                  caudal1_14_18 = c(14, 18), 
                  caudal2_15_17 = c(15, 17), 
                  body_length = c(1, 16),
                  row.names = c('start', 
                                'end'))

A = F2_whole_body_gpa$coords
F2_univariate_traits = interlmkdist(A, 
                                  lmks)

# arrayspecs(F2_univariate_traits, 
#            4, 
#            3)



