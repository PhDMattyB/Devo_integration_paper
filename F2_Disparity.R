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

# Metadata ----------------------------------------------------------------
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

# Shape trait data --------------------------------------------------------


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


## operculum shape
F2_whole_body = readland.tps('F2_No_GT.TPS', 
                               specID = 'imageID')
F2_whole_body = gpagen(F2_whole_body)
F2_whole_body = geomorph.data.frame(coords = two.d.array(F2_whole_body$coords), 
                                   Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                   parent_temp = identifiers$Parent_temp, 
                                   offspring_temp = identifiers$Offspring_temp,
                                   grand_temp = identifiers$Grand_temp,
                                   morph = identifiers$Morph, 
                                   population = identifiers$Lake,
                                   lake_morph = identifiers$Lake_morph,
                                   lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)

# whole body ecotype Disparity analysis ------------------------------------------------------

F2_whole_body

# lake_ecotype_disparity = morphol.disparity(coords ~ 1,
#                                            groups = ~lake_morph,
#                                            data = F2_cranio_geo_df,
#                                            iter = 999)

F2_whole_body_disp = morphol.disparity(coords ~ 1,
                                  groups = ~lake_morph_full,
                                  data = F2_whole_body,
                                  iter = 999)

Whole_body_pval = F2_whole_body_disp$PV.dist.Pval
Whole_body_disp = F2_whole_body_disp$PV.dist
disp_proc_var = full_mod_disp$Procrustes.var


Whole_body_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('Full_LM_set_disparity_pvals.csv')

Whole_body_disp %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('Full_LM_set_disparity_estimates.csv')


F2_whole_body_disp = morphol.disparity(coords ~ 1,
                                       groups = ~lake_morph,
                                       data = F2_whole_body,
                                       iter = 999)

Whole_body_pval = F2_whole_body_disp$PV.dist.Pval
Whole_body_disp = F2_whole_body_disp$PV.dist
disp_proc_var = F2_whole_body_disp$Procrustes.var

# F1 corrected landmarks --------------------------------------------------


F1_effects = readland.tps('F1_Corrected_landmarks.tps', 
                             specID = 'imageID')
F1_effects = gpagen(F1_effects)
F1_effects = geomorph.data.frame(coords = two.d.array(F1_effects$coords), 
                                    Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                    parent_temp = identifiers$Parent_temp, 
                                    offspring_temp = identifiers$Offspring_temp,
                                    grand_temp = identifiers$Grand_temp,
                                    morph = identifiers$Morph, 
                                    population = identifiers$Lake,
                                    lake_morph = identifiers$Lake_morph,
                                    lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)





F1_effects_mod_disp = morphol.disparity(coords ~ 1,
                                  groups = ~lake_morph_full,
                                  data = F1_effects,
                                  iter = 999)
F1_disp_pval = F1_effects_mod_disp$PV.dist.Pval
F1_disp_vals = F1_effects_mod_disp$PV.dist
F1_disp_proc_var = F1_effects_mod_disp$Procrustes.var

F1_disp_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('F1_Effects_whole_body.csv')

F1_disp_vals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('F1_Effects_whole_body.csv')

F1_effects_morph_mod_disp = morphol.disparity(coords ~ 1,
                                        groups = ~lake_morph,
                                        data = F1_effects,
                                        iter = 999)

F1_disp_pval = F1_effects_morph_mod_disp$PV.dist.Pval
F1_disp_vals = F1_effects_morph_mod_disp$PV.dist
F1_disp_proc_var = F1_effects_morph_mod_disp$Procrustes.var

F1_disp_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('F1_morph_Effects_whole_body_pvals.csv')

F1_disp_vals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('F1_morph_Effects_whole_body_estimates.csv')

# F2 effects --------------------------------------------------------------
F2_effects = readland.tps('F2_Corrected_landmarks.tps', 
                          specID = 'imageID')
F2_effects = gpagen(F2_effects)
F2_effects = geomorph.data.frame(coords = two.d.array(F2_effects$coords), 
                                 Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                 parent_temp = identifiers$Parent_temp, 
                                 offspring_temp = identifiers$Offspring_temp,
                                 grand_temp = identifiers$Grand_temp,
                                 morph = identifiers$Morph, 
                                 population = identifiers$Lake,
                                 lake_morph = identifiers$Lake_morph,
                                 lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)


F2_effects_mod_disp = morphol.disparity(coords ~ 1,
                                        groups = ~lake_morph_full,
                                        data = F2_effects,
                                        iter = 999)

F2_disp_pval = F2_effects_mod_disp$PV.dist.Pval
F2_disp_vals = F2_effects_mod_disp$PV.dist
F2_disp_proc_var = F2_effects_mod_disp$Procrustes.var

F2_disp_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('F2_Effects_whole_body.csv')

F2_disp_vals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('F2_Effects_whole_body.csv')

F2_effects_morph_mod_disp = morphol.disparity(coords ~ 1,
                                        groups = ~lake_morph,
                                        data = F2_effects,
                                        iter = 999)

F2_disp_pval = F2_effects_morph_mod_disp$PV.dist.Pval
F2_disp_vals = F2_effects_morph_mod_disp$PV.dist
F2_disp_proc_var = F2_effects_morph_mod_disp$Procrustes.var

F2_disp_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('F2_morph_Effects_whole_body_pvals.csv')

F2_disp_vals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('F2_morph_Effects_whole_body_estimates.csv')


# craniofacial ecotype Disparity analysis ------------------------------------------------------

F2_cranio_geo_df

lake_ecotype_disparity = morphol.disparity(coords ~ 1,
                  groups = ~lake_morph,
                  data = F2_cranio_geo_df,
                  iter = 999)

full_mod_disp = morphol.disparity(coords ~ 1,
                                           groups = ~lake_morph_full,
                                           data = F2_cranio_geo_df,
                                           iter = 999)

disp_pval = full_mod_disp$PV.dist.Pval
disp_vals = full_mod_disp$PV.dist
disp_proc_var = full_mod_disp$Procrustes.var


disp_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('Craniofacial_disparity_pvals.csv')

disp_vals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('Craniofacial_disparity_estimates.csv')

# body shape disparity ----------------------------------------------------
# lake_ecotype_disparity = morphol.disparity(coords ~ 1,
#                                            groups = ~lake_morph,
#                                            data = F2_cranio_geo_df,
#                                            iter = 999)

full_mod_disp_body = morphol.disparity(coords ~ 1,
                                  groups = ~lake_morph_full,
                                  data = F2_body_geo_df,
                                  iter = 999)

body_disp_pval = full_mod_disp_body$PV.dist.Pval
body_disp_vals = full_mod_disp_body$PV.dist
body_disp_proc_var = full_mod_disp_body$Procrustes.var

body_disp_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('Body_disparity_pvals.csv')

body_disp_vals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('Body_disparity_estimates.csv')

# 4bar shape disparity ----------------------------------------------------
# lake_ecotype_disparity = morphol.disparity(coords ~ 1,
#                                            groups = ~lake_morph,
#                                            data = F2_cranio_geo_df,
#                                            iter = 999)

full_mod_disp_4bar = morphol.disparity(coords ~ 1,
                                       groups = ~lake_morph_full,
                                       data = F2_4bar_geo_df,
                                       iter = 999)

fbar_disp_pval = full_mod_disp_4bar$PV.dist.Pval
fbar_disp_vals = full_mod_disp_4bar$PV.dist
fbar_disp_proc_var = full_mod_disp_4bar$Procrustes.var


fbar_disp_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('fbar_disparity_pvals.csv')

fbar_disp_vals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('fbar_disparity_estimates.csv')

# multishape 345 shape disparity ----------------------------------------------------
# lake_ecotype_disparity = morphol.disparity(coords ~ 1,
#                                            groups = ~lake_morph,
#                                            data = F2_cranio_geo_df,
#                                            iter = 999)

full_mod_disp_eye = morphol.disparity(coords ~ 1,
                                       groups = ~lake_morph_full,
                                       data = F2_multi_345,
                                       iter = 999)

eye_disp_pval = full_mod_disp_eye$PV.dist.Pval
eye_disp_vals = full_mod_disp_eye$PV.dist
eye_disp_proc_var = full_mod_disp_eye$Procrustes.var

eye_disp_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('Eye_disparity_pvals.csv')

eye_disp_vals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('Eye_disparity_estimates.csv')


# multishape 789 shape disparity ----------------------------------------------------
# lake_ecotype_disparity = morphol.disparity(coords ~ 1,
#                                            groups = ~lake_morph,
#                                            data = F2_cranio_geo_df,
#                                            iter = 999)

full_mod_disp_operculum = morphol.disparity(coords ~ 1,
                                      groups = ~lake_morph_full,
                                      data = F2_multi_789,
                                      iter = 999)

operculum_disp_pval = full_mod_disp_operculum$PV.dist.Pval
operculum_disp_vals = full_mod_disp_operculum$PV.dist
operculum_disp_proc_var = full_mod_disp_operculum$Procrustes.var

operculum_disp_pval %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('Operculum_disparity_pvals.csv')

operculum_disp_vals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  write_csv('Operculum_disparity_estimates.csv')
