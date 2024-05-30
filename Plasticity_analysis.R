##############################
##  Plasticity of integrated traits
##
## Matt Brachmann (PhDMattyB)
##
##  29.05.2024
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')


library(geomorph)
library(RRPP)
library(MASS)
library(ppcor)
library(igraph)
library(tidyverse)
library(reshape2)
library(candisc)

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
                factor)) 
# %>% 
#   arrange(individualID)


# Body shape data ---------------------------------------------------------

F2_tps = readland.tps('F2_No_GT.TPS', 
                      specID = 'imageID')

## superimposition on the entire dataset
F2_gpa = gpagen(F2_tps, 
                print.progress = F)

allometry_model1 = procD.lm(F2_gpa$coords ~ log(F2_gpa$Csize), 
                            iter = 999, 
                            RRPP = T)
summary(allometry_model1)


F2_shape_resid = arrayspecs(allometry_model1$residuals, 
                            p = dim(F2_gpa$coords)[1], 
                            k = dim(F2_gpa$coords)[2])
F2_allometry_adj_shape = F2_shape_resid + array(F2_gpa$consensus, 
                                                dim(F2_shape_resid))

mean_shape = mshape(F2_gpa$coords)
matrix_mean_shape = as.matrix(mean_shape)
mean_shape_array = array(matrix_mean_shape, 
                         dim = c(27, 2, 1))
# univariate trait data ---------------------------------------------------


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
                  head_depth = c(1, 22), 
                  jaw_2_6 = c(2, 6), 
                  lm_6_12 = c(6, 12), 
                  lm_12_13 = c(12, 13), 
                  lm_13_14 = c(13, 14), 
                  lm_14_15 = c(14, 15), 
                  lm_6_21 = c(6, 21), 
                  lm_20_21 = c(20, 21), 
                  lm_21_13 = c(21, 13), 
                  lm_20_13 = c(20, 13), 
                  lm_12_19 = c(12, 19), 
                  lm_13_19 = c(13, 19), 
                  lm_19_18 = c(19, 18), 
                  lm_18_17 = c(18, 17), 
                  lm_1_23 = c(1, 23), 
                  lm_23_2 = c(23, 2),
                  row.names = c('start', 
                                'end'))

A = F2_gpa$coords
# A = F2_whole_body_gpa$coords
F2_univariate_traits = interlmkdist(A, 
                                    lmks)

# arrayspecs(F2_univariate_traits, 
#            4, 
#            3)

F2_univariate_traits = F2_univariate_traits %>% 
  as.data.frame() %>% 
  rownames_to_column() 
# %>% 
#   arrange(rowname)


F2_univariate_traits = bind_cols(F2_univariate_traits, 
                                 identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)



# Plasticity  shape --------------------------------------------------------------

F2_temp_mod = procD.lm(F2_gpa$coords ~ identifiers$Offspring_temp, 
                       iter = 999)

## All individuals have the same fitted values.
## pull individual from offspring temp of 12 degrees
F2_temp_fitted = F2_temp_mod$GM$fitted[,,1]
F2_temp_matrix_12deg = as.matrix(F2_temp_fitted)
F2_temp_12deg_array = array(F2_temp_matrix_12deg, dim = c(27, 2, 1))

F2_temp_fitted_18deg = F2_temp_mod$GM$fitted[,,31]
F2_temp_matrix_18deg = as.matrix(F2_temp_fitted_18deg)
F2_temp_18deg_array = array(F2_temp_matrix_18deg, dim = c(27,2, 1))

# identifiers %>% 
#   filter(Offspring_temp == '18') %>% 
#   View()

F2_12deg_range = c(1:30, 61:91, 182:211, 244:273, 304:333, 364:382, 
                   413:442, 474:503, 534:563, 594:623, 655:683,
                   714:743, 774:803, 834:857, 871:900)
F2_18deg_range = c(31:60, 92:121, 152:181, 212:243, 274:303, 
                   334:363, 383:412, 443:473, 504:533, 564:593, 
                   624:654, 684:713, 744:773, 804:833, 858:870, 
                   901:931)

F2_array = array(0, dim = c(27, 2, 931))
for(i in F2_12deg_range){
  F2_array[,,i] = F2_gpa$coords[,,i] - F2_temp_12deg_array[,,1]
}

for(i in F2_18deg_range){
  F2_array[,,i] = F2_gpa$coords[,,i] - F2_temp_18deg_array[,,1]
}

## This is the array to use too pull out the linear traits due
## to plasticity
F2_array_consensus = array(0, dim = c(27, 2, 931))
for(i in 1:931){
  F2_array_consensus[,,i] = F2_array[,,i] + mean_shape_array[,,1]
}


plasticity_gpa = gpagen(F2_array_consensus)
# test_lm = geomorph.data.frame(plasticity_gpa)


F2_plasticity_data = geomorph.data.frame(plasticity_gpa, 
                              Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                              parent_temp = identifiers$Parent_temp, 
                              offspring_temp = identifiers$Offspring_temp,
                              grand_temp = identifiers$Grand_temp,
                              morph = identifiers$Morph, 
                              population = identifiers$Lake,
                              lake_morph = identifiers$Lake_morph,
                              lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)


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

# WC_ecotype_residuals = arrayspecs(WC_ecotype_residuals, 
#                                   27, 
#                                   2)

C = F2_plasticity_data$coords
# A = F2_whole_body_gpa$coords
F2_off_plasticity_traits = interlmkdist(C, 
                                    lmks)


F2_off_plasticity_traits = F2_off_plasticity_traits %>% 
  as.data.frame() %>% 
  rownames_to_column() 
# %>% 
#   arrange(rowname)


F2_off_plasticity_traits = bind_cols(F2_off_plasticity_traits, 
                                 identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)


off_plasticity_traits = F2_off_plasticity_traits %>% 
  as_tibble() %>% 
  group_by(lake_morph_Pair_Full_Temp) %>% 
  select(jaw_length:body_length)

vars_keep = names(off_plasticity_traits)[c(2,3,4,5,6,7,8,9,10,11)]
off_plasticity_trait_cor = off_plasticity_traits %>% 
  ungroup() %>% 
  split(.$lake_morph_Pair_Full_Temp) %>% 
  # ungroup() %>% 
  map(select, vars_keep) %>% 
  map(cor)

off_plasticity_graph = off_plasticity_trait_cor %>% 
  reshape2::melt() %>% 
  rename(lake_morph_full = L1)

off_plasticity_trait_cor_graph = ggplot(off_plasticity_graph, 
                                    aes(x = Var1, 
                                        y = Var2, 
                                        fill = value))+
  geom_tile()+
  facet_wrap(~lake_morph_full, 
             ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1))

ggsave('Univariate_offtemp_plasticity_trait_integration.tiff', 
       plot = off_plasticity_trait_cor_graph, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 40)


# parental temp effects ---------------------------------------------------

F2_parent_temp_mod = procD.lm(F2_gpa$coords ~ identifiers$Parent_temp, 
                       iter = 999)

## All individuals have the same fitted values.
## pull individual from offspring temp of 12 degrees
F2_parent_temp_fitted = F2_parent_temp_mod$GM$fitted[,,1]
F2_parent_temp_matrix_12deg = as.matrix(F2_parent_temp_fitted)
F2_parent_temp_12deg_array = array(F2_parent_temp_matrix_12deg, dim = c(27, 2, 1))

F2_parent_temp_fitted_18deg = F2_parent_temp_mod$GM$fitted[,,61]
F2_parent_temp_matrix_18deg = as.matrix(F2_parent_temp_fitted_18deg)
F2_parent_temp_18deg_array = array(F2_parent_temp_matrix_18deg, dim = c(27,2, 1))

identifiers %>%
  filter(Parent_temp == '12') %>%
  View()

F2_parent_12deg_range = c(1:60, 122:181, 244:303, 364:442, 474:533, 594:654, 
                   714:773, 834:870)
F2_parent_18deg_range = c(61:121, 182:243, 304:363, 443:473, 534:593, 
                          655:713, 774:833, 871:931)

F2_parent_array = array(0, dim = c(27, 2, 931))
for(i in F2_parent_12deg_range){
  F2_parent_array[,,i] = F2_gpa$coords[,,i] - F2_parent_temp_12deg_array[,,1]
}

for(i in F2_parent_18deg_range){
  F2_parent_array[,,i] = F2_gpa$coords[,,i] - F2_parent_temp_18deg_array[,,1]
}

## This is the array to use too pull out the linear traits due
## to plasticity
F2_parent_array_consensus = array(0, dim = c(27, 2, 931))
for(i in 1:931){
  F2_parent_array_consensus[,,i] = F2_parent_array[,,i] + mean_shape_array[,,1]
}


parent_effect_gpa = gpagen(F2_parent_array_consensus)
# test_lm = geomorph.data.frame(plasticity_gpa)


F2_parent_effect_data = geomorph.data.frame(parent_effect_gpa, 
                                         Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                         parent_temp = identifiers$Parent_temp, 
                                         offspring_temp = identifiers$Offspring_temp,
                                         grand_temp = identifiers$Grand_temp,
                                         morph = identifiers$Morph, 
                                         population = identifiers$Lake,
                                         lake_morph = identifiers$Lake_morph,
                                         lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)


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

# WC_ecotype_residuals = arrayspecs(WC_ecotype_residuals, 
#                                   27, 
#                                   2)

D = F2_parent_effect_data$coords
# A = F2_whole_body_gpa$coords
F2_parent_effect_traits = interlmkdist(D, 
                                        lmks)


F2_parent_plasticity_traits = F2_parent_effect_traits %>% 
  as.data.frame() %>% 
  rownames_to_column() 
# %>% 
#   arrange(rowname)


F2_parent_plasticity_traits = bind_cols(F2_parent_plasticity_traits, 
                                     identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)


parent_plasticity_traits = F2_parent_plasticity_traits %>% 
  as_tibble() %>% 
  group_by(lake_morph_Pair_Full_Temp) %>% 
  select(jaw_length:body_length)

vars_keep = names(parent_plasticity_traits)[c(2,3,4,5,6,7,8,9,10,11)]
parent_plasticity_trait_cor = parent_plasticity_traits %>% 
  ungroup() %>% 
  split(.$lake_morph_Pair_Full_Temp) %>% 
  # ungroup() %>% 
  map(select, vars_keep) %>% 
  map(cor)

parent_plasticity_graph = parent_plasticity_trait_cor %>% 
  reshape2::melt() %>% 
  rename(lake_morph_full = L1)

parent_plasticity_trait_cor_graph = ggplot(parent_plasticity_graph, 
                                        aes(x = Var1, 
                                            y = Var2, 
                                            fill = value))+
  geom_tile()+
  facet_wrap(~lake_morph_full, 
             ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1))

ggsave('Univariate_parent_temp_plasticity_trait_integration.tiff', 
       plot = parent_plasticity_trait_cor_graph, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 40)


# ecotype effects ---------------------------------------------------------

F2_ecotype_temp_mod = procD.lm(F2_gpa$coords ~ identifiers$Morph, 
                              iter = 999)

## All individuals have the same fitted values.
## pull individual from offspring temp of 12 degrees
F2_ecotype_temp_fitted = F2_ecotype_temp_mod$GM$fitted[,,1]
F2_ecotype_temp_matrix_12deg = as.matrix(F2_ecotype_temp_fitted)
F2_ecotype_temp_12deg_array = array(F2_ecotype_temp_matrix_12deg, dim = c(27, 2, 1))

F2_ecotype_temp_fitted_18deg = F2_ecotype_temp_mod$GM$fitted[,,243]
F2_ecotype_temp_matrix_18deg = as.matrix(F2_ecotype_temp_fitted_18deg)
F2_ecotype_temp_18deg_array = array(F2_ecotype_temp_matrix_18deg, dim = c(27,2, 1))

# identifiers %>%
#   filter(Morph == 'Cold') %>%
#   View()

F2_ecotype_12deg_range = c(1:121, 244:363, 474:593, 714:833)
F2_ecotype_18deg_range = c(122:243, 364:473, 594:713, 834:931)

F2_ecotype_array = array(0, dim = c(27, 2, 931))
for(i in F2_ecotype_12deg_range){
  F2_ecotype_array[,,i] = F2_gpa$coords[,,i] - F2_ecotype_temp_12deg_array[,,1]
}

for(i in F2_ecotype_18deg_range){
  F2_ecotype_array[,,i] = F2_gpa$coords[,,i] - F2_ecotype_temp_18deg_array[,,1]
}

## This is the array to use too pull out the linear traits due
## to plasticity
F2_ecotype_array_consensus = array(0, dim = c(27, 2, 931))
for(i in 1:931){
  F2_ecotype_array_consensus[,,i] = F2_ecotype_array[,,i] + mean_shape_array[,,1]
}


ecotype_effect_gpa = gpagen(F2_ecotype_array_consensus)
# test_lm = geomorph.data.frame(plasticity_gpa)


F2_ecotype_effect_data = geomorph.data.frame(ecotype_effect_gpa, 
                                            Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                            parent_temp = identifiers$Parent_temp, 
                                            offspring_temp = identifiers$Offspring_temp,
                                            grand_temp = identifiers$Grand_temp,
                                            morph = identifiers$Morph, 
                                            population = identifiers$Lake,
                                            lake_morph = identifiers$Lake_morph,
                                            lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)


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

# WC_ecotype_residuals = arrayspecs(WC_ecotype_residuals, 
#                                   27, 
#                                   2)

E = F2_ecotype_effect_data$coords
# A = F2_whole_body_gpa$coords
F2_ecotype_effect_traits = interlmkdist(E, 
                                       lmks)


F2_ecotype_plasticity_traits = F2_ecotype_effect_traits %>% 
  as.data.frame() %>% 
  rownames_to_column() 
# %>% 
#   arrange(rowname)


F2_ecotype_plasticity_traits = bind_cols(F2_ecotype_plasticity_traits, 
                                        identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)


ecotype_plasticity_traits = F2_ecotype_plasticity_traits %>% 
  as_tibble() %>% 
  group_by(lake_morph_Pair_Full_Temp) %>% 
  select(jaw_length:body_length)

vars_keep = names(ecotype_plasticity_traits)[c(2,3,4,5,6,7,8,9,10,11)]
ecotype_plasticity_trait_cor = ecotype_plasticity_traits %>% 
  ungroup() %>% 
  split(.$lake_morph_Pair_Full_Temp) %>% 
  # ungroup() %>% 
  map(select, vars_keep) %>% 
  map(cor)

ecotype_plasticity_graph = ecotype_plasticity_trait_cor %>% 
  reshape2::melt() %>% 
  rename(lake_morph_full = L1)

ecotype_plasticity_trait_cor_graph = ggplot(ecotype_plasticity_graph, 
                                           aes(x = Var1, 
                                               y = Var2, 
                                               fill = value))+
  geom_tile()+
  facet_wrap(~lake_morph_full, 
             ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1))

ggsave('Univariate_ecotype_temp_plasticity_trait_integration.tiff', 
       plot = ecotype_plasticity_trait_cor_graph, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 40)



# Matrix comparisons ------------------------------------------------------

# install.packages('lineup')
library(lineup)

## Need the traits for each ecotype*F1*F2 pair
## Or maybe just for each ecotype? 

F2_univariate_traits = F2_univariate_traits %>% 
  select(-rowname)
F2_off_plasticity_traits = F2_off_plasticity_traits %>% 
  select(-rowname)
F2_parent_plasticity_traits = F2_parent_plasticity_traits %>% 
  select(-rowname)
F2_ecotype_plasticity_traits = F2_ecotype_plasticity_traits %>% 
  select(-rowname)



# ASHNC matrix compare ----------------------------------------------------


ASHNC_original = F2_univariate_traits %>% 
  # select(-rowname) %>% 
  filter(Lake_morph == 'ASHNC') %>% 
  select(jaw_length:body_length)
ASHNC_F2_temp = F2_off_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'ASHNC')%>% 
  select(jaw_length:body_length)
ASHNC_parent_temp = F2_parent_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'ASHNC')%>% 
  select(jaw_length:body_length)
ASHNC_ecotype = F2_ecotype_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'ASHNC')%>% 
  select(jaw_length:body_length)

corbetw2mat(ASHNC_original, 
            ASHNC_F2_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(ASHNC_original, 
            ASHNC_parent_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(ASHNC_original, 
            ASHNC_ecotype, 
            what = 'paired', 
            corthresh = 0.7)

ASHNC_F2_temp = corbetw2mat(ASHNC_original, 
            ASHNC_F2_temp, 
            what = 'all', 
            corthresh = 0.7)

ASHNC_F1_temp = corbetw2mat(ASHNC_original, 
            ASHNC_parent_temp, 
            what = 'all', 
            corthresh = 0.7)

ASHNC_ecotype = corbetw2mat(ASHNC_original, 
            ASHNC_ecotype, 
            what = 'all', 
            corthresh = 0.7)


ASHNC_F2_temp = ASHNC_F2_temp %>% 
      reshape2::melt() %>% 
  mutate(comparison = 'F2_original', 
         morph = 'ASHNC')
ASHNC_F1_temp = ASHNC_F1_temp %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'F1_original', 
         morph = 'ASHNC')
ASHNC_ecotype = ASHNC_ecotype %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'Ecotype_original', 
         morph = 'ASHNC')

ASHNC_df = bind_rows(ASHNC_F2_temp, 
                     ASHNC_F1_temp, 
                     ASHNC_ecotype)


# ASHNW matrix compare ----------------------------------------------------

ASHNW_original = F2_univariate_traits %>% 
  # select(-rowname) %>% 
  filter(Lake_morph == 'ASHNW') %>% 
  select(jaw_length:body_length)
ASHNW_F2_temp = F2_off_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'ASHNW')%>% 
  select(jaw_length:body_length)
ASHNW_parent_temp = F2_parent_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'ASHNW')%>% 
  select(jaw_length:body_length)
ASHNW_ecotype = F2_ecotype_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'ASHNW')%>% 
  select(jaw_length:body_length)

corbetw2mat(ASHNW_original, 
            ASHNW_F2_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(ASHNW_original, 
            ASHNW_parent_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(ASHNW_original, 
            ASHNW_ecotype, 
            what = 'paired', 
            corthresh = 0.7)

ASHNW_F2_temp = corbetw2mat(ASHNW_original, 
                            ASHNW_F2_temp, 
                            what = 'all', 
                            corthresh = 0.7)

ASHNW_F1_temp = corbetw2mat(ASHNW_original, 
                            ASHNW_parent_temp, 
                            what = 'all', 
                            corthresh = 0.7)

ASHNW_ecotype = corbetw2mat(ASHNW_original, 
                            ASHNW_ecotype, 
                            what = 'all', 
                            corthresh = 0.7)


ASHNW_F2_temp = ASHNW_F2_temp %>% 
  reshape2::melt() %>% 
  mutate(comparison = 'F2_original', 
         morph = 'ASHNW')
ASHNW_F1_temp = ASHNW_F1_temp %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'F1_original', 
         morph = 'ASHNW')
ASHNW_ecotype = ASHNW_ecotype %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'Ecotype_original', 
         morph = 'ASHNW')

ASHNW_df = bind_rows(ASHNW_F2_temp, 
                     ASHNW_F1_temp, 
                     ASHNW_ecotype)

# MYVC matrix compare ----------------------------------------------------

MYVC_original = F2_univariate_traits %>% 
  # select(-rowname) %>% 
  filter(Lake_morph == 'MYVC') %>% 
  select(jaw_length:body_length)
MYVC_F2_temp = F2_off_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'MYVC')%>% 
  select(jaw_length:body_length)
MYVC_parent_temp = F2_parent_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'MYVC')%>% 
  select(jaw_length:body_length)
MYVC_ecotype = F2_ecotype_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'MYVC')%>% 
  select(jaw_length:body_length)

corbetw2mat(MYVC_original, 
            MYVC_F2_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(MYVC_original, 
            MYVC_parent_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(MYVC_original, 
            MYVC_ecotype, 
            what = 'paired', 
            corthresh = 0.7)

MYVC_F2_temp = corbetw2mat(MYVC_original, 
                            MYVC_F2_temp, 
                            what = 'all', 
                            corthresh = 0.7)

MYVC_F1_temp = corbetw2mat(MYVC_original, 
                            MYVC_parent_temp, 
                            what = 'all', 
                            corthresh = 0.7)

MYVC_ecotype = corbetw2mat(MYVC_original, 
                            MYVC_ecotype, 
                            what = 'all', 
                            corthresh = 0.7)


MYVC_F2_temp = MYVC_F2_temp %>% 
  reshape2::melt() %>% 
  mutate(comparison = 'F2_original', 
         morph = 'MYVC')
MYVC_F1_temp = MYVC_F1_temp %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'F1_original', 
         morph = 'MYVC')
MYVC_ecotype = MYVC_ecotype %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'Ecotype_original', 
         morph = 'MYVC')

MYVC_df = bind_rows(MYVC_F2_temp, 
                     MYVC_F1_temp, 
                     MYVC_ecotype)

# MYVW matrix compare ----------------------------------------------------

MYVW_original = F2_univariate_traits %>% 
  # select(-rowname) %>% 
  filter(Lake_morph == 'MYVW') %>% 
  select(jaw_length:body_length)
MYVW_F2_temp = F2_off_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'MYVW')%>% 
  select(jaw_length:body_length)
MYVW_parent_temp = F2_parent_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'MYVW')%>% 
  select(jaw_length:body_length)
MYVW_ecotype = F2_ecotype_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'MYVW')%>% 
  select(jaw_length:body_length)

corbetw2mat(MYVW_original, 
            MYVW_F2_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(MYVW_original, 
            MYVW_parent_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(MYVW_original, 
            MYVW_ecotype, 
            what = 'paired', 
            corthresh = 0.7)
MYVW_F2_temp = corbetw2mat(MYVW_original, 
                            MYVW_F2_temp, 
                            what = 'all', 
                            corthresh = 0.7)

MYVW_F1_temp = corbetw2mat(MYVW_original, 
                            MYVW_parent_temp, 
                            what = 'all', 
                            corthresh = 0.7)

MYVW_ecotype = corbetw2mat(MYVW_original, 
                            MYVW_ecotype, 
                            what = 'all', 
                            corthresh = 0.7)


MYVW_F2_temp = MYVW_F2_temp %>% 
  reshape2::melt() %>% 
  mutate(comparison = 'F2_original', 
         morph = 'MYVW')
MYVW_F1_temp = MYVW_F1_temp %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'F1_original', 
         morph = 'MYVW')
MYVW_ecotype = MYVW_ecotype %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'Ecotype_original', 
         morph = 'MYVW')

MYVW_df = bind_rows(MYVW_F2_temp, 
                     MYVW_F1_temp, 
                     MYVW_ecotype)


# SKRC matrix compare ----------------------------------------------------

SKRC_original = F2_univariate_traits %>% 
  # select(-rowname) %>% 
  filter(Lake_morph == 'SKRC') %>% 
  select(jaw_length:body_length)
SKRC_F2_temp = F2_off_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'SKRC')%>% 
  select(jaw_length:body_length)
SKRC_parent_temp = F2_parent_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'SKRC')%>% 
  select(jaw_length:body_length)
SKRC_ecotype = F2_ecotype_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'SKRC')%>% 
  select(jaw_length:body_length)

corbetw2mat(SKRC_original, 
            SKRC_F2_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(SKRC_original, 
            SKRC_parent_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(SKRC_original, 
            SKRC_ecotype, 
            what = 'paired', 
            corthresh = 0.7)

SKRC_F2_temp = corbetw2mat(SKRC_original, 
                            SKRC_F2_temp, 
                            what = 'all', 
                            corthresh = 0.7)

SKRC_F1_temp = corbetw2mat(SKRC_original, 
                            SKRC_parent_temp, 
                            what = 'all', 
                            corthresh = 0.7)

SKRC_ecotype = corbetw2mat(SKRC_original, 
                            SKRC_ecotype, 
                            what = 'all', 
                            corthresh = 0.7)


SKRC_F2_temp = SKRC_F2_temp %>% 
  reshape2::melt() %>% 
  mutate(comparison = 'F2_original', 
         morph = 'SKRC')
SKRC_F1_temp = SKRC_F1_temp %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'F1_original', 
         morph = 'SKRC')
SKRC_ecotype = SKRC_ecotype %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'Ecotype_original', 
         morph = 'SKRC')

SKRC_df = bind_rows(SKRC_F2_temp, 
                     SKRC_F1_temp, 
                     SKRC_ecotype)

# SKRW matrix compare ----------------------------------------------------

SKRW_original = F2_univariate_traits %>% 
  # select(-rowname) %>% 
  filter(Lake_morph == 'SKRW') %>% 
  select(jaw_length:body_length)
SKRW_F2_temp = F2_off_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'SKRW')%>% 
  select(jaw_length:body_length)
SKRW_parent_temp = F2_parent_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'SKRW')%>% 
  select(jaw_length:body_length)
SKRW_ecotype = F2_ecotype_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'SKRW')%>% 
  select(jaw_length:body_length)

corbetw2mat(SKRW_original, 
            SKRW_F2_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(SKRW_original, 
            SKRW_parent_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(SKRW_original, 
            SKRW_ecotype, 
            what = 'paired', 
            corthresh = 0.7)

SKRW_F2_temp = corbetw2mat(SKRW_original, 
                            SKRW_F2_temp, 
                            what = 'all', 
                            corthresh = 0.7)

SKRW_F1_temp = corbetw2mat(SKRW_original, 
                            SKRW_parent_temp, 
                            what = 'all', 
                            corthresh = 0.7)

SKRW_ecotype = corbetw2mat(SKRW_original, 
                            SKRW_ecotype, 
                            what = 'all', 
                            corthresh = 0.7)


SKRW_F2_temp = SKRW_F2_temp %>% 
  reshape2::melt() %>% 
  mutate(comparison = 'F2_original', 
         morph = 'SKRW')
SKRW_F1_temp = SKRW_F1_temp %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'F1_original', 
         morph = 'SKRW')
SKRW_ecotype = SKRW_ecotype %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'Ecotype_original', 
         morph = 'SKRW')

SKRW_df = bind_rows(SKRW_F2_temp, 
                     SKRW_F1_temp, 
                     SKRW_ecotype)

# CSWYC matrix compare ----------------------------------------------------

CSWYC_original = F2_univariate_traits %>% 
  # select(-rowname) %>% 
  filter(Lake_morph == 'CSWYC') %>% 
  select(jaw_length:body_length)
CSWYC_F2_temp = F2_off_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'CSWYC')%>% 
  select(jaw_length:body_length)
CSWYC_parent_temp = F2_parent_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'CSWYC')%>% 
  select(jaw_length:body_length)
CSWYC_ecotype = F2_ecotype_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'CSWYC')%>% 
  select(jaw_length:body_length)

corbetw2mat(CSWYC_original, 
            CSWYC_F2_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(CSWYC_original, 
            CSWYC_parent_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(CSWYC_original, 
            CSWYC_ecotype, 
            what = 'paired', 
            corthresh = 0.7)
CSWYC_F2_temp = corbetw2mat(CSWYC_original, 
                            CSWYC_F2_temp, 
                            what = 'all', 
                            corthresh = 0.7)

CSWYC_F1_temp = corbetw2mat(CSWYC_original, 
                            CSWYC_parent_temp, 
                            what = 'all', 
                            corthresh = 0.7)

CSWYC_ecotype = corbetw2mat(CSWYC_original, 
                            CSWYC_ecotype, 
                            what = 'all', 
                            corthresh = 0.7)


CSWYC_F2_temp = CSWYC_F2_temp %>% 
  reshape2::melt() %>% 
  mutate(comparison = 'F2_original', 
         morph = 'CSWYC')
CSWYC_F1_temp = CSWYC_F1_temp %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'F1_original', 
         morph = 'CSWYC')
CSWYC_ecotype = CSWYC_ecotype %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'Ecotype_original', 
         morph = 'CSWYC')

CSWYC_df = bind_rows(CSWYC_F2_temp, 
                     CSWYC_F1_temp, 
                     CSWYC_ecotype)

# GTSW matrix compare ----------------------------------------------------

GTSW_original = F2_univariate_traits %>% 
  # select(-rowname) %>% 
  filter(Lake_morph == 'GTSW') %>% 
  select(jaw_length:body_length)
GTSW_F2_temp = F2_off_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'GTSW')%>% 
  select(jaw_length:body_length)
GTSW_parent_temp = F2_parent_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'GTSW')%>% 
  select(jaw_length:body_length)
GTSW_ecotype = F2_ecotype_plasticity_traits %>% 
  # select(-rowname)%>% 
  filter(Lake_morph == 'GTSW')%>% 
  select(jaw_length:body_length)

corbetw2mat(GTSW_original, 
            GTSW_F2_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(GTSW_original, 
            GTSW_parent_temp, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(GTSW_original, 
            GTSW_ecotype, 
            what = 'paired', 
            corthresh = 0.7)
GTSW_F2_temp = corbetw2mat(GTSW_original, 
                            GTSW_F2_temp, 
                            what = 'all', 
                            corthresh = 0.7)

GTSW_F1_temp = corbetw2mat(GTSW_original, 
                            GTSW_parent_temp, 
                            what = 'all', 
                            corthresh = 0.7)

GTSW_ecotype = corbetw2mat(GTSW_original, 
                            GTSW_ecotype, 
                            what = 'all', 
                            corthresh = 0.7)


GTSW_F2_temp = GTSW_F2_temp %>% 
  reshape2::melt() %>% 
  mutate(comparison = 'F2_original', 
         morph = 'GTSW')
GTSW_F1_temp = GTSW_F1_temp %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'F1_original', 
         morph = 'GTSW')
GTSW_ecotype = GTSW_ecotype %>% 
  reshape2::melt()%>% 
  mutate(comparison = 'Ecotype_original', 
         morph = 'GTSW')

GTSW_df = bind_rows(GTSW_F2_temp, 
                     GTSW_F1_temp, 
                     GTSW_ecotype)


# BIGASS DATAFRAME - FINAL FORM -------------------------------------------

big_ass_df = bind_rows(ASHNC_df, 
                       ASHNW_df, 
                       MYVC_df, 
                       MYVW_df, 
                       SKRC_df, 
                       SKRW_df, 
                       CSWYC_df, 
                       GTSW_df)

F2_effects = big_ass_df %>% 
  filter(comparison == 'F2_original')

F2_effect_graph = ggplot(F2_effects,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile()+
  facet_grid(. ~ morph + comparison)+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

ggsave('Effects_F2_temp_integration.tiff', 
       plot = F2_effect_graph, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 20)

F1_effects = big_ass_df %>% 
  filter(comparison == 'F1_original')

F1_effect_graph = ggplot(F1_effects,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile()+
  facet_grid(. ~ morph + comparison)+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

ggsave('Effects_F1_temp_integration.tiff', 
       plot = F1_effect_graph, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 20)


Ecotype_effects = big_ass_df %>% 
  filter(comparison == 'Ecotype_original')

ecotype_effect_graph = ggplot(Ecotype_effects,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile()+
  facet_grid(. ~ morph + comparison)+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

ggsave('Effects_ecotype_temp_integration.tiff', 
       plot = ecotype_effect_graph, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 20)


 integration_effects = ggplot(big_ass_df,
        aes(x = Var1,
            y = Var2,
            fill = value))+
  geom_tile()+
   facet_grid(. ~ morph + comparison)+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))


 ggsave('Effects_Factors_on_integration.tiff', 
        plot = integration_effects, 
        dpi = 'retina', 
        units = 'cm', 
        width = 70, 
        height = 20)
 