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
                  # lm_1_23 = c(1, 13),
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
        remove = F) %>% 
  mutate(ratio1 = lm_1_23/fbar_23_27, 
         ratio2 = lm_1_23/lm_23_2) %>% 
  select(rowname, 
         jaw_length:lm_23_2, 
         ratio1:ratio2, 
         everything())

orig_uni_traits = F2_univariate_traits %>%
  as_tibble() %>%
  group_by(Lake_morph) %>%
  select(jaw_length:ratio2)

vars_keep = names(orig_uni_traits)[c(2,3,4,5,6,7,8,9,10,11, 
                                           12,13,14,15,16,17,18, 
                                           19,20,21,22,23,24,25,26,
                                           27,28,29)]
orig_uni_trait_cor = orig_uni_traits %>%
  ungroup() %>%
  # split(.$lake_morph_Pair_Full_Temp) %>%
  split(.$Lake_morph) %>% 
  # ungroup() %>%
  map(select, vars_keep) %>%
  map(cor)

orig_uni_graph = orig_uni_trait_cor %>%
  reshape2::melt() %>%
  rename(lake_morph = L1)

orig_uni_trait_cor_graph = ggplot(orig_uni_graph,
                                        aes(x = Var1,
                                            y = Var2,
                                            fill = value))+
  geom_tile()+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  facet_wrap(~lake_morph,
             ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

ggsave('Univariate_original_plasticity_trait_ecotype_integration.tiff',
       plot = orig_uni_trait_cor_graph,
       dpi = 'retina',
       units = 'cm',
       width = 35,
       height = 20)



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
        remove = F) %>% 
  mutate(ratio1 = lm_1_23/fbar_23_27, 
         ratio2 = lm_1_23/lm_23_2) %>% 
  select(rowname, 
         jaw_length:lm_23_2, 
         ratio1:ratio2, 
         everything())


off_plasticity_traits = F2_off_plasticity_traits %>%
  as_tibble() %>%
  group_by(Lake_morph) %>%
  select(jaw_length:ratio2)

vars_keep = names(off_plasticity_traits)[c(2,3,4,5,6,7,8,9,10,11, 
                                           12,13,14,15,16,17,18, 
                                           19,20,21,22,23,24,25,26,
                                           27,28,29)]
off_plasticity_trait_cor = off_plasticity_traits %>%
  ungroup() %>%
  # split(.$lake_morph_Pair_Full_Temp) %>%
  split(.$Lake_morph) %>% 
  # ungroup() %>%
  map(select, vars_keep) %>%
  map(cor)

off_plasticity_graph = off_plasticity_trait_cor %>%
  reshape2::melt() %>%
  rename(lake_morph = L1)

off_plasticity_trait_cor_graph = ggplot(off_plasticity_graph,
                                    aes(x = Var1,
                                        y = Var2,
                                        fill = value))+
  geom_tile()+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  facet_wrap(~lake_morph,
             ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

ggsave('Univariate_offtemp_plasticity_trait_ecotype_integration.tiff',
       plot = off_plasticity_trait_cor_graph,
       dpi = 'retina',
       units = 'cm',
       width = 35,
       height = 20)

# ASHN F2 effect matrix compare -------------------------------------------

## ASHN
off_plasticity_trait_cor$ASHNC
off_plasticity_trait_cor$ASHNW

test_ASHNC = off_plasticity_trait_cor$ASHNC %>%
  reshape2::melt() 

ggplot(test_ASHNC,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

test_ASHNW = off_plasticity_trait_cor$ASHNW %>%
  reshape2::melt() 

ASHNC_cor_plot = ggplot(test_ASHNW,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))



ASHNW_cor_plot = ASHN_F2_temp_cor = corbetw2mat(off_plasticity_trait_cor$ASHNC, 
            off_plasticity_trait_cor$ASHNW, 
            what = 'all', 
            corthresh = 0.7)

ASHN_F2_effect = ASHN_F2_temp_cor %>%
  reshape2::melt() 

ASHN_F2_plast = ggplot(ASHN_F2_effect,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ade8f4",
                       mid = "#ff006e",
                       high = "#ade8f4") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())

# MYV F2 effect matrix compare -------------------------------------------
off_plasticity_trait_cor$MYVC
off_plasticity_trait_cor$MYVW

corbetw2mat(off_plasticity_trait_cor$MYVC, 
            off_plasticity_trait_cor$MYVW, 
            what = 'paired', 
            corthresh = 0.7)

MYV_F2_temp_cor = corbetw2mat(off_plasticity_trait_cor$MYVC, 
            off_plasticity_trait_cor$MYVW, 
            what = 'all', 
            corthresh = 0.7)

MYVC_F2_temp = off_plasticity_trait_cor$MYVC %>%
  reshape2::melt() 

MYVC_cor_plot = ggplot(MYVC_F2_temp,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

MYVW_F2_temp = off_plasticity_trait_cor$MYVW %>%
  reshape2::melt() 

MYVW_cor_plot = ggplot(MYVW_F2_temp,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))



MYV_F2_temp_cor = corbetw2mat(off_plasticity_trait_cor$MYVC, 
                               off_plasticity_trait_cor$MYVW, 
                               what = 'all', 
                               corthresh = 0.7)

MYV_F2_temp_cor = MYV_F2_temp_cor %>%
  reshape2::melt() 

MYV_F2_plasticity = ggplot(MYV_F2_temp_cor,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ade8f4",
                       mid = "#ff006e",
                       high = "#ade8f4") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())

# SKR F2 effect matrix compare -------------------------------------------
off_plasticity_trait_cor$SKRC
off_plasticity_trait_cor$SKRW

corbetw2mat(off_plasticity_trait_cor$SKRC, 
            off_plasticity_trait_cor$SKRW, 
            what = 'paired', 
            corthresh = 0.7)

SKR_F2_temp_cor = corbetw2mat(off_plasticity_trait_cor$SKRC, 
            off_plasticity_trait_cor$SKRW, 
            what = 'all', 
            corthresh = 0.7)

SKRC_F2_temp = off_plasticity_trait_cor$SKRC %>%
  reshape2::melt() 

SKRC_cor_plot = ggplot(SKRC_F2_temp,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

SKRW_F2_temp = off_plasticity_trait_cor$SKRW %>%
  reshape2::melt() 

SKRW_cor_plot = ggplot(SKRW_F2_temp,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))



SKR_F2_temp_cor = corbetw2mat(off_plasticity_trait_cor$SKRC, 
                              off_plasticity_trait_cor$SKRW, 
                              what = 'all', 
                              corthresh = 0.7)

SKR_F2_temp_cor = SKR_F2_temp_cor %>%
  reshape2::melt() 

SKR_F2_plasticity = ggplot(SKR_F2_temp_cor,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ade8f4",
                       mid = "#ff006e",
                       high = "#ade8f4") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())


# GTSCSWY F2 effect matrix compare -------------------------------------------
off_plasticity_trait_cor$CSWYC
off_plasticity_trait_cor$GTSW

corbetw2mat(off_plasticity_trait_cor$CSWYC, 
            off_plasticity_trait_cor$GTSW, 
            what = 'paired', 
            corthresh = 0.7)

GTS_CSWY_F2_temp_cor = corbetw2mat(off_plasticity_trait_cor$CSWYC, 
            off_plasticity_trait_cor$GTSW, 
            what = 'all', 
            corthresh = 0.7)

CSWYC_F2_temp = off_plasticity_trait_cor$CSWYC %>%
  reshape2::melt() 

CSWY_cor_plot = ggplot(CSWYC_F2_temp,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

GTSW_F2_temp = off_plasticity_trait_cor$GTSW %>%
  reshape2::melt() 

GTS_cor_plot = ggplot(GTSW_F2_temp,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))



GTSCSWY_F2_temp_cor = corbetw2mat(off_plasticity_trait_cor$CSWYC, 
                              off_plasticity_trait_cor$GTSW, 
                              what = 'all', 
                              corthresh = 0.7)

GTSCSWY_F2_temp_cor = GTSCSWY_F2_temp_cor %>%
  reshape2::melt() 

GTSCSWY_F2_plasticity = ggplot(GTSCSWY_F2_temp_cor,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ade8f4",
                       mid = "#ff006e",
                       high = "#ade8f4") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())

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

# identifiers %>%
#   filter(Parent_temp == '12') %>%
#   View()

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
        remove = F) %>% 
  mutate(ratio1 = lm_1_23/fbar_23_27, 
         ratio2 = lm_1_23/lm_23_2) %>% 
  select(rowname, 
         jaw_length:lm_23_2, 
         ratio1:ratio2, 
         everything())


parent_plasticity_traits = F2_parent_plasticity_traits %>%
  as_tibble() %>%
  group_by(Lake_morph) %>%
  select(jaw_length:ratio2)

vars_keep = names(off_plasticity_traits)[c(2,3,4,5,6,7,8,9,10,11, 
                                           12,13,14,15,16,17,18, 
                                           19,20,21,22,23,24,25,26,
                                           27,28,29)]
parent_plasticity_trait_cor = parent_plasticity_traits %>%
  ungroup() %>%
  split(.$Lake_morph) %>%
  # ungroup() %>%
  map(select, vars_keep) %>%
  map(cor)

parent_plasticity_graph = parent_plasticity_trait_cor %>%
  reshape2::melt() %>%
  rename(lake_morph = L1)

parent_plasticity_trait_cor_graph = ggplot(parent_plasticity_graph,
                                        aes(x = Var1,
                                            y = Var2,
                                            fill = value))+
  geom_tile()+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  facet_wrap(~lake_morph,
             ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

ggsave('Univariate_parent_temp_plasticity_trait_ecotype_integration.tiff',
       plot = parent_plasticity_trait_cor_graph,
       dpi = 'retina',
       units = 'cm',
       width = 35,
       height = 20)

# ASHN F1 effect matrix compare -------------------------------------------

## ASHN

ASHNC_F1_cor = parent_plasticity_trait_cor$ASHNC %>%
  reshape2::melt() 

ASHNC_F1_cor_plot = ggplot(ASHNC_F1_cor,
       aes(x = Var1,
           y = Var2,
           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

ASHNW_F1_cor = parent_plasticity_trait_cor$ASHNW %>%
  reshape2::melt() 

ASHNW_cor_plot = ggplot(ASHNW_F1_cor,
                        aes(x = Var1,
                            y = Var2,
                            fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))



ASHN_F1_temp_cor = corbetw2mat(parent_plasticity_trait_cor$ASHNC, 
                                                parent_plasticity_trait_cor$ASHNW, 
                                                what = 'all', 
                                                corthresh = 0.7)

ASHN_F1_effect = ASHN_F1_temp_cor %>%
  reshape2::melt() 

ASHN_F1_plast = ggplot(ASHN_F1_effect,
                       aes(x = Var1,
                           y = Var2,
                           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ade8f4",
                       mid = "#ff006e",
                       high = "#ade8f4") +
  labs(title = 'A) ASHN')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

# MYV F2 effect matrix compare -------------------------------------------
corbetw2mat(parent_plasticity_trait_cor$MYVC, 
            parent_plasticity_trait_cor$MYVW, 
            what = 'paired', 
            corthresh = 0.7)

MYV_F1_temp_cor = corbetw2mat(parent_plasticity_trait_cor$MYVC, 
                              parent_plasticity_trait_cor$MYVW, 
                              what = 'all', 
                              corthresh = 0.7)

MYVC_F1_temp = parent_plasticity_trait_cor$MYVC %>%
  reshape2::melt() 

MYVC_cor_plot = ggplot(MYVC_F1_temp,
                       aes(x = Var1,
                           y = Var2,
                           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

MYVW_F1_temp = parent_plasticity_trait_cor$MYVW %>%
  reshape2::melt() 

MYVW_cor_plot = ggplot(MYVW_F1_temp,
                       aes(x = Var1,
                           y = Var2,
                           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))



MYV_F1_temp_cor = corbetw2mat(parent_plasticity_trait_cor$MYVC, 
                              parent_plasticity_trait_cor$MYVW, 
                              what = 'all', 
                              corthresh = 0.7)

MYV_F1_temp_cor = MYV_F1_temp_cor %>%
  reshape2::melt() 

MYV_F1_plasticity = ggplot(MYV_F1_temp_cor,
                           aes(x = Var1,
                               y = Var2,
                               fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ade8f4",
                       mid = "#ff006e",
                       high = "#ade8f4") +
  labs(title = 'B) MYV')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1), 
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

# SKR F2 effect matrix compare -------------------------------------------

corbetw2mat(parent_plasticity_trait_cor$SKRC, 
            parent_plasticity_trait_cor$SKRW, 
            what = 'paired', 
            corthresh = 0.7)

SKR_F1_temp_cor = corbetw2mat(parent_plasticity_trait_cor$SKRC, 
                              parent_plasticity_trait_cor$SKRW, 
                              what = 'all', 
                              corthresh = 0.7)

SKRC_F1_temp = parent_plasticity_trait_cor$SKRC %>%
  reshape2::melt() 

SKRC_cor_plot = ggplot(SKRC_F1_temp,
                       aes(x = Var1,
                           y = Var2,
                           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

SKRW_F1_temp = parent_plasticity_trait_cor$SKRW %>%
  reshape2::melt() 

SKRW_cor_plot = ggplot(SKRW_F1_temp,
                       aes(x = Var1,
                           y = Var2,
                           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))



SKR_F1_temp_cor = corbetw2mat(parent_plasticity_trait_cor$SKRC, 
                              parent_plasticity_trait_cor$SKRW, 
                              what = 'all', 
                              corthresh = 0.7)

SKR_F1_temp_cor = SKR_F1_temp_cor %>%
  reshape2::melt() 

SKR_F1_plasticity = ggplot(SKR_F1_temp_cor,
                           aes(x = Var1,
                               y = Var2,
                               fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ade8f4",
                       mid = "#ff006e",
                       high = "#ade8f4") +
  labs(title = 'C) SKR')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1), 
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())


# GTSCSWY F2 effect matrix compare -------------------------------------------
corbetw2mat(parent_plasticity_trait_cor$CSWYC, 
            parent_plasticity_trait_cor$GTSW, 
            what = 'paired', 
            corthresh = 0.7)

GTS_CSWY_F1_temp_cor = corbetw2mat(parent_plasticity_trait_cor$CSWYC, 
                                   parent_plasticity_trait_cor$GTSW, 
                                   what = 'all', 
                                   corthresh = 0.7)

CSWYC_F1_temp = parent_plasticity_trait_cor$CSWYC %>%
  reshape2::melt() 

CSWY_cor_plot = ggplot(CSWYC_F1_temp,
                       aes(x = Var1,
                           y = Var2,
                           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

GTSW_F1_temp = parent_plasticity_trait_cor$GTSW %>%
  reshape2::melt() 

GTS_cor_plot = ggplot(GTSW_F1_temp,
                      aes(x = Var1,
                          y = Var2,
                          fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#003049",
                       mid = "#fefae0",
                       high = "#d62828") +
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))



GTSCSWY_F1_temp_cor = corbetw2mat(parent_plasticity_trait_cor$CSWYC, 
                                  parent_plasticity_trait_cor$GTSW, 
                                  what = 'all', 
                                  corthresh = 0.7)

GTSCSWY_F1_temp_cor = GTSCSWY_F1_temp_cor %>%
  reshape2::melt() 

GTSCSWY_F1_plasticity = ggplot(GTSCSWY_F1_temp_cor,
                               aes(x = Var1,
                                   y = Var2,
                                   fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ade8f4",
                       mid = "#ff006e",
                       high = "#ade8f4") +
  labs(title = 'D) GTS-CSWY')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1), 
        legend.position = 'none')



# F1 + F2 combo graphs ----------------------------------------------------
library(patchwork)

ASHN_combo_graphs = ASHN_F1_plast | ASHN_F2_plast
MYV_combo_graphs = MYV_F1_plasticity | MYV_F2_plasticity
SKR_combo_graphs = SKR_F1_plasticity | SKR_F2_plasticity
GTSCSWY_combo_graphs = GTSCSWY_F1_plasticity | GTSCSWY_F2_plasticity

Big_graph = ASHN_combo_graphs/MYV_combo_graphs/SKR_combo_graphs/GTSCSWY_combo_graphs

ggsave('Figure1_Effects_F1_F2_on_Integration.tiff',
       plot = Big_graph,
       dpi = 'retina',
       units = 'cm',
       width = 40,
       height = 40)
##




# Compare F1 and F2 effects to original data ------------------------------

orig_uni_graph = orig_uni_trait_cor %>%
  reshape2::melt() %>%
  rename(lake_morph = L1)

off_plasticity_graph = off_plasticity_trait_cor %>%
  reshape2::melt() %>%
  rename(lake_morph = L1)

parent_plasticity_graph = parent_plasticity_trait_cor %>%
  reshape2::melt() %>%
  rename(lake_morph = L1)



# ASHNC F1+F2 = Original --------------------------------------------------
ASHNC_original = orig_uni_trait_cor$ASHNC 

ASHNC_F2 = off_plasticity_trait_cor$ASHNC 

ASHNC_F1 = parent_plasticity_trait_cor$ASHNC


ASHNC_F1_original = corbetw2mat(ASHNC_original, 
                                ASHNC_F1, 
                                what = 'all', 
                                corthresh = 0.7) 

ASHNC_F2_original = corbetw2mat(ASHNC_original, 
                                ASHNC_F2, 
                                what = 'all', 
                                corthresh = 0.7) 


ASHNC_F1_original = ASHNC_F1_original %>%
  reshape2::melt() 

ASHNC_F2_original = ASHNC_F2_original %>%
  reshape2::melt() 


ASHNC_F1_original_plot = ggplot(ASHNC_F1_original,
                       aes(x = Var1,
                           y = Var2,
                           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'A) ASHNC F1 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

ASHNC_F2_original_plot = ggplot(ASHNC_F2_original,
                                aes(x = Var1,
                                    y = Var2,
                                    fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'B) ASHNC F2 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        # legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

ASHNC_compare = ASHNC_F1_original_plot|ASHNC_F2_original_plot

# ASHNW F1+F2 = Original --------------------------------------------------
ASHNW_original = orig_uni_trait_cor$ASHNW 

ASHNW_F2 = off_plasticity_trait_cor$ASHNW 

ASHNW_F1 = parent_plasticity_trait_cor$ASHNW


ASHNW_F1_original = corbetw2mat(ASHNW_original, 
                                ASHNW_F1, 
                                what = 'all', 
                                corthresh = 0.7) 

ASHNW_F2_original = corbetw2mat(ASHNW_original, 
                                ASHNW_F2, 
                                what = 'all', 
                                corthresh = 0.7) 


ASHNW_F1_original = ASHNW_F1_original %>%
  reshape2::melt() 

ASHNW_F2_original = ASHNW_F2_original %>%
  reshape2::melt() 


ASHNW_F1_original_plot = ggplot(ASHNW_F1_original,
                                aes(x = Var1,
                                    y = Var2,
                                    fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'C) ASHNW F1 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

ASHNW_F2_original_plot = ggplot(ASHNW_F2_original,
                                aes(x = Var1,
                                    y = Var2,
                                    fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'D) ASHNW F2 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        # legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

ASHNW_compare = ASHNW_F1_original_plot|ASHNW_F2_original_plot

ASHN_Matrix_evo = ASHNC_compare/ASHNW_compare


# MYVC F1+F2 = Original --------------------------------------------------
MYVC_original = orig_uni_trait_cor$MYVC 

MYVC_F2 = off_plasticity_trait_cor$MYVC 

MYVC_F1 = parent_plasticity_trait_cor$MYVC


MYVC_F1_original = corbetw2mat(MYVC_original, 
                                MYVC_F1, 
                                what = 'all', 
                                corthresh = 0.7) 

MYVC_F2_original = corbetw2mat(MYVC_original, 
                                MYVC_F2, 
                                what = 'all', 
                                corthresh = 0.7) 


MYVC_F1_original = MYVC_F1_original %>%
  reshape2::melt() 

MYVC_F2_original = MYVC_F2_original %>%
  reshape2::melt() 


MYVC_F1_original_plot = ggplot(MYVC_F1_original,
                                aes(x = Var1,
                                    y = Var2,
                                    fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'A) MYVC F1 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

MYVC_F2_original_plot = ggplot(MYVC_F2_original,
                                aes(x = Var1,
                                    y = Var2,
                                    fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'B) MYVC F2 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        # legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

MYVC_compare = MYVC_F1_original_plot|MYVC_F2_original_plot

# MYVW F1+F2 = Original --------------------------------------------------
MYVW_original = orig_uni_trait_cor$MYVW 

MYVW_F2 = off_plasticity_trait_cor$MYVW 

MYVW_F1 = parent_plasticity_trait_cor$MYVW


MYVW_F1_original = corbetw2mat(MYVW_original, 
                                MYVW_F1, 
                                what = 'all', 
                                corthresh = 0.7) 

MYVW_F2_original = corbetw2mat(MYVW_original, 
                                MYVW_F2, 
                                what = 'all', 
                                corthresh = 0.7) 


MYVW_F1_original = MYVW_F1_original %>%
  reshape2::melt() 

MYVW_F2_original = MYVW_F2_original %>%
  reshape2::melt() 


MYVW_F1_original_plot = ggplot(MYVW_F1_original,
                                aes(x = Var1,
                                    y = Var2,
                                    fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'C) MYVW F1 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

MYVW_F2_original_plot = ggplot(MYVW_F2_original,
                                aes(x = Var1,
                                    y = Var2,
                                    fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'D) MYVW F2 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        # legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

MYVW_compare = MYVW_F1_original_plot|MYVW_F2_original_plot

MYV_Matrix_evo = MYVC_compare/MYVW_compare

# SKRC F1+F2 = Original --------------------------------------------------
SKRC_original = orig_uni_trait_cor$SKRC 

SKRC_F2 = off_plasticity_trait_cor$SKRC 

SKRC_F1 = parent_plasticity_trait_cor$SKRC


SKRC_F1_original = corbetw2mat(SKRC_original, 
                               SKRC_F1, 
                               what = 'all', 
                               corthresh = 0.7) 

SKRC_F2_original = corbetw2mat(SKRC_original, 
                               SKRC_F2, 
                               what = 'all', 
                               corthresh = 0.7) 


SKRC_F1_original = SKRC_F1_original %>%
  reshape2::melt() 

SKRC_F2_original = SKRC_F2_original %>%
  reshape2::melt() 


SKRC_F1_original_plot = ggplot(SKRC_F1_original,
                               aes(x = Var1,
                                   y = Var2,
                                   fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'A) SKRC F1 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

SKRC_F2_original_plot = ggplot(SKRC_F2_original,
                               aes(x = Var1,
                                   y = Var2,
                                   fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'B) SKRC F2 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        # legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

SKRC_compare = SKRC_F1_original_plot|SKRC_F2_original_plot

# SKRW F1+F2 = Original --------------------------------------------------
SKRW_original = orig_uni_trait_cor$SKRW 

SKRW_F2 = off_plasticity_trait_cor$SKRW 

SKRW_F1 = parent_plasticity_trait_cor$SKRW


SKRW_F1_original = corbetw2mat(SKRW_original, 
                               SKRW_F1, 
                               what = 'all', 
                               corthresh = 0.7) 

SKRW_F2_original = corbetw2mat(SKRW_original, 
                               SKRW_F2, 
                               what = 'all', 
                               corthresh = 0.7) 


SKRW_F1_original = SKRW_F1_original %>%
  reshape2::melt() 

SKRW_F2_original = SKRW_F2_original %>%
  reshape2::melt() 


SKRW_F1_original_plot = ggplot(SKRW_F1_original,
                               aes(x = Var1,
                                   y = Var2,
                                   fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'C) SKRW F1 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

SKRW_F2_original_plot = ggplot(SKRW_F2_original,
                               aes(x = Var1,
                                   y = Var2,
                                   fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'D) SKRW F2 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        # legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

SKRW_compare = SKRW_F1_original_plot|SKRW_F2_original_plot

SKR_Matrix_evo = SKRC_compare/SKRW_compare

# CSWYC F1+F2 = Original --------------------------------------------------
CSWYC_original = orig_uni_trait_cor$CSWYC 

CSWYC_F2 = off_plasticity_trait_cor$CSWYC 

CSWYC_F1 = parent_plasticity_trait_cor$CSWYC


CSWYC_F1_original = corbetw2mat(CSWYC_original, 
                               CSWYC_F1, 
                               what = 'all', 
                               corthresh = 0.7) 

CSWYC_F2_original = corbetw2mat(CSWYC_original, 
                               CSWYC_F2, 
                               what = 'all', 
                               corthresh = 0.7) 


CSWYC_F1_original = CSWYC_F1_original %>%
  reshape2::melt() 

CSWYC_F2_original = CSWYC_F2_original %>%
  reshape2::melt() 


CSWYC_F1_original_plot = ggplot(CSWYC_F1_original,
                               aes(x = Var1,
                                   y = Var2,
                                   fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'A) CSWYC F1 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

CSWYC_F2_original_plot = ggplot(CSWYC_F2_original,
                               aes(x = Var1,
                                   y = Var2,
                                   fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'B) CSWYC F2 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        # legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

CSWYC_compare = CSWYC_F1_original_plot|CSWYC_F2_original_plot

# GTSW F1+F2 = Original --------------------------------------------------
GTSW_original = orig_uni_trait_cor$GTSW 

GTSW_F2 = off_plasticity_trait_cor$GTSW 

GTSW_F1 = parent_plasticity_trait_cor$GTSW


GTSW_F1_original = corbetw2mat(GTSW_original, 
                               GTSW_F1, 
                               what = 'all', 
                               corthresh = 0.7) 

GTSW_F2_original = corbetw2mat(GTSW_original, 
                               GTSW_F2, 
                               what = 'all', 
                               corthresh = 0.7) 


GTSW_F1_original = GTSW_F1_original %>%
  reshape2::melt() 

GTSW_F2_original = GTSW_F2_original %>%
  reshape2::melt() 


GTSW_F1_original_plot = ggplot(GTSW_F1_original,
                               aes(x = Var1,
                                   y = Var2,
                                   fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'C) GTSW F1 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

GTSW_F2_original_plot = ggplot(GTSW_F2_original,
                               aes(x = Var1,
                                   y = Var2,
                                   fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#ffff3f",
                       mid = "#007f5f",
                       high = "#ffff3f") +
  labs(title = 'D) GTSW F2 vs Original')+
  # facet_wrap(~lake_morph_full,
  #            ncol = 4)+
  # facet_wrap(~lake_morph,
  #            ncol = 4)+
  theme_bw()+
  theme(strip.background = element_rect(fill = 'white'),
        strip.text = element_text(face = 'bold'),
        axis.title = element_blank(),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        # legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

GTSW_compare = GTSW_F1_original_plot|GTSW_F2_original_plot

GTSCSWY_Matrix_evo = CSWYC_compare/GTSW_compare



# nmds between trait matrices ---------------------------------------------

ASHNC_original = orig_uni_trait_cor$ASHNC %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  mutate(matrix_group = 'Original')

ASHNC_F2 = off_plasticity_trait_cor$ASHNC %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  mutate(matrix_group = 'F2 effect')

ASHNC_F1 = parent_plasticity_trait_cor$ASHNC%>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  mutate(matrix_group = 'F1 effect')




# Code I dont need but scared to delete -----------------------------------

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


# ASHN ecotype comparisons F2 temp effect ---------------------------------
ASHNC_F2_temp_effects = F2_off_plasticity_traits %>% 
  filter(Lake_morph == 'ASHNC')  %>% 
  select(jaw_length:ratio2)
  # View()
  # distinct(jaw_length)
# slice(1)

ASHNW_F2_temp_effects = F2_off_plasticity_traits %>% 
  filter(Lake_morph == 'ASHNW') %>% 
  slice(-122) %>% 
  select(jaw_length:ratio2)


corbetw2mat(ASHNC_F2_temp_effects, 
            ASHNW_F2_temp_effects, 
            what = 'paired', 
            corthresh = 0.7)

corbetw2mat(ASHNC_F2_temp_effects, 
            ASHNW_F2_temp_effects, 
            what = 'all', 
            corthresh = 0.7)



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
 