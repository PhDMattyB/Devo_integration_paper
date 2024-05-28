##############################
## Univaritate trait integration
##
## Matt Brachmann (PhDMattyB)
##
## 28.05.2024
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
                factor)) %>% 
  arrange(individualID)


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




# univariate traits - interlandmark distances -------------------------------------------------------

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

A = F2_gpa$coords
# A = F2_whole_body_gpa$coords
F2_univariate_traits = interlmkdist(A, 
                                    lmks)

# arrayspecs(F2_univariate_traits, 
#            4, 
#            3)

F2_univariate_traits = F2_univariate_traits %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  arrange(rowname)


F2_univariate_traits = bind_cols(F2_univariate_traits, 
                                 identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)


# trait correlations - with allometry ------------------------------------------------------

traits = F2_univariate_traits %>% 
  as_tibble() %>% 
  group_by(lake_morph_Pair_Full_Temp) %>% 
  select(jaw_length:body_length)

vars_keep = names(traits)[c(2,3,4,5,6,7,8,9,10,11)]
trait_cor = traits %>% 
  ungroup() %>% 
  split(.$lake_morph_Pair_Full_Temp) %>% 
  # ungroup() %>% 
  map(select, vars_keep) %>% 
  map(cor)

graph = trait_cor %>% 
  reshape2::melt() %>% 
  rename(lake_morph_full = L1)
  
trait_cor_graph = ggplot(graph, 
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

ggsave('Univariate_trait_integration.tiff', 
       plot = trait_cor_graph, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 40)

# univariate traits - interlandmark distances no allometry -------------------------------------------------------

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

B = F2_allometry_adj_shape
# A = F2_whole_body_gpa$coords
F2_noallo_traits = interlmkdist(B, 
                                    lmks)

# arrayspecs(F2_univariate_traits, 
#            4, 
#            3)

F2_noallo_traits = F2_noallo_traits %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  arrange(rowname)


F2_noallo_traits = bind_cols(F2_noallo_traits, 
                                 identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)

# trait correlations - NO allometry ------------------------------------------------------

noallo_traits = F2_noallo_traits %>% 
  as_tibble() %>% 
  group_by(lake_morph_Pair_Full_Temp) %>% 
  select(jaw_length:body_length)

vars_keep = names(noallo_traits)[c(2,3,4,5,6,7,8,9,10,11)]
noallo_trait_cor = noallo_traits %>% 
  ungroup() %>% 
  split(.$lake_morph_Pair_Full_Temp) %>% 
  # ungroup() %>% 
  map(select, vars_keep) %>% 
  map(cor)

noallo_graph = noallo_trait_cor %>% 
  reshape2::melt() %>% 
  rename(lake_morph_full = L1)

noallo_trait_cor_graph = ggplot(noallo_graph, 
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

ggsave('Univariate_NOallo_trait_integration.tiff', 
       plot = noallo_trait_cor_graph, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 40)



# Ecotype effects ----------------------------------------------

F2_tps = readland.tps('F2_No_GT.TPS', 
                      specID = 'imageID')

## superimposition on the entire dataset
F2_gpa = gpagen(F2_tps, 
                print.progress = F)

F2_data = geomorph.data.frame(coords = two.d.array(F2_gpa$coords), 
                                       Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                                       parent_temp = identifiers$Parent_temp, 
                                       offspring_temp = identifiers$Offspring_temp,
                                       grand_temp = identifiers$Grand_temp,
                                       morph = identifiers$Morph, 
                                       population = identifiers$Lake,
                                       lake_morph = identifiers$Lake_morph,
                                       lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)


# F2_landmark_data = two.d.array(F2_gpa$coords) %>%
#   as.data.frame() %>% 
#   rownames_to_column() %>% 
#   as_tibble() %>% 
#   arrange(rowname)
# 
# mod_data = bind_cols(F2_landmark_data, 
#                      identifiers)

### candisc won't work for landmark coordinates. 
### Here's the work around
WC_ecotype_mod = procD.lm(coords ~ offspring_temp*parent_temp,
                          data = F2_data)
# mod = lm(fixed_coords ~ Offspring_temp + Parent_temp, 
#                data = identifiers)
summary(WC_ecotype_mod)

WC_ecotype_residuals = WC_ecotype_mod$residuals

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

WC_ecotype_residuals = arrayspecs(WC_ecotype_residuals, 
           27, 
           2)

C = WC_ecotype_residuals
# A = F2_whole_body_gpa$coords
F2_WC_ecotype_traits = interlmkdist(C, 
                                lmks)


F2_WC_ecotype_traits = F2_WC_ecotype_traits %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  arrange(rowname)


F2_WC_ecotype_traits = bind_cols(F2_WC_ecotype_traits, 
                             identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)


WC_ecotype_traits = F2_WC_ecotype_traits %>% 
  as_tibble() %>% 
  group_by(lake_morph_Pair_Full_Temp) %>% 
  select(jaw_length:body_length)

vars_keep = names(WC_ecotype_traits)[c(2,3,4,5,6,7,8,9,10,11)]
WC_ecotype_trait_cor = WC_ecotype_traits %>% 
  ungroup() %>% 
  split(.$lake_morph_Pair_Full_Temp) %>% 
  # ungroup() %>% 
  map(select, vars_keep) %>% 
  map(cor)

WC_ecotype_graph = WC_ecotype_trait_cor %>% 
  reshape2::melt() %>% 
  rename(lake_morph_full = L1)

WC_ecotype_trait_cor_graph = ggplot(WC_ecotype_graph, 
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

ggsave('Univariate_WC_Ecotype_trait_integration.tiff', 
       plot = noallo_trait_cor_graph, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 40)

# Offspring temp effects ----------------------------------------------

F2_tps = readland.tps('F2_No_GT.TPS', 
                      specID = 'imageID')

## superimposition on the entire dataset
F2_gpa = gpagen(F2_tps, 
                print.progress = F)

F2_data = geomorph.data.frame(coords = two.d.array(F2_gpa$coords), 
                              Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                              parent_temp = identifiers$Parent_temp, 
                              offspring_temp = identifiers$Offspring_temp,
                              grand_temp = identifiers$Grand_temp,
                              morph = identifiers$Morph, 
                              population = identifiers$Lake,
                              lake_morph = identifiers$Lake_morph,
                              lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)


# F2_landmark_data = two.d.array(F2_gpa$coords) %>%
#   as.data.frame() %>% 
#   rownames_to_column() %>% 
#   as_tibble() %>% 
#   arrange(rowname)
# 
# mod_data = bind_cols(F2_landmark_data, 
#                      identifiers)

### candisc won't work for landmark coordinates. 
### Here's the work around
offspring_temp_mod = procD.lm(coords ~ morph*parent_temp,
                          data = F2_data)
# mod = lm(fixed_coords ~ Offspring_temp + Parent_temp, 
#                data = identifiers)
summary(offspring_temp_mod)

offspring_temp_residuals = offspring_temp_mod$residuals

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

offspring_temp_residuals = arrayspecs(offspring_temp_residuals, 
                                  27, 
                                  2)

D = offspring_temp_residuals
# A = F2_whole_body_gpa$coords
F2_offspring_temp_traits = interlmkdist(D, 
                                    lmks)


F2_offspring_temp_traits = F2_offspring_temp_traits %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  arrange(rowname)


F2_offspring_temp_traits = bind_cols(F2_offspring_temp_traits, 
                                 identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)


offspring_temp_traits = F2_offspring_temp_traits %>% 
  as_tibble() %>% 
  group_by(lake_morph_Pair_Full_Temp) %>% 
  select(jaw_length:body_length)

vars_keep = names(offspring_temp_traits)[c(2,3,4,5,6,7,8,9,10,11)]
offspring_temp_trait_cor = offspring_temp_traits %>% 
  ungroup() %>% 
  split(.$lake_morph_Pair_Full_Temp) %>% 
  # ungroup() %>% 
  map(select, vars_keep) %>% 
  map(cor)

offspring_temp_graph = offspring_temp_trait_cor %>% 
  reshape2::melt() %>% 
  rename(lake_morph_full = L1)

offspring_temp_trait_cor_graph = ggplot(offspring_temp_graph, 
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

ggsave('Univariate_offspring_temp_trait_integration.tiff', 
       plot = noallo_trait_cor_graph, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 40)

# Parent temp effects ----------------------------------------------

F2_tps = readland.tps('F2_No_GT.TPS', 
                      specID = 'imageID')

## superimposition on the entire dataset
F2_gpa = gpagen(F2_tps, 
                print.progress = F)

F2_data = geomorph.data.frame(coords = two.d.array(F2_gpa$coords), 
                              Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                              parent_temp = identifiers$Parent_temp, 
                              offspring_temp = identifiers$Offspring_temp,
                              grand_temp = identifiers$Grand_temp,
                              morph = identifiers$Morph, 
                              population = identifiers$Lake,
                              lake_morph = identifiers$Lake_morph,
                              lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)


# F2_landmark_data = two.d.array(F2_gpa$coords) %>%
#   as.data.frame() %>% 
#   rownames_to_column() %>% 
#   as_tibble() %>% 
#   arrange(rowname)
# 
# mod_data = bind_cols(F2_landmark_data, 
#                      identifiers)

### candisc won't work for landmark coordinates. 
### Here's the work around
parent_temp_mod = procD.lm(coords ~ morph*offspring_temp,
                              data = F2_data)
# mod = lm(fixed_coords ~ Offspring_temp + Parent_temp, 
#                data = identifiers)
summary(parent_temp_mod)

parent_temp_residuals = parent_temp_mod$residuals

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

parent_temp_residuals = arrayspecs(parent_temp_residuals, 
                                      27, 
                                      2)

E = parent_temp_residuals
# A = F2_whole_body_gpa$coords
F2_parent_temp_traits = interlmkdist(E, 
                                        lmks)


F2_parent_temp_traits = F2_parent_temp_traits %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  arrange(rowname)


F2_parent_temp_traits = bind_cols(F2_parent_temp_traits, 
                                     identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)


parent_temp_traits = F2_parent_temp_traits %>% 
  as_tibble() %>% 
  group_by(lake_morph_Pair_Full_Temp) %>% 
  select(jaw_length:body_length)

vars_keep = names(parent_temp_traits)[c(2,3,4,5,6,7,8,9,10,11)]
parent_temp_trait_cor = parent_temp_traits %>% 
  ungroup() %>% 
  split(.$lake_morph_Pair_Full_Temp) %>% 
  # ungroup() %>% 
  map(select, vars_keep) %>% 
  map(cor)

parent_temp_graph = parent_temp_trait_cor %>% 
  reshape2::melt() %>% 
  rename(lake_morph_full = L1)

parent_temp_trait_cor_graph = ggplot(parent_temp_graph, 
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

ggsave('Univariate_parent_temp_trait_integration.tiff', 
       plot = noallo_trait_cor_graph, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 40)


# Matrix correlations -----------------------------------------------------

## How do the morph, offspring temp, and parent temp data sets
## correlate back to the original data

##Going to use the RV coefficient

library(MatrixCorrelation)

F2_univariate_traits
F2_noallo_traits
F2_WC_ecotype_traits
F2_offspring_temp_traits
F2_parent_temp_traits

## Ecotype to orginal
RV(F2_univariate_traits, 
   F2_WC_ecotype_traits)
RV(F2_univariate_traits, 
   F2_offspring_temp_traits)
RV(F2_univariate_traits, 
   F2_parent_temp_traits)

RV(F2_noallo_traits, 
   F2_WC_ecotype_traits)
RV(F2_noallo_traits, 
   F2_offspring_temp_traits)
RV(F2_noallo_traits, 
   F2_parent_temp_traits)
