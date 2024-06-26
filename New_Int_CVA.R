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

# multi_shape_1011 = readland.tps('multi_shape_1011.TPS', 
#                                specID = 'imageID')
# multi_shape_1011 = gpagen(multi_shape_1011)
# F2_multi_1011 = geomorph.data.frame(coords = two.d.array(multi_shape_1011$coords), 
#                                    Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
#                                    parent_temp = identifiers$Parent_temp, 
#                                    offspring_temp = identifiers$Offspring_temp,
#                                    grand_temp = identifiers$Grand_temp,
#                                    morph = identifiers$Morph, 
#                                    population = identifiers$Lake,
#                                    lake_morph = identifiers$Lake_morph,
#                                    lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)
# 
# multi_shape_12 = readland.tps('multi_shape_12.TPS', 
#                                specID = 'imageID')
# multi_shape_12 = gpagen(multi_shape_12)
# F2_multi_12 = geomorph.data.frame(coords = two.d.array(multi_shape_12$coords), 
#                                    Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
#                                    parent_temp = identifiers$Parent_temp, 
#                                    offspring_temp = identifiers$Offspring_temp,
#                                    grand_temp = identifiers$Grand_temp,
#                                    morph = identifiers$Morph, 
#                                    population = identifiers$Lake,
#                                    lake_morph = identifiers$Lake_morph,
#                                    lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)
# 

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

# operculum shape 789 -----------------------------------------------------

multi_shape789_mod = procD.lm(coords ~ offspring_temp * lake_morph, 
                              data = F2_multi_789)
summary(multi_shape789_mod)
prep.lda(multi_shape789_mod, 
         inherent.groups = TRUE) # see groups available

lda.args = prep.lda(multi_shape789_mod) 
multi_shape789_CVA = do.call(lda, lda.args)

## CVS scores for each individual
multi_shape789_cva_scores = predict(multi_shape789_CVA)

multi_shape789_cva_scores = multi_shape789_cva_scores$x %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(individualID = rowname) %>% 
  arrange(individualID)

multi_shape789_cva = bind_cols(multi_shape789_cva_scores, 
                               identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)

multi_shape789_cva %>% 
  write_csv('F2_operculum_shape_cva_per_pop.csv')


# body shape ------------------------------------------------------

multi_body_mod = procD.lm(coords ~ offspring_temp * lake_morph, 
                              data = F2_body_geo_df)
summary(multi_body_mod)
prep.lda(multi_body_mod, 
         inherent.groups = TRUE) # see groups available

lda.args = prep.lda(multi_body_mod) 
multi_body_CVA = do.call(lda, lda.args)

## CVS scores for each individual
multi_body_cva_scores = predict(multi_body_CVA)

multi_body_cva_scores = multi_body_cva_scores$x %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(individualID = rowname) %>% 
  arrange(individualID)

multi_body_cva = bind_cols(multi_body_cva_scores, 
                               identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)

multi_body_cva %>% 
  write_csv('F2_body_shape_cva_per_pop.csv')


# craniofacial shape ------------------------------------------------------


multi_cranio_mod = procD.lm(coords ~ offspring_temp * lake_morph, 
                          data = F2_cranio_geo_df)
summary(multi_cranio_mod)
prep.lda(multi_cranio_mod, 
         inherent.groups = TRUE) # see groups available

lda.args = prep.lda(multi_cranio_mod) 
multi_cranio_CVA = do.call(lda, lda.args)

## CVS scores for each individual
multi_cranio_cva_scores = predict(multi_cranio_CVA)

multi_cranio_cva_scores = multi_cranio_cva_scores$x %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(individualID = rowname) %>% 
  arrange(individualID)

multi_cranio_cva = bind_cols(multi_cranio_cva_scores, 
                           identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)

multi_cranio_cva %>% 
  write_csv('F2_cranio_shape_cva_per_pop.csv')


# 4bar shape lda ----------------------------------------------------------


multi_4bar_mod = procD.lm(coords ~ offspring_temp * lake_morph, 
                            data = F2_4bar_geo_df)
summary(multi_4bar_mod)
prep.lda(multi_4bar_mod, 
         inherent.groups = TRUE) # see groups available

lda.args = prep.lda(multi_4bar_mod) 
multi_4bar_CVA = do.call(lda, lda.args)

## CVS scores for each individual
multi_4bar_cva_scores = predict(multi_4bar_CVA)

multi_4bar_cva_scores = multi_4bar_cva_scores$x %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(individualID = rowname) %>% 
  arrange(individualID)

multi_4bar_cva = bind_cols(multi_4bar_cva_scores, 
                             identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)

multi_4bar_cva %>% 
  write_csv('F2_4bar_shape_cva_per_pop.csv')




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


# univariate jaw length ---------------------------------------------------

jaw_length_mod = procD.lm(jaw_length ~ Offspring_temp * Lake_morph, 
                          data = F2_univariate_traits)

plot(jaw_length_mod, type = "diagnostics") 

summary(jaw_length_mod)

jaw_length_fitted = jaw_length_mod$fitted

# jaw_length_mod$residuals
jaw_length_mod$LM$fitted

jaw_length_fitted %>% 
  as.data.frame() %>% 
  bind_cols(identifiers)%>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F) %>% 
  rename(jaw_length = V1) %>% 
  write_csv('F2_jaw_length_fitted_per_pop.csv')


# univariate fbANOVA# univariate fbar 23-24 ---------------------------------------------------

fbar2324_mod = procD.lm(fbar_23_24 ~ Offspring_temp * Lake_morph, 
                          data = F2_univariate_traits)

summary(fbar2324_mod)

fbar2324_fitted = fbar2324_mod$fitted

fbar2324_fitted %>% 
  as.data.frame() %>% 
  bind_cols(identifiers)%>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F) %>% 
  rename(fbar_2324 = V1) %>% 
  write_csv('F2_fbar2324_fitted_per_pop.csv')

# univariate fbar 8-24 ---------------------------------------------------

fbar824_mod = procD.lm(fbar_8_24 ~ Offspring_temp * Lake_morph, 
                        data = F2_univariate_traits)

summary(fbar824_mod)

fbar824_fitted = fbar824_mod$fitted

fbar824_fitted %>% 
  as.data.frame() %>% 
  bind_cols(identifiers)%>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F) %>% 
  rename(fbar_824 = V1) %>% 
  write_csv('F2_fbar824_fitted_per_pop.csv')

# univariate fbar 8-27 ---------------------------------------------------

fbar827_mod = procD.lm(fbar_8_27 ~ Offspring_temp * Lake_morph, 
                        data = F2_univariate_traits)

summary(fbar827_mod)

fbar827_fitted = fbar827_mod$fitted

fbar827_fitted %>% 
  as.data.frame() %>% 
  bind_cols(identifiers)%>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F) %>% 
  rename(fbar_827 = V1) %>% 
  write_csv('F2_fbar827_fitted_per_pop.csv')


# univariate fbar 23-27 ---------------------------------------------------

fbar2327_mod = procD.lm(fbar_23_27 ~ Offspring_temp * Lake_morph, 
                       data = F2_univariate_traits)

summary(fbar2327_mod)

fbar2327_fitted = fbar2327_mod$fitted

fbar2327_fitted %>% 
  as.data.frame() %>% 
  bind_cols(identifiers)%>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F) %>% 
  rename(fbar_2327 = V1) %>% 
  write_csv('F2_fbar2327_fitted_per_pop.csv')

# univariate fbar 25-26 ---------------------------------------------------

fbar2526_mod = procD.lm(fbar_25_26 ~ Offspring_temp * Lake_morph, 
                        data = F2_univariate_traits)

summary(fbar2526_mod)

fbar2526_fitted = fbar2526_mod$fitted

fbar2526_fitted %>% 
  as.data.frame() %>% 
  bind_cols(identifiers)%>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F) %>% 
  rename(fbar_2526 = V1) %>% 
  write_csv('F2_fbar2526_fitted_per_pop.csv')


# univariate body depth ---------------------------------------------------

body_depth_mod = procD.lm(body_width ~ Offspring_temp * Lake_morph, 
                        data = F2_univariate_traits)

summary(body_depth_mod)

body_depth_fitted = body_depth_mod$fitted

body_depth_fitted %>% 
  as.data.frame() %>% 
  bind_cols(identifiers)%>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F) %>% 
  rename(body_depth = V1) %>% 
  write_csv('F2_body_depth_fitted_per_pop.csv')

# univariate caudal depth 1 ---------------------------------------------------

caudal_depth1_mod = procD.lm(caudal1_14_18 ~ Offspring_temp * Lake_morph, 
                          data = F2_univariate_traits)

summary(caudal_depth1_mod)

caudal_depth1_fitted = caudal_depth1_mod$fitted

caudal_depth1_fitted %>% 
  as.data.frame() %>% 
  bind_cols(identifiers)%>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F) %>% 
  rename(caudal_depth1 = V1) %>% 
  write_csv('F2_caudal_depth1_fitted_per_pop.csv')
# univariate caudal depth 2 ---------------------------------------------------

caudal_depth2_mod = procD.lm(caudal2_15_17 ~ Offspring_temp * Lake_morph, 
                             data = F2_univariate_traits)

summary(caudal_depth2_mod)

caudal_depth2_fitted = caudal_depth2_mod$fitted

caudal_depth2_fitted %>% 
  as.data.frame() %>% 
  bind_cols(identifiers)%>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F) %>% 
  rename(caudal_depth2 = V1) %>% 
  write_csv('F2_caudal_depth2_fitted_per_pop.csv')


# univariate body length ---------------------------------------------------

body_length_mod = procD.lm(body_length ~ Offspring_temp * Lake_morph, 
                             data = F2_univariate_traits)

summary(body_length_mod)

body_length_fitted = body_length_mod$fitted

body_length_fitted %>% 
  as.data.frame() %>% 
  bind_cols(identifiers)%>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F) %>% 
  rename(body_length = V1) %>% 
  write_csv('F2_body_length_fitted_per_pop.csv')

# ASHN compare integration between lm sets -------------------------------------

## multivariate traits
cranio_ashn = read_csv('F2_cranio_shape_cva_per_pop.csv') %>% 
  # filter(Lake == 'ASHN')
  filter(Lake_morph == 'ASHNC')
body_shape_ashn = read_csv('F2_body_shape_cva_per_pop.csv') %>% 
  # filter(Lake == 'ASHN')
  filter(Lake_morph == 'ASHNC')
fbar_ashn = read_csv('F2_4bar_shape_cva_per_pop.csv') %>% 
  # filter(Lake == 'ASHN')
  filter(Lake_morph == 'ASHNC')
eye_shape_ashn = read_csv('F2_eye_shape_cva_per_pop.csv') %>% 
  # filter(Lake == 'ASHN')
  filter(Lake_morph == 'ASHNC')
operculum_shape_ashn = read_csv('F2_operculum_shape_cva_per_pop.csv') %>% 
  # filter(Lake == 'ASHN')
  filter(Lake_morph == 'ASHNC')


##
# ASHN ld axes data viz ---------------------------------------------------

ASHN_means = operculum_shape_ashn %>% 
  group_by(Ecotype_off_temp) %>% 
  summarize(mean_LD1 = mean(LD1), 
            mean_LD2 = mean(LD2)) %>% 
  filter(Ecotype_off_temp %in% c('ASHNC_12', 
                                 'ASHNC_18', 
                                 'ASHNW_12', 
                                 'ASHNW_18'))
ASHN_data = operculum_shape_ashn %>% 
  filter(Ecotype_off_temp %in% c('ASHNC_12', 
                                 'ASHNC_18', 
                                 'ASHNW_12', 
                                 'ASHNW_18'))
ASHN_data$Offspring_temp = as.factor(ASHN_data$Offspring_temp)

ggplot(data = ASHN_data, 
       aes(x = LD1, 
           y = LD2)) + 
  geom_point(aes(col = Ecotype_off_temp, 
                 shape = Offspring_temp))+
  geom_point(data = ASHN_means, 
             # col = 'Black', 
             size = 4,
             aes(col = Ecotype_off_temp, 
                 x = mean_LD1, 
                 y = mean_LD2))






# ASHNC integration -------------------------------------------------------

ASHNC_lm_integration_ecotype = bind_cols(cranio_ashn$LD1, 
                                        body_shape_ashn$LD1, 
                                        fbar_ashn$LD2, 
                                        eye_shape_ashn$LD1, 
                                        operculum_shape_ashn$LD2) %>% 
  rename(cranio_shape = 1, 
         body_shape = 2, 
         fbar_shape = 3, 
         eye_shape = 4, 
         operculum_shape = 5)

## integration of multivariate traits
ASHNC_multivariate = ASHNC_lm_integration_ecotype %>% 
  dplyr::select(cranio_shape, 
         body_shape, 
         fbar_shape, 
         eye_shape, 
         operculum_shape)

ASHNC_multi_int_ecotype_pcor = pcor(x = ASHNC_multivariate, 
                                     method = 'pearson')

ASHN_multi_int_ecotype_pcor$p.value
ASHN_multi_int_ecotype_pcor$estimate

ASHN_multi_int_temp = bind_cols(cranio_ashn$LD2, 
          body_shape_ashn$LD2, 
          fbar_ashn$LD1, 
          eye_shape_ashn$LD2, 
          operculum_shape_ashn$LD1) %>% 
  rename(cranio_shape = 1, 
         body_shape = 2, 
         fbar_shape = 3, 
         eye_shape = 4, 
         operculum_shape = 5)

ASHN_multi_int_temp = pcor(x = ASHN_multi_int_temp, 
                                    method = 'pearson')

ASHN_multi_int_temp$p.value
ASHN_multi_int_temp$estimate

