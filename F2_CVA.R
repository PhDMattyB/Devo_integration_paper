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
                                lake_morph = identifiers$Lake_morph,
                                lake_morph_full = identifiers$lake_morph_Pair_Full_Temp)

# F2_off_temp = lm.rrpp(coords ~ offspring_temp * lake_morph,
#                      data = F2_geo_df)
# 


# F2 off temp per lake morph CVA per morph -----------------------------------------------


F2_off_temp = procD.lm(coords ~ offspring_temp * lake_morph, 
                data = F2_geo_df)

summary(F2_off_temp)

prep.lda(F2_off_temp, 
         inherent.groups = TRUE) # see groups available

lda.args = prep.lda(F2_off_temp) 

CVA = do.call(lda, lda.args)
CVA # 3 CVs produced

## Axis scaling
CVA$scaling # will return all CVs

CVA$means
## divide axis 1 and sum of all axes to get var explained
CVA$svd

## CVS scores for each individual
F2_off_temp_cva_scores = predict(CVA)

F2_off_temp_cva_scores = F2_off_temp_cva_scores$x %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(individualID = rowname) %>% 
  arrange(individualID)

F2_off_temp_cva = bind_cols(F2_off_temp_cva_scores, 
          identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)


# F2 off temp per lake morph data viz ----------------------------------------------------

## graphing them all together is hard to see
F2_off_temp_means = F2_off_temp_cva %>% 
  group_by(Ecotype_off_temp) %>% 
  summarize(mean_LD1 = mean(LD1), 
            mean_LD2 = mean(LD2))

ggplot(data = F2_off_temp_cva, 
       aes(x = LD1, 
           y = LD2)) + 
  geom_point(aes(col = Ecotype_off_temp, 
                 shape = Offspring_temp))+
  geom_point(data = F2_off_temp_means, 
             # col = 'Black', 
             size = 4,
             aes(col = Ecotype_off_temp, 
                 x = mean_LD1, 
                 y = mean_LD2))

## Graph for each Lake

F2_off_temp_means = F2_off_temp_cva %>% 
  group_by(Ecotype_off_temp) %>% 
  summarize(mean_LD1 = mean(LD1), 
            mean_LD2 = mean(LD2)) %>% 
  filter(Ecotype_off_temp %in% c('ASHNC_12', 
                                 'ASHNC_18', 
                                 'ASHNW_12', 
                                 'ASHNW_18'))

ASHN_off_temp = F2_off_temp_cva %>% 
  filter(Ecotype_off_temp %in% c('ASHNC_12', 
                                 'ASHNC_18', 
                                 'ASHNW_12', 
                                 'ASHNW_18'))

ggplot(data = ASHN_off_temp, 
       aes(x = LD1, 
           y = LD2)) + 
  geom_point(aes(col = Ecotype_off_temp, 
                 shape = Offspring_temp))+
  geom_point(data = F2_off_temp_means, 
             # col = 'Black', 
             size = 4,
             aes(col = Ecotype_off_temp, 
                 x = mean_LD1, 
                 y = mean_LD2))


F2_off_temp_means = F2_off_temp_cva %>% 
  group_by(Ecotype_off_temp) %>% 
  summarize(mean_LD1 = mean(LD1), 
            mean_LD2 = mean(LD2)) %>% 
  filter(Ecotype_off_temp %in% c('MYVC_12', 
                                 'MYVC_18', 
                                 'MYVW_12', 
                                 'MYVW_18'))

MYV_off_temp = F2_off_temp_cva %>% 
  filter(Ecotype_off_temp %in% c('MYVC_12', 
                                 'MYVC_18', 
                                 'MYVW_12', 
                                 'MYVW_18'))

ggplot(data = MYV_off_temp, 
       aes(x = LD1, 
           y = LD2)) + 
  geom_point(aes(col = Ecotype_off_temp, 
                 shape = Offspring_temp))+
  geom_point(data = F2_off_temp_means, 
             # col = 'Black', 
             size = 4,
             aes(col = Ecotype_off_temp, 
                 x = mean_LD1, 
                 y = mean_LD2))


F2_off_temp_means = F2_off_temp_cva %>% 
  group_by(Ecotype_off_temp) %>% 
  summarize(mean_LD1 = mean(LD1), 
            mean_LD2 = mean(LD2)) %>% 
  filter(Ecotype_off_temp %in% c('SKRC_12', 
                                 'SKRC_18', 
                                 'SKRW_12', 
                                 'SKRW_18'))

SKR_off_temp = F2_off_temp_cva %>% 
  filter(Ecotype_off_temp %in% c('SKRC_12', 
                                 'SKRC_18', 
                                 'SKRW_12', 
                                 'SKRW_18'))

ggplot(data = SKR_off_temp, 
       aes(x = LD1, 
           y = LD2)) + 
  geom_point(aes(col = Ecotype_off_temp, 
                 shape = Offspring_temp))+
  geom_point(data = F2_off_temp_means, 
             # col = 'Black', 
             size = 4,
             aes(col = Ecotype_off_temp, 
                 x = mean_LD1, 
                 y = mean_LD2))


F2_off_temp_means = F2_off_temp_cva %>% 
  group_by(Ecotype_off_temp) %>% 
  summarize(mean_LD1 = mean(LD1), 
            mean_LD2 = mean(LD2)) %>% 
  filter(Ecotype_off_temp %in% c('CSWYC_12', 
                                 'CSWYC_18', 
                                 'GTSW_12', 
                                 'GTSW_18'))

GTS_CSWY_off_temp = F2_off_temp_cva %>% 
  filter(Ecotype_off_temp %in% c('CSWYC_12', 
                                 'CSWYC_18', 
                                 'GTSW_12', 
                                 'GTSW_18'))

ggplot(data = GTS_CSWY_off_temp, 
       aes(x = LD1, 
           y = LD2)) + 
  geom_point(aes(col = Ecotype_off_temp, 
                 shape = Offspring_temp))+
  geom_point(data = F2_off_temp_means, 
             # col = 'Black', 
             size = 4,
             aes(col = Ecotype_off_temp, 
                 x = mean_LD1, 
                 y = mean_LD2))

# F2 off temp CVA temp treatment -----------------------------------------------


F2_off_temp_only = procD.lm(coords ~ offspring_temp, 
                       data = F2_geo_df)

summary(F2_off_temp_only)

prep.lda(F2_off_temp_only, 
         inherent.groups = TRUE) # see groups available

lda.args = prep.lda(F2_off_temp_only) 

CVA = do.call(lda, lda.args)
CVA # 3 CVs produced

## Axis scaling
CVA$scaling # will return all CVs

CVA$means
## divide axis 1 and sum of all axes to get var explained
CVA$svd

## CVS scores for each individual
F2_off_temp_only_cva_scores = predict(CVA)

F2_off_temp_only_cva_scores = F2_off_temp_only_cva_scores$x %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(individualID = rowname) %>% 
  arrange(individualID)

F2_off_temp_only_cva = bind_cols(F2_off_temp_only_cva_scores, 
                            identifiers) %>% 
  unite('Ecotype_off_temp', 
        Lake_morph, 
        Offspring_temp, 
        sep = '_', 
        remove = F)

F2_off_temp_only_means = F2_off_temp_only_cva %>% 
  group_by(Ecotype_off_temp) %>% 
  summarize(mean_LD1 = mean(LD1))

ggplot(data = F2_off_temp_only_cva, 
       aes(x = LD1)) + 
  geom_density(aes(col = Ecotype_off_temp, 
                   fill = Ecotype_off_temp))




# ASHN_off_temp_means = F2_off_temp_only_cva %>% 
#   group_by(Ecotype_off_temp) %>% 
#   summarize(mean_LD1 = mean(LD1)) %>% 
#   filter(Ecotype_off_temp %in% c('ASHNC_12', 
#                                  'ASHNC_18', 
#                                  'ASHNW_12', 
#                                  'ASHNW_18'))

ASHN_off_temp_only = F2_off_temp_only_cva %>% 
  filter(Ecotype_off_temp %in% c('ASHNC_12', 
                                 'ASHNC_18', 
                                 'ASHNW_12', 
                                 'ASHNW_18'))

ggplot(data = ASHN_off_temp_only, 
       aes(x = LD1)) +
  geom_density(aes(col = Ecotype_off_temp, 
                   fill = Ecotype_off_temp))
  

# F2_off_temp_means = F2_off_temp_cva %>% 
#   group_by(Ecotype_off_temp) %>% 
#   summarize(mean_LD1 = mean(LD1), 
#             mean_LD2 = mean(LD2)) %>% 
#   filter(Ecotype_off_temp %in% c('MYVC_12', 
#                                  'MYVC_18', 
#                                  'MYVW_12', 
#                                  'MYVW_18'))

MYV_off_temp_only = F2_off_temp_only_cva %>% 
  filter(Ecotype_off_temp %in% c('MYVC_12', 
                                 'MYVC_18', 
                                 'MYVW_12', 
                                 'MYVW_18'))

ggplot(data = MYV_off_temp_only, 
       aes(x = LD1)) +
  geom_density(aes(col = Ecotype_off_temp, 
                   fill = Ecotype_off_temp))


# F2_off_temp_means = F2_off_temp_cva %>% 
#   group_by(Ecotype_off_temp) %>% 
#   summarize(mean_LD1 = mean(LD1), 
#             mean_LD2 = mean(LD2)) %>% 
#   filter(Ecotype_off_temp %in% c('SKRC_12', 
#                                  'SKRC_18', 
#                                  'SKRW_12', 
#                                  'SKRW_18'))

SKR_off_temp_only = F2_off_temp_only_cva %>% 
  filter(Ecotype_off_temp %in% c('SKRC_12', 
                                 'SKRC_18', 
                                 'SKRW_12', 
                                 'SKRW_18'))

ggplot(data = SKR_off_temp_only, 
       aes(x = LD1)) +
  geom_density(aes(col = Ecotype_off_temp, 
                   fill = Ecotype_off_temp))



# F2_off_temp_means = F2_off_temp_cva %>% 
#   group_by(Ecotype_off_temp) %>% 
#   summarize(mean_LD1 = mean(LD1), 
#             mean_LD2 = mean(LD2)) %>% 
#   filter(Ecotype_off_temp %in% c('CSWYC_12', 
#                                  'CSWYC_18', 
#                                  'GTSW_12', 
#                                  'GTSW_18'))

GTS_CSWY_off_temp_only = F2_off_temp_only_cva %>% 
  filter(Ecotype_off_temp %in% c('CSWYC_12', 
                                 'CSWYC_18', 
                                 'GTSW_12', 
                                 'GTSW_18'))

ggplot(data = GTS_CSWY_off_temp_only, 
       aes(x = LD1)) +
  geom_density(aes(col = Ecotype_off_temp, 
                   fill = Ecotype_off_temp))

##
