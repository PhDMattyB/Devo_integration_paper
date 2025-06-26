##############################
## Integration pattern cor
##
## Matt Brachmann (PhDMattyB)
##
## 10.01.2025
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')


library(geomorph)
library(RRPP)
library(MASS)
library(ppcor)
library(igraph)
library(reshape2)
library(candisc)
library(tidyverse)
library(lineup)

theme_set(theme_bw())

wc_cols = c('#003566', 
            '#ff006e')

# Wild data ---------------------------------------------------------------

wild_univariate_traits = read_csv("Wild_univar_traits_nokinetics.csv")
wild_kinetics = read_csv('Wild_Jaw_kinetic_traits.csv') %>% 
  dplyr::select(-PreMax_Rotation, 
                -Opercular_Rotation)

lake_morph = wild_univariate_traits %>% 
  dplyr::select(Lake_morph)

wild_uni_traits = wild_univariate_traits %>%
  as_tibble() %>%
  group_by(Lake_morph) %>%
  dplyr::select(jaw_length:ratio2)

wild_uni_traits = bind_cols(wild_uni_traits, 
                            wild_kinetics) 
wild_traits_scaled = wild_uni_traits %>% 
  ungroup() %>% 
  dplyr::select(-Lake_morph) %>% 
  scale(., center = T, scale = T) %>% 
  as_tibble() %>% 
  bind_cols(lake_morph, 
            .) %>% 
  rename(CMA = ratio1, 
         OMA = ratio2) %>% 
  dplyr::select('Lake_morph',
                'jaw_length', 
                'head_depth', 
                'Opercular_KT', 
                'PreMax_KT', 
                'CMA',
                'OMA', 
                'jaw_2_6', 
                'fbar_23_24', 
                'fbar_8_24', 
                'fbar_8_27', 
                'fbar_23_27', 
                'fbar_25_26', 
                'max_27_3', 
                'max_3_28', 
                'max_28_27', 
                'body_length', 
                'body_width', 
                'lm_6_12', 
                'lm_12_13', 
                'lm_13_14',
                'lm_14_15', 
                'lm_6_21', 
                'lm_20_21', 
                'lm_21_13', 
                'lm_20_13', 
                'lm_12_19', 
                'lm_13_19', 
                'lm_19_18', 
                'lm_18_17', 
                'lm_1_23', 
                'lm_23_2', 
                'caudal1_14_18', 
                'caudal2_15_17')

# vars_keep = names(wild_traits_scaled)[c(2,3,4,5,6,7,8,9,10,11,
#                                         12,13,14,15,16,17,18,
#                                         19,20,21,22,23,24,25,26,
#                                         27,28,29, 30, 31, 32,
#                                         33, 34)]
# wild_uni_trait_cor = wild_traits_scaled %>%
#   ungroup() %>%
#   # split(.$lake_morph_Pair_Full_Temp) %>%
#   split(.$Lake_morph) %>% 
#   # ungroup() %>%
#   map(dplyr::select, vars_keep) %>%
#   map(cor)
# 
# ASHN_wild_cor = corbetw2mat(wild_uni_trait_cor$ASHNC, 
#                             wild_uni_trait_cor$ASHNW, 
#                             what = 'all', 
#                             corthresh = 0.5)
# 

ASHN_cor_traits = wild_traits_scaled %>% 
  # dplyr::select(Lake_morph, 
  #               Opercular_KT, 
  #               PreMax_KT, 
  #               CMA, 
  #               OMA) %>% 
  filter(Lake_morph %in% c('ASHNC', 
                           'ASHNW'))

# MYV_cor_traits = wild_traits_scaled %>% 
#   dplyr::select(Lake_morph, 
#                 Opercular_KT, 
#                 PreMax_KT, 
#                 CMA, 
#                 OMA) %>% 
#   filter(Lake_morph %in% c('MYVC', 
#                            'MYVW'))


# positive correlation wild ----------------------------------------


positive_correlation_wild = ggplot(data = ASHN_cor_traits)+
  geom_point(aes(x = Opercular_KT, 
                 y = OMA, 
                 group = Lake_morph, 
                 col = Lake_morph), 
             size = 2)+
  geom_smooth(aes(x = Opercular_KT, 
                  y = OMA, 
                  group = Lake_morph, 
                  col = Lake_morph), 
              method = 'lm', 
              position = 'identity', 
              se = F)+
  labs(x = 'Opercular kinetics', 
       y = 'OMA', 
       title = 'Positive correlation')+
  scale_color_manual(values = wc_cols)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        title = element_text(size = 14), 
        legend.position = 'none')


ggplot(data = ASHN_cor_traits)+
  geom_point(aes(x = lm_12_13, 
                 y = lm_21_13, 
                 group = Lake_morph, 
                 col = Lake_morph), 
             size = 2)+
  geom_smooth(aes(x = lm_12_13, 
                  y = lm_21_13, 
                  group = Lake_morph, 
                  col = Lake_morph), 
              method = 'lm', 
              position = 'identity', 
              se = F)+
  labs(x = 'lm 12 13',
       y = 'lm 21 13',
       title = 'Positive correlation')+
  scale_color_manual(values = wc_cols)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        title = element_text(size = 14), 
        legend.position = 'none')

# Negative correlation wild ----------------------------------------


negative_correlation_wild = ggplot(data = ASHN_cor_traits)+
  geom_point(aes(x = Opercular_KT, 
                 y = PreMax_KT, 
                 group = Lake_morph, 
                 col = Lake_morph), 
             size = 2)+
  geom_smooth(aes(x = Opercular_KT, 
                  y = PreMax_KT, 
                  group = Lake_morph, 
                  col = Lake_morph), 
              method = 'lm', 
              position = 'identity', 
              se = F)+
  labs(x = 'Opercular kinetics', 
       y = 'Pre-maxilla kinetics', 
       title = 'Negative correlation')+
  scale_color_manual(values = wc_cols)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        title = element_text(size = 14), 
        legend.position = 'none')


# Neutral correlation wild -----------------------------------------------------

neutral_correlation_wild = ggplot(data = ASHN_cor_traits)+
  geom_point(aes(x = Opercular_KT, 
                 y = CMA, 
                 group = Lake_morph, 
                 col = Lake_morph), 
             size = 2)+
  geom_smooth(aes(x = Opercular_KT, 
                  y = CMA, 
                  group = Lake_morph, 
                  col = Lake_morph), 
              method = 'lm', 
              position = 'identity', 
              se = F)+
  labs(x = 'Opercular kinetics', 
       y = 'CMA', 
       title = 'Neutral correlation')+
  scale_color_manual(values = wc_cols)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        title = element_text(size = 14), 
        legend.position = 'none')



# F2_data -----------------------------------------------------------------

F2_traits_scaled = F2_orig_traits %>% 
  ungroup() %>% 
  dplyr::select(-Lake_morph) %>% 
  scale(., 
        center = T, 
        scale = T) %>% 
  as_tibble() %>% 
  bind_cols(lake_morph, 
            .) %>% 
  dplyr::rename(OMA = ratio1, 
                CMA = ratio2)%>% 
  dplyr::select('Lake_morph',
                'jaw_length', 
                'head_depth', 
                'Opercular_KT', 
                'PreMax_KT', 
                'CMA',
                'OMA', 
                'jaw_2_6', 
                'fbar_23_24', 
                'fbar_8_24', 
                'fbar_8_27', 
                'fbar_23_27', 
                'fbar_25_26', 
                'max_27_3', 
                'max_3_28', 
                'max_28_27', 
                'body_length', 
                'body_width', 
                'lm_6_12', 
                'lm_12_13', 
                'lm_13_14',
                'lm_14_15', 
                'lm_6_21', 
                'lm_20_21', 
                'lm_21_13', 
                'lm_20_13', 
                'lm_12_19', 
                'lm_13_19', 
                'lm_19_18', 
                'lm_18_17', 
                'lm_1_23', 
                'lm_23_2', 
                'caudal1_14_18', 
                'caudal2_15_17')

ASHN_F2 = F2_traits_scaled %>% 
  filter(Lake_morph %in% c('ASHNC', 
                           'ASHNW'))
# F2 - what pink color means ----------------------------------------

test = lm(PreMax_KT ~ max_28_27*Lake_morph, 
   data = ASHN_F2)
summary(test)

ggplot(data = ASHN_F2)+
  geom_point(aes(x = PreMax_KT, 
                 y = max_28_27, 
                 group = Lake_morph, 
                 col = Lake_morph), 
             size = 2)+
  geom_smooth(aes(x = PreMax_KT, 
                  y = max_28_27,
                  group = Lake_morph, 
                  col = Lake_morph), 
              method = 'lm', 
              position = 'identity', 
              se = F)+
  # labs(x = 'OMA',
  #      y = 'Opercular kinetics',
  #      title = 'Positive correlation')+
  scale_color_manual(values = wc_cols)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        title = element_text(size = 14), 
        legend.position = 'none')

ASHN_orig_cor %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>%
  # select(2:34) %>% 
  filter_all(any_vars(. < 0.01)) %>% 
  View()

# F2 - what light blue color means ----------------------------------------
ASHN_orig_cor %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>%
  # select(2:34) %>% 
  filter_all(any_vars(. < 0.01)) %>% 
  View()

test = lm(jaw_length ~ max_27_3*Lake_morph, 
          data = ASHN_F2)
summary(test)

ggplot(data = ASHN_F2)+
  geom_point(aes(x = jaw_length, 
                 y = max_27_3, 
                 group = Lake_morph, 
                 col = Lake_morph), 
             size = 2)+
  geom_smooth(aes(x = jaw_length, 
                  y = max_27_3,
                  group = Lake_morph, 
                  col = Lake_morph), 
              method = 'lm', 
              position = 'identity', 
              se = F)+
  # labs(x = 'OMA',
  #      y = 'Opercular kinetics',
  #      title = 'Positive correlation')+
  scale_color_manual(values = wc_cols)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        title = element_text(size = 14), 
        legend.position = 'none')



# F2 - what Dark blue color means ----------------------------------------
ASHN_orig_cor %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  as_tibble() %>%
  # select(2:34) %>% 
  filter_all(any_vars(. < 0.01)) %>% 
  View()

test = lm(PreMax_KT ~ OMA*Lake_morph, 
          data = ASHN_F2)
summary(test)

ggplot(data = ASHN_F2)+
  geom_point(aes(x = PreMax_KT, 
                 y = OMA, 
                 group = Lake_morph, 
                 col = Lake_morph), 
             size = 2)+
  geom_smooth(aes(x = PreMax_KT, 
                  y = OMA,
                  group = Lake_morph, 
                  col = Lake_morph), 
              method = 'lm', 
              position = 'identity', 
              se = F)+
  # labs(x = 'OMA',
  #      y = 'Opercular kinetics',
  #      title = 'Positive correlation')+
  scale_color_manual(values = wc_cols)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        title = element_text(size = 14), 
        legend.position = 'none')


test = lm(fbar_23_27 ~ lm_20_21*Lake_morph, 
          data = ASHN_F2)
summary(test)

ggplot(data = ASHN_F2)+
  geom_point(aes(x = fbar_23_27, 
                 y = lm_20_21, 
                 group = Lake_morph, 
                 col = Lake_morph), 
             size = 2)+
  geom_smooth(aes(x = fbar_23_27, 
                  y = lm_20_21,
                  group = Lake_morph, 
                  col = Lake_morph), 
              method = 'lm', 
              position = 'identity', 
              se = F)+
  # labs(x = 'OMA',
  #      y = 'Opercular kinetics',
  #      title = 'Positive correlation')+
  scale_color_manual(values = wc_cols)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        title = element_text(size = 14), 
        legend.position = 'none')

test = lm(OMA ~ fbar_8_27*Lake_morph, 
          data = ASHN_F2)
summary(test)

ggplot(data = ASHN_F2)+
  geom_point(aes(x = OMA, 
                 y = fbar_8_27, 
                 group = Lake_morph, 
                 col = Lake_morph), 
             size = 2)+
  geom_smooth(aes(x = OMA, 
                  y = fbar_8_27,
                  group = Lake_morph, 
                  col = Lake_morph), 
              method = 'lm', 
              position = 'identity', 
              se = F)+
  # labs(x = 'OMA',
  #      y = 'Opercular kinetics',
  #      title = 'Positive correlation')+
  scale_color_manual(values = wc_cols)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        title = element_text(size = 14), 
        legend.position = 'none')
