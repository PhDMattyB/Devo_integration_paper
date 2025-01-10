##############################
## Pattern integration - wild fish
##
## Matt Brachmann (PhDMattyB)
##
## 06.07.2024
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


# read multiple tps files -------------------------------------------------

wild_f2_data = readmulti.tps(c('Wild_28LM_Final.tps', 
                       'F2_No_GT_ALL_LMs.TPS'), 
                       specID = 'imageID')

# writeland.tps(A = wild_f2_data, 
#               file = 'wild_F2_combo.tps')

## common gpa
common_gpa = gpagen(wild_f2_data, 
                    print.progress = T)

id = read_csv('wild_F2_id.csv')
common_df = geomorph.data.frame(coords = two.d.array(common_gpa$coords),
                                split = id$generation, 
                                id = id$rowname)

coord_sub = coords.subset(common_gpa$coords, 
                                   id$generation)

wild_lmk_coords = coord_sub$wild
F2_lmk_coords = coord_sub$F2
 

# wild fish data  -----------------------------------------------------------

wild_identifiers = read_csv('TPS_Wild_metadata.csv') 

# wild_tps = readland.tps('Wild_Final.TPS', 
#                         specID = 'imageID')
# 
# ## superimposition on the entire dataset
# wild_gpa = gpagen(wild_tps, 
#                   print.progress = F)
# 

# wild disparity ----------------------------------------------------------
# wild_df = geomorph.data.frame(coords = two.d.array(wild_gpa$coords), 
#                                  Lake = wild_identifiers$Lake, 
#                                  Morph = wild_identifiers$Morph, 
#                                  Ecotypes = wild_identifiers$Lake_morph)
# 


# wild univariate traits --------------------------------------------------

lmks = data.frame(jaw_length = c(1, 2), 
                  fbar_23_24 = c(23, 24), 
                  fbar_8_24 = c(8, 24), 
                  fbar_8_27 = c(8, 27), 
                  fbar_23_27 = c(23, 27), 
                  fbar_25_26 = c(25, 26),
                  max_27_3 = c(27, 3), 
                  max_3_28 = c(3, 28), 
                  max_28_27 = c(28, 17),
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
# common_traits = interlmkdist(common_coords, 
#                              lmks) %>%
#   as.data.frame() %>% 
#   rownames_to_column() %>% 
#   as_tibble() %>% 
#   mutate(ratio1 = lm_1_23/fbar_23_27, 
#          ratio2 = lm_1_23/lm_23_2) %>% 
#   dplyr::select(rowname, 
#          jaw_length:lm_23_2, 
#          ratio1:ratio2, 
#          everything()) 
# 
# 
#   write_csv('univariate_traits_common_gpa_wild_F2.csv')

# wild_coords = wild_gpa$coords
# A = F2_whole_body_gpa$coords
wild_univariate_traits = interlmkdist(wild_lmk_coords, 
                                    lmks)

# arrayspecs(F2_univariate_traits, 
#            4, 
#            3)

wild_univariate_traits = wild_univariate_traits %>% 
  as.data.frame() %>% 
  rownames_to_column() 
# %>% 
#   arrange(rowname)

wild_univariate_traits = bind_cols(wild_univariate_traits, 
          wild_identifiers)%>% 
  mutate(ratio1 = lm_1_23/fbar_23_27, 
         ratio2 = lm_1_23/lm_23_2) %>% 
  dplyr::select(rowname, 
         jaw_length:lm_23_2, 
         ratio1:ratio2, 
         everything())

# wild_univariate_traits %>%
#   write_csv('Wild_univar_traits_nokinetics.csv')

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
  rename(OMA = ratio1, 
         CMA = ratio2) %>% 
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

# vars_keep = names(wild_uni_traits)[c(2,3,4,5,6,7,8,9,10,11, 
#                                      12,13,14,15,16,17,18, 
#                                      19,20,21,22,23,24,25,26,
#                                      27,28,29)]
vars_keep = names(wild_traits_scaled)[c(2,3,4,5,6,7,8,9,10,11, 
                                     12,13,14,15,16,17,18, 
                                     19,20,21,22,23,24,25,26,
                                     27,28,29, 30, 31, 32, 
                                     33, 34)]
wild_uni_trait_cor = wild_traits_scaled %>%
  ungroup() %>%
  # split(.$lake_morph_Pair_Full_Temp) %>%
  split(.$Lake_morph) %>% 
  # ungroup() %>%
  map(dplyr::select, vars_keep) %>%
  map(cor)

wild_uni_graph = wild_uni_trait_cor %>%
  reshape2::melt() %>%
  rename(lake_morph = L1)

wild_uni_trait_cor_graph = ggplot(wild_uni_graph,
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

ggsave('Wild_Univariate__trait_ecotype_integration.tiff',
       plot = wild_uni_trait_cor_graph,
       dpi = 'retina',
       units = 'cm',
       width = 35,
       height = 20)



# ASHN trait correlations -------------------------------------------------

ASHN_wild_cor = corbetw2mat(wild_uni_trait_cor$ASHNC, 
                               wild_uni_trait_cor$ASHNW, 
                               what = 'all', 
                               corthresh = 0.5)

ASHN_wild_cor = ASHN_wild_cor %>% 
  reshape2::melt()

ASHN_wild_cor_graph = ggplot(ASHN_wild_cor,
                       aes(x = Var1,
                           y = Var2,
                           fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#219ebc",
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
        plot.title = element_text(size = 22),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())


# MYV trait correlations -------------------------------------------------

MYV_wild_cor = corbetw2mat(wild_uni_trait_cor$MYVC, 
                            wild_uni_trait_cor$MYVW, 
                            what = 'all', 
                            corthresh = 0.5)

MYV_wild_cor = MYV_wild_cor %>%
  reshape2::melt() 

MYV_wild_cor_graph = ggplot(MYV_wild_cor,
                             aes(x = Var1,
                                 y = Var2,
                                 fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#219ebc",
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
        plot.title = element_text(size = 22),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())


# SKR trait correlations -------------------------------------------------

SKR_wild_cor = corbetw2mat(wild_uni_trait_cor$SKRC, 
                           wild_uni_trait_cor$SKRW, 
                           what = 'all', 
                           corthresh = 0.5)

SKR_wild_cor = SKR_wild_cor %>%
  reshape2::melt() 

SKR_wild_cor_graph = ggplot(SKR_wild_cor,
                            aes(x = Var1,
                                y = Var2,
                                fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#219ebc",
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
        plot.title = element_text(size = 22),
        # axis.text.x = element_text(angle = 90,
        #                            vjust = 0.5,
        #                            hjust=1),
        legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())
# GTS_CSWY trait correlations -------------------------------------------------

GTS_CSWY_wild_cor = corbetw2mat(wild_uni_trait_cor$CSWY, 
                                wild_uni_trait_cor$GTS, 
                                what = 'all', 
                                corthresh = 0.5)

GTS_CSWY_wild_cor = GTS_CSWY_wild_cor %>%
  reshape2::melt() 

GTS_CSWY_wild_cor_graph = ggplot(GTS_CSWY_wild_cor,
                                 aes(x = Var1,
                                     y = Var2,
                                     fill = value))+
  geom_tile(col = 'black')+
  # scale_fill_gradient2(low = "#003049",
  #                      mid = "#eae2b7",
  #                      high = "#d62828") +
  scale_fill_gradient2(low = "#219ebc",
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
        plot.title = element_text(size = 22),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1),
        legend.position = 'none') 
        # axis.text.x = element_blank(), 
        # axis.ticks.x = element_blank())


##
# RKLT trait correlations -------------------------------------------------

RKLT_wild_cor = corbetw2mat(wild_uni_trait_cor$RKLTC, 
                           wild_uni_trait_cor$RKLTW, 
                           what = 'all', 
                           corthresh = 0.5)

RKLT_wild_cor = RKLT_wild_cor %>%
  reshape2::melt() 

RKLT_wild_cor_graph = ggplot(RKLT_wild_cor,
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
  labs(title = 'A) RKLT')+
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

# STN trait correlations -------------------------------------------------

STN_wild_cor = corbetw2mat(wild_uni_trait_cor$STNC, 
                           wild_uni_trait_cor$STNW, 
                           what = 'all', 
                           corthresh = 0.5)

STN_wild_cor = STN_wild_cor %>%
  reshape2::melt() 

STN_wild_cor_graph = ggplot(STN_wild_cor,
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
  labs(title = 'A) STN')+
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





# Combine correlation graphs ----------------------------------------------
wild_cor_graphs = ASHN_wild_cor_graph/MYV_wild_cor_graph/SKR_wild_cor_graph/GTS_CSWY_wild_cor_graph

ASHN_combo_graphs = ASHN_wild_cor_graph|ASHN_orig_cor|ASHN_F2_plast| ASHN_F1_plast
MYV_combo_graphs = MYV_wild_cor_graph|MYV_orig_cor_plot|MYV_F2_plasticity|MYV_F1_plasticity
SKR_combo_graphs = SKR_wild_cor_graph|SKR_orig_cor_plot|SKR_F2_plasticity|SKR_F1_plasticity
GTSCSWY_combo_graphs = GTS_CSWY_wild_cor_graph|GTSCSWY_orig_cor_plot|GTSCSWY_F2_plasticity|GTSCSWY_F1_plasticity

Big_graph = ASHN_combo_graphs/MYV_combo_graphs/SKR_combo_graphs/GTSCSWY_combo_graphs

# Big_graph2 = Big_graph+
#   theme(panel.background = element_rect(fill='transparent'), 
#   plot.background = element_rect(fill = 'transparent', 
#                                        color = NA))

# ggsave <- function(..., bg = 'transparent') ggplot2::ggsave(..., bg = bg)

ggsave('transback_08.01.2025_SCALED_Figure1_Effects_wild_F1_F2_on_Integration_version2.tiff',
       plot = Big_graph,
       bg = 'transparent',
       dpi = 'retina',
       units = 'cm',
       width = 60,
       height = 60)
##
