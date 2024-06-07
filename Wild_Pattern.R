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
library(tidyverse)
library(reshape2)
library(candisc)



# wild fish data  -----------------------------------------------------------

wild_identifiers = read_csv('TPS_Wild_metadata.csv') 

wild_tps = readland.tps('Wild_Final.TPS', 
                        specID = 'imageID')

## superimposition on the entire dataset
wild_gpa = gpagen(wild_tps, 
                  print.progress = F)



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

wild_coords = wild_gpa$coords
# A = F2_whole_body_gpa$coords
wild_univariate_traits = interlmkdist(wild_coords, 
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
  select(rowname, 
         jaw_length:lm_23_2, 
         ratio1:ratio2, 
         everything())


wild_uni_traits = wild_univariate_traits %>%
  as_tibble() %>%
  group_by(Lake_morph) %>%
  select(jaw_length:ratio2)

vars_keep = names(wild_uni_traits)[c(2,3,4,5,6,7,8,9,10,11, 
                                     12,13,14,15,16,17,18, 
                                     19,20,21,22,23,24,25,26,
                                     27,28,29)]
wild_uni_trait_cor = wild_uni_traits %>%
  ungroup() %>%
  # split(.$lake_morph_Pair_Full_Temp) %>%
  split(.$Lake_morph) %>% 
  # ungroup() %>%
  map(select, vars_keep) %>%
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

