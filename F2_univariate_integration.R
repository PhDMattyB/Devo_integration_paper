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


# trait correlations ------------------------------------------------------

df = F2_univariate_traits %>% 
  as_tibble() %>% 
  group_by(lake_morph_Pair_Full_Temp) %>% 
  select(rowname, 
         jaw_length:body_length)

trait_matrix = df %>% 
  ungroup() %>% 
  select(jaw_length:body_length) %>% 
  as.matrix()


cor(x = trait_matrix, 
    use = 'everything', 
    method = 'pearson')



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
  facet_wrap(~lake_morph_full)+
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
       width = 70, 
       height = 30)
