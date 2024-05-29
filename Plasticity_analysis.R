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

identifiers %>% 
  filter(Offspring_temp == '18') %>% 
  View()

F2_12deg_range = c(1:30, 61:91, 182:211, 244:273, 304:333, 364:382, 
                   413:442, 474:503, 534:563, 594:623, 655:683,
                   714:743, 774:803, 834:857, 871:900)
F2_18deg_range = c(31:60, 92:121, 152:181, 212:243, 274:303, 
                   334:363, 383:412, 443:473, 504:533, 564:593, 
                   624:654, 684:713, 744:773, 804:833, 858:870)

F2_array = array(0, dim = c(27, 2, 900))
for(i in F2_12deg_range){
  F2_array[,,i] = F2_gpa$coords[,,i] - F2_temp_12deg_array[,,1]
}

for(i in F2_18deg_range){
  F2_array[,,i] = F2_gpa$coords[,,i] - F2_temp_18deg_array[,,1]
}

## This is the array to use too pull out the linear traits due
## to plasticity
F2_array_consensus = array(0, dim = c(27, 2, 900))
for(i in 1:900){
  F2_array_consensus[,,i] = F2_array[,,i] + mean_shape_array[,,1]
}


test_gpa = gpagen(F2_array_consensus)
test_lm = geomorph.data.frame(test_gpa)

