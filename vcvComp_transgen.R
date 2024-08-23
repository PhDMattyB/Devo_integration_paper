##############################
## vcvComp transgen integration
##
## Matt Brachmann (PhDMattyB)
##
## 23.08.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(geomorph)
library(vcvComp)


# Metadata ----------------------------------------------------------------
F2_identifiers = read_csv('F2_metadata.csv') %>% 
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
  select(individualID, 
         Lake, 
         Ecotype_pair, 
         Morph, 
         Lake_morph)


# Univariate trait data ---------------------------------------------------


wild_univariate = read_csv('Wild_Univariate_traits.csv') %>% 
  select(-Order, 
         -rowname) %>% 
  select(ImageID, 
         Lake, 
         Morph, 
         Lake_morph, 
         everything()) %>% 
  rename(individualID = ImageID) 

F2_parental_effects = read_csv('F1_Plasticity_Corrected.csv') %>% 
  select(-Order, 
         -rowname, 
         -Lake, 
         -Ecotype_Pair_Full_Temp, 
         -lake_morph_Pair_Full_Temp, 
         -Grand_temp, 
         -Parent_temp, 
         -Offspring_temp, 
         -Full_temp) %>% 
  select(individualID, 
         Ecotype_pair, 
         Morph, 
         Lake_morph,
         everything()) %>% 
  rename(Lake = Ecotype_pair)

# F2_parental_effects = bind_cols(F2_identifiers, 
#                                 F2_parental_effects)

F2_offspring_effects = read_csv('F2_Corrected_F2_temp_only.csv')%>% 
  select(-Order, 
         -rowname, 
         -Lake, 
         -Ecotype_Pair_Full_Temp, 
         -lake_morph_Pair_Full_Temp, 
         -Grand_temp, 
         -Parent_temp, 
         -Offspring_temp, 
         -Full_temp) %>% 
  select(individualID, 
         Ecotype_pair, 
         Morph, 
         Lake_morph,
         everything())

