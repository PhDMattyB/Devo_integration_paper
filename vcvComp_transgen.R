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
  mutate(Group = 'Wild') %>% 
  select(-Order, 
         -rowname) %>% 
  select(Group, 
         ImageID, 
         Lake, 
         Morph, 
         Lake_morph, 
         everything()) %>% 
  rename(individualID = ImageID) 

F2_parental_effects = read_csv('F1_Plasticity_Corrected.csv') %>% 
  mutate(Group = 'Transgen') %>% 
  select(-Order, 
         -rowname, 
         -Lake, 
         -Ecotype_Pair_Full_Temp, 
         -lake_morph_Pair_Full_Temp, 
         -Grand_temp, 
         -Parent_temp, 
         -Offspring_temp, 
         -Full_temp,
         -Ecotype_off_temp) %>% 
  select(Group, 
         individualID, 
         Ecotype_pair, 
         Morph, 
         Lake_morph,
         everything()) %>% 
  rename(Lake = Ecotype_pair)

# F2_parental_effects = bind_cols(F2_identifiers, 
#                                 F2_parental_effects)

F2_offspring_effects = read_csv('F2_Corrected_F2_temp_only.csv')%>%
  mutate(Group = 'Withingen') %>% 
  select(-Order, 
         -rowname, 
         -Lake, 
         -Ecotype_Pair_Full_Temp, 
         -lake_morph_Pair_Full_Temp, 
         -Grand_temp, 
         -Parent_temp, 
         -Offspring_temp, 
         -Full_temp, 
         -Ecotype_off_temp) %>% 
  select(Group, 
         individualID, 
         Ecotype_pair, 
         Morph, 
         Lake_morph,
         everything()) %>% 
  rename(Lake = Ecotype_pair)


Full_data = bind_rows(wild_univariate, 
                      F2_parental_effects, 
                      F2_offspring_effects)


# vcvComp analysis --------------------------------------------------------

Traits = Full_data %>% 
  select(-Group, 
         -individualID, 
         -Lake, 
         -Morph, 
         -Lake_morph)

Trait_ID = Full_data %>% 
  select(Group, 
         individualID, 
         Lake, 
         Morph, 
         Lake_morph) %>% 
  unite(col = 'Big_group', 
        c('Lake_morph', 
          'Group'),
        sep = '_')


PCA = prcomp(Traits, 
                       rank. = 5, 
                       tol = sqrt(.Machine$double.eps))
pca_scores = PCA$x

## scree plot to determiine number of pc axes to use in model
factoextra::fviz_eig(PCA)


phenotypes_pooled_var = cov.group(pca_scores, 
                                  groups = Trait_ID$Big_group)

phenotype_eigen_vals = mat.sq.dist(phenotypes_pooled_var, 
                                   dist. = 'Riemannian')

prcoa = pr.coord(phenotype_eigen_vals)
prcoa$Variance

