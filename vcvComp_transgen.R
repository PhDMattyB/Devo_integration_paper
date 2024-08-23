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
  rename(individualID = ImageID) %>% 
  filter(Lake %in% c('ASHN', 
                     'MYV', 
                     'SKR', 
                     'GTS',  
                     'CSWY'))

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

pooled_pc_coords = prcoa$PCoords %>% 
  as.data.frame() %>% 
  rownames_to_column('GROUP') %>% 
  as.tibble()


# var_exp = prcoa$Variance$exVar %>% 
#   as.tibble() %>% 
#   rename(var_explained = value)
# var_dim = 1:nrow(prcoa$Variance) %>% 
#   as.tibble() %>% 
#   rename(dimensions = value)
# 
# var_plot_data = bind_cols(var_dim,
#                           var_exp) %>%
#   mutate_if(is.integer, as.character)
# var_plot_data$dimensions = factor(var_plot_data$dimensions, 
#                                   levels=c('1', 
#                                            '2', 
#                                            '3', 
#                                            '4', 
#                                            '5', 
#                                            '6', 
#                                            '7', 
#                                            '8', 
#                                            '9', 
#                                            '10', 
#                                            '11'))
# # as.character(dimensions) %>% 
# ggplot(data = var_plot_data, 
#        aes(x = dimensions, 
#            y = var_explained))+
#   geom_col(aes(fill = (as.numeric(dimensions) %% 2 == 0)))+
#   scale_fill_manual(values = c('#2a9d8f',
#                                '#e76f51'))+
#   labs(x = 'Dimensions', 
#        y = 'Variance explained')+
#   theme(legend.position = 'none', 
#         panel.grid = element_blank(), 
#         axis.title = element_text(size = 14), 
#         axis.text = element_text(size = 12))
# 


## Good enough for now

pooled_pc_coords = pooled_pc_coords %>% 
  separate(col = GROUP, 
           into = c('Lake_morph', 
                    'Effect'), 
           sep = '_', 
           remove = FALSE)

pooled_pc_coords = mutate(.data = pooled_pc_coords, 
                          POP = as.factor(case_when(
                            Lake_morph == 'ASHNC' ~ 'ASHN',
                            Lake_morph == 'ASHNW' ~ 'ASHN',
                            Lake_morph == 'CSWY' ~ 'CSWY',
                            Lake_morph == 'GTS' ~ 'GTS',
                            Lake_morph == 'MYVC' ~ 'MYV',
                            Lake_morph == 'MYVW' ~ 'MYV',
                            Lake_morph == 'SKRC' ~ 'SKR',
                            Lake_morph == 'SKRW' ~ 'SKR')))

# pooled_pc_coords %>% 
#   write_csv('vcvComp_generation_effect.csv')

pooled_pc_coords = read_csv('vcvComp_generation_effect.csv')

WC_colour_palette = c('#023047', 
                      '#ffb703', 
                      '#219ebc', 
                      '#fb8500', 
                      '#8ecae6', 
                      '#c1121f', 
                      '#2a9d8f', 
                      '#ffafcc')

principal_coord_analysis = ggplot(data = pooled_pc_coords, 
                                  aes(x = PCo1, 
                                      y = PCo2, 
                                      col = Lake_morph, 
                                      group = Effect))+
  # geom_line(col = 'black')+
  geom_point(size = 3,
             aes(shape = Effect))+
  scale_color_manual(values = WC_colour_palette)+
  guides(col=guide_legend(title = 'Population'))+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

ggsave('~/Parsons_Postdoc/Stickleback_Morphometric_data/Principal_coord_analysis_pooled_covariance.tiff', 
       plot = principal_coord_analysis, 
       dpi = 'retina',
       units = 'cm',
       height = 10, 
       width = 15)