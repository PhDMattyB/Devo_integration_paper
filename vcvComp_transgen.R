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
library(effectsize)
library(MASS)
library(tidyverse)
library(patchwork)

theme_set(theme_bw())
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
                factor))

wild_identifiers = read_csv('TPS_Wild_metadata.csv') 

# Univariate trait data ---------------------------------------------------


wild_univariate = read_csv('Wild_scaled_kinetic_traits.csv') %>% 
  mutate(Group = 'Wild') %>%
  bind_cols(wild_identifiers, 
            .) %>% 
  dplyr::select(-Order, 
         # -rowname, 
         -Lake_morph...6) %>%
  rename(Lake_morph = Lake_morph...5) %>% 
  dplyr::select(Group, 
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

F2_orig = read_csv('F2_traits_scaled_kinetic.csv') %>% 
  mutate(Group = 'F2 generation') %>% 
  bind_cols(F2_identifiers, 
            .) %>% 
  dplyr::select(-Order, 
                # -rowname, 
                -Lake, 
                -Ecotype_Pair_Full_Temp, 
                -lake_morph_Pair_Full_Temp, 
                -Grand_temp, 
                -Parent_temp, 
                -Offspring_temp, 
                -Full_temp, 
                -Lake_morph...13) %>% 
  dplyr::select(Group, 
                individualID, 
                Ecotype_pair, 
                Morph, 
                Lake_morph...8,
                everything()) %>% 
  rename(Lake = Ecotype_pair, 
         Lake_morph = Lake_morph...8, 
         OMA = ratio1, 
         CMA = ratio2)



F2_parental_effects = read_csv('TGP_traits_scaled_kinetic.csv') %>% 
  mutate(Group = 'Transgen') %>% 
  bind_cols(F2_identifiers, 
            .) %>% 
  dplyr::select(-Order, 
         # -rowname, 
         -Lake, 
         -Ecotype_Pair_Full_Temp, 
         -lake_morph_Pair_Full_Temp, 
         -Grand_temp, 
         -Parent_temp, 
         -Offspring_temp, 
         -Full_temp, 
         -Lake_morph...13) %>% 
  dplyr::select(Group, 
         individualID, 
         Ecotype_pair, 
         Morph, 
         Lake_morph...8,
         everything()) %>% 
  rename(Lake = Ecotype_pair, 
         Lake_morph = Lake_morph...8)

# F2_parental_effects = bind_cols(F2_identifiers, 
#                                 F2_parental_effects)

F2_offspring_effects = read_csv('WGP_traits_scaled_kinetic.csv')%>%
  mutate(Group = 'Withingen') %>% 
  bind_cols(F2_identifiers, 
            .) %>% 
  dplyr::select(-Order, 
         # -rowname, 
         -Lake, 
         -Ecotype_Pair_Full_Temp, 
         -lake_morph_Pair_Full_Temp, 
         -Grand_temp, 
         -Parent_temp, 
         -Offspring_temp, 
         -Full_temp, 
         -Lake_morph...13) %>% 
  dplyr::select(Group, 
         individualID, 
         Ecotype_pair, 
         Morph, 
         Lake_morph...8,
         everything()) %>% 
  rename(Lake = Ecotype_pair, 
         Lake_morph = Lake_morph...8)


Full_data = bind_rows(wild_univariate, 
                      F2_orig,
                      F2_parental_effects, 
                      F2_offspring_effects)


# vcvComp analysis --------------------------------------------------------

Traits = Full_data %>% 
  dplyr::select(-Group, 
         -individualID, 
         -Lake, 
         -Morph, 
         -Lake_morph)

# scaled_traits = scale(Traits, 
#                       center = T, 
#                       scale = T) %>% 
#   as_tibble() 

Trait_ID = Full_data %>% 
  dplyr::select(Group, 
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


var_exp = prcoa$Variance$exVar %>%
  as.tibble() %>%
  rename(var_explained = value)
var_dim = 1:nrow(prcoa$Variance) %>%
  as.tibble() %>%
  rename(dimensions = value)

var_plot_data = bind_cols(var_dim,
                          var_exp) %>%
  mutate_if(is.integer, as.character) %>% 
  na.omit()
var_plot_data$dimensions = factor(var_plot_data$dimensions,
                                  levels=c('1',
                                           '2',
                                           '3',
                                           '4',
                                           '5',
                                           '6',
                                           '7',
                                           '8',
                                           '9',
                                           '10'))
# as.character(dimensions) %>%
Variance_explained = ggplot(data = var_plot_data,
       aes(x = dimensions,
           y = var_explained))+
  geom_col(aes(fill = (as.numeric(dimensions) %% 2 == 0)))+
  scale_fill_manual(values = c('#2a9d8f',
                               '#e76f51'))+
  labs(x = 'Dimensions',
       y = 'Variance explained')+
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

ggsave('~/Parsons_Postdoc/Stickleback_Morphometric_data/scaled_PCoA_Variance_axes_Euclidean.tiff', 
       plot = Variance_explained, 
       dpi = 'retina',
       units = 'cm',
       height = 10, 
       width = 15)


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


pooled_pc_coords %>%
  write_csv('scaled_vcvComp_generation_effect_Euclidean.csv')

# pooled_pc_coords = read_csv('scaled_vcvComp_generation_effect.csv')
# pooled_pc_coords = read_csv('vcvComp_generation_effect_F2orig.csv')
pooled_pc_coords = read_csv('scaled_vcvComp_generation_effect_Euclidean.csv')

pooled_pc_coords = mutate(.data = pooled_pc_coords, 
                          Effect = as.factor(case_when(
                            Effect == 'Wild' ~ 'Wild',
                            Effect == 'F2 generation' ~ 'F2 generation',
                            Effect == 'Transgen' ~ 'Trans-generational',
                            Effect == 'Withingen' ~ 'Within-generational')))

pooled_pc_coords = mutate(.data = pooled_pc_coords, 
                          Lake_morph = as.factor(case_when(
                            Lake_morph == 'ASHNC' ~ 'ASHNC',
                            Lake_morph == 'ASHNW' ~ 'ASHNW',
                            Lake_morph == 'CSWY' ~ 'CSWYC',
                            Lake_morph == 'GTS' ~ 'GTSW',
                            Lake_morph == 'CSWYC' ~ 'CSWYC',
                            Lake_morph == 'GTSW' ~ 'GTSW',
                            Lake_morph == 'MYVC' ~ 'MYVC',
                            Lake_morph == 'MYVW' ~ 'MYVW',
                            Lake_morph == 'SKRC' ~ 'SKRC',
                            Lake_morph == 'SKRW' ~ 'SKRW')))
WC_colour_palette = c('#22577a', 
                      '#f94144', 
                      '#38a3a5', 
                      '#f3722c', 
                      '#57cc99', 
                      '#f8961e', 
                      '#80ed99', 
                      '#f9c74f')

pooled_pc_coords %>% 
  group_by(Effect) %>% 
  summarize(variance_PCo1 = var(PCo1), 
            variance_PCo2 = var(PCo2))

pooled_pc_coords$Effect = factor(pooled_pc_coords$Effect, 
                                 levels = c('Wild',
                                            'F2 generation',
                                            'Within-generational',
                                            'Trans-generational'))

principal_coord_analysis = ggplot(data = pooled_pc_coords, 
                                  aes(x = PCo1, 
                                      y = PCo2, 
                                      col = Lake_morph, 
                                      group = Lake_morph))+
  geom_point(size = 3,
             aes(shape = Effect, 
                 col = Lake_morph))+
  stat_ellipse(aes(group = Effect), 
               col = 'black')+
  geom_line(col = 'black')+
  scale_color_manual(values = WC_colour_palette)+
  guides(col=guide_legend(title = 'Population'))+
  labs(x = 'PCoA1 (28.0%)',
       y = 'PCoA2 (19.9%)')+
  # labs(x = 'PCoA1 (81.0%)', 
  #      y = 'PCoA2 (12.6%)')+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

ggsave('~/Parsons_Postdoc/Stickleback_Morphometric_data/NEWNEW_Scaled_Principal_coord_analysis_pooled_covariance_FINAL.tiff', 
       plot = principal_coord_analysis, 
       dpi = 'retina',
       units = 'cm',
       height = 15, 
       width = 20)


# Plotting shape for each ecotype -----------------------------------------

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





# Model disparity between datasets ----------------------------------------

## wild disparity
wild_univariate

wild_traits = wild_univariate %>% 
  select(-Group, 
         -individualID, 
         -Lake, 
         -Morph, 
         -Lake_morph)


F2_whole_body_disp = morphol.disparity(coords ~ 1,
                                       groups = ~lake_morph_full,
                                       data = F2_whole_body,
                                       iter = 999)

Whole_body_pval = F2_whole_body_disp$PV.dist.Pval
Whole_body_disp = F2_whole_body_disp$PV.dist
disp_proc_var = F2_whole_body_disp$Procrustes.var





# Variable vs non-variable environments -----------------------------------


## F2 original data doesn't add anything to the plot
F2_original = read_csv('F2_Original_univariate_traits.csv') %>%
  bind_cols(F2_identifiers) %>%
  mutate(Group = 'F2') %>%
  # select(-Lake_morph...34) %>%
  select(-Lake_morph...37) %>%
  # select(Group,
  #        individualID,
  #        Lake,
  #        Morph,
  #        Lake_morph...1,
  #        everything()) %>%
  rename(Lake_morph = Lake_morph...1) %>%
  filter(Lake %in% c('ASHN',
                     'MYV',
                     'SKR',
                     'GTS',
                     'CSWY'))


F2_data = F2_original %>% 
  select(Group, 
         Lake, 
         Morph, 
         Lake_morph, 
         Full_temp, 
         Grand_temp,
         Parent_temp, 
         Offspring_temp,
         jaw_length:lm_23_2)

Traits = cbind(F2_data$jaw_length, 
               F2_data$fbar_23_24, 
               F2_data$fbar_8_24, 
               F2_data$fbar_8_27, 
               F2_data$fbar_23_27, 
               F2_data$fbar_25_26, 
               F2_data$body_width, 
               F2_data$caudal1_14_18, 
               F2_data$caudal2_15_17, 
               F2_data$body_length,
               F2_data$head_depth, 
               F2_data$jaw_2_6, 
               F2_data$lm_6_12, 
               F2_data$lm_12_13, 
               F2_data$lm_13_14, 
               F2_data$lm_14_15, 
               F2_data$lm_6_21, 
               F2_data$lm_20_21, 
               F2_data$lm_21_13, 
               F2_data$lm_20_13, 
               F2_data$lm_12_19, 
               F2_data$lm_13_19, 
               F2_data$lm_19_18, 
               F2_data$lm_18_17, 
               F2_data$lm_1_23, 
               F2_data$lm_23_2)

temp_manova = manova(Traits ~ Lake_morph*Parent_temp*Offspring_temp, 
       data = F2_data)


summary(temp_manova)


full_temp_manova = manova(Traits ~ Full_temp, 
                          data = F2_data)

summary(full_temp_manova)


Morph_temp_manova = manova(Traits ~ Morph*Full_temp, 
                          data = F2_data)

summary(Morph_temp_manova)


effectsize::eta_squared(temp_manova)
effectsize::eta_squared(full_temp_manova)
effectsize::eta_squared(Morph_temp_manova)

post_hoc_1 = lda(F2_data$Lake_morph ~ Traits, 
                 CV=F)
post_hoc_1

post_hoc_2 = lda(F2_data$Parent_temp ~ Traits, 
                 CV=F)
post_hoc_2

post_hoc_3 = lda(F2_data$Offspring_temp ~ Traits, 
                 CV=F)
post_hoc_3

post_hoc_4 = lda(F2_data$Full_temp ~ Traits, 
                 CV = F)
post_hoc_4

post_hoc_5 = lda(F2_data$Morph ~ Traits, 
                 CV = F)
post_hoc_5


temp_col_pal = c('#023047', 
                 '#e9c46a', 
                 '#0077b6', 
                 '#e76f51')



plot_lda = data.frame(F2_data[, "Full_temp"], 
                      lda = predict(post_hoc_4)$x)
lda_full_temp = ggplot(plot_lda) + 
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 colour = Full_temp), 
             size = 4)+
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  scale_color_manual(values = temp_col_pal)+
  labs(x = 'LD1 (88.35)', 
       y = 'LD2 (9.0%)')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))
ggsave('~/Parsons_Postdoc/Stickleback_Morphometric_data/WC_LAKES_COMBINED.tiff', 
       plot = lda_full_temp, 
       dpi = 'retina',
       units = 'cm',
       height = 20, 
       width = 20)


Morph_col_pal = c('#eb5e28', 
                  '#457b9d')

plot_lda_morph = data.frame(F2_data[, "Morph"], 
                      lda = predict(post_hoc_5)$x)
lda_morph = ggplot(plot_lda_morph) + 
  geom_density(aes(x = LD1, 
                 colour = Morph, 
                 fill = Morph), 
             size = 4)+
  # geom_hline(yintercept = 0, 
             # col = 'black')+
  # geom_vline(xintercept = 0, 
             # col = 'black')+
  scale_color_manual(values = Morph_col_pal)+
  labs(x = 'LD1', 
       y = 'Density')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')


WC_colour_palette = c('#22577a', 
                      '#f94144', 
                      '#38a3a5', 
                      '#f3722c', 
                      '#57cc99', 
                      '#f8961e', 
                      '#80ed99', 
                      '#f9c74f')
plot_lda_lake_morph = data.frame(F2_data[, "Lake_morph"], 
                            lda = predict(post_hoc_1)$x)

lda_lake_morph = ggplot(plot_lda_lake_morph) + 
  geom_point(aes(x = lda.LD1,
                   y = lda.LD2,
                   colour = Lake_morph, 
                   fill = Lake_morph), 
               size = 4)+
  geom_hline(yintercept = 0,
  col = 'black')+
  geom_vline(xintercept = 0,
  col = 'black')+
  scale_color_manual(values = WC_colour_palette)+
  labs(x = 'LD1 (31.2%)', 
       y = 'LD2 (28.2)')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

lda_lake_morph|lda_morph|lda_full_temp



# Divergent selection + variation due to temp -----------------------------


# ASHN Discriminant plots -------------------------------------------------
## ASHN
F2_ASHN_data = F2_data %>% 
  filter(Lake == 'ASHN')

ASHN_traits = cbind(F2_ASHN_data$jaw_length, 
                    F2_ASHN_data$fbar_23_24, 
                    F2_ASHN_data$fbar_8_24, 
                    F2_ASHN_data$fbar_8_27, 
                    F2_ASHN_data$fbar_23_27, 
                    F2_ASHN_data$fbar_25_26, 
                    F2_ASHN_data$body_width, 
                    F2_ASHN_data$caudal1_14_18, 
                    F2_ASHN_data$caudal2_15_17, 
                    F2_ASHN_data$body_length,
                    F2_ASHN_data$head_depth, 
                    F2_ASHN_data$jaw_2_6, 
                    F2_ASHN_data$lm_6_12, 
                    F2_ASHN_data$lm_12_13, 
                    F2_ASHN_data$lm_13_14, 
                    F2_ASHN_data$lm_14_15, 
                    F2_ASHN_data$lm_6_21, 
                    F2_ASHN_data$lm_20_21, 
                    F2_ASHN_data$lm_21_13, 
                    F2_ASHN_data$lm_20_13, 
                    F2_ASHN_data$lm_12_19, 
                    F2_ASHN_data$lm_13_19, 
                    F2_ASHN_data$lm_19_18, 
                    F2_ASHN_data$lm_18_17, 
                    F2_ASHN_data$lm_1_23, 
                    F2_ASHN_data$lm_23_2)

ASHN_manova = manova(ASHN_traits ~ Morph*Full_temp, 
                           data = F2_ASHN_data)
summary(ASHN_manova)

effectsize::eta_squared(ASHN_manova)

ASHN_posthoc_1 = lda(F2_ASHN_data$Morph ~ ASHN_traits, 
                 CV=F)

## tells us how well the analysis did without a cross validation
table = table(F2_ASHN_data$Morph, 
              predict(ASHN_posthoc_1)$class)
sum(table[row(table) == col(table)])/sum(table)

#leave one out cross validation
## NEED TO ADD CV = TRUE IN THE LDFA
CV = table(F2_ASHN_data$Morph, 
             ASHN_posthoc_1$class)
##Calculate the re-substitution error for the cross validation
sum(CV[row(CV) == col(CV)])/sum(CV)

ASHN_morph_lda = data.frame(F2_ASHN_data[, "Morph"], 
                            lda = predict(ASHN_posthoc_1)$x)
ASHN_lda_morph = ggplot(ASHN_morph_lda) + 
  geom_density(aes(x = LD1, 
                   colour = Morph, 
                   fill = Morph), 
               size = 4)+
  xlim(min = -5, 
       max = 5)+
  # geom_hline(yintercept = 0, 
  # col = 'black')+
  # geom_vline(xintercept = 0, 
  # col = 'black')+
  scale_color_manual(values = Morph_col_pal)+
  labs(x = 'LD1', 
       y = 'Density', 
       title = 'A)')+
  # annotate("text", 
  #          x = -3, 
  #          y = 0.4, 
  #          label= "79.4% discrimination")+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

ASHN_posthoc_2 = lda(F2_ASHN_data$Full_temp ~ ASHN_traits, 
                     CV=F)

## tells us how well the analysis did without a cross validation
table = table(F2_ASHN_data$Full_temp, 
              predict(ASHN_posthoc_2)$class)
sum(table[row(table) == col(table)])/sum(table)

#leave one out cross validation
## NEED TO ADD CV = TRUE IN THE LDFA
CV = lda(F2_ASHN_data$Full_temp ~ ASHN_traits, 
    CV=T)
CV = table(F2_ASHN_data$Full_temp, 
           CV$class)
##Calculate the re-substitution error for the cross validation
sum(CV[row(CV) == col(CV)])/sum(CV)



ASHN_ld_temp = data.frame(F2_ASHN_data[, "Full_temp"], 
                          F2_ASHN_data[,"Morph"],
                      lda = predict(ASHN_posthoc_2)$x)
ASHN_lda_temp = ggplot(ASHN_ld_temp) + 
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 colour = Full_temp), 
             size = 3)+
  xlim(min = -5, 
       max = 5)+
  ylim(min = -5, 
       max = 5)+
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  scale_color_manual(values = temp_col_pal)+
  labs(x = 'LD1 (78.2%)', 
       y = 'LD2 (17.8%)', 
       title = 'A)')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

ASHN_LDA_PLOT = ASHN_lda_morph|ASHN_lda_temp


# MYV Discriminant plots -------------------------------------------------
## MYV
F2_MYV_data = F2_data %>% 
  filter(Lake == 'MYV')

MYV_traits = cbind(F2_MYV_data$jaw_length, 
                    F2_MYV_data$fbar_23_24, 
                    F2_MYV_data$fbar_8_24, 
                    F2_MYV_data$fbar_8_27, 
                    F2_MYV_data$fbar_23_27, 
                    F2_MYV_data$fbar_25_26, 
                    F2_MYV_data$body_width, 
                    F2_MYV_data$caudal1_14_18, 
                    F2_MYV_data$caudal2_15_17, 
                    F2_MYV_data$body_length,
                    F2_MYV_data$head_depth, 
                    F2_MYV_data$jaw_2_6, 
                    F2_MYV_data$lm_6_12, 
                    F2_MYV_data$lm_12_13, 
                    F2_MYV_data$lm_13_14, 
                    F2_MYV_data$lm_14_15, 
                    F2_MYV_data$lm_6_21, 
                    F2_MYV_data$lm_20_21, 
                    F2_MYV_data$lm_21_13, 
                    F2_MYV_data$lm_20_13, 
                    F2_MYV_data$lm_12_19, 
                    F2_MYV_data$lm_13_19, 
                    F2_MYV_data$lm_19_18, 
                    F2_MYV_data$lm_18_17, 
                    F2_MYV_data$lm_1_23, 
                    F2_MYV_data$lm_23_2)

MYV_manova = manova(MYV_traits ~ Morph*Full_temp, 
                     data = F2_MYV_data)
summary(MYV_manova)

effectsize::eta_squared(MYV_manova)

MYV_posthoc_1 = lda(F2_MYV_data$Morph ~ MYV_traits, 
                     CV=F)

## tells us how well the analysis did without a cross validation
table = table(F2_MYV_data$Morph, 
              predict(MYV_posthoc_1)$class)
sum(table[row(table) == col(table)])/sum(table)

#leave one out cross validation
## NEED TO ADD CV = TRUE IN THE LDFA
CV = lda(F2_MYV_data$Morph ~ MYV_traits, 
    CV=T)
CV = table(F2_MYV_data$Morph, 
           CV$class)
##Calculate the re-substitution error for the cross validation
sum(CV[row(CV) == col(CV)])/sum(CV)

MYV_morph_lda = data.frame(F2_MYV_data[, "Morph"], 
                            lda = predict(MYV_posthoc_1)$x)
MYV_lda_morph = ggplot(MYV_morph_lda) + 
  geom_density(aes(x = LD1, 
                   colour = Morph, 
                   fill = Morph), 
               size = 4)+
  xlim(min = -5, 
       max = 5)+
  # geom_hline(yintercept = 0, 
  # col = 'black')+
  # geom_vline(xintercept = 0, 
  # col = 'black')+
  scale_color_manual(values = Morph_col_pal)+
  labs(x = 'LD1', 
       y = 'Density', 
       title = 'B)')+
  # annotate("text", 
  #          x = -2, 
  #          y = 0.4, 
  #          label= "62.9% discrimination")+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

MYV_posthoc_2 = lda(F2_MYV_data$Full_temp ~ MYV_traits, 
                     CV=F)

## tells us how well the analysis did without a cross validation
table = table(F2_MYV_data$Full_temp, 
              predict(MYV_posthoc_2)$class)
sum(table[row(table) == col(table)])/sum(table)

#leave one out cross validation
## NEED TO ADD CV = TRUE IN THE LDFA
CV = lda(F2_MYV_data$Full_temp ~ MYV_traits, 
         CV=T)
CV = table(F2_MYV_data$Full_temp, 
           CV$class)
##Calculate the re-substitution error for the cross validation
sum(CV[row(CV) == col(CV)])/sum(CV)



MYV_ld_temp = data.frame(F2_MYV_data[, "Full_temp"],
                         F2_MYV_data[, 'Morph'],
                          lda = predict(MYV_posthoc_2)$x)
MYV_lda_temp = ggplot(MYV_ld_temp) + 
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 colour = Full_temp), 
             size = 3,)+
  xlim(min = -5, 
       max = 5)+
  ylim(min = -5, 
       max = 5)+
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  scale_color_manual(values = temp_col_pal)+
  labs(x = 'LD1 (64.9%)', 
       y = 'LD2 (27.0%)', 
       title = 'B)')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

MYV_LDA_PLOT = MYV_lda_morph|MYV_lda_temp

# SKR Discriminant plots -------------------------------------------------
## SKR
F2_SKR_data = F2_data %>% 
  filter(Lake == 'SKR')

SKR_traits = cbind(F2_SKR_data$jaw_length, 
                   F2_SKR_data$fbar_23_24, 
                   F2_SKR_data$fbar_8_24, 
                   F2_SKR_data$fbar_8_27, 
                   F2_SKR_data$fbar_23_27, 
                   F2_SKR_data$fbar_25_26, 
                   F2_SKR_data$body_width, 
                   F2_SKR_data$caudal1_14_18, 
                   F2_SKR_data$caudal2_15_17, 
                   F2_SKR_data$body_length,
                   F2_SKR_data$head_depth, 
                   F2_SKR_data$jaw_2_6, 
                   F2_SKR_data$lm_6_12, 
                   F2_SKR_data$lm_12_13, 
                   F2_SKR_data$lm_13_14, 
                   F2_SKR_data$lm_14_15, 
                   F2_SKR_data$lm_6_21, 
                   F2_SKR_data$lm_20_21, 
                   F2_SKR_data$lm_21_13, 
                   F2_SKR_data$lm_20_13, 
                   F2_SKR_data$lm_12_19, 
                   F2_SKR_data$lm_13_19, 
                   F2_SKR_data$lm_19_18, 
                   F2_SKR_data$lm_18_17, 
                   F2_SKR_data$lm_1_23, 
                   F2_SKR_data$lm_23_2)

SKR_manova = manova(SKR_traits ~ Morph*Full_temp, 
                    data = F2_SKR_data)
summary(SKR_manova)

effectsize::eta_squared(SKR_manova)

SKR_posthoc_1 = lda(F2_SKR_data$Morph ~ SKR_traits, 
                    CV=F)

## tells us how well the analysis did without a cross validation
table = table(F2_SKR_data$Morph, 
              predict(SKR_posthoc_1)$class)
sum(table[row(table) == col(table)])/sum(table)

#leave one out cross validation
## NEED TO ADD CV = TRUE IN THE LDFA
CV = lda(F2_SKR_data$Morph ~ SKR_traits, 
         CV=T)
CV = table(F2_SKR_data$Morph, 
           CV$class)
##Calculate the re-substitution error for the cross validation
sum(CV[row(CV) == col(CV)])/sum(CV)

SKR_morph_lda = data.frame(F2_SKR_data[, "Morph"], 
                           lda = predict(SKR_posthoc_1)$x)
SKR_lda_morph = ggplot(SKR_morph_lda) + 
  geom_density(aes(x = LD1, 
                   colour = Morph, 
                   fill = Morph), 
               size = 4)+
  xlim(min = -5, 
       max = 5)+
  # geom_hline(yintercept = 0, 
  # col = 'black')+
  # geom_vline(xintercept = 0, 
  # col = 'black')+
  scale_color_manual(values = Morph_col_pal)+
  labs(x = 'LD1', 
       y = 'Density', 
       title = 'C)')+
  # annotate("text", 
  #          x = -2.5, 
  #          y = 0.4, 
  #          label= "78.0% discrimination")+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

SKR_posthoc_2 = lda(F2_SKR_data$Full_temp ~ SKR_traits, 
                    CV=F)

## tells us how well the analysis did without a cross validation
table = table(F2_SKR_data$Full_temp, 
              predict(SKR_posthoc_2)$class)
sum(table[row(table) == col(table)])/sum(table)

#leave one out cross validation
## NEED TO ADD CV = TRUE IN THE LDFA
CV = lda(F2_SKR_data$Full_temp ~ SKR_traits, 
         CV=T)
CV = table(F2_SKR_data$Full_temp, 
           CV$class)
##Calculate the re-substitution error for the cross validation
sum(CV[row(CV) == col(CV)])/sum(CV)



SKR_ld_temp = data.frame(F2_SKR_data[, "Full_temp"],
                         F2_SKR_data[, 'Morph'],
                         lda = predict(SKR_posthoc_2)$x)
SKR_lda_temp = ggplot(SKR_ld_temp) + 
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 colour = Full_temp), 
             size = 3)+
  xlim(min = -5, 
       max = 5)+
  ylim(min = -5, 
       max = 5)+
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  scale_color_manual(values = temp_col_pal)+
  labs(x = 'LD1 (81.9%)', 
       y = 'LD2 (10.3%)', 
       title = 'C)')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

SKR_LDA_PLOT = SKR_lda_morph|SKR_lda_temp


# GTS_CSWY Discriminant plots -------------------------------------------------
## GTS_CSWY
F2_GTS_CSWY_data = F2_data %>% 
  filter(Lake %in% c('GTS', 
                     'CSWY'))

GTS_CSWY_traits = cbind(F2_GTS_CSWY_data$jaw_length, 
                   F2_GTS_CSWY_data$fbar_23_24, 
                   F2_GTS_CSWY_data$fbar_8_24, 
                   F2_GTS_CSWY_data$fbar_8_27, 
                   F2_GTS_CSWY_data$fbar_23_27, 
                   F2_GTS_CSWY_data$fbar_25_26, 
                   F2_GTS_CSWY_data$body_width, 
                   F2_GTS_CSWY_data$caudal1_14_18, 
                   F2_GTS_CSWY_data$caudal2_15_17, 
                   F2_GTS_CSWY_data$body_length,
                   F2_GTS_CSWY_data$head_depth, 
                   F2_GTS_CSWY_data$jaw_2_6, 
                   F2_GTS_CSWY_data$lm_6_12, 
                   F2_GTS_CSWY_data$lm_12_13, 
                   F2_GTS_CSWY_data$lm_13_14, 
                   F2_GTS_CSWY_data$lm_14_15, 
                   F2_GTS_CSWY_data$lm_6_21, 
                   F2_GTS_CSWY_data$lm_20_21, 
                   F2_GTS_CSWY_data$lm_21_13, 
                   F2_GTS_CSWY_data$lm_20_13, 
                   F2_GTS_CSWY_data$lm_12_19, 
                   F2_GTS_CSWY_data$lm_13_19, 
                   F2_GTS_CSWY_data$lm_19_18, 
                   F2_GTS_CSWY_data$lm_18_17, 
                   F2_GTS_CSWY_data$lm_1_23, 
                   F2_GTS_CSWY_data$lm_23_2)

GTS_CSWY_manova = manova(GTS_CSWY_traits ~ Morph*Full_temp, 
                    data = F2_GTS_CSWY_data)
summary(GTS_CSWY_manova)

effectsize::eta_squared(GTS_CSWY_manova)

GTS_CSWY_posthoc_1 = lda(F2_GTS_CSWY_data$Morph ~ GTS_CSWY_traits, 
                    CV=F)

## tells us how well the analysis did without a cross validation
table = table(F2_GTS_CSWY_data$Morph, 
              predict(GTS_CSWY_posthoc_1)$class)
sum(table[row(table) == col(table)])/sum(table)

#leave one out cross validation
## NEED TO ADD CV = TRUE IN THE LDFA
CV = lda(F2_GTS_CSWY_data$Morph ~ GTS_CSWY_traits, 
         CV=T)
CV = table(F2_GTS_CSWY_data$Morph, 
           CV$class)
##Calculate the re-substitution error for the cross validation
sum(CV[row(CV) == col(CV)])/sum(CV)

GTS_CSWY_morph_lda = data.frame(F2_GTS_CSWY_data[, "Morph"], 
                           lda = predict(GTS_CSWY_posthoc_1)$x)
GTS_CSWY_lda_morph = ggplot(GTS_CSWY_morph_lda) + 
  geom_density(aes(x = LD1, 
                   colour = Morph, 
                   fill = Morph), 
               size = 4)+
  xlim(min = -5, 
       max = 5)+
  # geom_hline(yintercept = 0, 
  # col = 'black')+
  # geom_vline(xintercept = 0, 
  # col = 'black')+
  scale_color_manual(values = Morph_col_pal)+
  labs(x = 'LD1', 
       y = 'Density', 
       title = 'D)')+
  # annotate("text", 
  #          x = -2.5, 
  #          y = 0.4, 
  #          label= "84.8% discrimination")+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

GTS_CSWY_posthoc_2 = lda(F2_GTS_CSWY_data$Full_temp ~ GTS_CSWY_traits, 
                    CV=F)

## tells us how well the analysis did without a cross validation
table = table(F2_GTS_CSWY_data$Full_temp, 
              predict(GTS_CSWY_posthoc_2)$class)
sum(table[row(table) == col(table)])/sum(table)

#leave one out cross validation
## NEED TO ADD CV = TRUE IN THE LDFA
CV = lda(F2_GTS_CSWY_data$Full_temp ~ GTS_CSWY_traits, 
         CV=T)
CV = table(F2_GTS_CSWY_data$Full_temp, 
           CV$class)
##Calculate the re-substitution error for the cross validation
sum(CV[row(CV) == col(CV)])/sum(CV)



GTS_CSWY_ld_temp = data.frame(F2_GTS_CSWY_data[, "Full_temp"], 
                              F2_GTS_CSWY_data[, 'Morph'],
                         lda = predict(GTS_CSWY_posthoc_2)$x)
GTS_CSWY_lda_temp = ggplot(GTS_CSWY_ld_temp) + 
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 colour = Full_temp), 
             size = 3)+
  xlim(min = -5, 
       max = 5)+
  ylim(min = -5, 
       max = 5)+
  geom_hline(yintercept = 0, 
             col = 'black')+
  geom_vline(xintercept = 0, 
             col = 'black')+
  scale_color_manual(values = temp_col_pal)+
  labs(x = 'LD1 (83.8%)', 
       y = 'LD2 (10.8%)', 
       title = 'D)')+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

GTS_CSWY_LDA_PLOT = GTS_CSWY_lda_morph|GTS_CSWY_lda_temp


# LDA PLOT COMBINED -------------------------------------------------------

Big_Plot_pappi = ASHN_LDA_PLOT/MYV_LDA_PLOT/SKR_LDA_PLOT/GTS_CSWY_LDA_PLOT

Morph_divergence = (ASHN_lda_morph|MYV_lda_morph)/(SKR_lda_morph|GTS_CSWY_lda_morph)

ggsave('~/Parsons_Postdoc/Stickleback_Morphometric_data/Multivariate_trait_divergence.tiff', 
       plot = Morph_divergence, 
       dpi = 'retina',
       units = 'cm',
       height = 10, 
       width = 15)


Temp_divergence = (ASHN_lda_temp|MYV_lda_temp)/(SKR_lda_temp|GTS_CSWY_lda_temp)

ggsave('~/Parsons_Postdoc/Stickleback_Morphometric_data/Transgenerational_Temp_Effects.tiff', 
       plot = Temp_divergence, 
       dpi = 'retina',
       units = 'cm',
       height = 15, 
       width = 20)
