##############################
## PC loading matrices per ecotype pair
##
## Matt Brachmann (PhDMattyB)
##
## 27.05.2026
##
##############################

setwd("~/Parsons_postdoc/Stickleback_Morphometric_data/Updated Landmarks/")

library(tidyverse)
library(lme4)
library(lmerTest)
library(brms)
library(paran)
library(broom.mixed)
library(MatrixCorrelation)


# functions ---------------------------------------------------------------

cor_fun <- function(data, pc_cols, trait_cols) {
  
  expand_grid(
    pc = pc_cols,
    trait = trait_cols) %>%
    mutate(test = map2(pc,trait,
                       \(x, y)
                       cor.test(
                         x = data[[x]],
                         y = data[[y]],
                         method = "pearson") %>%
                         tidy())) %>%
    unnest(test) %>%
    transmute(
      pc,
      trait,
      r = estimate,
      p = p.value)
}

# wild data ----------------------------------------------------------------
wild_data = read_csv("WILD_SCALED_FIXED_27.05.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                rowname, 
                Lake, 
                Morph, 
                5:39)
ashn_wild_data = wild_data %>% 
  filter(Lake_morph %in% c('ASHNC', 
                           'ASHNW'))
ashn_wild_trait_mat = ashn_wild_data %>% 
  dplyr::select(5:39) %>% 
  as.matrix()

myv_wild_data = wild_data %>% 
  filter(Lake_morph %in% c('MYVC', 
                           'MYVW'))
myv_wild_trait_mat = myv_wild_data %>% 
  dplyr::select(5:39) %>% 
  as.matrix()

skr_wild_data = wild_data %>% 
  filter(Lake_morph %in% c('SKRC', 
                           'SKRW'))
skr_wild_trait_mat = skr_wild_data %>% 
  dplyr::select(5:39) %>% 
  as.matrix()

gts_gar_wild_data = wild_data %>% 
  filter(Lake_morph %in% c('CSWY', 
                           'GTS'))
gts_gar_wild_trait_mat = gts_gar_wild_data %>% 
  dplyr::select(5:39) %>% 
  as.matrix()


# PCA analysis on each wild ecotype pair ----------------------------------

ashn_wild_pca = prcomp(ashn_wild_trait_mat, 
                       center = F, 
                       scale. = F)
paran(ashn_wild_trait_mat, 
      iterations = 1000, 
      graph = T)
## 5 pc axes retained

ashn_wild_scores = as.data.frame(ashn_wild_pca$x[,1:5])
ashn_wild_PCA = bind_cols(ashn_wild_data, 
                          ashn_wild_scores)


myv_wild_pca = prcomp(myv_wild_trait_mat, 
                       center = F, 
                       scale. = F)
paran(myv_wild_trait_mat, 
      iterations = 1000, 
      graph = T)
## 5 pc axes retained

myv_wild_scores = as.data.frame(myv_wild_pca$x[,1:5])
myv_wild_PCA = bind_cols(myv_wild_data, 
                          myv_wild_scores)


skr_wild_pca = prcomp(skr_wild_trait_mat, 
                       center = F, 
                       scale. = F)
paran(skr_wild_trait_mat, 
      iterations = 1000, 
      graph = T)
## 5 pc axes retained

skr_wild_scores = as.data.frame(skr_wild_pca$x[,1:5])
skr_wild_PCA = bind_cols(skr_wild_data, 
                          skr_wild_scores)


gts_gar_wild_pca = prcomp(gts_gar_wild_trait_mat, 
                       center = F, 
                       scale. = F)
paran(gts_gar_wild_trait_mat, 
      iterations = 1000, 
      graph = T)
## 6 pc axes retained

gts_gar_wild_scores = as.data.frame(gts_gar_wild_pca$x[,1:6])
gts_gar_wild_PCA = bind_cols(gts_gar_wild_data, 
                          gts_gar_wild_scores)



# wild data - trait/pc corelations -------------------------------------------

# ASHN --------------------------------------------------------------------

ashn_trait_cols = ashn_wild_PCA %>% 
  dplyr::select(5:39) %>% 
  names()

ashn_pca_cols = ashn_wild_PCA %>% 
  dplyr::select(40:44) %>% 
  names()

ashn_wild_WC_loads = ashn_wild_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = ashn_pca_cols, 
    trait_cols = ashn_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

ashn_wild_WC_loads %>% 
  ggplot(aes(pc, 
             trait, 
             fill = r))+
  geom_tile()+
  facet_wrap(~Morph) +
  scale_fill_gradient2(
    low = "#277da1",
    mid = "white",
    high = "#f94144",
    limits = c(-1, 1)) +
  theme(axis.text.y = element_text(size = 7), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.title = element_blank())

ashn_wild_cold = ashn_wild_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
ashn_wild_warm = ashn_wild_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# MYV ---------------------------------------------------------------------
myv_trait_cols = myv_wild_PCA %>% 
  dplyr::select(5:39) %>% 
  names()

myv_pca_cols = myv_wild_PCA %>% 
  dplyr::select(40:44) %>% 
  names()

myv_wild_WC_loads = myv_wild_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = myv_pca_cols, 
    trait_cols = myv_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

myv_wild_WC_loads %>% 
  ggplot(aes(pc, 
             trait, 
             fill = r))+
  geom_tile()+
  facet_wrap(~Morph) +
  scale_fill_gradient2(
    low = "#277da1",
    mid = "white",
    high = "#f94144",
    limits = c(-1, 1)) +
  theme(axis.text.y = element_text(size = 7), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.title = element_blank())

myv_wild_cold = myv_wild_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
myv_wild_warm = myv_wild_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# SKR ---------------------------------------------------------------------
skr_trait_cols = skr_wild_PCA %>% 
  dplyr::select(5:39) %>% 
  names()

skr_pca_cols = skr_wild_PCA %>% 
  dplyr::select(40:44) %>% 
  names()

skr_wild_WC_loads = skr_wild_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = skr_pca_cols, 
    trait_cols = skr_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

skr_wild_WC_loads %>% 
  ggplot(aes(pc, 
             trait, 
             fill = r))+
  geom_tile()+
  facet_wrap(~Morph) +
  scale_fill_gradient2(
    low = "#277da1",
    mid = "white",
    high = "#f94144",
    limits = c(-1, 1)) +
  theme(axis.text.y = element_text(size = 7), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.title = element_blank())

skr_wild_cold = skr_wild_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
skr_wild_warm = skr_wild_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# GTS/GAR -----------------------------------------------------------------
gts_gar_trait_cols = gts_gar_wild_PCA %>% 
  dplyr::select(5:39) %>% 
  names()

gts_gar_pca_cols = gts_gar_wild_PCA %>% 
  dplyr::select(40:45) %>% 
  names()

gts_gar_wild_WC_loads = gts_gar_wild_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = gts_gar_pca_cols, 
    trait_cols = gts_gar_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

gts_gar_wild_WC_loads %>% 
  ggplot(aes(pc, 
             trait, 
             fill = r))+
  geom_tile()+
  facet_wrap(~Morph) +
  scale_fill_gradient2(
    low = "#277da1",
    mid = "white",
    high = "#f94144",
    limits = c(-1, 1)) +
  theme(axis.text.y = element_text(size = 7), 
        strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        axis.title = element_blank())

gts_gar_wild_cold = gts_gar_wild_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
gts_gar_wild_warm = gts_gar_wild_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# wild data - matrix differences -------------------------------------------

# ASHN --------------------------------------------------------------------

ashn_wild_cold_mat = ashn_wild_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

ashn_wild_warm_mat = ashn_wild_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

ashn_wild_diff_mat = ashn_wild_warm_mat - ashn_wild_cold_mat

ashn_wild_diff_mat_long = ashn_wild_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

ashn_wild_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

ashn_wild_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

ashn_wild_diff_mat_plot = 
  ggplot(ashn_wild_diff_mat_long, 
         aes(x = PC, 
             y = trait, 
             fill = diff))+
  geom_tile() +
  scale_fill_gradient2(
    low = "#277da1",
    mid = "white",
    high = "#f94144") +
  labs(
    x = "PC axis",
    y = "Trait",
    fill = "Δ loading",
    title = "A) ASHN warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# MYV ---------------------------------------------------------------------
myv_wild_cold_mat = myv_wild_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

myv_wild_warm_mat = myv_wild_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

myv_wild_diff_mat = myv_wild_warm_mat - myv_wild_cold_mat

myv_wild_diff_mat_long = myv_wild_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

myv_wild_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

myv_wild_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

myv_wild_diff_mat_plot = 
  ggplot(myv_wild_diff_mat_long, 
         aes(x = PC, 
             y = trait, 
             fill = diff))+
  geom_tile() +
  scale_fill_gradient2(
    low = "#277da1",
    mid = "white",
    high = "#f94144") +
  labs(
    x = "PC axis",
    y = "Trait",
    fill = "Δ loading",
    title = "B) MYV warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# SKR ---------------------------------------------------------------------
skr_wild_cold_mat = skr_wild_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

skr_wild_warm_mat = skr_wild_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

skr_wild_diff_mat = skr_wild_warm_mat - skr_wild_cold_mat

skr_wild_diff_mat_long = skr_wild_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

skr_wild_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

skr_wild_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

skr_wild_diff_mat_plot = 
  ggplot(skr_wild_diff_mat_long, 
         aes(x = PC, 
             y = trait, 
             fill = diff))+
  geom_tile() +
  scale_fill_gradient2(
    low = "#277da1",
    mid = "white",
    high = "#f94144") +
  labs(
    x = "PC axis",
    y = "Trait",
    fill = "Δ loading",
    title = "C) SKR warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# GTS/GAR -----------------------------------------------------------------
gts_gar_wild_cold_mat = gts_gar_wild_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

gts_gar_wild_warm_mat = gts_gar_wild_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

gts_gar_wild_diff_mat = gts_gar_wild_warm_mat - gts_gar_wild_cold_mat

gts_gar_wild_diff_mat_long = gts_gar_wild_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

gts_gar_wild_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

gts_gar_wild_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

gts_gar_wild_diff_mat_plot = 
  ggplot(gts_gar_wild_diff_mat_long, 
         aes(x = PC, 
             y = trait, 
             fill = diff))+
  geom_tile() +
  scale_fill_gradient2(
    low = "#277da1",
    mid = "white",
    high = "#f94144") +
  labs(
    x = "PC axis",
    y = "Trait",
    fill = "Δ loading",
    title = "D) GTS/GAR warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# wild data - Effect size --------------------------------------------------

# ASHN --------------------------------------------------------------------
## Taking the difference matrix and finding per trait effect size

ashn_wild_trait_effect = ashn_wild_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

ashn_trait_effect_plot = 
  ggplot(ashn_wild_trait_effect, 
         aes(reorder(trait, 
                     effect_size), 
             effect_size))+
  geom_col(col = 'black',
           fill = '#480355') +
  coord_flip()+
  ylim(0.0, 3)+
  labs(y = 'Absolute sum of r-values across axes', 
       title = 'A) ASHN')+
  theme(axis.title.y = element_blank())


# MYV ---------------------------------------------------------------------
myv_wild_trait_effect = myv_wild_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

myv_trait_effect_plot = 
  ggplot(myv_wild_trait_effect, 
         aes(reorder(trait, 
                     effect_size), 
             effect_size))+
  geom_col(col = 'black',
           fill = '#480355') +
  coord_flip()+
  ylim(0.0, 3)+
  labs(y = 'Absolute sum of r-values across axes', 
       title = 'A) MYV')+
  theme(axis.title.y = element_blank())

# SKR ---------------------------------------------------------------------
skr_wild_trait_effect = skr_wild_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

skr_trait_effect_plot = 
  ggplot(skr_wild_trait_effect, 
         aes(reorder(trait, 
                     effect_size), 
             effect_size))+
  geom_col(col = 'black',
           fill = '#480355') +
  coord_flip()+
  ylim(0.0, 3)+
  labs(y = 'Absolute sum of r-values across axes', 
       title = 'A) SKR')+
  theme(axis.title.y = element_blank())

# GTS/GAR -----------------------------------------------------------------

gts_gar_wild_trait_effect = gts_gar_wild_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

gts_gar_trait_effect_plot = 
  ggplot(gts_gar_wild_trait_effect, 
         aes(reorder(trait, 
                     effect_size), 
             effect_size))+
  geom_col(col = 'black',
           fill = '#480355') +
  coord_flip()+
  ylim(0.0, 3)+
  labs(y = 'Absolute sum of r-values across axes', 
       title = 'A) GTS/GAR')+
  theme(axis.title.y = element_blank())

# wild data - Permutation test RV coefficient ------------------------------

# ASHN --------------------------------------------------------------------

ashn_wild_PCA

ashn_trait_cols = ashn_wild_PCA %>% 
  dplyr::select(5:39) %>% 
  names()

ashn_pca_cols = ashn_wild_PCA %>% 
  dplyr::select(40:44) %>% 
  names()

ashn_wild_perm_test_data = ashn_wild_PCA %>% 
  dplyr::select(Morph, 
                5:39, 
                40:44)

ashn_wild_obs = RV(ashn_wild_cold_mat, 
                   ashn_wild_warm_mat)
set.seed(123)

ashn_wild_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- ashn_wild_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = ashn_pca_cols,
        trait_cols = ashn_trait_cols)) %>% 
    # dplyr::select(-data) %>%
    unnest(correlations)
  
  # 3. build matrices
  cold_perm <- perm_results %>%
    filter(Morph == "Cold") %>%
    ungroup() %>% 
    dplyr::select(trait, pc, r) %>%
    pivot_wider(names_from = pc, values_from = r) %>%
    column_to_rownames("trait") %>%
    as.matrix()
  
  warm_perm <- perm_results %>%
    filter(Morph == "Warm") %>%
    ungroup() %>% 
    dplyr::select(trait, pc, r) %>%
    pivot_wider(names_from = pc, values_from = r) %>%
    column_to_rownames("trait") %>%
    as.matrix()
  
  # 4. RV statistic
  RV(cold_perm, warm_perm)
})

ashn_wild_delta_RV = mean(ashn_wild_null_rv) - ashn_wild_obs

ashn_wild_z = (ashn_wild_obs - mean(ashn_wild_null_rv)) / sd(ashn_wild_null_rv)

ashn_wild_rv_df = data.frame(ashn_wild_null_rv = ashn_wild_null_rv)

ashn_wild_RV_perm_plot = 
  ggplot(ashn_wild_rv_df, 
         aes(x = ashn_wild_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = ashn_wild_obs), 
             color = "#480355", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "A) ASHN",
    subtitle = paste0("Observed RV = ", round(ashn_wild_obs, 3)))

