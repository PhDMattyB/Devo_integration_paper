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


# MYV ---------------------------------------------------------------------
myv_wild_PCA

myv_trait_cols = myv_wild_PCA %>% 
  dplyr::select(5:39) %>% 
  names()

myv_pca_cols = myv_wild_PCA %>% 
  dplyr::select(40:44) %>% 
  names()

myv_wild_perm_test_data = myv_wild_PCA %>% 
  dplyr::select(Morph, 
                5:39, 
                40:44)

myv_wild_obs = RV(myv_wild_cold_mat, 
                   myv_wild_warm_mat)
set.seed(123)

myv_wild_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- myv_wild_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = myv_pca_cols,
        trait_cols = myv_trait_cols)) %>% 
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

myv_wild_delta_RV = mean(myv_wild_null_rv) - myv_wild_obs

myv_wild_z = (myv_wild_obs - mean(myv_wild_null_rv)) / sd(myv_wild_null_rv)

myv_wild_rv_df = data.frame(myv_wild_null_rv = myv_wild_null_rv)

myv_wild_RV_perm_plot = 
  ggplot(myv_wild_rv_df, 
         aes(x = myv_wild_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = myv_wild_obs), 
             color = "#480355", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "B) MYV",
    subtitle = paste0("Observed RV = ", round(myv_wild_obs, 3)))


# SKR ---------------------------------------------------------------------
skr_wild_PCA

skr_trait_cols = skr_wild_PCA %>% 
  dplyr::select(5:39) %>% 
  names()

skr_pca_cols = skr_wild_PCA %>% 
  dplyr::select(40:44) %>% 
  names()

skr_wild_perm_test_data = skr_wild_PCA %>% 
  dplyr::select(Morph, 
                5:39, 
                40:44)

skr_wild_obs = RV(skr_wild_cold_mat, 
                  skr_wild_warm_mat)
set.seed(123)

skr_wild_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- skr_wild_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = skr_pca_cols,
        trait_cols = skr_trait_cols)) %>% 
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

skr_wild_delta_RV = mean(skr_wild_null_rv) - skr_wild_obs

skr_wild_z = (skr_wild_obs - mean(skr_wild_null_rv)) / sd(skr_wild_null_rv)

skr_wild_rv_df = data.frame(skr_wild_null_rv = skr_wild_null_rv)

skr_wild_RV_perm_plot = 
  ggplot(skr_wild_rv_df, 
         aes(x = skr_wild_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = skr_wild_obs), 
             color = "#480355", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "C) SKR",
    subtitle = paste0("Observed RV = ", round(skr_wild_obs, 3)))



# GTS/GAR -----------------------------------------------------------------
gts_gar_wild_PCA

gts_gar_trait_cols = gts_gar_wild_PCA %>% 
  dplyr::select(5:39) %>% 
  names()

gts_gar_pca_cols = gts_gar_wild_PCA %>% 
  dplyr::select(40:45) %>% 
  names()

gts_gar_wild_perm_test_data = gts_gar_wild_PCA %>% 
  dplyr::select(Morph, 
                5:39, 
                40:45)

gts_gar_wild_obs = RV(gts_gar_wild_cold_mat, 
                  gts_gar_wild_warm_mat)
set.seed(123)

gts_gar_wild_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- gts_gar_wild_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = gts_gar_pca_cols,
        trait_cols = gts_gar_trait_cols)) %>% 
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

gts_gar_wild_delta_RV = mean(gts_gar_wild_null_rv) - gts_gar_wild_obs

gts_gar_wild_z = (gts_gar_wild_obs - mean(gts_gar_wild_null_rv)) / sd(gts_gar_wild_null_rv)

gts_gar_wild_rv_df = data.frame(gts_gar_wild_null_rv = gts_gar_wild_null_rv)

gts_gar_wild_RV_perm_plot = 
  ggplot(gts_gar_wild_rv_df, 
         aes(x = gts_gar_wild_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = gts_gar_wild_obs), 
             color = "#480355", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "D) GTS/GAR",
    subtitle = paste0("Observed RV = ", round(gts_gar_wild_obs, 3)))



# F2 DATA -----------------------------------------------------------------


F2_data = read_csv("F2_Original_univariate_traits_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17)

ashn_F2_orig_data = F2_data %>% 
  filter(Lake_morph %in% c('ASHNC', 
                           'ASHNW'))

ashn_F2_orig_trait_mat = ashn_F2_orig_data %>% 
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  as.matrix()

myv_F2_orig_data = F2_data %>% 
  filter(Lake_morph %in% c('MYVC', 
                           'MYVW'))

myv_F2_orig_trait_mat = myv_F2_orig_data %>% 
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  as.matrix()

skr_F2_orig_data = F2_data %>% 
  filter(Lake_morph %in% c('SKRC', 
                           'SKRW'))
skr_F2_orig_trait_mat = skr_F2_orig_data %>% 
  dplyr::select(jaw_length:caudal2_15_17) %>%  
  as.matrix()

gts_gar_F2_orig_data = F2_data %>% 
  filter(Lake_morph %in% c('CSWYC', 
                           'GTSW'))
gts_gar_F2_orig_trait_mat = gts_gar_F2_orig_data %>% 
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  as.matrix()


# F2 Orig PCA -------------------------------------------------------------
ashn_F2_orig_pca = prcomp(ashn_F2_orig_trait_mat, 
                       center = F, 
                       scale. = F)
paran(ashn_F2_orig_trait_mat, 
      iterations = 1000, 
      graph = T)
## 6 pc axes retained

ashn_F2_orig_scores = as.data.frame(ashn_F2_orig_pca$x[,1:6])
ashn_F2_orig_PCA = bind_cols(ashn_F2_orig_data, 
                          ashn_F2_orig_scores)


myv_F2_orig_pca = prcomp(myv_F2_orig_trait_mat, 
                      center = F, 
                      scale. = F)
paran(myv_F2_orig_trait_mat, 
      iterations = 1000, 
      graph = T)
## 7 pc axes retained

myv_F2_orig_scores = as.data.frame(myv_F2_orig_pca$x[,1:7])
myv_F2_orig_PCA = bind_cols(myv_F2_orig_data, 
                         myv_F2_orig_scores)


skr_F2_orig_pca = prcomp(skr_F2_orig_trait_mat, 
                      center = F, 
                      scale. = F)
paran(skr_F2_orig_trait_mat, 
      iterations = 1000, 
      graph = T)
## 7 pc axes retained

skr_F2_orig_scores = as.data.frame(skr_F2_orig_pca$x[,1:7])
skr_F2_orig_PCA = bind_cols(skr_F2_orig_data, 
                         skr_F2_orig_scores)


gts_gar_F2_orig_pca = prcomp(gts_gar_F2_orig_trait_mat, 
                          center = F, 
                          scale. = F)
paran(gts_gar_F2_orig_trait_mat, 
      iterations = 1000, 
      graph = T)
## 7 pc axes retained

gts_gar_F2_orig_scores = as.data.frame(gts_gar_F2_orig_pca$x[,1:7])
gts_gar_F2_orig_PCA = bind_cols(gts_gar_F2_orig_data, 
                             gts_gar_F2_orig_scores)

# F2 ORIG data - trait/pc corelations -------------------------------------------

# ASHN --------------------------------------------------------------------

ashn_F2_orig_trait_cols = ashn_F2_orig_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

ashn_F2_orig_pca_cols = ashn_F2_orig_PCA %>% 
  dplyr::select(40:45) %>% 
  names()

ashn_F2_orig_WC_loads = ashn_F2_orig_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = ashn_F2_orig_pca_cols, 
    trait_cols = ashn_F2_orig_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

ashn_F2_orig_WC_loads %>% 
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

ashn_F2_orig_cold = ashn_F2_orig_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
ashn_F2_orig_warm = ashn_F2_orig_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# MYV ---------------------------------------------------------------------
myv_F2_orig_trait_cols = myv_F2_orig_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

myv_F2_orig_pca_cols = myv_F2_orig_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

myv_F2_orig_WC_loads = myv_F2_orig_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = myv_F2_orig_pca_cols, 
    trait_cols = myv_F2_orig_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

myv_F2_orig_WC_loads %>% 
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

myv_F2_orig_cold = myv_F2_orig_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
myv_F2_orig_warm = myv_F2_orig_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# SKR ---------------------------------------------------------------------
skr_F2_orig_trait_cols = skr_F2_orig_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

skr_F2_orig_pca_cols = skr_F2_orig_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

skr_F2_orig_WC_loads = skr_F2_orig_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = skr_F2_orig_pca_cols, 
    trait_cols = skr_F2_orig_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

skr_F2_orig_WC_loads %>% 
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

skr_F2_orig_cold = skr_F2_orig_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
skr_F2_orig_warm = skr_F2_orig_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# GTS/GAR -----------------------------------------------------------------
gts_gar_F2_orig_trait_cols = gts_gar_F2_orig_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

gts_gar_F2_orig_pca_cols = gts_gar_F2_orig_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

gts_gar_F2_orig_WC_loads = gts_gar_F2_orig_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = gts_gar_F2_orig_pca_cols, 
    trait_cols = gts_gar_F2_orig_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

gts_gar_F2_orig_WC_loads %>% 
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

gts_gar_F2_orig_cold = gts_gar_F2_orig_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
gts_gar_F2_orig_warm = gts_gar_F2_orig_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# F2 ORIG data - matrix differences -------------------------------------------

# ASHN --------------------------------------------------------------------

ashn_F2_orig_cold_mat = ashn_F2_orig_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

ashn_F2_orig_warm_mat = ashn_F2_orig_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

ashn_F2_orig_diff_mat = ashn_F2_orig_warm_mat - ashn_F2_orig_cold_mat

ashn_F2_orig_diff_mat_long = ashn_F2_orig_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

ashn_F2_orig_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

ashn_F2_orig_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

ashn_F2_orig_diff_mat_plot = 
  ggplot(ashn_F2_orig_diff_mat_long, 
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
    title = "A) ASHN F2 warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# MYV ---------------------------------------------------------------------
myv_F2_orig_cold_mat = myv_F2_orig_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

myv_F2_orig_warm_mat = myv_F2_orig_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

myv_F2_orig_diff_mat = myv_F2_orig_warm_mat - myv_F2_orig_cold_mat

myv_F2_orig_diff_mat_long = myv_F2_orig_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

myv_F2_orig_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

myv_F2_orig_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

myv_F2_orig_diff_mat_plot = 
  ggplot(myv_F2_orig_diff_mat_long, 
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
    title = "B) MYV F2 warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# SKR ---------------------------------------------------------------------
skr_F2_orig_cold_mat = skr_F2_orig_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

skr_F2_orig_warm_mat = skr_F2_orig_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

skr_F2_orig_diff_mat = skr_F2_orig_warm_mat - skr_F2_orig_cold_mat

skr_F2_orig_diff_mat_long = skr_F2_orig_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

skr_F2_orig_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

skr_F2_orig_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

skr_F2_orig_diff_mat_plot = 
  ggplot(skr_F2_orig_diff_mat_long, 
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
    title = "C) SKR F2 warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# GTS/GAR -----------------------------------------------------------------
gts_gar_F2_orig_cold_mat = gts_gar_F2_orig_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

gts_gar_F2_orig_warm_mat = gts_gar_F2_orig_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

gts_gar_F2_orig_diff_mat = gts_gar_F2_orig_warm_mat - gts_gar_F2_orig_cold_mat

gts_gar_F2_orig_diff_mat_long = gts_gar_F2_orig_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

gts_gar_F2_orig_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

gts_gar_F2_orig_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

gts_gar_F2_orig_diff_mat_plot = 
  ggplot(gts_gar_F2_orig_diff_mat_long, 
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
    title = "D) GTS/GAR F2 warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')







# F2 ORIG data - Effect size --------------------------------------------------

# ASHN --------------------------------------------------------------------
## Taking the difference matrix and finding per trait effect size

ashn_F2_orig_trait_effect = ashn_F2_orig_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

ashn_trait_effect_plot = 
  ggplot(ashn_F2_orig_trait_effect, 
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
myv_F2_orig_trait_effect = myv_F2_orig_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

myv_trait_effect_plot = 
  ggplot(myv_F2_orig_trait_effect, 
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
skr_F2_orig_trait_effect = skr_F2_orig_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

skr_trait_effect_plot = 
  ggplot(skr_F2_orig_trait_effect, 
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

gts_gar_F2_orig_trait_effect = gts_gar_F2_orig_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

gts_gar_trait_effect_plot = 
  ggplot(gts_gar_F2_orig_trait_effect, 
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

# F2 ORIG data - Permutation test RV coefficient ------------------------------

# ASHN --------------------------------------------------------------------

ashn_F2_orig_PCA

# ashn_trait_cols = ashn_F2_orig_PCA %>% 
#   dplyr::select(7:39) %>% 
#   names()
# 
# ashn_pca_cols = ashn_F2_orig_PCA %>% 
#   dplyr::select(40:45) %>% 
#   names()

ashn_F2_orig_perm_test_data = ashn_F2_orig_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:45)

ashn_F2_orig_obs = RV(ashn_F2_orig_cold_mat, 
                   ashn_F2_orig_warm_mat)
set.seed(123)

ashn_F2_orig_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- ashn_F2_orig_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = ashn_F2_orig_pca_cols,
        trait_cols = ashn_F2_orig_trait_cols)) %>% 
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

ashn_F2_orig_delta_RV = mean(ashn_F2_orig_null_rv) - ashn_F2_orig_obs

ashn_F2_orig_z = (ashn_F2_orig_obs - mean(ashn_F2_orig_null_rv)) / sd(ashn_F2_orig_null_rv)

ashn_F2_orig_rv_df = data.frame(ashn_F2_orig_null_rv = ashn_F2_orig_null_rv)

ashn_F2_orig_RV_perm_plot = 
  ggplot(ashn_F2_orig_rv_df, 
         aes(x = ashn_F2_orig_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = ashn_F2_orig_obs), 
             color = "#480355", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "A) ASHN F2 unmodified",
    subtitle = paste0("Observed RV = ", round(ashn_F2_orig_obs, 3)))


# MYV ---------------------------------------------------------------------
myv_F2_orig_PCA

myv_F2_orig_trait_cols = myv_F2_orig_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

myv_F2_orig_pca_cols = myv_F2_orig_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

myv_F2_orig_perm_test_data = myv_F2_orig_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)

myv_F2_orig_obs = RV(myv_F2_orig_cold_mat, 
                  myv_F2_orig_warm_mat)
set.seed(123)

myv_F2_orig_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- myv_F2_orig_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = myv_F2_orig_pca_cols,
        trait_cols = myv_F2_orig_trait_cols)) %>% 
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

myv_F2_orig_delta_RV = mean(myv_F2_orig_null_rv) - myv_F2_orig_obs

myv_F2_orig_z = (myv_F2_orig_obs - mean(myv_F2_orig_null_rv)) / sd(myv_F2_orig_null_rv)

myv_F2_orig_rv_df = data.frame(myv_F2_orig_null_rv = myv_F2_orig_null_rv)

myv_F2_orig_RV_perm_plot = 
  ggplot(myv_F2_orig_rv_df, 
         aes(x = myv_F2_orig_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = myv_F2_orig_obs), 
             color = "#480355", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "B) MYV F2 Unmodified",
    subtitle = paste0("Observed RV = ", round(myv_F2_orig_obs, 3)))


# SKR ---------------------------------------------------------------------
skr_F2_orig_PCA

skr_F2_orig_trait_cols = skr_F2_orig_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

skr_F2_orig_pca_cols = skr_F2_orig_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

skr_F2_orig_perm_test_data = skr_F2_orig_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)

skr_F2_orig_obs = RV(skr_F2_orig_cold_mat, 
                  skr_F2_orig_warm_mat)
set.seed(123)

skr_F2_orig_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- skr_F2_orig_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = skr_F2_orig_pca_cols,
        trait_cols = skr_F2_orig_trait_cols)) %>% 
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

skr_F2_orig_delta_RV = mean(skr_F2_orig_null_rv) - skr_F2_orig_obs

skr_F2_orig_z = (skr_F2_orig_obs - mean(skr_F2_orig_null_rv)) / sd(skr_F2_orig_null_rv)

skr_F2_orig_rv_df = data.frame(skr_F2_orig_null_rv = skr_F2_orig_null_rv)

skr_F2_orig_RV_perm_plot = 
  ggplot(skr_F2_orig_rv_df, 
         aes(x = skr_F2_orig_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = skr_F2_orig_obs), 
             color = "#480355", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "C) SKR F2 Unmodified",
    subtitle = paste0("Observed RV = ", round(skr_F2_orig_obs, 3)))



# GTS/GAR -----------------------------------------------------------------
gts_gar_F2_orig_PCA

gts_gar_F2_orig_trait_cols = gts_gar_F2_orig_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

gts_gar_F2_orig_pca_cols = gts_gar_F2_orig_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

gts_gar_F2_orig_perm_test_data = gts_gar_F2_orig_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)

gts_gar_F2_orig_obs = RV(gts_gar_F2_orig_cold_mat, 
                      gts_gar_F2_orig_warm_mat)
set.seed(123)

gts_gar_F2_orig_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- gts_gar_F2_orig_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = gts_gar_F2_orig_pca_cols,
        trait_cols = gts_gar_F2_orig_trait_cols)) %>% 
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

gts_gar_F2_orig_delta_RV = mean(gts_gar_F2_orig_null_rv) - gts_gar_F2_orig_obs

gts_gar_F2_orig_z = (gts_gar_F2_orig_obs - mean(gts_gar_F2_orig_null_rv)) / sd(gts_gar_F2_orig_null_rv)

gts_gar_F2_orig_rv_df = data.frame(gts_gar_F2_orig_null_rv = gts_gar_F2_orig_null_rv)

gts_gar_F2_orig_RV_perm_plot = 
  ggplot(gts_gar_F2_orig_rv_df, 
         aes(x = gts_gar_F2_orig_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = gts_gar_F2_orig_obs), 
             color = "#480355", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "D) GTS/GAR F2 Unmodified",
    subtitle = paste0("Observed RV = ", round(gts_gar_F2_orig_obs, 3)))




# WGP data ----------------------------------------------------------------

WGP_data = read_csv("WGP_TRAITS_SCALED_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17)


# WGP_trait_mat <- WGP_data %>%
#   dplyr::select(jaw_length:caudal2_15_17) %>% 
#   # dplyr::select(starts_with('value'))
#   # dplyr::select(all_of(trait_names)) %>%
#   as.matrix()
# 
# 
ashn_WGP_data = WGP_data %>% 
  filter(Lake_morph %in% c('ASHNC', 
                           'ASHNW'))

ashn_WGP_trait_mat = ashn_WGP_data %>% 
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  as.matrix()

myv_WGP_data = F2_data %>% 
  filter(Lake_morph %in% c('MYVC', 
                           'MYVW'))

myv_WGP_trait_mat = myv_WGP_data %>% 
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  as.matrix()

skr_WGP_data = F2_data %>% 
  filter(Lake_morph %in% c('SKRC', 
                           'SKRW'))
skr_WGP_trait_mat = skr_WGP_data %>% 
  dplyr::select(jaw_length:caudal2_15_17) %>%  
  as.matrix()

gts_gar_WGP_data = F2_data %>% 
  filter(Lake_morph %in% c('CSWYC', 
                           'GTSW'))
gts_gar_WGP_trait_mat = gts_gar_WGP_data %>% 
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  as.matrix()


# WGP PCA -------------------------------------------------------------
ashn_WGP_pca = prcomp(ashn_WGP_trait_mat, 
                          center = F, 
                          scale. = F)
paran(ashn_WGP_trait_mat, 
      iterations = 1000, 
      graph = T)
## 7 pc axes retained

ashn_WGP_scores = as.data.frame(ashn_WGP_pca$x[,1:7])
ashn_WGP_PCA = bind_cols(ashn_WGP_data, 
                             ashn_WGP_scores)


myv_WGP_pca = prcomp(myv_WGP_trait_mat, 
                         center = F, 
                         scale. = F)
paran(myv_WGP_trait_mat, 
      iterations = 1000, 
      graph = T)
## 7 pc axes retained

myv_WGP_scores = as.data.frame(myv_WGP_pca$x[,1:7])
myv_WGP_PCA = bind_cols(myv_WGP_data, 
                            myv_WGP_scores)


skr_WGP_pca = prcomp(skr_WGP_trait_mat, 
                         center = F, 
                         scale. = F)
paran(skr_WGP_trait_mat, 
      iterations = 1000, 
      graph = T)
## 7 pc axes retained

skr_WGP_scores = as.data.frame(skr_WGP_pca$x[,1:7])
skr_WGP_PCA = bind_cols(skr_WGP_data, 
                            skr_WGP_scores)


gts_gar_WGP_pca = prcomp(gts_gar_WGP_trait_mat, 
                             center = F, 
                             scale. = F)
paran(gts_gar_WGP_trait_mat, 
      iterations = 1000, 
      graph = T)
## 7 pc axes retained

gts_gar_WGP_scores = as.data.frame(gts_gar_WGP_pca$x[,1:7])
gts_gar_WGP_PCA = bind_cols(gts_gar_WGP_data, 
                                gts_gar_WGP_scores)

# WGP data - trait/pc corelations -------------------------------------------

# ASHN --------------------------------------------------------------------

ashn_WGP_trait_cols = ashn_WGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

ashn_WGP_pca_cols = ashn_WGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

ashn_WGP_WC_loads = ashn_WGP_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = ashn_WGP_pca_cols, 
    trait_cols = ashn_WGP_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

ashn_WGP_WC_loads %>% 
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

ashn_WGP_cold = ashn_WGP_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
ashn_WGP_warm = ashn_WGP_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# MYV ---------------------------------------------------------------------
myv_WGP_trait_cols = myv_WGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

myv_WGP_pca_cols = myv_WGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

myv_WGP_WC_loads = myv_WGP_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = myv_WGP_pca_cols, 
    trait_cols = myv_WGP_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

myv_WGP_WC_loads %>% 
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

myv_WGP_cold = myv_WGP_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
myv_WGP_warm = myv_WGP_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# SKR ---------------------------------------------------------------------
skr_WGP_trait_cols = skr_WGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

skr_WGP_pca_cols = skr_WGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

skr_WGP_WC_loads = skr_WGP_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = skr_WGP_pca_cols, 
    trait_cols = skr_WGP_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

skr_WGP_WC_loads %>% 
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

skr_WGP_cold = skr_WGP_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
skr_WGP_warm = skr_WGP_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# GTS/GAR -----------------------------------------------------------------
gts_gar_WGP_trait_cols = gts_gar_WGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

gts_gar_WGP_pca_cols = gts_gar_WGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

gts_gar_WGP_WC_loads = gts_gar_WGP_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = gts_gar_WGP_pca_cols, 
    trait_cols = gts_gar_WGP_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

gts_gar_WGP_WC_loads %>% 
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

gts_gar_WGP_cold = gts_gar_WGP_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
gts_gar_WGP_warm = gts_gar_WGP_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# WGP data - matrix differences -------------------------------------------

# ASHN --------------------------------------------------------------------

ashn_WGP_cold_mat = ashn_WGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

ashn_WGP_warm_mat = ashn_WGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

ashn_WGP_diff_mat = ashn_WGP_warm_mat - ashn_WGP_cold_mat

ashn_WGP_diff_mat_long = ashn_WGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

ashn_WGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

ashn_WGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

ashn_WGP_diff_mat_plot = 
  ggplot(ashn_WGP_diff_mat_long, 
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
    title = "A) ASHN WGP warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# MYV ---------------------------------------------------------------------
myv_WGP_cold_mat = myv_WGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

myv_WGP_warm_mat = myv_WGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

myv_WGP_diff_mat = myv_WGP_warm_mat - myv_WGP_cold_mat

myv_WGP_diff_mat_long = myv_WGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

myv_WGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

myv_WGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

myv_WGP_diff_mat_plot = 
  ggplot(myv_WGP_diff_mat_long, 
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
    title = "B) MYV WGP warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# SKR ---------------------------------------------------------------------
skr_WGP_cold_mat = skr_WGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

skr_WGP_warm_mat = skr_WGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

skr_WGP_diff_mat = skr_WGP_warm_mat - skr_WGP_cold_mat

skr_WGP_diff_mat_long = skr_WGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

skr_WGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

skr_WGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

skr_WGP_diff_mat_plot = 
  ggplot(skr_WGP_diff_mat_long, 
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
    title = "C) SKR WGP warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# GTS/GAR -----------------------------------------------------------------
gts_gar_WGP_cold_mat = gts_gar_WGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

gts_gar_WGP_warm_mat = gts_gar_WGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

gts_gar_WGP_diff_mat = gts_gar_WGP_warm_mat - gts_gar_WGP_cold_mat

gts_gar_WGP_diff_mat_long = gts_gar_WGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

gts_gar_WGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

gts_gar_WGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

gts_gar_WGP_diff_mat_plot = 
  ggplot(gts_gar_WGP_diff_mat_long, 
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
    title = "D) GTS/GAR WGP warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# WGP data - Effect size --------------------------------------------------

# ASHN --------------------------------------------------------------------
## Taking the difference matrix and finding per trait effect size

ashn_WGP_trait_effect = ashn_WGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

ashn_trait_effect_plot = 
  ggplot(ashn_WGP_trait_effect, 
         aes(reorder(trait, 
                     effect_size), 
             effect_size))+
  geom_col(col = 'black',
           fill = '#92b4a7') +
  coord_flip()+
  ylim(0.0, 3)+
  labs(y = 'Absolute sum of r-values across axes', 
       title = 'A) ASHN')+
  theme(axis.title.y = element_blank())


# MYV ---------------------------------------------------------------------
myv_WGP_trait_effect = myv_WGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

myv_trait_effect_plot = 
  ggplot(myv_WGP_trait_effect, 
         aes(reorder(trait, 
                     effect_size), 
             effect_size))+
  geom_col(col = 'black',
           fill = '#92b4a7') +
  coord_flip()+
  ylim(0.0, 3)+
  labs(y = 'Absolute sum of r-values across axes', 
       title = 'B) MYV')+
  theme(axis.title.y = element_blank())

# SKR ---------------------------------------------------------------------
skr_WGP_trait_effect = skr_WGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

skr_trait_effect_plot = 
  ggplot(skr_WGP_trait_effect, 
         aes(reorder(trait, 
                     effect_size), 
             effect_size))+
  geom_col(col = 'black',
           fill = '#92b4a7') +
  coord_flip()+
  ylim(0.0, 3)+
  labs(y = 'Absolute sum of r-values across axes', 
       title = 'C) SKR')+
  theme(axis.title.y = element_blank())

# GTS/GAR -----------------------------------------------------------------

gts_gar_WGP_trait_effect = gts_gar_WGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

gts_gar_trait_effect_plot = 
  ggplot(gts_gar_WGP_trait_effect, 
         aes(reorder(trait, 
                     effect_size), 
             effect_size))+
  geom_col(col = 'black',
           fill = '#92b4a7') +
  coord_flip()+
  ylim(0.0, 3)+
  labs(y = 'Absolute sum of r-values across axes', 
       title = 'D) GTS/GAR')+
  theme(axis.title.y = element_blank())

# WGP data - Permutation test RV coefficient ------------------------------

# ASHN --------------------------------------------------------------------

ashn_WGP_PCA

# ashn_trait_cols = ashn_WGP_PCA %>% 
#   dplyr::select(7:39) %>% 
#   names()
# 
# ashn_pca_cols = ashn_WGP_PCA %>% 
#   dplyr::select(40:45) %>% 
#   names()

ashn_WGP_perm_test_data = ashn_WGP_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)

ashn_WGP_obs = RV(ashn_WGP_cold_mat, 
                      ashn_WGP_warm_mat)
set.seed(123)

ashn_WGP_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- ashn_WGP_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = ashn_WGP_pca_cols,
        trait_cols = ashn_WGP_trait_cols)) %>% 
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

ashn_WGP_delta_RV = mean(ashn_WGP_null_rv) - ashn_WGP_obs

ashn_WGP_z = (ashn_WGP_obs - mean(ashn_WGP_null_rv)) / sd(ashn_WGP_null_rv)

ashn_WGP_rv_df = data.frame(ashn_WGP_null_rv = ashn_WGP_null_rv)

ashn_WGP_RV_perm_plot = 
  ggplot(ashn_WGP_rv_df, 
         aes(x = ashn_WGP_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = ashn_WGP_obs), 
             color = "#92b4a7", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "A) ASHN WGP",
    subtitle = paste0("Observed RV = ", round(ashn_WGP_obs, 3)))


# MYV ---------------------------------------------------------------------
myv_WGP_PCA

myv_WGP_trait_cols = myv_WGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

myv_WGP_pca_cols = myv_WGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

myv_WGP_perm_test_data = myv_WGP_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)

myv_WGP_obs = RV(myv_WGP_cold_mat, 
                     myv_WGP_warm_mat)
set.seed(123)

myv_WGP_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- myv_WGP_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = myv_WGP_pca_cols,
        trait_cols = myv_WGP_trait_cols)) %>% 
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

myv_WGP_delta_RV = mean(myv_WGP_null_rv) - myv_WGP_obs

myv_WGP_z = (myv_WGP_obs - mean(myv_WGP_null_rv)) / sd(myv_WGP_null_rv)

myv_WGP_rv_df = data.frame(myv_WGP_null_rv = myv_WGP_null_rv)

myv_WGP_RV_perm_plot = 
  ggplot(myv_WGP_rv_df, 
         aes(x = myv_WGP_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = myv_WGP_obs), 
             color = "#92b4a7", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "B) MYV WGP",
    subtitle = paste0("Observed RV = ", round(myv_WGP_obs, 3)))


# SKR ---------------------------------------------------------------------
skr_WGP_PCA

skr_WGP_trait_cols = skr_WGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

skr_WGP_pca_cols = skr_WGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

skr_WGP_perm_test_data = skr_WGP_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)

skr_WGP_obs = RV(skr_WGP_cold_mat, 
                     skr_WGP_warm_mat)
set.seed(123)

skr_WGP_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- skr_WGP_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = skr_WGP_pca_cols,
        trait_cols = skr_WGP_trait_cols)) %>% 
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

skr_WGP_delta_RV = mean(skr_WGP_null_rv) - skr_WGP_obs

skr_WGP_z = (skr_WGP_obs - mean(skr_WGP_null_rv)) / sd(skr_WGP_null_rv)

skr_WGP_rv_df = data.frame(skr_WGP_null_rv = skr_WGP_null_rv)

skr_WGP_RV_perm_plot = 
  ggplot(skr_WGP_rv_df, 
         aes(x = skr_WGP_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = skr_WGP_obs), 
             color = "#92b4a7", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "C) SKR WGP",
    subtitle = paste0("Observed RV = ", round(skr_WGP_obs, 3)))



# GTS/GAR -----------------------------------------------------------------
gts_gar_WGP_PCA

gts_gar_WGP_trait_cols = gts_gar_WGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

gts_gar_WGP_pca_cols = gts_gar_WGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

gts_gar_WGP_perm_test_data = gts_gar_WGP_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)

gts_gar_WGP_obs = RV(gts_gar_WGP_cold_mat, 
                         gts_gar_WGP_warm_mat)
set.seed(123)

gts_gar_WGP_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- gts_gar_WGP_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = gts_gar_WGP_pca_cols,
        trait_cols = gts_gar_WGP_trait_cols)) %>% 
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

gts_gar_WGP_delta_RV = mean(gts_gar_WGP_null_rv) - gts_gar_WGP_obs

gts_gar_WGP_z = (gts_gar_WGP_obs - mean(gts_gar_WGP_null_rv)) / sd(gts_gar_WGP_null_rv)

gts_gar_WGP_rv_df = data.frame(gts_gar_WGP_null_rv = gts_gar_WGP_null_rv)

gts_gar_WGP_RV_perm_plot = 
  ggplot(gts_gar_WGP_rv_df, 
         aes(x = gts_gar_WGP_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = gts_gar_WGP_obs), 
             color = "#92b4a7", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "D) GTS/GAR WGP",
    subtitle = paste0("Observed RV = ", round(gts_gar_WGP_obs, 3)))






# TGP data ----------------------------------------------------------------

TGP_data = read_csv("TGP_TRAITS_SCALED_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17)


# TGP_trait_mat <- TGP_data %>%
#   dplyr::select(jaw_length:caudal2_15_17) %>% 
#   # dplyr::select(starts_with('value'))
#   # dplyr::select(all_of(trait_names)) %>%
#   as.matrix()
# 
# 
ashn_TGP_data = TGP_data %>% 
  filter(Lake_morph %in% c('ASHNC', 
                           'ASHNW'))

ashn_TGP_trait_mat = ashn_TGP_data %>% 
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  as.matrix()

myv_TGP_data = F2_data %>% 
  filter(Lake_morph %in% c('MYVC', 
                           'MYVW'))

myv_TGP_trait_mat = myv_TGP_data %>% 
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  as.matrix()

skr_TGP_data = F2_data %>% 
  filter(Lake_morph %in% c('SKRC', 
                           'SKRW'))
skr_TGP_trait_mat = skr_TGP_data %>% 
  dplyr::select(jaw_length:caudal2_15_17) %>%  
  as.matrix()

gts_gar_TGP_data = F2_data %>% 
  filter(Lake_morph %in% c('CSWYC', 
                           'GTSW'))
gts_gar_TGP_trait_mat = gts_gar_TGP_data %>% 
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  as.matrix()


# TGP PCA -------------------------------------------------------------
ashn_TGP_pca = prcomp(ashn_TGP_trait_mat, 
                      center = F, 
                      scale. = F)
paran(ashn_TGP_trait_mat, 
      iterations = 1000, 
      graph = T)
## 6 pc axes retained

ashn_TGP_scores = as.data.frame(ashn_TGP_pca$x[,1:6])
ashn_TGP_PCA = bind_cols(ashn_TGP_data, 
                         ashn_TGP_scores)


myv_TGP_pca = prcomp(myv_TGP_trait_mat, 
                     center = F, 
                     scale. = F)
paran(myv_TGP_trait_mat, 
      iterations = 1000, 
      graph = T)
## 7 pc axes retained

myv_TGP_scores = as.data.frame(myv_TGP_pca$x[,1:7])
myv_TGP_PCA = bind_cols(myv_TGP_data, 
                        myv_TGP_scores)


skr_TGP_pca = prcomp(skr_TGP_trait_mat, 
                     center = F, 
                     scale. = F)
paran(skr_TGP_trait_mat, 
      iterations = 1000, 
      graph = T)
## 7 pc axes retained

skr_TGP_scores = as.data.frame(skr_TGP_pca$x[,1:7])
skr_TGP_PCA = bind_cols(skr_TGP_data, 
                        skr_TGP_scores)


gts_gar_TGP_pca = prcomp(gts_gar_TGP_trait_mat, 
                         center = F, 
                         scale. = F)
paran(gts_gar_TGP_trait_mat, 
      iterations = 1000, 
      graph = T)
## 7 pc axes retained

gts_gar_TGP_scores = as.data.frame(gts_gar_TGP_pca$x[,1:7])
gts_gar_TGP_PCA = bind_cols(gts_gar_TGP_data, 
                            gts_gar_TGP_scores)

# TGP data - trait/pc corelations -------------------------------------------

# ASHN --------------------------------------------------------------------

ashn_TGP_trait_cols = ashn_TGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

ashn_TGP_pca_cols = ashn_TGP_PCA %>% 
  dplyr::select(40:45) %>% 
  names()

ashn_TGP_WC_loads = ashn_TGP_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = ashn_TGP_pca_cols, 
    trait_cols = ashn_TGP_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

ashn_TGP_WC_loads %>% 
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

ashn_TGP_cold = ashn_TGP_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
ashn_TGP_warm = ashn_TGP_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# MYV ---------------------------------------------------------------------
myv_TGP_trait_cols = myv_TGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

myv_TGP_pca_cols = myv_TGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

myv_TGP_WC_loads = myv_TGP_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = myv_TGP_pca_cols, 
    trait_cols = myv_TGP_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

myv_TGP_WC_loads %>% 
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

myv_TGP_cold = myv_TGP_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
myv_TGP_warm = myv_TGP_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# SKR ---------------------------------------------------------------------
skr_TGP_trait_cols = skr_TGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

skr_TGP_pca_cols = skr_TGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

skr_TGP_WC_loads = skr_TGP_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = skr_TGP_pca_cols, 
    trait_cols = skr_TGP_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

skr_TGP_WC_loads %>% 
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

skr_TGP_cold = skr_TGP_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
skr_TGP_warm = skr_TGP_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# GTS/GAR -----------------------------------------------------------------
gts_gar_TGP_trait_cols = gts_gar_TGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

gts_gar_TGP_pca_cols = gts_gar_TGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

gts_gar_TGP_WC_loads = gts_gar_TGP_PCA %>% 
  group_by(Morph) %>% 
  nest() %>% 
  mutate(correlations = map(
    data, 
    cor_fun, 
    pc_cols = gts_gar_TGP_pca_cols, 
    trait_cols = gts_gar_TGP_trait_cols)) %>% 
  dplyr::select(-data) %>% 
  unnest(correlations)

## Graphical check

gts_gar_TGP_WC_loads %>% 
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

gts_gar_TGP_cold = gts_gar_TGP_WC_loads %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(pc, trait, r_cold = r)
gts_gar_TGP_warm = gts_gar_TGP_WC_loads %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(pc, trait, r_warm = r)


# TGP data - matrix differences -------------------------------------------

# ASHN --------------------------------------------------------------------

ashn_TGP_cold_mat = ashn_TGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

ashn_TGP_warm_mat = ashn_TGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

ashn_TGP_diff_mat = ashn_TGP_warm_mat - ashn_TGP_cold_mat

ashn_TGP_diff_mat_long = ashn_TGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

ashn_TGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

ashn_TGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

ashn_TGP_diff_mat_plot = 
  ggplot(ashn_TGP_diff_mat_long, 
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
    title = "A) ASHN TGP warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# MYV ---------------------------------------------------------------------
myv_TGP_cold_mat = myv_TGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

myv_TGP_warm_mat = myv_TGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

myv_TGP_diff_mat = myv_TGP_warm_mat - myv_TGP_cold_mat

myv_TGP_diff_mat_long = myv_TGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

myv_TGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

myv_TGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

myv_TGP_diff_mat_plot = 
  ggplot(myv_TGP_diff_mat_long, 
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
    title = "B) MYV TGP warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# SKR ---------------------------------------------------------------------
skr_TGP_cold_mat = skr_TGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

skr_TGP_warm_mat = skr_TGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

skr_TGP_diff_mat = skr_TGP_warm_mat - skr_TGP_cold_mat

skr_TGP_diff_mat_long = skr_TGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

skr_TGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

skr_TGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

skr_TGP_diff_mat_plot = 
  ggplot(skr_TGP_diff_mat_long, 
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
    title = "C) SKR TGP warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# GTS/GAR -----------------------------------------------------------------
gts_gar_TGP_cold_mat = gts_gar_TGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Cold') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

gts_gar_TGP_warm_mat = gts_gar_TGP_WC_loads %>% 
  ungroup() %>% 
  filter(Morph == 'Warm') %>% 
  dplyr::select(trait, pc, r) %>% 
  pivot_wider(names_from = pc, 
              values_from = r) %>% 
  column_to_rownames('trait') %>% 
  as.matrix()

gts_gar_TGP_diff_mat = gts_gar_TGP_warm_mat - gts_gar_TGP_cold_mat

gts_gar_TGP_diff_mat_long = gts_gar_TGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  pivot_longer(cols = -trait, 
               names_to = 'PC', 
               values_to = 'diff')

gts_gar_TGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff > 0.3) %>% 
  arrange(PC, 
          diff)

gts_gar_TGP_diff_mat_long %>% 
  group_by(PC) %>% 
  filter(diff < -0.3) %>% 
  arrange(PC, 
          diff)

gts_gar_TGP_diff_mat_plot = 
  ggplot(gts_gar_TGP_diff_mat_long, 
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
    title = "D) GTS/GAR TGP warm-cold trait correlations") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# TGP data - Effect size --------------------------------------------------

# ASHN --------------------------------------------------------------------
## Taking the difference matrix and finding per trait effect size

ashn_TGP_trait_effect = ashn_TGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

ashn_trait_effect_plot = 
  ggplot(ashn_TGP_trait_effect, 
         aes(reorder(trait, 
                     effect_size), 
             effect_size))+
  geom_col(col = 'black',
           fill = '#d1f0b1') +
  coord_flip()+
  ylim(0.0, 3)+
  labs(y = 'Absolute sum of r-values across axes', 
       title = 'A) ASHN')+
  theme(axis.title.y = element_blank())


# MYV ---------------------------------------------------------------------
myv_TGP_trait_effect = myv_TGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

myv_trait_effect_plot = 
  ggplot(myv_TGP_trait_effect, 
         aes(reorder(trait, 
                     effect_size), 
             effect_size))+
  geom_col(col = 'black',
           fill = '#d1f0b1') +
  coord_flip()+
  ylim(0.0, 3)+
  labs(y = 'Absolute sum of r-values across axes', 
       title = 'B) MYV')+
  theme(axis.title.y = element_blank())

# SKR ---------------------------------------------------------------------
skr_TGP_trait_effect = skr_TGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

skr_trait_effect_plot = 
  ggplot(skr_TGP_trait_effect, 
         aes(reorder(trait, 
                     effect_size), 
             effect_size))+
  geom_col(col = 'black',
           fill = '#d1f0b1') +
  coord_flip()+
  ylim(0.0, 3)+
  labs(y = 'Absolute sum of r-values across axes', 
       title = 'C) SKR')+
  theme(axis.title.y = element_blank())

# GTS/GAR -----------------------------------------------------------------

gts_gar_TGP_trait_effect = gts_gar_TGP_diff_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('trait') %>% 
  mutate(effect_size = rowSums(abs(across(where(is.numeric))))) %>% 
  arrange(desc(effect_size))

gts_gar_trait_effect_plot = 
  ggplot(gts_gar_TGP_trait_effect, 
         aes(reorder(trait, 
                     effect_size), 
             effect_size))+
  geom_col(col = 'black',
           fill = '#d1f0b1') +
  coord_flip()+
  ylim(0.0, 3)+
  labs(y = 'Absolute sum of r-values across axes', 
       title = 'D) GTS/GAR')+
  theme(axis.title.y = element_blank())

# TGP data - Permutation test RV coefficient ------------------------------

# ASHN --------------------------------------------------------------------

ashn_TGP_PCA

# ashn_trait_cols = ashn_TGP_PCA %>% 
#   dplyr::select(7:39) %>% 
#   names()
# 
# ashn_pca_cols = ashn_TGP_PCA %>% 
#   dplyr::select(40:45) %>% 
#   names()

ashn_TGP_perm_test_data = ashn_TGP_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:45)

ashn_TGP_obs = RV(ashn_TGP_cold_mat, 
                  ashn_TGP_warm_mat)
set.seed(123)

ashn_TGP_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- ashn_TGP_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = ashn_TGP_pca_cols,
        trait_cols = ashn_TGP_trait_cols)) %>% 
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

ashn_TGP_delta_RV = mean(ashn_TGP_null_rv) - ashn_TGP_obs

ashn_TGP_z = (ashn_TGP_obs - mean(ashn_TGP_null_rv)) / sd(ashn_TGP_null_rv)

ashn_TGP_rv_df = data.frame(ashn_TGP_null_rv = ashn_TGP_null_rv)

ashn_TGP_RV_perm_plot = 
  ggplot(ashn_TGP_rv_df, 
         aes(x = ashn_TGP_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = ashn_TGP_obs), 
             color = "#d1f0b1", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "A) ASHN TGP",
    subtitle = paste0("Observed RV = ", round(ashn_TGP_obs, 3)))


# MYV ---------------------------------------------------------------------
myv_TGP_PCA

myv_TGP_trait_cols = myv_TGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

myv_TGP_pca_cols = myv_TGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

myv_TGP_perm_test_data = myv_TGP_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)

myv_TGP_obs = RV(myv_TGP_cold_mat, 
                 myv_TGP_warm_mat)
set.seed(123)

myv_TGP_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- myv_TGP_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = myv_TGP_pca_cols,
        trait_cols = myv_TGP_trait_cols)) %>% 
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

myv_TGP_delta_RV = mean(myv_TGP_null_rv) - myv_TGP_obs

myv_TGP_z = (myv_TGP_obs - mean(myv_TGP_null_rv)) / sd(myv_TGP_null_rv)

myv_TGP_rv_df = data.frame(myv_TGP_null_rv = myv_TGP_null_rv)

myv_TGP_RV_perm_plot = 
  ggplot(myv_TGP_rv_df, 
         aes(x = myv_TGP_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = myv_TGP_obs), 
             color = "#d1f0b1", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "B) MYV TGP",
    subtitle = paste0("Observed RV = ", round(myv_TGP_obs, 3)))


# SKR ---------------------------------------------------------------------
skr_TGP_PCA

skr_TGP_trait_cols = skr_TGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

skr_TGP_pca_cols = skr_TGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

skr_TGP_perm_test_data = skr_TGP_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)

skr_TGP_obs = RV(skr_TGP_cold_mat, 
                 skr_TGP_warm_mat)
set.seed(123)

skr_TGP_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- skr_TGP_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = skr_TGP_pca_cols,
        trait_cols = skr_TGP_trait_cols)) %>% 
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

skr_TGP_delta_RV = mean(skr_TGP_null_rv) - skr_TGP_obs

skr_TGP_z = (skr_TGP_obs - mean(skr_TGP_null_rv)) / sd(skr_TGP_null_rv)

skr_TGP_rv_df = data.frame(skr_TGP_null_rv = skr_TGP_null_rv)

skr_TGP_RV_perm_plot = 
  ggplot(skr_TGP_rv_df, 
         aes(x = skr_TGP_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = skr_TGP_obs), 
             color = "#d1f0b1", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "C) SKR TGP",
    subtitle = paste0("Observed RV = ", round(skr_TGP_obs, 3)))



# GTS/GAR -----------------------------------------------------------------
gts_gar_TGP_PCA

gts_gar_TGP_trait_cols = gts_gar_TGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

gts_gar_TGP_pca_cols = gts_gar_TGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

gts_gar_TGP_perm_test_data = gts_gar_TGP_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)

gts_gar_TGP_obs = RV(gts_gar_TGP_cold_mat, 
                     gts_gar_TGP_warm_mat)
set.seed(123)

gts_gar_TGP_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- gts_gar_TGP_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = gts_gar_TGP_pca_cols,
        trait_cols = gts_gar_TGP_trait_cols)) %>% 
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

gts_gar_TGP_delta_RV = mean(gts_gar_TGP_null_rv) - gts_gar_TGP_obs

gts_gar_TGP_z = (gts_gar_TGP_obs - mean(gts_gar_TGP_null_rv)) / sd(gts_gar_TGP_null_rv)

gts_gar_TGP_rv_df = data.frame(gts_gar_TGP_null_rv = gts_gar_TGP_null_rv)

gts_gar_TGP_RV_perm_plot = 
  ggplot(gts_gar_TGP_rv_df, 
         aes(x = gts_gar_TGP_null_rv))+
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = gts_gar_TGP_obs), 
             color = "#d1f0b1", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "D) GTS/GAR TGP",
    subtitle = paste0("Observed RV = ", round(gts_gar_TGP_obs, 3)))



