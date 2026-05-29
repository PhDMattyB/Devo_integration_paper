##############################
## PC loading matrices
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
# F2 data -----------------------------------------------------------------

F2_data = read_csv("F2_Original_univariate_traits_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17)


F2_orig_trait_mat <- F2_data %>%
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  # dplyr::select(starts_with('value'))
  # dplyr::select(all_of(trait_names)) %>%
  as.matrix()


# F2 data - PCA -----------------------------------------------------------

pca = prcomp(F2_orig_trait_mat, 
              center = F, 
              scale. = F)

paran(F2_orig_trait_mat, iterations = 1000, graph = TRUE)

scores = as.data.frame(pca$x[, 1:7])
F2_PCA = bind_cols(F2_data, scores)



# F2 data - trait/PC correlations -----------------------------------------


trait_cols = F2_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

pc_cols = F2_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

F2_WC_loads = F2_PCA %>%
  group_by(Morph) %>%
  nest() %>%
  mutate(
    correlations = map(
      data,
      cor_fun,
      pc_cols = pc_cols,
      trait_cols = trait_cols)) %>%
  dplyr::select(-data) %>%
  unnest(correlations)

## check graphically
F2_WC_loads %>%
  ggplot(aes(pc, 
             trait, 
             fill = r)) +
  geom_tile() +
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

F2_cold = F2_WC_loads %>%
  filter(Morph == "Cold") %>%
  dplyr::select(pc, trait, r_cold = r)

F2_warm = F2_WC_loads %>%
  filter(Morph == "Warm") %>%
  dplyr::select(pc, trait, r_warm = r)


# F2 data - matrix differences --------------------------------------------


F2_cold_mat <- F2_WC_loads %>%
  ungroup() %>% 
  filter(Morph == "Cold") %>%
  dplyr::select(trait, pc, r) %>%
  pivot_wider(names_from = pc,
    values_from = r) %>%
  column_to_rownames("trait") %>%
  as.matrix()


F2_warm_mat <- F2_WC_loads %>%
  ungroup() %>% 
  filter(Morph == "Warm") %>%
  dplyr::select(trait, pc, r) %>%
  ungroup() %>% 
  pivot_wider(names_from = pc,
              values_from = r) %>%
  column_to_rownames("trait") %>%
  as.matrix()
F2_diff_mat <- F2_warm_mat - F2_cold_mat


F2_diff_mat_long <- F2_diff_mat %>%
  as.data.frame() %>%
  rownames_to_column("trait") %>%
  pivot_longer(
    cols = -trait,
    names_to = "PC",
    values_to = "diff"
  )

F2_Ecotype_diff_mat_plot = ggplot(F2_diff_mat_long, 
       aes(x = PC, 
           y = trait, 
           fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#277da1",
    mid = "white",
    high = "#f94144") +
  labs(
    x = "PC axis",
    y = "Trait",
    fill = "Δ loading",
    title = "A) F2 phenotype") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')



# F2 data - Trait effect size ---------------------------------------------

## Taking the difference matrix and finding per trait effect size
F2_trait_effect <- F2_diff_mat %>%
  as.data.frame() %>%
  rownames_to_column("trait") %>%
  mutate(
    effect_size = rowSums(abs(across(where(is.numeric))))
  ) %>%
  arrange(desc(effect_size))

F2_trait_effect_plot = ggplot(F2_trait_effect, 
       aes(reorder(trait, 
                   effect_size), 
           effect_size)) +
  geom_col(col = 'black',
    fill = '#81667a') +
  ylim(0.0,1.2)+
  labs(y = 'Absolute sum of r-values across axes')+
  coord_flip()+
  theme(axis.title.y = element_blank())


# allCorrelations(
#   F2_cold_mat,
#   F2_warm_mat,
#   ncomp1 = 7,
#   ncomp2 = 7)



# F2 data - Permutation testing RV coefficient ----------------------------

perm_test_data = F2_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)


obs <- RV(F2_cold_mat, F2_warm_mat)
set.seed(123)

null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = pc_cols,
        trait_cols = trait_cols)) %>% 
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


delta_RV = mean(null_rv) - obs

z = (obs - mean(null_rv)) / sd(null_rv)

F2_rv_df <- data.frame(null_rv = null_rv)

F2_RV_Perm_plot = ggplot(F2_rv_df, aes(x = null_rv)) +
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = obs), 
             color = "#81667a", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "Permutation test of matrix similarity (RV)",
    subtitle = paste0("Observed RV = ", round(obs, 3))
  )



# WGP data ----------------------------------------------------------------
WGP_data = read_csv("WGP_TRAITS_SCALED_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17)


WGP_trait_mat <- WGP_data %>%
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  # dplyr::select(starts_with('value'))
  # dplyr::select(all_of(trait_names)) %>%
  as.matrix()


# WGP data - PCA ----------------------------------------------------------
WGP_pca = prcomp(WGP_trait_mat, 
             center = F, 
             scale. = F)

paran(WGP_trait_mat, iterations = 1000, graph = TRUE)

WGP_scores = as.data.frame(WGP_pca$x[, 1:7])
WGP_PCA = bind_cols(WGP_data, WGP_scores)


# WGP data - trait/pc cor -------------------------------------------------
WGP_trait_cols = WGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

WGP_pc_cols = WGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

WGP_WC_loads = WGP_PCA %>%
  group_by(Morph) %>%
  nest() %>%
  mutate(
    correlations = map(
      data,
      cor_fun,
      pc_cols = WGP_pc_cols,
      trait_cols = WGP_trait_cols)) %>%
  dplyr::select(-data) %>%
  unnest(correlations)

## check graphically
WGP_WC_loads %>%
  ggplot(aes(pc, 
             trait, 
             fill = r)) +
  geom_tile() +
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

WGP_cold = WGP_WC_loads %>%
  filter(Morph == "Cold") %>%
  dplyr::select(pc, trait, r_cold = r)

WGP_warm = WGP_WC_loads %>%
  filter(Morph == "Warm") %>%
  dplyr::select(pc, trait, r_warm = r)


# WGP data - matrix differences -------------------------------------------
WGP_cold_mat <- WGP_WC_loads %>%
  ungroup() %>% 
  filter(Morph == "Cold") %>%
  dplyr::select(trait, pc, r) %>%
  pivot_wider(names_from = pc,
              values_from = r) %>%
  column_to_rownames("trait") %>%
  as.matrix()


WGP_warm_mat <- WGP_WC_loads %>%
  ungroup() %>% 
  filter(Morph == "Warm") %>%
  dplyr::select(trait, pc, r) %>%
  ungroup() %>% 
  pivot_wider(names_from = pc,
              values_from = r) %>%
  column_to_rownames("trait") %>%
  as.matrix()
WGP_diff_mat <- WGP_warm_mat - WGP_cold_mat


WGP_diff_mat_long <- WGP_diff_mat %>%
  as.data.frame() %>%
  rownames_to_column("trait") %>%
  pivot_longer(
    cols = -trait,
    names_to = "PC",
    values_to = "diff"
  )
WGP_Ecotype_diff_mat_plot = ggplot(WGP_diff_mat_long, 
                                  aes(x = PC, 
                                      y = trait, 
                                      fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#277da1",
    mid = "white",
    high = "#f94144") +
  labs(
    x = "PC axis",
    y = "Trait",
    fill = "Δ loading",
    title = "B) Within-generational plasticity") +
  theme(
    axis.text.y = element_text(size = 6),
    axis.title = element_blank(), 
    legend.position = 'none')


# WGP data - Effect size --------------------------------------------------
## Taking the difference matrix and finding per trait effect size
WGP_trait_effect <- WGP_diff_mat %>%
  as.data.frame() %>%
  rownames_to_column("trait") %>%
  mutate(
    effect_size = rowSums(abs(across(where(is.numeric))))) %>%
  arrange(desc(effect_size))

WGP_trait_effect_plot = ggplot(WGP_trait_effect, 
                              aes(reorder(trait, 
                                          effect_size), 
                                  effect_size)) +
  geom_col(col = 'black',
           fill = '#92b4a7') +
  coord_flip()+
  ylim(0.0, 1.2)+
  labs(y = 'Absolute sum of r-values across axes')+
  theme(axis.title.y = element_blank())


# WGP data - Permutation test RV coefficient ------------------------------

WGP_perm_test_data = WGP_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)


WGP_obs <- RV(WGP_cold_mat, WGP_warm_mat)
set.seed(123)

WGP_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- WGP_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = pc_cols,
        trait_cols = trait_cols)) %>% 
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


WGP_delta_RV = mean(WGP_null_rv) - WGP_obs

WGP_z = (WGP_obs - mean(WGP_null_rv)) / sd(WGP_null_rv)

WGP_rv_df <- data.frame(WGP_null_rv = WGP_null_rv)

WGP_RV_Perm_plot = ggplot(WGP_rv_df, 
                          aes(x = WGP_null_rv)) +
  geom_histogram(bins = 30, 
                 fill = "#adb5bd", 
                 color = "black") +
  geom_vline(aes(xintercept = WGP_obs), 
             color = "#92b4a7", 
             linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "Permutation test of matrix similarity (RV)",
    subtitle = paste0("Observed RV = ", round(WGP_obs, 3))
  )


# TGP data ----------------------------------------------------------------
TGP_data = read_csv("TGP_TRAITS_SCALED_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17)


TGP_trait_mat <- TGP_data %>%
  dplyr::select(jaw_length:caudal2_15_17) %>% 
  # dplyr::select(starts_with('value'))
  # dplyr::select(all_of(trait_names)) %>%
  as.matrix()


# TGP data - PCA ----------------------------------------------------------
TGP_pca = prcomp(TGP_trait_mat, 
                 center = F, 
                 scale. = F)

paran(TGP_trait_mat, iterations = 1000, graph = TRUE)

TGP_scores = as.data.frame(TGP_pca$x[, 1:7])
TGP_PCA = bind_cols(TGP_data, TGP_scores)


# TGP data - trait/pc cor -------------------------------------------------
TGP_trait_cols = TGP_PCA %>% 
  dplyr::select(7:39) %>% 
  names()

TGP_pc_cols = TGP_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

TGP_WC_loads = TGP_PCA %>%
  group_by(Morph) %>%
  nest() %>%
  mutate(
    correlations = map(
      data,
      cor_fun,
      pc_cols = TGP_pc_cols,
      trait_cols = TGP_trait_cols)) %>%
  dplyr::select(-data) %>%
  unnest(correlations)

## check graphically
TGP_WC_loads %>%
  ggplot(aes(pc, 
             trait, 
             fill = r)) +
  geom_tile() +
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

TGP_cold = TGP_WC_loads %>%
  filter(Morph == "Cold") %>%
  dplyr::select(pc, trait, r_cold = r)

TGP_warm = TGP_WC_loads %>%
  filter(Morph == "Warm") %>%
  dplyr::select(pc, trait, r_warm = r)


# TGP data - matrix differences -------------------------------------------
TGP_cold_mat <- TGP_WC_loads %>%
  ungroup() %>% 
  filter(Morph == "Cold") %>%
  dplyr::select(trait, pc, r) %>%
  pivot_wider(names_from = pc,
              values_from = r) %>%
  column_to_rownames("trait") %>%
  as.matrix()


TGP_warm_mat <- TGP_WC_loads %>%
  ungroup() %>% 
  filter(Morph == "Warm") %>%
  dplyr::select(trait, pc, r) %>%
  ungroup() %>% 
  pivot_wider(names_from = pc,
              values_from = r) %>%
  column_to_rownames("trait") %>%
  as.matrix()
TGP_diff_mat <- TGP_warm_mat - TGP_cold_mat


TGP_diff_mat_long <- TGP_diff_mat %>%
  as.data.frame() %>%
  rownames_to_column("trait") %>%
  pivot_longer(
    cols = -trait,
    names_to = "PC",
    values_to = "diff"
  )
TGP_Ecotype_diff_mat_plot = ggplot(TGP_diff_mat_long, 
                                   aes(x = PC, 
                                       y = trait, 
                                       fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#277da1",
    mid = "white",
    high = "#f94144") +
  labs(
    x = "PC axis",
    y = "Trait",
    fill = "Δ loading",
    title = "C) Trans-generational plasticity") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# TGP data - Effect size --------------------------------------------------
## Taking the difference matrix and finding per trait effect size
TGP_trait_effect <- TGP_diff_mat %>%
  as.data.frame() %>%
  rownames_to_column("trait") %>%
  mutate(
    effect_size = rowSums(abs(across(where(is.numeric))))) %>%
  arrange(desc(effect_size))

TGP_trait_effect_plot = ggplot(TGP_trait_effect, 
                               aes(reorder(trait, 
                                           effect_size), 
                                   effect_size)) +
  geom_col(col = 'black',
           fill = '#d1f0b1') +
  coord_flip()+
  ylim(0.0, 1.2)+
  labs(y = 'Absolute sum of r-values across axes')+
  theme(axis.title.y = element_blank())


# TGP data - Permutation test RV coefficient ------------------------------

TGP_perm_test_data = TGP_PCA %>% 
  dplyr::select(Morph, 
                7:39, 
                40:46)


TGP_obs <- RV(TGP_cold_mat, TGP_warm_mat)
set.seed(123)

TGP_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- TGP_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = pc_cols,
        trait_cols = trait_cols)) %>% 
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


TGP_delta_RV = mean(TGP_null_rv) - TGP_obs

TGP_z = (TGP_obs - mean(TGP_null_rv)) / sd(TGP_null_rv)

TGP_rv_df <- data.frame(TGP_null_rv = TGP_null_rv)

TGP_RV_Perm_plot = ggplot(TGP_rv_df, aes(x = TGP_null_rv)) +
  geom_histogram(bins = 30, fill = "#adb5bd", color = "black") +
  geom_vline(aes(xintercept = TGP_obs), color = "#d1f0b1", linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "Permutation test of matrix similarity (RV)",
    subtitle = paste0("Observed RV = ", round(TGP_obs, 3))
  )





# Compare effect sizes to diff mats ---------------------------------------

F2_trait_effect_0.75 = F2_trait_effect %>% 
  as_tibble() %>% 
  filter(effect_size >= 0.75) %>% 
  dplyr::select(trait, 
                effect_size)
WGP_trait_effect_0.75 = WGP_trait_effect %>% 
  as_tibble() %>% 
  filter(effect_size >= 0.75) %>% 
  dplyr::select(trait, 
                effect_size)
TGP_trait_effect_0.75 = TGP_trait_effect %>% 
  as_tibble() %>% 
  filter(effect_size >= 0.75) %>% 
  dplyr::select(trait, 
                effect_size)


log_rats = inner_join(F2_trait_effect, 
           WGP_trait_effect, 
           by = 'trait') %>% 
  inner_join(., 
             TGP_trait_effect, 
             by = 'trait') %>% 
  dplyr::select(-starts_with('PC')) %>% 
  rename(F2_eff_size = 2, 
         WGP_eff_size = 3, 
         TGP_eff_size = 4)%>%
  mutate(
    log_rat_WGP_vs_F2 = log2(WGP_eff_size / F2_eff_size),
    log_rat_TGP_vs_F2 = log2(TGP_eff_size / F2_eff_size), 
    log_rat_WGP_vs_TGP = log2(WGP_eff_size/TGP_eff_size)) %>% 
  dplyr::select(trait, contains('log_rat'))

log_rats %>% 
  filter(log_rat_WGP_vs_TGP > 0.5)

log_rats %>% 
  filter(log_rat_WGP_vs_TGP < 0.5, 
         log_rat_WGP_vs_TGP > 0.1)

log_rats %>% 
  filter(log_rat_WGP_vs_TGP < -0.5)

log_rats %>% 
  filter(log_rat_WGP_vs_TGP < -0.1, 
         log_rat_WGP_vs_TGP > -0.5)

log_rats %>% 
  filter(log_rat_WGP_vs_TGP < 0.1, 
         log_rat_WGP_vs_TGP > -0.1)


log_rats_plot = log_rats %>%
  dplyr::select(trait,
                log_rat_WGP_vs_TGP) %>%
  pivot_longer(
    -trait,
    names_to = "comparison",
    values_to = "pct_diff") %>%
  ggplot(aes(comparison, trait, fill = pct_diff)) +
  geom_tile(col = 'black') +
  scale_fill_gradient2(
    low = "#47126b",
    mid = "white",
    high = "#ea698b"
  ) 

ggsave('log_ratio_WGP_pink_TGP_purple.tiff', 
       plot = log_rats_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 10, 
       height = 15)

## intersections between data sets

inner_join(F2_trait_effect_0.75, 
          WGP_trait_effect_0.75, 
          by = 'trait')

inner_join(F2_trait_effect_0.75, 
           TGP_trait_effect_0.75, 
           by = 'trait')

inner_join(WGP_trait_effect_0.75, 
           TGP_trait_effect_0.75, 
           by = 'trait')

anti_join(F2_trait_effect_0.75, 
          WGP_trait_effect_0.75, 
          by = 'trait')

anti_join(WGP_trait_effect_0.75,
          F2_trait_effect_0.75, 
          by = 'trait')







# wild plots --------------------------------------------------------------

# wild data ----------------------------------------------------------------
wild_data = read_csv("WILD_SCALED_FIXED_27.05.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                rowname, 
                Lake, 
                Morph, 
                5:39)


wild_trait_mat <- wild_data %>%
  dplyr::select(5:39) %>% 
  # dplyr::select(starts_with('value'))
  # dplyr::select(all_of(trait_names)) %>%
  as.matrix()


# wild data - PCA ----------------------------------------------------------
wild_pca = prcomp(wild_trait_mat, 
                 center = F, 
                 scale. = F)

paran(wild_trait_mat, iterations = 1000, graph = TRUE)

wild_scores = as.data.frame(wild_pca$x[, 1:7])
wild_PCA = bind_cols(wild_data, wild_scores)


# wild data - trait/pc cor -------------------------------------------------
wild_trait_cols = wild_PCA %>% 
  dplyr::select(5:39) %>% 
  names()

wild_pc_cols = wild_PCA %>% 
  dplyr::select(40:46) %>% 
  names()

wild_WC_loads = wild_PCA %>%
  group_by(Morph) %>%
  nest() %>%
  mutate(
    correlations = map(
      data,
      cor_fun,
      pc_cols = wild_pc_cols,
      trait_cols = wild_trait_cols)) %>%
  dplyr::select(-data) %>%
  unnest(correlations)

## check graphically
wild_WC_loads %>%
  ggplot(aes(pc, 
             trait, 
             fill = r)) +
  geom_tile() +
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

wild_cold = wild_WC_loads %>%
  filter(Morph == "Cold") %>%
  dplyr::select(pc, trait, r_cold = r)

wild_warm = wild_WC_loads %>%
  filter(Morph == "Warm") %>%
  dplyr::select(pc, trait, r_warm = r)


# wild data - matrix differences -------------------------------------------
wild_cold_mat <- wild_WC_loads %>%
  ungroup() %>% 
  filter(Morph == "Cold") %>%
  dplyr::select(trait, pc, r) %>%
  pivot_wider(names_from = pc,
              values_from = r) %>%
  column_to_rownames("trait") %>%
  as.matrix()


wild_warm_mat <- wild_WC_loads %>%
  ungroup() %>% 
  filter(Morph == "Warm") %>%
  dplyr::select(trait, pc, r) %>%
  ungroup() %>% 
  pivot_wider(names_from = pc,
              values_from = r) %>%
  column_to_rownames("trait") %>%
  as.matrix()
wild_diff_mat <- wild_warm_mat - wild_cold_mat


wild_diff_mat_long <- wild_diff_mat %>%
  as.data.frame() %>%
  rownames_to_column("trait") %>%
  pivot_longer(
    cols = -trait,
    names_to = "PC",
    values_to = "diff"
  )
wild_Ecotype_diff_mat_plot = ggplot(wild_diff_mat_long, 
                                   aes(x = PC, 
                                       y = trait, 
                                       fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#277da1",
    mid = "white",
    high = "#f94144") +
  labs(
    x = "PC axis",
    y = "Trait",
    fill = "Δ loading",
    title = "A) wild - grandparental generation") +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.title = element_blank(), 
    legend.position = 'none')


# wild data - Effect size --------------------------------------------------
## Taking the difference matrix and finding per trait effect size
wild_trait_effect <- wild_diff_mat %>%
  as.data.frame() %>%
  rownames_to_column("trait") %>%
  mutate(
    effect_size = rowSums(abs(across(where(is.numeric))))) %>%
  arrange(desc(effect_size))

wild_trait_effect_plot = ggplot(wild_trait_effect, 
                               aes(reorder(trait, 
                                           effect_size), 
                                   effect_size)) +
  geom_col(col = 'black',
           fill = '#480355') +
  coord_flip()+
  ylim(0.0, 1.8)+
  labs(y = 'Absolute sum of r-values across axes')+
  theme(axis.title.y = element_blank())


# wild data - Permutation test RV coefficient ------------------------------

wild_perm_test_data = wild_PCA %>% 
  dplyr::select(Morph, 
                5:39, 
                40:46)


wild_obs <- RV(wild_cold_mat, wild_warm_mat)
set.seed(123)

wild_null_rv <- replicate(1000, {
  
  # 1. permute group labels in RAW data
  permuted <- wild_perm_test_data %>%
    mutate(Morph = sample(Morph))
  
  # 2. recompute correlations
  perm_results <- permuted %>%
    group_by(Morph) %>%
    nest() %>%
    mutate(
      correlations = map(
        data,
        cor_fun,
        pc_cols = pc_cols,
        trait_cols = trait_cols)) %>% 
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


wild_delta_RV = mean(wild_null_rv) - wild_obs

wild_z = (wild_obs - mean(wild_null_rv)) / sd(wild_null_rv)

wild_rv_df <- data.frame(wild_null_rv = wild_null_rv)

wild_RV_Perm_plot = ggplot(wild_rv_df, aes(x = wild_null_rv)) +
  geom_histogram(bins = 30, fill = "#adb5bd", color = "black") +
  geom_vline(aes(xintercept = wild_obs), color = "#d1f0b1", linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "Permutation test of matrix similarity (RV)",
    subtitle = paste0("Observed RV = ", round(wild_obs, 3))
  )




# Combo plots -------------------------------------------------------------

eco_diff_combo = F2_Ecotype_diff_mat_plot+WGP_Ecotype_diff_mat_plot+TGP_Ecotype_diff_mat_plot


trait_effect_combo = F2_trait_effect_plot+WGP_trait_effect_plot+TGP_trait_effect_plot


RV_perm_plot_combo = F2_RV_Perm_plot+WGP_RV_Perm_plot+TGP_RV_Perm_plot

combo_of_kings = eco_diff_combo/trait_effect_combo/RV_perm_plot_combo

ggsave('Patterns_Integration_Fixed.svg', 
       plot = combo_of_kings, 
       dpi = 'retina', 
       units = 'cm', 
       width = 40, 
       height = 30)  

