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
    title = "Difference in trait–PC loadings (Warm - Cold)") +
  theme(
    axis.text.y = element_text(size = 6))



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
    fill = '#54478c') +
  coord_flip()+
  labs(x = 'Effect size', 
       y = 'Trait')


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
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(aes(xintercept = obs), color = "red", linewidth = 1.2) +
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
    title = "Difference in trait–PC loadings (Warm - Cold)") +
  theme(
    axis.text.y = element_text(size = 6))


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
           fill = '#7785ac') +
  coord_flip()+
  labs(x = 'Effect size', 
       y = 'Trait')


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

WGP_RV_Perm_plot = ggplot(WGP_rv_df, aes(x = WGP_null_rv)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(aes(xintercept = WGP_obs), color = "red", linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "Permutation test of matrix similarity (RV)",
    subtitle = paste0("Observed RV = ", round(WGP_obs, 3))
  )


# WGP data ----------------------------------------------------------------
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


# WGP data - PCA ----------------------------------------------------------
TGP_pca = prcomp(TGP_trait_mat, 
                 center = F, 
                 scale. = F)

paran(TGP_trait_mat, iterations = 1000, graph = TRUE)

TGP_scores = as.data.frame(TGP_pca$x[, 1:7])
TGP_PCA = bind_cols(TGP_data, TGP_scores)


# WGP data - trait/pc cor -------------------------------------------------
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


# WGP data - matrix differences -------------------------------------------
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
    title = "Difference in trait–PC loadings (Warm - Cold)") +
  theme(
    axis.text.y = element_text(size = 6))


# WGP data - Effect size --------------------------------------------------
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
           fill = '#7785ac') +
  coord_flip()+
  labs(x = 'Effect size', 
       y = 'Trait')


# WGP data - Permutation test RV coefficient ------------------------------

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
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(aes(xintercept = TGP_obs), color = "red", linewidth = 1.2) +
  labs(
    x = "RV (null distribution)",
    y = "Frequency",
    title = "Permutation test of matrix similarity (RV)",
    subtitle = paste0("Observed RV = ", round(TGP_obs, 3))
  )


