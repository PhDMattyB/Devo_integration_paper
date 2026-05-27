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

pca = prcomp(F2_orig_trait_mat, 
              center = F, 
              scale. = F)

paran(F2_orig_trait_mat, iterations = 1000, graph = TRUE)

scores = as.data.frame(pca$x[, 1:7])
F2_PCA = bind_cols(F2_data, scores)

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


F2_WC_loads %>% 
  filter(pc == 'PC1')

## check graphically
F2_WC_loads %>%
  ggplot(aes(pc, 
             trait, 
             fill = r)) +
  geom_tile() +
  facet_wrap(~Morph) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    limits = c(-1, 1)) +
  theme(
    axis.text.y = element_text(size = 7))

F2_cold = F2_WC_loads %>%
  filter(Morph == "Cold") %>%
  dplyr::select(pc, trait, r_cold = r)

F2_warm = F2_WC_loads %>%
  filter(Morph == "Warm") %>%
  dplyr::select(pc, trait, r_warm = r)

F2_diff_df <- F2_cold %>%
  inner_join(F2_warm,
    by = c("pc", "trait")) %>%
  mutate(diff = r_warm - r_cold)


ggplot(diff_df, 
       aes(pc,
           trait,
           fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red") 


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

pheatmap(
  diff_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)


F2_trait_effect <- diff_mat %>%
  as.data.frame() %>%
  rownames_to_column("trait") %>%
  mutate(
    effect_size = rowSums(abs(across(where(is.numeric))))
  ) %>%
  arrange(desc(effect_size))

ggplot(F2_trait_effect, 
       aes(reorder(trait, 
                   effect_size), 
           effect_size)) +
  geom_col() +
  coord_flip()

pc_effect <- diff_mat %>%
  as.data.frame() %>%
  summarise(across(everything(), ~ sum(abs(.)))) %>%
  pivot_longer(
    everything(),
    names_to = "PC",
    values_to = "effect"
  )

ggplot(pc_effect, aes(PC, effect)) +
  geom_col()

##how much a trait’s pattern of associations with the 7 PCs differs between ecotypes
trait_contrib <- (F2_warm_mat - F2_cold_mat)^2 %>%
  rowSums() %>%
  sort(decreasing = TRUE)

trait_df <- data.frame(
  trait = names(trait_contrib),
  contribution = as.numeric(trait_contrib)
)

ggplot(trait_df, aes(x = reorder(trait, contribution), y = contribution)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Trait",
    y = "Contribution to ecotype difference",
    title = "Trait-level contribution to divergence in trait–PC structure"
  ) 


allCorrelations(
  F2_cold_mat,
  F2_warm_mat,
  ncomp1 = 7,
  ncomp2 = 7)

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

hist(null_rv)
abline(v = obs, col = "red", lwd = 3)
