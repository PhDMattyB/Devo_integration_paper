##############################
## RRPP trajectory analysis 
##
## Matt Brachmann (PhDMattyB)
##
## 27.05.2026
##
##############################


setwd("~/Parsons_postdoc/Stickleback_Morphometric_data/Updated Landmarks/")

library(tidyverse)
library(paran)
library(RRPP)



# STRAIGHT TRAITS HOMIE ---------------------------------------------------

F2_data = read_csv("F2_Original_univariate_traits_FIXED_11.02.2026.csv") %>% 
  dplyr::select(Lake_morph, 
                Ecotype_pair, 
                lake_morph_Pair_Full_Temp, 
                Morph, 
                Offspring_temp, 
                Parent_temp, 
                jaw_length:caudal2_15_17) %>% 
  unite(col = off_parent_temp, 
        c('Offspring_temp', 
          'Parent_temp'), 
        sep = '_', 
        remove = F)

F2_traits = F2_data %>% 
  dplyr::select(8:40) %>% 
  as.matrix()

F2_RRPP_mod <- lm.rrpp(F2_traits ~ Morph*Offspring_temp*Parent_temp , 
                       SS.type = "I", 
                       data = F2_data, 
                       print.progress = FALSE, 
                       iter = 999,
                       verbose = TRUE) 
summary(F2_RRPP_mod, 
        formula = FALSE)
anova(F2_RRPP_mod) 

coef(F2_RRPP_mod, test = TRUE)

betaTest(F2_RRPP_mod)

F2_off_traj = trajectory.analysis(F2_RRPP_mod, 
                                  fit.null = NULL, 
                                  groups = F2_data$Morph, 
                                  traj.pts = F2_data$Offspring_temp, 
                                  pca = TRUE, 
                                  print.progress = FALSE)

summary(F2_off_traj, attribute = "MD", angle.type = "deg", 
        show.trajectories = T)
summary(F2_off_traj, attribute = "TC", angle.type = "deg")

F2_parent_traj = trajectory.analysis(F2_RRPP_mod, 
                                     fit.null = NULL, 
                                     groups = F2_data$Morph, 
                                     traj.pts = F2_data$Parent_temp, 
                                     pca = TRUE, 
                                     print.progress = FALSE)
summary(F2_parent_traj, attribute = "MD", angle.type = "deg")
summary(F2_parent_traj, attribute = "TC", angle.type = "deg")


F2_off_parent_traj = trajectory.analysis(F2_RRPP_mod, 
                                         fit.null = NULL, 
                                         groups = F2_data$Morph, 
                                         traj.pts = F2_data$off_parent_temp, 
                                         pca = TRUE, 
                                         print.progress = FALSE)

summary(F2_off_parent_traj, 
        attribute = "MD", 
        angle.type = "deg", 
        show.trajectories = T)
summary(F2_off_parent_traj, 
        attribute = "TC", 
        angle.type = "deg", 
        show.trajectories = T)

F2_groups <- interaction(F2_data$Morph, F2_data$Offspring_temp, F2_data$Parent_temp)
F2_groups <- interaction(F2_data$Morph, F2_data$Offspring_temp)
F2_groups <- interaction(F2_data$Morph, F2_data$Parent_temp)

pw_traj <- pairwise(F2_RRPP_mod, 
                    groups = F2_groups)
summary(pw_traj)



# PC axes -----------------------------------------------------------------



F2_pca_data = read_csv('F2_orig_PCA_sig_axes.csv') %>% 
  unite(col = off_parent_temp, 
        c('Offspring_temp', 
          'Parent_temp'), 
        sep = '_', 
        remove = F)

WGP_pca_data = read_csv('WGP_PCA_sig_axes.csv')

TGP_pca_data = read_csv('TGP_PCA_sig_axes.csv')

wild_pca_data = read_csv('wild_PCA_sig_axes.csv')

F2_PC_axes = F2_pca_data %>% dplyr::select(PC1:PC7) %>% 
  as.matrix()

F2_RRPP_mod <- lm.rrpp(F2_PC_axes ~ Morph*Offspring_temp*Parent_temp , 
                       SS.type = "I", 
               data = F2_pca_data, 
               print.progress = FALSE, 
               iter = 999,
               verbose = TRUE) 
summary(F2_RRPP_mod, 
        formula = FALSE)
anova(F2_RRPP_mod) 

coef(F2_RRPP_mod, test = TRUE)


F2_off_traj = trajectory.analysis(F2_RRPP_mod, 
                    fit.null = NULL, 
                    groups = F2_pca_data$Morph, 
                    traj.pts = F2_pca_data$Offspring_temp, 
                    pca = TRUE, 
                    print.progress = FALSE)

summary(F2_off_traj, attribute = "MD", angle.type = "deg", 
        show.trajectories = T)
summary(F2_off_traj, attribute = "TC", angle.type = "deg")


F2_parent_traj = trajectory.analysis(F2_RRPP_mod, 
                                  fit.null = NULL, 
                                  groups = F2_pca_data$Morph, 
                                  traj.pts = F2_pca_data$Parent_temp, 
                                  pca = TRUE, 
                                  print.progress = FALSE)
summary(F2_parent_traj, attribute = "MD", angle.type = "deg")


F2_off_parent_traj = trajectory.analysis(F2_RRPP_mod, 
                                     fit.null = NULL, 
                                     groups = F2_pca_data$Morph, 
                                     traj.pts = F2_pca_data$off_parent_temp, 
                                     pca = TRUE, 
                                     print.progress = FALSE)

summary(F2_off_parent_traj, 
        attribute = "MD", 
        angle.type = "deg", 
        show.trajectories = T)
