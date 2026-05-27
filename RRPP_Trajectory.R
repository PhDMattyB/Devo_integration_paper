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

plot(F2_off_traj, 
     pch = as.numeric(F2_pca_data$Morph) + 20,
     bg = as.numeric(F2_pca_data$Offspring_temp),
          cex = 0.7, col = "gray")
add.trajectories(F2_off_traj, 
                 traj.pch = c(21, 22, 23, 24), 
                 start.bg = 1, 
                 end.bg = 2)
legend("topright", levels(Pupfish$Pop), pch =  c(21, 22), pt.bg = 1)


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
