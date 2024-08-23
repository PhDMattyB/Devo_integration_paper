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

wild_univariate = read_csv('Wild_Univariate_traits.csv')

F2_parental_effects = read_csv('F1_Plasticity_Corrected.csv')

F2_offspring_effects = read_csv('F2_Corrected_F2_temp_only.csv')
