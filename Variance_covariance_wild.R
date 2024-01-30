##############################
## Variance-Covariance wild stickleback
##
## Matt Brachmann (PhDMattyB)
##
## 25.01.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/')

# install.packages('geomorph')
# install.packages('vcvComp')

library(geomorph)
library(vcvComp)
library(tidyverse)

## Example data from cichlids
data("Tropheus")

landmarks = read_csv('allometry minimised data (XY) with ID (6 population pairs).csv')

identifiers = landmarks %>% 
  select(ID, 
         POP, 
         Morph, 
         CS)

LM_data = landmarks %>% 
  select(-ID, 
         -POP, 
         -Morph, 
         -CS) %>% 
  select(-starts_with('LMS'))

phenotypes = as.matrix(LM_data[which(names(LM_data) == 'LM1X'):
                    which(names(LM_data) == 'LM22Y')])

# dim(phenotypes)

phenotype_array = arrayspecs(phenotypes, p = 22, k = 2)


phenotype_gpa = gpagen(phenotype_array, 
                       print.progress = F)
