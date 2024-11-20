##############################
## Family effects - integration paper
##
## Matt Brachmann (PhDMattyB)
##
## 20.11.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(tidyverse)

F1_lmk_dist = read_csv('F1_Plasticity_Corrected.csv')

F1_lmk_dist %>% 
  select(individualID) %>% 
  separate(col = individualID, 
           into = c('Ecotype', 
                    'Temp', 
                    'Exp_unit', 
                    'Individual_num'), 
           sep = '_') %>% 
  group_by(Ecotype, 
           Temp, 
           Exp_unit) %>% 
  summarize(n = n()) %>% 
  arrange(Ecotype, 
          Temp, 
          Exp_unit) %>% 
  View()

