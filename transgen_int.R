##############################
## Transgenerational integration models
##
## Matt Brachmann (PhDMattyB)
##
## 06.05.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(geomorph)
library(vcvComp)
library(factoextra)
library(tidyverse)



# Transgenerational plasticity all lakes combo ----------------------------


F2_tps = readland.tps('F2_No_GT.TPS', 
                      specID = 'imageID')

identifiers = read_csv('F2_Metadata.CSV', 
                       col_names = T) 

F2_craniofacial_gpa = gpagen(F2_craniofacial,
                             print.progress = F)


subset_F2_craniofacial_coords = coords.subset(F2_craniofacial_gpa$coords,
                                              identifiers$Full_temp)




# Transgenerational plasticity per lake -----------------------------------


F2_craniofacial = readland.tps('F2_Craniofacial_LM.TPS',
                               specID = 'imageID')

identifiers = read_csv('F2_metadata.csv') %>% 
  unite('Ecotype_Pair_Full_Temp', 
        Ecotype_pair, 
        Full_temp, 
        sep = '_', 
        remove = F)


F2_craniofacial_gpa = gpagen(F2_craniofacial,
                             print.progress = F)

subset_F2_craniofacial_coords = coords.subset(F2_craniofacial_gpa$coords,
                                              identifiers$Ecotype_Pair_Full_Temp)


