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
                       col_names = T) %>% 
  unite('Ecotype_Pair_Full_Temp', 
        Ecotype_pair, 
        Full_temp, 
        sep = '_', 
        remove = F)

F2_craniofacial_gpa = gpagen(F2_craniofacial,
                             print.progress = F)

# subset_F2_craniofacial_coords = coords.subset(F2_craniofacial_gpa$coords,
#                                               identifiers$Full_temp)

F2_geo_df = geomorph.data.frame(F2_craniofacial_gpa, 
                    Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                    parent_temp = identifiers$Parent_temp, 
                    offspring_temp = identifiers$Offspring_temp, 
                    morph = identifiers$Morph, 
                    population = identifiers$Lake)

lm.rrpp(coords ~ population * morph, 
        data = F2_geo_df)



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

subset_F2_craniofacial_coords = coords.subset(F2_craniofacial_gpa$coords,
                                              identifiers$Lake)


F2_data_frame = data.frame(Y = two.d.array(F2_craniofacial_gpa$coords), 
           Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
           parent_temp = identifiers$Parent_temp, 
           offspring_temp = identifiers$Offspring_temp, 
           morph = identifiers$Morph, 
           population = identifiers$Lake)

ASHN_df = F2_data_frame[F2_data_frame$population == 'ASHN',]

str(ASHN_df)

ASHN_df = ASHN_df %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  rename(individual = rowname)

ASHN_coords = ASHN_df %>% 
  select(Y.1.X:Y.15.Y) %>% 
  as.data.frame()


lm.rrpp(ASHN_coords ~ ASHN_df$parent_temp + ASHN_df$offspring_temp)

ASHN_df[[3]]


fit <- lm.rrpp(coords ~ logSize + Sex*Pop, SS.type = "I", 
               data = Pupfish, print.progress = FALSE,
               turbo = FALSE, verbose = TRUE) 
summary(fit, formula = FALSE)