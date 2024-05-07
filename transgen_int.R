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

F2_geo_df = geomorph.data.frame(coords = two.d.array(F2_craniofacial_gpa$coords), 
                    Full_factor = identifiers$Ecotype_Pair_Full_Temp, 
                    parent_temp = identifiers$Parent_temp, 
                    offspring_temp = identifiers$Offspring_temp,
                    grand_temp = identifiers$Grand_temp,
                    morph = identifiers$Morph, 
                    population = identifiers$Lake)

F2_off_mod = lm.rrpp(coords ~ offspring_temp, 
                     data = F2_geo_df, 
                     turbo = F, 
                     verbose = T)

summary(F2_off_mod)
anova(F2_off_mod)
coef(F2_off_mod, 
     test = T)


F2_par_mod = lm.rrpp(coords ~ parent_temp, 
                     data = F2_geo_df, 
                     turbo = F, 
                     verbose = T)

summary(F2_par_mod)
anova(F2_par_mod)
coef(F2_par_mod, 
     test = T)
F2_par_fitted = F2_par_mod$LM$fitted

F2_grandpar_mod = lm.rrpp(coords ~ grand_temp, 
                     data = F2_geo_df, 
                     turbo = F, 
                     verbose = T)

summary(F2_grandpar_mod)
anova(F2_grandpar_mod)
coef(F2_grandpar_mod, 
     test = T)
F2_grandpar_fitted = F2_grandpar_mod$LM$fitted

F2_par_off_mod = lm.rrpp(coords ~ parent_temp * offspring_temp, 
                     data = F2_geo_df, 
                     turbo = F, 
                     verbose = T)

summary(F2_par_off_mod)
anova(F2_par_off_mod)
coef(F2_par_off_mod, 
     test = T)
F2_par_off_fitted = F2_par_off_mod$LM$fitted

F2_grandpar_off_mod = lm.rrpp(coords ~ grand_temp * offspring_temp, 
                         data = F2_geo_df, 
                         turbo = F, 
                         verbose = T)

summary(F2_grandpar_off_mod)
anova(F2_grandpar_off_mod)
coef(F2_grandpar_off_mod, 
     test = T)
F2_grandpar_off_fitted = F2_grandpar_off_mod$LM$fitted

F2_grandpar_par_mod = lm.rrpp(coords ~ grand_temp * parent_temp, 
                              data = F2_geo_df, 
                              turbo = F, 
                              verbose = T)

summary(F2_grandpar_par_mod)
anova(F2_grandpar_par_mod)
coef(F2_grandpar_par_mod, 
     test = T)
F2_grandpar_par_fitted = F2_grandpar_par_mod$LM$fitted


F2_Full_mod = lm.rrpp(coords ~ grand_temp * parent_temp * offspring_temp, 
                      data = F2_geo_df, 
                      turbo = F, 
                      verbose = T)

summary(F2_Full_mod)
anova(F2_Full_mod)
coef(F2_Full_mod, 
     test = T)
F2_full_fitted = F2_Full_mod$LM$fitted


# Integration test using fitted values ------------------------------------


# Offspring fitted integration --------------------------------------------

## So this is a test of integration only due to offspring temperature
## So if offspring were raised at either 12 or 18 and does
## not account for any transgenerational differences

F2_off_fitted = F2_off_mod$LM$fitted
## need to convert this to a 3d array


identifiers = read_csv('F2_Metadata.CSV', 
                       col_names = T) %>% 
  unite('Ecotype_Pair_Full_Temp', 
        Ecotype_pair, 
        Full_temp, 
        sep = '_', 
        remove = F)

F2_off_fit_sub = coords.subset(F2_off_fitted,
                               identifiers$Ecotype_Pair_Full_Temp)


## integration analysis code
subset_F2_craniofacial_coords = coords.subset(F2_craniofacial_gpa$coords,
                                              identifiers$Ecotype_Pair_Full_Temp)

vrel_F2_craniofacial = Map(function(x) integration.Vrel(x),
                           subset_F2_craniofacial_coords)

## ASHN shows a high degree of plasticity in the F2 generation
## and no effect of the partental F1 generation on integration
ASHN_12_plasticity = compare.ZVrel(vrel_F2_craniofacial$`ASHN_12@12`,
                                   vrel_F2_craniofacial$`ASHN_12@18`)

ASHN_18_plasticity = compare.ZVrel(vrel_F2_craniofacial$`ASHN_18@12`,
                                   vrel_F2_craniofacial$`ASHN_18@18`)

ASHN_12_Trans = compare.ZVrel(vrel_F2_craniofacial$`ASHN_12@12`,
                              vrel_F2_craniofacial$`ASHN_18@12`)

ASHN_18_trans = compare.ZVrel(vrel_F2_craniofacial$`ASHN_12@18`,
                              vrel_F2_craniofacial$`ASHN_18@18`)

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
