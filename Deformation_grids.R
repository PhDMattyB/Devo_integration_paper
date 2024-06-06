##############################
## Integration Deformation grids
##
## Matt Brachmann (PhDMattyB)
##
##  06.06.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')


library(geomorph)
library(RRPP)
library(tidyverse)
library(car)
library(Morpho)


# metadata ----------------------------------------------------------------

identifiers = read_csv('F2_metadata.csv') %>% 
  rename(individualID = Names) %>% 
  unite('lake_morph_Pair_Full_Temp', 
        Lake_morph, 
        Full_temp, 
        sep = '_', 
        remove = F) %>% 
  unite('Ecotype_Pair_Full_Temp', 
        Ecotype_pair, 
        Full_temp, 
        sep = '_', 
        remove = F) %>% 
  mutate(across(c('Lake_morph',
                  'Offspring_temp',
                  'Parent_temp',
                  'Grand_temp'),
                factor))

# Original F2 traits ------------------------------------------------------

F2_tps = readland.tps('F2_No_GT.TPS', 
                      specID = 'imageID')

## superimposition on the entire dataset
F2_gpa = gpagen(F2_tps, 
                print.progress = F)


F2_original_data = two.d.array(F2_gpa$coords) # get coords into an array
F2_original_data = cbind(identifiers, F2_original_data) # bind with the variables; all in the correct order

test_lm = procD.lm(F2_gpa$coords ~ identifiers$Parent_temp*identifiers$Offspring_temp*identifiers$Lake_morph, 
         iter = 999)
summary(test_lm)

F2_orig_ordered = F2_oringial_data[order(F2_original_data$Lake_morph),]  # order to make it easier to split the dataset

rownames(F2_orig_ordered) <- as.integer(1:931) #name rows to make it easier to organise


# ASHN original trait deformation grids -----------------------------------

ASHNC_orig_shape = F2_orig_ordered[1:121,]
ASHNW_orig_shape = F2_orig_ordered[122:243,]

# test = ASHNC_orig_shape[,-c(1:3)]
ASHNC_orig_shape = arrayspecs(ASHNC_orig_shape, 27, 2)
ASHNC_orig_mean = mshape(ASHNC_orig_shape)

ASHNW_orig_shape = arrayspecs(ASHNW_orig_shape, 27, 2)
ASHNW_orig_mean = mshape(ASHNW_orig_shape)

####plotRefToTarget (reference, target species)

links <- cbind(c(1,2,6,12,13,14,15,16,17,18,19,20,21,22,2,23,23,24,8,25,10,8,7,7), 
               c(2,6,12,13,14,15,16,17,18,19,20,21,22,1,23,27,24,8,27,26,11,9,8,9))
ASHNC_def <- gridPar(tar.pt.size=0.8, 
                     tar.pt.bg = "#ade8f4", 
                     tar.link.col="black", 
                     tar.link.lwd = 2.5,
                     grid.col='black', 
                     grid.lty = 10, 
                     n.col.cell = 30)
plotRefToTarget(ASHNW_orig_mean, 
                ASHNC_orig_mean, 
                mag=4, 
                links=links, 
                gridPars=ASHNC_def)

ASHNW_def <- gridPar(tar.pt.size=0.8, 
                     tar.pt.bg = "#ff006e", 
                     tar.link.col="black", 
                     tar.link.lwd = 2.5,
                     grid.col='black', 
                     grid.lty = 10, 
                     n.col.cell = 30)
plotRefToTarget(ASHNC_orig_mean, 
                ASHNW_orig_mean, 
                mag=4, 
                links=links, 
                gridPars=ASHNW_def)



# ASHN F1 effect deformation grids ----------------------------------------
F1_effect_tps = readland.tps('F1_Corrected_landmarks.tps', 
                             specID = 'imageID')

F1_effect_gpa = gpagen(F1_effect_tps, 
                print.progress = F)


F1_effect_data = two.d.array(F1_effect_gpa$coords) # get coords into an array
F1_effect_data = cbind(identifiers, F1_effect_data) # bind with the variables; all in the correct order




F1_effect_ordered = F1_effect_data[order(F1_effect_data$Lake_morph),]  # order to make it easier to split the dataset

rownames(F1_effect_ordered) <- as.integer(1:931) #name rows to make it easier to organise


ASHNC_F1_shape = F1_effect_ordered[1:121,]
ASHNW_F1_shape = F1_effect_ordered[122:243,]

ASHNC_F1_shape = ASHNC_F1_shape[,-c(1:12)]
ASHNC_F1_shape = arrayspecs(ASHNC_F1_shape, 27, 2)
ASHNC_F1_mean = mshape(ASHNC_F1_shape)

ASHNW_F1_shape = ASHNW_F1_shape[,-c(1:12)]
ASHNW_F1_shape = arrayspecs(ASHNW_F1_shape, 27, 2)
ASHNW_F1_mean = mshape(ASHNW_F1_shape)

####plotRefToTarget (reference, target species)

links <- cbind(c(1,2,6,12,13,14,15,16,17,18,19,20,21,22,2,23,23,24,8,25,10,8,7,7), 
               c(2,6,12,13,14,15,16,17,18,19,20,21,22,1,23,27,24,8,27,26,11,9,8,9))
ASHNC_F1_def <- gridPar(tar.pt.size=0.8, 
                        tar.pt.bg = "#ade8f4", 
                        tar.link.col="black", 
                        tar.link.lwd = 2.5,
                        grid.col='black', 
                        grid.lty = 10, 
                        n.col.cell = 30)
plotRefToTarget(ASHNW_F1_mean, 
                ASHNC_F1_mean, 
                mag=4, 
                links=links, 
                gridPars=ASHNC_F1_def)

ASHNW_F1_def <- gridPar(tar.pt.size=0.8, 
                        tar.pt.bg = "#ff006e", 
                        tar.link.col="black", 
                        tar.link.lwd = 2.5,
                        grid.col='black', 
                        grid.lty = 10, 
                        n.col.cell = 30)
plotRefToTarget(ASHNC_F1_mean, 
                ASHNW_F1_mean, 
                mag=4, 
                links=links, 
                gridPars=ASHNW_F1_def)


# ASHN F2 effect deformation grids ----------------------------------------
F2_effect_tps = readland.tps('F2_Corrected_landmarks.tps', 
                             specID = 'imageID')

F2_effect_gpa = gpagen(F2_effect_tps, 
                       print.progress = F)


F2_effect_data = two.d.array(F2_effect_gpa$coords) # get coords into an array
F2_effect_data = cbind(identifiers, F2_effect_data) # bind with the variables; all in the correct order

F2_effect_ordered = F2_effect_data[order(F2_effect_data$Lake_morph),]  # order to make it easier to split the dataset

rownames(F2_effect_ordered) <- as.integer(1:931) #name rows to make it easier to organise


ASHNC_F2_shape = F2_effect_ordered[1:121,]
ASHNW_F2_shape = F2_effect_ordered[122:243,]

ASHNC_F2_shape = ASHNC_F2_shape[,-c(1:12)]
ASHNC_F2_shape = arrayspecs(ASHNC_F2_shape, 27, 2)
ASHNC_F2_mean = mshape(ASHNC_F2_shape)

ASHNW_F2_shape = ASHNW_F2_shape[,-c(1:12)]
ASHNW_F2_shape = arrayspecs(ASHNW_F2_shape, 27, 2)
ASHNW_F2_mean = mshape(ASHNW_F2_shape)

####plotRefToTarget (reference, target species)

links <- cbind(c(1,2,6,12,13,14,15,16,17,18,19,20,21,22,2,23,23,24,8,25,10,8,7,7), 
               c(2,6,12,13,14,15,16,17,18,19,20,21,22,1,23,27,24,8,27,26,11,9,8,9))
ASHNC_F2_def <- gridPar(tar.pt.size=0.8, 
                        tar.pt.bg = "#ade8f4", 
                        tar.link.col="black", 
                        tar.link.lwd = 2.5,
                        grid.col='black', 
                        grid.lty = 10, 
                        n.col.cell = 30)
plotRefToTarget(ASHNW_F2_mean, 
                ASHNC_F2_mean, 
                mag=4, 
                links=links, 
                gridPars=ASHNC_F2_def)

ASHNW_F2_def <- gridPar(tar.pt.size=0.8, 
                        tar.pt.bg = "#ff006e", 
                        tar.link.col="black", 
                        tar.link.lwd = 2.5,
                        grid.col='black', 
                        grid.lty = 10, 
                        n.col.cell = 30)
plotRefToTarget(ASHNC_F2_mean, 
                ASHNW_F2_mean, 
                mag=4, 
                links=links, 
                gridPars=ASHNW_F2_def)



# ASHN Original vs F1 and F2 effect ----------------------------------------------------


ASHNC_origin_f1_def <- gridPar(tar.pt.size=0.8,
                              tar.pt.bg = "#ade8f4", 
                              tar.link.col="black", 
                              tar.link.lwd = 2.5,
                              grid.col='black', 
                              grid.lty = 10, 
                              n.col.cell = 30)
plotRefToTarget(ASHNC_orig_mean, 
                ASHNC_F1_mean, 
                mag=20, 
                links=links, 
                gridPars=ASHNC_origin_f1_def)


ASHNC_origin_f2_def <- gridPar(tar.pt.size=0.8, 
                               tar.pt.bg = "#ade8f4", 
                               tar.link.col="black", 
                               tar.link.lwd = 2.5,
                               grid.col='black', 
                               grid.lty = 10, 
                               n.col.cell = 30)
plotRefToTarget(ASHNC_orig_mean, 
                ASHNC_F2_mean, 
                mag=20, 
                links=links, 
                gridPars=ASHNC_origin_f2_def)


ASHNW_origin_f1_def <- gridPar(tar.pt.size=0.8, 
                              tar.pt.bg = "#ff006e", 
                              tar.link.col="black", 
                              tar.link.lwd = 2.5,
                              grid.col='black', 
                              grid.lty = 10, 
                              n.col.cell = 30)
plotRefToTarget(ASHNW_orig_mean, 
                ASHNW_F1_mean, 
                mag=20, 
                links=links, 
                gridPars=ASHNW_origin_f1_def)


ASHNW_origin_f2_def <- gridPar(tar.pt.size=0.8, 
                               tar.pt.bg = "#ff006e", 
                               tar.link.col="black", 
                               tar.link.lwd = 2.5,
                               grid.col='black', 
                               grid.lty = 10, 
                               n.col.cell = 30)
plotRefToTarget(ASHNW_orig_mean, 
                ASHNW_F2_mean, 
                mag=20, 
                links=links, 
                gridPars=ASHNW_origin_f2_def)



# ASHN wild deformation grids ---------------------------------------------

wild_tps = readland.tps('Wild_Final.TPS', 
                        specID = 'imageID')

wild_identifiers = read_csv('TPS_Wild_metadata.csv') 

## superimposition on the entire dataset
wild_gpa = gpagen(wild_tps, 
                  print.progress = F)


wild_effect_data = two.d.array(wild_gpa$coords) # get coords into an array
wild_effect_data = cbind(wild_identifiers, wild_effect_data) # bind with the variables; all in the correct order

wild_effect_ordered = wild_effect_data[order(wild_effect_data$Lake_morph),]  # order to make it easier to split the dataset

rownames(wild_effect_ordered) <- as.integer(1:331) #name rows to make it easier to organise


ASHNC_wild_shape = wild_effect_ordered[1:30,]
ASHNW_wild_shape = wild_effect_ordered[31:60,]

ASHNC_wild_shape = ASHNC_wild_shape[,-c(1:5)]
ASHNC_wild_shape = arrayspecs(ASHNC_wild_shape, 27, 2)
ASHNC_wild_mean = mshape(ASHNC_wild_shape)

ASHNW_wild_shape = ASHNW_wild_shape[,-c(1:5)]
ASHNW_wild_shape = arrayspecs(ASHNW_wild_shape, 27, 2)
ASHNW_wild_mean = mshape(ASHNW_wild_shape)

####plotRefToTarget (reference, target species)

links <- cbind(c(1,2,6,12,13,14,15,16,17,18,19,20,21,22,2,23,23,24,8,25,10,8,7,7), 
               c(2,6,12,13,14,15,16,17,18,19,20,21,22,1,23,27,24,8,27,26,11,9,8,9))
ASHNC_wild_def <- gridPar(tar.pt.size=0.8, 
                          tar.pt.bg = "#ade8f4", 
                          tar.link.col="black", 
                          tar.link.lwd = 2.5,
                          grid.col='black', 
                          grid.lty = 10, 
                          n.col.cell = 30)
plotRefToTarget(ASHNW_wild_mean, 
                ASHNC_wild_mean, 
                mag=3, 
                links=links, 
                gridPars=ASHNC_wild_def)

ASHNW_wild_def <- gridPar(tar.pt.size=0.8, 
                          tar.pt.bg = "#ff006e", 
                          tar.link.col="black", 
                          tar.link.lwd = 2.5,
                          grid.col='black', 
                          grid.lty = 10, 
                          n.col.cell = 30)
plotRefToTarget(ASHNC_wild_mean, 
                ASHNW_wild_mean, 
                mag=4, 
                links=links, 
                gridPars=ASHNW_wild_def)
