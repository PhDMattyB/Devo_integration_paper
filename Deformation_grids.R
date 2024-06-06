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

F2_orig_ordered = F2_oringial_data[order(F2_original_data$Lake_morph),]  # order to make it easier to split the dataset

rownames(F2_orig_ordered) <- as.integer(1:931) #name rows to make it easier to organise

ASHNC_orig_shape = F2_orig_ordered[1:121,]
ASHNW_orig_shape = F2_orig_ordered[122:243,]

# test = ASHNC_orig_shape[,-c(1:3)]
ASHNC_orig_shape = arrayspecs(ASHNC_orig_shape, 27, 2)
ASHNC_orig_mean = mshape(ASHNC_orig_shape)

ASHNW_orig_shape = arrayspecs(ASHNW_orig_shape, 27, 2)
ASHNW_orig_mean = mshape(ASHNW_orig_shape)

links <- cbind(c(1,2,7,8,8,9), c(2,7,8,9,13,10))
LFBen <- gridPar(tar.pt.size=0.8, tar.pt.bg = "green", tar.link.col="black", tar.link.lwd = 2.5,grid.col='darkgreen', grid.lty = 10, n.col.cell = 30)
plotRefToTarget(LFLimneticMean, LFBenthicMean, mag=4, links=links, gridPars=LFBen)
