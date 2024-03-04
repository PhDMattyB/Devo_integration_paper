##############################
## Wild fish integration analysis
##
## Matt Brachmann (PhDMattyB)
##
## 04.03.2024
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(geomorph)
library(vcvComp)
library(factoextra)
library(tidyverse)

wild_tps = readland.tps('Wild_Final.TPS', 
                        specID = 'imageID')

identifiers = read_csv('TPS_Wild_metadata.csv') 

## superimposition on the entire dataset
wild_gpa = gpagen(wild_tps, 
                       print.progress = F)

wild_gpa$coords
## plot the landmarks to see the updated ones on the fish
plotAllSpecimens(wild_gpa$coords)

## integration test without removing allometric scaling


## Need to standardize for allometric variation
procD.lm()
