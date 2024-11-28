##############################
## Jaw kinetics (Corins code)
##
## Matt Brachmann (PhDMattyB)
##
## 28.11.2024
##
##############################

setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

library(geomorph)
library(tidyverse)
library(linkR)


## landmarks read during the wild_pattern script
## TGP and WPG effects isolated in plasticity_analysis script
## All LMs have been put through a common gpa
## they should be good to go straight into this script


# Wild kinetics -----------------------------------------------------------

###Calculate transmission coefficients for Anterior/Pre-Maxillary linkage (23,27,37,3)
PM_LMs <- c(23,27,28,3)
wild_coords = wild_lmk_coords[PM_LMs,,]


#Create empty table for data
KT_PreMax <- numeric(length=331)
PreMax_Rotation <- numeric(length=331)
for (i in 1:331) {
  #Define linkages, joints etc.
  PM_ind <- wild_coords[,,i]
  lms <- PM_ind
  joint.types <- c('R','S','S','R')
  joint.conn <- rbind(c(0,1),c(1,2),c(2,3),c(3,0))
  joint.cons <- list(NA,NA,NA,NA)
  
  linkage <- defineLinkage(joint.coor = lms, joint.types = joint.types,
                           joint.cons = joint.cons, joint.conn=joint.conn)
  #Overcome error message
  linkage2 <- linkage
  linkage2$joint.cons <- list(c(0,0,1), NA, NA, c(0,0,1))  
  
  #Animate linkage, calculate transmission coefficients, add to table
  drawLinkage(linkage2, animate=FALSE)
  animate <- animateLinkage(linkage2, input.param=c(0,pi/1800), input.joint=4)
  kine <- linkageKinematics(animate)
  KT <- abs(kine$links.rdis.d['Link2',] / kine$links.rdis.d['Link1',])
  KT_PreMax[i] <- KT[2]
  KT_PM_df <- as.data.frame(KT_PreMax)
  
  #Calculate Output rotation
  PM_Output <- print(kine$links.rdis.d['Link2',])
  PreMax_Rotation[i] <- PM_Output[2]
  PM_Rotation_DF <- as.data.frame(PreMax_Rotation)
}

Pheno_4Bar <- cbind(KT_PM_df, PM_Rotation_DF)

