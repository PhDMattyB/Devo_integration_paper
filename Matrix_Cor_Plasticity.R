##############################
## matrix correlation - plasticity
##
## Matt Brachmann (PhDMattyB)
##
## 26.08.2024
##
##############################


setwd('~/Parsons_Postdoc/Stickleback_Morphometric_data/Updated Landmarks/')

# install.packages('MatrixCorrelation')
library(tidyverse)
library(MatrixCorrelation)


# Mat cor F2 offspring effect -------------------------------------------------------
## Run the plasticity_analysis.R script to get the 
## corrected trait correlations

## if this next object doesn't run, re run the script
## Plasticity due to offspring temperature 
off_plasticity_trait_cor

off_plast_ASHNW = off_plasticity_trait_cor$ASHNW
off_plast_ASHNC = off_plasticity_trait_cor$ASHNC

allCorrelations(X1 = off_plast_ASHNW, 
                X2 = off_plast_ASHNC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

off_plast_MYVW = off_plasticity_trait_cor$MYVW
off_plast_MYVC = off_plasticity_trait_cor$MYVC


allCorrelations(X1 = off_plast_MYVW, 
                X2 = off_plast_MYVC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

off_plast_SKRW = off_plasticity_trait_cor$SKRW
off_plast_SKRC = off_plasticity_trait_cor$SKRC

allCorrelations(X1 = off_plast_SKRW, 
                X2 = off_plast_SKRC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

off_plast_GTSW = off_plasticity_trait_cor$GTSW
off_plast_CSWYC = off_plasticity_trait_cor$CSWYC

allCorrelations(X1 = off_plast_GTSW, 
                X2 = off_plast_CSWYC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)


# Mat cor F1 parental effect -------------------------------------------------

parent_plasticity_trait_cor

parent_plast_ASHNW = parent_plasticity_trait_cor$ASHNW
parent_plast_ASHNC = parent_plasticity_trait_cor$ASHNC

allCorrelations(X1 = parent_plast_ASHNW, 
                X2 = parent_plast_ASHNC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

parent_plast_MYVW = parent_plasticity_trait_cor$MYVW
parent_plast_MYVC = parent_plasticity_trait_cor$MYVC


allCorrelations(X1 = parent_plast_MYVW, 
                X2 = parent_plast_MYVC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

parent_plast_SKRW = parent_plasticity_trait_cor$SKRW
parent_plast_SKRC = parent_plasticity_trait_cor$SKRC

allCorrelations(X1 = parent_plast_SKRW, 
                X2 = parent_plast_SKRC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

parent_plast_GTSW = parent_plasticity_trait_cor$GTSW
parent_plast_CSWYC = parent_plasticity_trait_cor$CSWYC

allCorrelations(X1 = parent_plast_GTSW, 
                X2 = parent_plast_CSWYC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)



# Mat cor wild data -------------------------------------------------------

wild_uni_trait_cor

wild_ASHNW = wild_uni_trait_cor$ASHNW
wild_ASHNC = wild_uni_trait_cor$ASHNC

allCorrelations(X1 = wild_ASHNW, 
                X2 = wild_ASHNC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

wild_MYVW = wild_uni_trait_cor$MYVW
wild_MYVC = wild_uni_trait_cor$MYVC


allCorrelations(X1 = wild_MYVW, 
                X2 = wild_MYVC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

wild_SKRW = wild_uni_trait_cor$SKRW
wild_SKRC = wild_uni_trait_cor$SKRC

allCorrelations(X1 = wild_SKRW, 
                X2 = wild_SKRC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

wild_GTSW = wild_uni_trait_cor$GTS
wild_CSWYC = wild_uni_trait_cor$CSWY

allCorrelations(X1 = wild_GTSW, 
                X2 = wild_CSWYC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)


# F2 orig matrix correlations ---------------------------------------------

orig_uni_trait_cor

orig_ASHNW = orig_uni_trait_cor$ASHNW
orig_ASHNC = orig_uni_trait_cor$ASHNC

allCorrelations(X1 = orig_ASHNW, 
                X2 = orig_ASHNC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

orig_MYVW = orig_uni_trait_cor$MYVW
orig_MYVC = orig_uni_trait_cor$MYVC


allCorrelations(X1 = orig_MYVW, 
                X2 = orig_MYVC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

orig_SKRW = orig_uni_trait_cor$SKRW
orig_SKRC = orig_uni_trait_cor$SKRC

allCorrelations(X1 = orig_SKRW, 
                X2 = orig_SKRC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

orig_GTSW = orig_uni_trait_cor$GTS
orig_CSWYC = orig_uni_trait_cor$CSWY

allCorrelations(X1 = orig_GTSW, 
                X2 = orig_CSWYC, 
                ncomp1 = 28,
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)




# within vs tran plasticity -----------------------------------------------

allCorrelations(X1 = off_plasticity_trait_cor$ASHNW, 
                X2 = parent_plasticity_trait_cor$ASHNW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = off_plasticity_trait_cor$ASHNC, 
                X2 = parent_plasticity_trait_cor$ASHNC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = off_plasticity_trait_cor$MYVW, 
                X2 = parent_plasticity_trait_cor$MYVW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = off_plasticity_trait_cor$MYVC, 
                X2 = parent_plasticity_trait_cor$MYVC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = off_plasticity_trait_cor$SKRW, 
                X2 = parent_plasticity_trait_cor$SKRW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)
allCorrelations(X1 = off_plasticity_trait_cor$SKRC, 
                X2 = parent_plasticity_trait_cor$SKRC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = off_plasticity_trait_cor$GTSW, 
                X2 = parent_plasticity_trait_cor$GTSW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = off_plasticity_trait_cor$CSWYC, 
                X2 = parent_plasticity_trait_cor$CSWYC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)



# F2 generation vs TGP matrix ---------------------------------------------


allCorrelations(X1 = orig_uni_trait_cor$ASHNW, 
                X2 = parent_plasticity_trait_cor$ASHNW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = orig_uni_trait_cor$ASHNC, 
                X2 = parent_plasticity_trait_cor$ASHNC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = orig_uni_trait_cor$MYVW, 
                X2 = parent_plasticity_trait_cor$MYVW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = orig_uni_trait_cor$MYVC, 
                X2 = parent_plasticity_trait_cor$MYVC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = orig_uni_trait_cor$SKRW, 
                X2 = parent_plasticity_trait_cor$SKRW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = orig_uni_trait_cor$SKRC, 
                X2 = parent_plasticity_trait_cor$SKRC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = orig_uni_trait_cor$GTSW, 
                X2 = parent_plasticity_trait_cor$GTSW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = orig_uni_trait_cor$CSWYC, 
                X2 = parent_plasticity_trait_cor$CSWYC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)



# F2 generation vs WGP matrix ---------------------------------------------

allCorrelations(X1 = off_plasticity_trait_cor$ASHNW, 
                X2 = orig_uni_trait_cor$ASHNW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = off_plasticity_trait_cor$ASHNC, 
                X2 = orig_uni_trait_cor$ASHNC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = off_plasticity_trait_cor$MYVW, 
                X2 = orig_uni_trait_cor$MYVW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = off_plasticity_trait_cor$MYVC, 
                X2 = orig_uni_trait_cor$MYVC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = off_plasticity_trait_cor$SKRW, 
                X2 = orig_uni_trait_cor$SKRW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)
allCorrelations(X1 = off_plasticity_trait_cor$SKRC, 
                X2 = orig_uni_trait_cor$SKRC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = off_plasticity_trait_cor$GTSW, 
                X2 = orig_uni_trait_cor$GTSW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = off_plasticity_trait_cor$CSWYC, 
                X2 = orig_uni_trait_cor$CSWYC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

