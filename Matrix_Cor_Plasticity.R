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

off_plast_ASHN_RV = mrv(y = off_plast_ASHNW, 
                   x = off_plast_ASHNC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_ASHNW <- off_plast_ASHNW[sample(nrow(off_plast_ASHNW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = off_plast_ASHNC, y = permuted_ASHNW)
}
off_plast_ASHN_RV_sig = (sum(perm_rv_vals >= off_plast_ASHN_RV) + 1) / (n_perm + 1)

# allCorrelations(X1 = off_plast_ASHNW, 
#                 X2 = off_plast_ASHNC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)

off_plast_MYVW = off_plasticity_trait_cor$MYVW
off_plast_MYVC = off_plasticity_trait_cor$MYVC

off_plast_MYV_RV = mrv(y = off_plast_MYVW, 
                        x = off_plast_MYVC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_MYVW <- off_plast_MYVW[sample(nrow(off_plast_MYVW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = off_plast_MYVC, y = permuted_MYVW)
}
off_plast_MYV_RV_sig = (sum(perm_rv_vals >= off_plast_MYV_RV) + 1) / (n_perm + 1)


# allCorrelations(X1 = off_plast_MYVW, 
#                 X2 = off_plast_MYVC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)

off_plast_SKRW = off_plasticity_trait_cor$SKRW
off_plast_SKRC = off_plasticity_trait_cor$SKRC

off_plast_SKR_RV = mrv(y = off_plast_SKRW, 
                        x = off_plast_SKRC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_SKRW <- off_plast_SKRW[sample(nrow(off_plast_SKRW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = off_plast_SKRC, y = permuted_SKRW)
}
off_plast_SKR_RV_sig = (sum(perm_rv_vals >= off_plast_SKR_RV) + 1) / (n_perm + 1)


# allCorrelations(X1 = off_plast_SKRW, 
#                 X2 = off_plast_SKRC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)

off_plast_GTSW = off_plasticity_trait_cor$GTSW
off_plast_CSWYC = off_plasticity_trait_cor$CSWYC

off_plast_GTS_CSWY_RV = mrv(y = off_plast_GTSW, 
                        x = off_plast_CSWYC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_GTSW <- off_plast_GTSW[sample(nrow(off_plast_GTSW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = off_plast_CSWYC, y = permuted_GTSW)
}
off_plast_GTS_CSWY_RV_sig = (sum(perm_rv_vals >= off_plast_GTS_CSWY_RV) + 1) / (n_perm + 1)


# allCorrelations(X1 = off_plast_GTSW, 
#                 X2 = off_plast_CSWYC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)


# Mat cor F1 parental effect -------------------------------------------------

parent_plasticity_trait_cor

parent_plast_ASHNW = parent_plasticity_trait_cor$ASHNW
parent_plast_ASHNC = parent_plasticity_trait_cor$ASHNC

parent_plast_ASHN_RV = mrv(y = parent_plast_ASHNW, 
                        x = parent_plast_ASHNC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_ASHNW <- parent_plast_ASHNW[sample(nrow(parent_plast_ASHNW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = parent_plast_ASHNC, y = permuted_ASHNW)
}
parent_plast_ASHN_RV_sig = (sum(perm_rv_vals >= parent_plast_ASHN_RV) + 1) / (n_perm + 1)

# allCorrelations(X1 = parent_plast_ASHNW, 
#                 X2 = parent_plast_ASHNC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)

parent_plast_MYVW = parent_plasticity_trait_cor$MYVW
parent_plast_MYVC = parent_plasticity_trait_cor$MYVC

parent_plast_MYV_RV = mrv(y = parent_plast_MYVW, 
                           x = parent_plast_MYVC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_MYVW <- parent_plast_MYVW[sample(nrow(parent_plast_MYVW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = parent_plast_MYVC, y = permuted_MYVW)
}
parent_plast_MYV_RV_sig = (sum(perm_rv_vals >= parent_plast_MYV_RV) + 1) / (n_perm + 1)

# allCorrelations(X1 = parent_plast_MYVW, 
#                 X2 = parent_plast_MYVC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)

parent_plast_SKRW = parent_plasticity_trait_cor$SKRW
parent_plast_SKRC = parent_plasticity_trait_cor$SKRC

parent_plast_SKR_RV = mrv(y = parent_plast_SKRW, 
                           x = parent_plast_SKRC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_SKRW <- parent_plast_SKRW[sample(nrow(parent_plast_SKRW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = parent_plast_SKRC, y = permuted_SKRW)
}
parent_plast_SKR_RV_sig = (sum(perm_rv_vals >= parent_plast_SKR_RV) + 1) / (n_perm + 1)

# allCorrelations(X1 = parent_plast_SKRW, 
#                 X2 = parent_plast_SKRC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)

parent_plast_GTSW = parent_plasticity_trait_cor$GTSW
parent_plast_CSWYC = parent_plasticity_trait_cor$CSWYC

parent_plast_GTS_CSWY_RV = mrv(y = parent_plast_GTSW, 
                           x = parent_plast_CSWYC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_GTSW <- parent_plast_GTSW[sample(nrow(parent_plast_GTSW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = parent_plast_CSWYC, y = permuted_GTSW)
}
parent_plast_GTS_CSWY_RV_sig = (sum(perm_rv_vals >= parent_plast_GTS_CSWY_RV) + 1) / (n_perm + 1)


# allCorrelations(X1 = parent_plast_GTSW, 
#                 X2 = parent_plast_CSWYC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)



# Mat cor wild data -------------------------------------------------------

wild_uni_trait_cor

wild_ASHNW = wild_uni_trait_cor$ASHNW
wild_ASHNC = wild_uni_trait_cor$ASHNC

Wild_ASHN_RV = mrv(y = wild_ASHNW, 
                   x = wild_ASHNC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_ASHNW <- wild_ASHNW[sample(nrow(wild_ASHNW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = wild_ASHNC, y = permuted_ASHNW)
}
Wild_ASHN_RV_sig = (sum(perm_rv_vals >= Wild_ASHN_RV) + 1) / (n_perm + 1)


# 
# allCorrelations(X1 = wild_ASHNW, 
#                 X2 = wild_ASHNC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)

wild_MYVW = wild_uni_trait_cor$MYVW
wild_MYVC = wild_uni_trait_cor$MYVC

Wild_MYV_RV = mrv(y = wild_MYVW, 
                   x = wild_MYVC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_MYVW <- wild_MYVW[sample(nrow(wild_MYVW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = wild_MYVC, y = permuted_MYVW)
}
Wild_MYV_RV_sig = (sum(perm_rv_vals >= Wild_MYV_RV) + 1) / (n_perm + 1)

# 
# allCorrelations(X1 = wild_MYVW, 
#                 X2 = wild_MYVC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)

wild_SKRW = wild_uni_trait_cor$SKRW
wild_SKRC = wild_uni_trait_cor$SKRC

Wild_SKR_RV = mrv(y = wild_SKRW, 
                   x = wild_SKRC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_SKRW <- wild_SKRW[sample(nrow(wild_SKRW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = wild_SKRC, y = permuted_SKRW)
}
Wild_SKR_RV_sig = (sum(perm_rv_vals >= Wild_SKR_RV) + 1) / (n_perm + 1)


# allCorrelations(X1 = wild_SKRW, 
#                 X2 = wild_SKRC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)

wild_GTSW = wild_uni_trait_cor$GTS
wild_CSWYC = wild_uni_trait_cor$CSWY

Wild_GTS_CSWY_RV = mrv(y = wild_GTSW, 
                   x = wild_CSWYC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_GTS_CSWY <- wild_GTSW[sample(nrow(wild_GTSW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = wild_CSWYC, y = permuted_GTS_CSWY)
}
Wild_GTS_CSWY_RV_sig = (sum(perm_rv_vals >= Wild_GTS_CSWY_RV) + 1) / (n_perm + 1)


# allCorrelations(X1 = wild_GTSW, 
#                 X2 = wild_CSWYC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)


# F2 orig matrix correlations ---------------------------------------------

orig_uni_trait_cor

orig_ASHNW = orig_uni_trait_cor$ASHNW
orig_ASHNC = orig_uni_trait_cor$ASHNC

orig_ASHN_RV = mrv(y = orig_ASHNW, 
                   x = orig_ASHNC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_ASHNW <- orig_ASHNW[sample(nrow(orig_ASHNW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = orig_ASHNC, y = permuted_ASHNW)
}
orig_ASHN_RV_sig = (sum(perm_rv_vals >= orig_ASHN_RV) + 1) / (n_perm + 1)


# allCorrelations(X1 = orig_ASHNW, 
#                 X2 = orig_ASHNC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)

orig_MYVW = orig_uni_trait_cor$MYVW
orig_MYVC = orig_uni_trait_cor$MYVC

orig_MYV_RV = mrv(y = orig_MYVW, 
                   x = orig_MYVC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_MYVW <- orig_MYVW[sample(nrow(orig_MYVW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = orig_MYVC, y = permuted_MYVW)
}
orig_MYV_RV_sig = (sum(perm_rv_vals >= orig_MYV_RV) + 1) / (n_perm + 1)


# allCorrelations(X1 = orig_MYVW, 
#                 X2 = orig_MYVC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)

orig_SKRW = orig_uni_trait_cor$SKRW
orig_SKRC = orig_uni_trait_cor$SKRC

orig_SKR_RV = mrv(y = orig_SKRW, 
                   x = orig_SKRC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_SKRW <- orig_SKRW[sample(nrow(orig_SKRW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = orig_SKRC, y = permuted_SKRW)
}
orig_SKR_RV_sig = (sum(perm_rv_vals >= orig_SKR_RV) + 1) / (n_perm + 1)

# allCorrelations(X1 = orig_SKRW, 
#                 X2 = orig_SKRC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)

orig_GTSW = orig_uni_trait_cor$GTS
orig_CSWYC = orig_uni_trait_cor$CSWY

orig_GTS_CSWY_RV = mrv(y = orig_GTSW, 
                   x = orig_CSWYC)

n_perm = 9999
perm_rv_vals = numeric(n_perm)
for (i in 1:n_perm) {
  # Shuffle the rows of Y
  permuted_GTSW <- orig_GTSW[sample(nrow(orig_GTSW)), ]
  
  # Calculate the mRV for the permuted data and store it
  perm_rv_vals[i] <- mrv(x = orig_CSWYC, y = permuted_GTSW)
}
orig_GTS_CSWY_RV_sig = (sum(perm_rv_vals >= orig_GTS_CSWY_RV) + 1) / (n_perm + 1)


# allCorrelations(X1 = orig_GTSW, 
#                 X2 = orig_CSWYC, 
#                 ncomp1 = 28,
#                 ncomp2 = 28, 
#                 method = 'RVadj', 
#                 plot = F)




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


# Wild to orig F2 generation ----------------------------------------------

allCorrelations(X1 = wild_uni_trait_cor$ASHNW, 
                X2 = orig_uni_trait_cor$ASHNW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$ASHNC, 
                X2 = orig_uni_trait_cor$ASHNC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$MYVW, 
                X2 = orig_uni_trait_cor$MYVW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$MYVC, 
                X2 = orig_uni_trait_cor$MYVC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$SKRW, 
                X2 = orig_uni_trait_cor$SKRW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)
allCorrelations(X1 = wild_uni_trait_cor$SKRC, 
                X2 = orig_uni_trait_cor$SKRC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$GTS, 
                X2 = orig_uni_trait_cor$GTSW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$CSWY, 
                X2 = orig_uni_trait_cor$CSWYC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

# RV wild to F1 effects ---------------------------------------------------

allCorrelations(X1 = wild_uni_trait_cor$ASHNW, 
                X2 = parent_plasticity_trait_cor$ASHNW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$ASHNC, 
                X2 = parent_plasticity_trait_cor$ASHNC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$MYVW, 
                X2 = parent_plasticity_trait_cor$MYVW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$MYVC, 
                X2 = parent_plasticity_trait_cor$MYVC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$SKRW, 
                X2 = parent_plasticity_trait_cor$SKRW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)
allCorrelations(X1 = wild_uni_trait_cor$SKRC, 
                X2 = parent_plasticity_trait_cor$SKRC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$GTS, 
                X2 = parent_plasticity_trait_cor$GTSW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$CSWY, 
                X2 = parent_plasticity_trait_cor$CSWYC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)


# RV wild to F2 effects ---------------------------------------------------

allCorrelations(X1 = wild_uni_trait_cor$ASHNW, 
                X2 = off_plasticity_trait_cor$ASHNW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$ASHNC, 
                X2 = off_plasticity_trait_cor$ASHNC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$MYVW, 
                X2 = off_plasticity_trait_cor$MYVW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$MYVC, 
                X2 = off_plasticity_trait_cor$MYVC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$SKRW, 
                X2 = off_plasticity_trait_cor$SKRW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)
allCorrelations(X1 = wild_uni_trait_cor$SKRC, 
                X2 = off_plasticity_trait_cor$SKRC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$GTS, 
                X2 = off_plasticity_trait_cor$GTSW, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

allCorrelations(X1 = wild_uni_trait_cor$CSWY, 
                X2 = off_plasticity_trait_cor$CSWYC, 
                ncomp1 = 28, 
                ncomp2 = 28, 
                method = 'RVadj', 
                plot = F)

