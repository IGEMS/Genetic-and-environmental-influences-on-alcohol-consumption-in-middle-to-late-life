# ------------------------------------#
#    Author: Daniel Gustavson         #
#        Date: 11/15/2023             #
#                                     #
# Read model output for bivariate     #
# age moderation and meta analyze     #
# across samples                      #
# ------------------------------------#

# Loading Required Libraries
# -----------------------------------------------------------------------
library(tidyverse)
library(metafor)
require(OpenMx)   #Loads OpenMx
require(psych)    #Loads Psych package
source("GenEpiHelperFunctions.R")


# #######################################
# ###       Reformat FTC Output       ###
# #######################################

## Load FTC and compute some extra parameters that weren't in there to begin with 
load(file="output/FTC_univModACEFit_5_23_2023.Rdata")

(FTC_AF <- mxEval(expression=ACE.aF %*% t(ACE.aF), model=univModACEFit))
(FTC_CF <- mxEval(expression=ACE.cF %*% t(ACE.cF), model=univModACEFit))
(FTC_EF <- mxEval(expression=ACE.eF %*% t(ACE.eF), model=univModACEFit))
(FTC_AM <- mxEval(expression=ACE.aM %*% t(ACE.aM), model=univModACEFit))
(FTC_CM <- mxEval(expression=ACE.cM %*% t(ACE.cM), model=univModACEFit))
(FTC_EM <- mxEval(expression=ACE.eM %*% t(ACE.eM), model=univModACEFit))
(FTC_VF <- FTC_AF+FTC_CF+FTC_EF)
(FTC_VM <- FTC_AM+FTC_CM+FTC_EM)
(FTC_stndAF <- mxEval(expression=diag2vec((ACE.aF %*% t(ACE.aF))/(ACE.aF %*% t(ACE.aF)+ACE.cF %*% t(ACE.cF)+ACE.eF %*% t(ACE.eF))), model=univModACEFit)[2])
(FTC_stndCF <- mxEval(expression=diag2vec((ACE.cF %*% t(ACE.cF))/(ACE.aF %*% t(ACE.aF)+ACE.cF %*% t(ACE.cF)+ACE.eF %*% t(ACE.eF))), model=univModACEFit)[2])
(FTC_stndEF <- mxEval(expression=diag2vec((ACE.eF %*% t(ACE.eF))/(ACE.aF %*% t(ACE.aF)+ACE.cF %*% t(ACE.cF)+ACE.eF %*% t(ACE.eF))), model=univModACEFit)[2])
(FTC_stndAM <- mxEval(expression=diag2vec((ACE.aM %*% t(ACE.aM))/(ACE.aM %*% t(ACE.aM)+ACE.cM %*% t(ACE.cM)+ACE.eM %*% t(ACE.eM))), model=univModACEFit)[2])
(FTC_stndCM <- mxEval(expression=diag2vec((ACE.cM %*% t(ACE.cM))/(ACE.aM %*% t(ACE.aM)+ACE.cM %*% t(ACE.cM)+ACE.eM %*% t(ACE.eM))), model=univModACEFit)[2])
(FTC_stndEM <- mxEval(expression=diag2vec((ACE.eM %*% t(ACE.eM))/(ACE.aM %*% t(ACE.aM)+ACE.cM %*% t(ACE.cM)+ACE.eM %*% t(ACE.eM))), model=univModACEFit)[2])

# Now get SEs for these computations
(se_stndAF <- mxSE(diag2vec((ACE.aF %*% t(ACE.aF))/(ACE.aF %*% t(ACE.aF)+ACE.cF %*% t(ACE.cF)+ACE.eF %*% t(ACE.eF))),model= univModACEFit)[2])
(se_stndCF <- mxSE(diag2vec((ACE.cF %*% t(ACE.cF))/(ACE.aF %*% t(ACE.aF)+ACE.cF %*% t(ACE.cF)+ACE.eF %*% t(ACE.eF))),model= univModACEFit)[2])
(se_stndEF <- mxSE(diag2vec((ACE.eF %*% t(ACE.eF))/(ACE.aF %*% t(ACE.aF)+ACE.cF %*% t(ACE.cF)+ACE.eF %*% t(ACE.eF))),model= univModACEFit)[2])
(se_stndAM <- mxSE(diag2vec((ACE.aM %*% t(ACE.aM))/(ACE.aM %*% t(ACE.aM)+ACE.cM %*% t(ACE.cM)+ACE.eM %*% t(ACE.eM))),model= univModACEFit)[2])
(se_stndCM <- mxSE(diag2vec((ACE.cM %*% t(ACE.cM))/(ACE.aM %*% t(ACE.aM)+ACE.cM %*% t(ACE.cM)+ACE.eM %*% t(ACE.eM))),model= univModACEFit)[2])
(se_stndEM <- mxSE(diag2vec((ACE.eM %*% t(ACE.eM))/(ACE.aM %*% t(ACE.aM)+ACE.cM %*% t(ACE.cM)+ACE.eM %*% t(ACE.eM))),model= univModACEFit)[2])

# Save for later
FTC_std <- c(FTC_stndAF, FTC_stndCF, FTC_stndEF, FTC_stndAM, FTC_stndCM, FTC_stndEM)
FTC_se  <- c(se_stndAF,  se_stndCF,  se_stndEF,  se_stndAM,  se_stndCM,  se_stndEM)


#######################################
###     Reading in Model Output     ###
#######################################

## Note: This version excludes samples n<1000 (which almost all had convergence problems except MTSATA)
##   GENDER also excluded (OS only) and SATSA also excluded (convergence problem & n=1185)
# N=67125 with pheno data

#loop_full <- c("FTC", "LSADT", "MADT","MIDT","MIDUS",
#                "OVER50", "SALT") # final 3 have no OSDZ
loop_full <- c("LSADT", "MADT","MIDT","MIDUS",
               "OVER50", "SALT") # final 3 have no OSDZ
loop_no_female <- c("NASNRC", "VETSA")

estimates_full <- list()
SEs_full <- list()
n_obs_full <- list()
STD_full <-list()
STDse_full <- list()
for (i in 1:length(loop_full)){
  print(loop_full[i])
  currsample <- loop_full[i]
   
  # Figure out filename and read in file
  infilename1 <- paste("output/modelout_univModACE_Fit_", paste(currsample, "2024", sep="_"), sep="")
  infilename2 <- paste("output/modelout_univModACE_Summ_", paste(currsample, "2024", sep="_"), sep="")
  #infilename1 <- paste("output/modelout_univModACE_Fit_", currsample, sep="")
  #infilename2 <- paste("output/modelout_univModACE_Summ_", currsample, sep="")
  indatname1 <- paste("Fit_", currsample, sep="")
  indatname2 <- paste("Summ_", currsample, sep="")
  load(file=infilename1)
  load(file=infilename2)
  
  # Rename Fit and Summ objects for later
  assign(indatname1, univModACEFit)
  assign(indatname2, univModACESumm)
  
  # Pull estimates and SEs
  est <- univModACESumm$parameters$Estimate
  SE  <- univModACESumm$parameters$Std.Error
  param_names <- univModACESumm$parameters$name

  # Pull number of observations
  n_obs_full[i] <- univModACESumm$numObs
  
  # add them to data frame
  estimates_full[[i]] <- est
  SEs_full[[i]] <- SE 
  
  
  
  # Standard errors for standardized ACEs
  std_AFse <- mxSE(ACE.stndAF,model= univModACEFit)[2,1]
  std_CFse <- mxSE(ACE.stndCF,model= univModACEFit)[2,1]
  std_EFse <- mxSE(ACE.stndEF,model= univModACEFit)[2,1]
  std_AMse <- mxSE(ACE.stndAM,model= univModACEFit)[2,1]
  std_CMse <- mxSE(ACE.stndCM,model= univModACEFit)[2,1]
  std_EMse <- mxSE(ACE.stndEM,model= univModACEFit)[2,1]
  std_se <- c(std_AFse, std_CFse, std_EFse, std_AMse, std_CMse, std_EMse)
  STDse_full[[i]] <- std_se
  
  std_AF <- univModACEFit$ACE$stndAF$result[2,1]
  std_CF <- univModACEFit$ACE$stndCF$result[2,1]
  std_EF <- univModACEFit$ACE$stndEF$result[2,1]
  std_AM <- univModACEFit$ACE$stndAM$result[2,1]
  std_CM <- univModACEFit$ACE$stndCM$result[2,1]
  std_EM <- univModACEFit$ACE$stndEM$result[2,1]
  std <- c(std_AF, std_CF, std_EF, std_AM, std_CM, std_EM)
  STD_full[[i]] <- std
  
  # Create separate dataframe for est and SE from standardized
  
}
# Convert to list to matrix and add row/column labels
estdat_full <- do.call("cbind",estimates_full)
rownames(estdat_full) <- param_names
colnames(estdat_full) <- loop_full
estdat_full

SEdat_full <- do.call("cbind",SEs_full)
rownames(SEdat_full) <- param_names
colnames(SEdat_full) <- loop_full
SEdat_full

STDdat_full <- do.call("cbind",STD_full)
rownames(STDdat_full) <- c("std_AF", "std_CF", "std_EF", "std_AM", "std_CM", "std_EM")
colnames(STDdat_full) <- loop_full
STDdat_full

STDsedat_full <- do.call("cbind", STDse_full)
rownames(STDsedat_full) <- c("std_AF", "std_CF", "std_EF", "std_AM", "std_CM", "std_EM")
colnames(STDsedat_full) <- loop_full
STDsedat_full



## FTC ONLY
loop_FTC <- c("FTC")
estimates_ftc <- list()
SEs_ftc <- list()
n_obs_ftc <- list()
for (i in 1:length(loop_FTC)){
  print(loop_FTC[i])
  currsample <- loop_FTC[i]
  
  # Figure out filename and read in file
  infilename1 <- paste("output/modelout_univModACE_Fit_", paste(currsample, "2024", sep="_"), sep="")
  infilename2 <- paste("output/modelout_univModACE_Summ_", paste(currsample, "2024", sep="_"), sep="")
  indatname1 <- paste("Fit_", currsample, sep="")
  indatname2 <- paste("Summ_", currsample, sep="")
  load(file=infilename1)
  load(file=infilename2)
  
  # Rename Fit and Summ objects for later
  assign(indatname1, univModACEFit)
  assign(indatname2, univModACESumm)
  
  # Pull estimates and SEs
  est <- univModACESumm$parameters$Estimate
  SE  <- univModACESumm$parameters$Std.Error
  param_names <- univModACESumm$parameters$name
  
  # Pull number of observations
  n_obs_ftc[i] <- univModACESumm$numObs
  
  # add them to data frame
  estimates_ftc[[i]] <- est
  SEs_ftc[[i]] <- SE 
  
}
# Convert to list to matrix and add row/column labels
estdat_ftc <- do.call("cbind",estimates_ftc)
rownames(estdat_ftc) <- param_names
colnames(estdat_ftc) <- loop_FTC
estdat_ftc

SEdat_ftc <- do.call("cbind",SEs_ftc)
rownames(SEdat_ftc) <- param_names
colnames(SEdat_ftc) <- loop_FTC
SEdat_ftc

n_obs_full2 <- c(n_obs_full, n_obs_ftc)

## Combine parameter ests and SEs together 
estdat_full2 <- cbind(estdat_full, estdat_ftc)
sedat_full2 <- cbind(SEdat_full, SEdat_ftc)

## Combine additional ACEs and SEs together
STDdat_full2 <- cbind(STDdat_full, FTC_std)
colnames(STDdat_full2)[7] <- "FTC"
STDsedat_full2 <- cbind(STDsedat_full, FTC_se)
colnames(STDsedat_full2)[7] <- "FTC"



## Now add the groups without females

estimates_no_female <- list()
SEs_no_female <- list()
n_obs_no_female <- list()
STD_no_female <-list()
STDse_no_female <- list()
for (i in 1:length(loop_no_female)){
  print(loop_no_female[i])
  currsample <- loop_no_female[i]
  
  infilename1 <- paste("output/modelout_univModACE_Fit_", paste(currsample, "2024", sep="_"), sep="")
  infilename2 <- paste("output/modelout_univModACE_Summ_", paste(currsample, "2024", sep="_"), sep="")
  indatname1 <- paste("Fit_", currsample, sep="")
  indatname2 <- paste("Summ_", currsample, sep="")
  load(file=infilename1)
  load(file=infilename2)
  assign(indatname1, univModACEFit)
  assign(indatname2, univModACESumm)
  
  # Pull estimates and SEs
  est <- univModACESumm$parameters$Estimate
  SE  <- univModACESumm$parameters$Std.Error
  param_names_no_female <- univModACESumm$parameters$name
  
  # add them to data frame
  estimates_no_female[[i]] <- est
  SEs_no_female[[i]] <- SE 
  
  # Pull number of observations
  n_obs_no_female[i] <- univModACESumm$numObs
  
  # Standard errors for standardized ACEs
  std_AMse <- mxSE(ACE.stndAM,model= univModACEFit)[2,1]
  std_CMse <- mxSE(ACE.stndCM,model= univModACEFit)[2,1]
  std_EMse <- mxSE(ACE.stndEM,model= univModACEFit)[2,1]
  std_se <- c(std_AMse, std_CMse, std_EMse)
  STDse_no_female[[i]] <- std_se
  
  std_AM <- univModACEFit$ACE$stndAM$result[2,1]
  std_CM <- univModACEFit$ACE$stndCM$result[2,1]
  std_EM <- univModACEFit$ACE$stndEM$result[2,1]
  std <- c(std_AM, std_CM, std_EM)
  STD_no_female[[i]] <- std
  
  
}
# Convert to list to matrix and add row/column labels
estdat_no_female <- do.call("cbind",estimates_no_female)
rownames(estdat_no_female) <- param_names_no_female
colnames(estdat_no_female) <- loop_no_female
estdat_no_female

SEdat_no_female <- do.call("cbind",SEs_no_female)
rownames(SEdat_no_female) <- param_names_no_female
colnames(SEdat_no_female) <- loop_no_female
SEdat_no_female

STDdat_no_female <- do.call("cbind",STD_no_female)
rownames(STDdat_no_female) <- c("std_AM", "std_CM", "std_EM")
colnames(STDdat_no_female) <- loop_no_female
STDdat_no_female

STDsedat_no_female <- do.call("cbind", STDse_no_female)
rownames(STDsedat_no_female) <- c("std_AM", "std_CM", "std_EM")
colnames(STDsedat_no_female) <- loop_no_female
STDsedat_no_female



#############################################################
###     Combine into single estimates and SE for meta     ###
#############################################################

estdat <- merge(estdat_full2, estdat_no_female, by='row.names', all=T)
row.names(estdat) <- estdat$Row.names
estdat <- select(estdat, -Row.names)
SEdat  <- merge(sedat_full2,  SEdat_no_female,  by='row.names', all=T)
row.names(SEdat) <- SEdat$Row.names
SEdat  <- select(SEdat, -Row.names)
n_obs <- unlist(c(n_obs_full2, n_obs_no_female))

STDdat <- merge(STDdat_full2, STDdat_no_female, by='row.names', all=T)
row.names(STDdat) <- STDdat$Row.names
STDdat <- select(STDdat, -Row.names)
STDsedat <- merge(STDsedat_full2, STDsedat_no_female, by='row.names', all=T)
row.names(STDsedat) <- STDsedat$Row.names
STDsedat <- select(STDsedat, -Row.names)

# save(estdat,  file="data/Est_for_meta.RData")
# save(SEdat,   file="data/SE_for_meta.RData")
# save(n_obs,   file="data/n_obs_for_meta.RData")
# save(STDdat,  file="data/STD_for_meta.RData")
# save(STDsedat,file="data/STDse_for_meta.RData")

# Can start here
# load(file= "data/Est_for_meta.RData")
# load(file= "data/SE_for_meta.RData")
# load(file= "data/n_obs_for_meta.RData")
# load(file= "data/STD_for_meta.RData")
# load(file= "data/STDse_for_meta.RData")

###########################################################
###     Loop through meta analyses and forest plots     ###
###########################################################

# need to check N - for now it's based on number of observations which should be families

pdf(file="output/output_Forest_Plots_May2024.pdf")  
slab <- colnames(estdat)

for (i in 1:length(estdat[,1])){
  #print(estdat[i,])  
  print(row.names(estdat)[i])
  param <- row.names(estdat)[i]
  
  # Run random effect model
  yi <- as.numeric(estdat[i,])
  vi <- as.numeric(SEdat[i,])*as.numeric(SEdat[i,])
  rma <- rma.uni(yi,vi,weights=n_obs,method="FE")
  
  # Create forest plot
  forest(rma, header=paste("Parameter = ", row.names(estdat)[i], sep=""),
         slab=slab, xlim=c(-2,3))
}

dev.off()



pdf(file="output/output_Forest_Plots_June2024_STD_ACE.pdf")  
slab <- colnames(STDdat)

for (i in 1:length(STDdat[,1])){
  print(row.names(STDdat)[i])
  param <- row.names(STDdat)[i]
  
  # Run random effect model
  yi <- as.numeric(STDdat[i,])
  vi <- as.numeric(STDsedat[i,])*as.numeric(STDsedat[i,])
  rma <- rma.uni(yi,vi,weights=n_obs,method="FE")
  
  # Create forest plot
  forest(rma, header=paste("Parameter = ", row.names(STDdat)[i], sep=""),
         slab=slab, xlim=c(-1,2))
}

dev.off()

