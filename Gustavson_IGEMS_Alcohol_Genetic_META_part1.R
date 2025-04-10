# ------------------------------------
#    Authors: Daniel Gustavson & Matt Panizzon
#
# Univariate Moderated Twin Analysis 
# Moderator can differ between twins
# Sex differences
# ------------------------------------

# Loading Required Libraries
# -----------------------------------------------------------------------
library(tidyverse)
require(OpenMx)   #Loads OpenMx
require(psych)    #Loads Psych package
source("GenEpiHelperFunctions.R")

########################################
# Reading in Data and create variables #
########################################
#data <- read.csv("Alcohol_GENETIC_data_v2.csv",header=T)
load("Twin_File_for_OpenMX.RData")
names(newtwins)

### Create age spline moderator variables
# Create centered age variable
newtwins <- newtwins %>% 
  rowwise() %>% 
  mutate(age_c=mean(c(AGE_1stassessed_1, AGE_1stassessed_2), na.rm=T))
describe(newtwins$age_c) # M = 56.1 SD = 11.27
(75-56.1)/11.27 # compute SD theshold for age 75 = 1.67
newtwins$zage_c <- scale(newtwins$age_c)
describe(newtwins$zage_c)
newtwins$zage_spline75_younger <- (pmin(0,(newtwins$zage_c-1.67)))
newtwins$zage_spline75_older <- (pmax(0,(newtwins$zage_c-1.67)))
hist(newtwins$zage_spline75_younger)
hist(newtwins$zage_spline75_older)

### Create orthogonal contrasts for country
table(newtwins$COUNTRY)

# Sweden first (most families)
(newtwins <- newtwins %>% mutate(Country_x1 = case_when(
  (COUNTRY=="SWE" ~ 1),
  ((COUNTRY=="DEN" | COUNTRY=="AUS" | COUNTRY=="USA")  ~ -(1/3))
  )))

# Denmark second (2nd most families)
(newtwins <- newtwins %>% mutate(Country_x2 = case_when(
  (COUNTRY=="SWE" ~ 0),
  (COUNTRY=="DEN" ~ 1),
  ((COUNTRY=="AUS" | COUNTRY=="USA")  ~ -.5)
)))

# Contrast 3 (USA vs. AUS)
(newtwins <- newtwins %>% mutate(Country_x3 = case_when(
  (COUNTRY=="SWE" ~ 0),
  (COUNTRY=="DEN" ~ 0),
  (COUNTRY=="AUS" ~ -1),
  (COUNTRY=="USA" ~ 1)
)))
table(newtwins$COUNTRY)
table(newtwins$Country_x1)
table(newtwins$Country_x2)
table(newtwins$Country_x3)

## Compute Zygosity / Sex Variable
(newtwins <- newtwins %>% mutate(zygSex = case_when(
  ((ZYGOS==1 & (sexF_1==0 | sexF_2==0)) ~ 1),   #MZ Male
  ((ZYGOS==1 & (sexF_1==1 | sexF_2==1)) ~ 2),   #MZ Female
  ((ZYGOS==2 & (sexF_1==0 | sexF_2==0)) ~ 3),   #DZ Male
  ((ZYGOS==2 & (sexF_1==1 | sexF_2==1)) ~ 4),   #DZ Female
  (ZYGOS==3 ~ 5)                  #DZ Opposite Sex
)))
table(newtwins$ZYGOS)
table(newtwins$zygSex)

cor.test(newtwins$AGE_1stassessed_1, newtwins$AGE_1stassessed_2)

#save(newtwins, file="Genetic_Data_Prepped.RData")
load("Genetic_Data_Prepped.RData")

### #IF RUNNING ON THE CLUSTER - DO THIS ####
# sinteractive --time=8:00:00 --ntasks=8 --part=ainteractive
# ml anaconda
# conda activate /projects/lessem/software/anaconda/envs/R-latest
# R
# library(OpenMx)


########################################
#        Create for loop to go         #
#  through all studies and run model   #
########################################

# Full list of studies - but we cannot loop through all because
#                        some are missing certain groups (e.g., OSDZ)
(loop <- unique(newtwins$SAMPLE))


# Subset studies based on data available
loop_full <- c("OVER50", "SALT","MIDUS", "MIDT", "LSADT", "MADT")
#loop_no_OSDZ <- c("MTSADA","OCTOTWIN","SATSA")
loop_no_female <- c("VETSA","NASNRC")

# Studies excluded entirely: GENDER (OSDZ is the only group)
# loop_full <- c("OVER50", "HARMONY", "OATS", "SALT","MIDUS", "MIDT", "LSADT", "MADT")
# loop_no_OSDZ <- c("MTSADA","OCTOTWIN","SATSA")
# loop_no_female <- c("VETSA","NASNRC")

#mxOption(NULL, "Default optimizer", "NPSOL") 
mxOption(NULL, "Default optimizer", "SLSQP") 


#loop_full <- c("HARMONY")
for (i in loop_full){
  print(i)
  currsample <- i
  currdat <- filter(newtwins, SAMPLE==i)
  #print(length(currdat$PAIRID)) 
  #print(describe(currdat$zISCED_trim_1))  
  #print(describe(currdat$zage_c))

  twinvar <- c('zISCED_trim', 'zgrams_sqrt_windsor') #variable of interest
  #modvar <- c('zage_c','Country_x1','Country_x2','Country_x3') #moderator variable
  modvar <- c('zage_c') #moderator variable
  
  
  # Defining Variables for OpenMx
  nv <- 2
  ntv <- nv*2
  (selvars <- paste(twinvar,c(rep("_1",nv),rep("_2",nv)),sep=""))
  (modvars <- modvar)
  #(modvars <- paste(modvar,c(rep("_1",4)),sep="tr"))
  #modvars <- paste(modvar,c(rep("_1",4),rep("_2",4)),sep="")
  
  # Create Separate Data Sets for MZ / DZ Twins by sex
  MZMdata <- as.data.frame(subset(currdat, zygSex==1,c(selvars,modvars))) # MZ Male Twins
  DZMdata <- as.data.frame(subset(currdat, zygSex==3,c(selvars,modvars))) # DZ Male Twins
  MZFdata <- as.data.frame(subset(currdat, zygSex==2,c(selvars,modvars))) # MZ Female Twins
  DZFdata <- as.data.frame(subset(currdat, zygSex==4,c(selvars,modvars))) # DZ Female Twins
  DZOdata <- as.data.frame(subset(currdat, zygSex==5,c(selvars,modvars))) # DZ Opp Sex Twins
  
  
  # Removing Twins with Missing Moderator Values   #
  ##################################################
  # only variable 1 is needed to check
  # var 2 is also based on age & everyone has country data
  MZMdata <- MZMdata[is.na(MZMdata[,modvars[1]])==0,] 
  DZMdata <- DZMdata[is.na(DZMdata[,modvars[1]])==0,]
  MZFdata <- MZFdata[is.na(MZFdata[,modvars[1]])==0,] 
  DZFdata <- DZFdata[is.na(DZFdata[,modvars[1]])==0,]
  DZOdata <- DZOdata[is.na(DZOdata[,modvars[1]])==0,]
  
  # Remove twins with missing education variable
  MZMdata <- MZMdata[is.na(MZMdata[,selvars[3]])==0,] 
  MZMdata <- MZMdata[is.na(MZMdata[,selvars[1]])==0,] 
  DZMdata <- DZMdata[is.na(DZMdata[,selvars[3]])==0,]
  DZMdata <- DZMdata[is.na(DZMdata[,selvars[1]])==0,] 
  MZFdata <- MZFdata[is.na(MZFdata[,selvars[3]])==0,] 
  MZFdata <- MZFdata[is.na(MZFdata[,selvars[1]])==0,] 
  DZFdata <- DZFdata[is.na(DZFdata[,selvars[3]])==0,]
  DZFdata <- DZFdata[is.na(DZFdata[,selvars[1]])==0,] 	
  DZOdata <- DZOdata[is.na(DZOdata[,selvars[3]])==0,]
  DZOdata <- DZOdata[is.na(DZOdata[,selvars[1]])==0,] 
  
  ########################################
  ### recalculate full MZ and DZ pairs ###
  ### (commented out when not needed)  ###
  ########################################
  # MZMdatafull <- subset(MZMdata, is.na(MZMdata$zgrams_sqrt_windsor_1)==F & is.na(MZMdata$zgrams_sqrt_windsor_2)==F)
  # MZFdatafull <- subset(MZFdata, is.na(MZFdata$zgrams_sqrt_windsor_1)==F & is.na(MZFdata$zgrams_sqrt_windsor_2)==F)
  # DZMdatafull <- subset(DZMdata, is.na(DZMdata$zgrams_sqrt_windsor_1)==F & is.na(DZMdata$zgrams_sqrt_windsor_2)==F)
  # DZFdatafull <- subset(DZFdata, is.na(DZFdata$zgrams_sqrt_windsor_1)==F & is.na(DZFdata$zgrams_sqrt_windsor_2)==F)
  # DZOdatafull <- subset(DZOdata, is.na(DZOdata$zgrams_sqrt_windsor_1)==F & is.na(DZOdata$zgrams_sqrt_windsor_2)==F)
  # (FullMZM <- describe(MZMdatafull$zgrams_sqrt_windsor_1)) #4505
  # (FullMZF <- describe(MZFdatafull$zgrams_sqrt_windsor_1)) #4968
  # (FullDZM <- describe(DZMdatafull$zgrams_sqrt_windsor_1)) #3484
  # (FullDZF <- describe(DZFdatafull$zgrams_sqrt_windsor_1)) #4080
  # (FullDZO <- describe(DZOdatafull$zgrams_sqrt_windsor_1)) #6744
  # print(full_mz_pairs <- FullMZM$n +FullMZF$n) #9473 mz pairs
  # print(full_dz_pairs <- FullDZM$n +FullDZF$n +FullDZO$n) #14308 dz pairs
  # print(full_pairs <- full_mz_pairs+full_dz_pairs) #23656 full pairs
  # 

  ######################################
  #    Print Descriptive Statistics    #
  ######################################
  # Male Twins
  # colMeans(MZMdata,na.rm=TRUE)
  # colMeans(DZMdata,na.rm=TRUE)
  # cor(MZMdata,use="complete")
  # cor(DZMdata,use="complete")
  # 
  # # Female Twins
  # colMeans(MZFdata,na.rm=TRUE)
  # colMeans(DZFdata,na.rm=TRUE)
  # cor(MZFdata,use="complete")
  # cor(DZFdata,use="complete")
  # 
  # # OS Twins
  # colMeans(DZOdata,na.rm=TRUE)
  # cor(DZOdata,use="complete")
  # 
  
  
univModACEModel <- mxModel("univModACE",
 mxModel("ACE",
         # Matrices a, c, and e to store a, c, and e path coefficients
         # Males
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.2,.1,.7), name="aM", lbound=-20, ubound=20),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.2,.1,.2), name="cM", lbound=-30, ubound=30),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.7,-.41,.7), name="eM", lbound=-20, ubound=20),
         # Females
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.2,.1,.7), name="aF", lbound=-20, ubound=20),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.2,.1,.2), name="cF", lbound=-20, ubound=20),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.7,.1,.7), name="eF", lbound=-20, ubound=20),
         
         # Compute heritability at mean
         mxAlgebra( expression=aM %*% t(aM), name="AM" ),
         mxAlgebra( expression=cM %*% t(cM), name="CM" ),
         mxAlgebra( expression=eM %*% t(eM), name="EM" ),
         mxAlgebra( expression=aF %*% t(aF), name="AF" ),
         mxAlgebra( expression=cF %*% t(cF), name="CF" ),
         mxAlgebra( expression=eF %*% t(eF), name="EF" ),
         # Algebra to compute total variances and standard deviations (diagonal only)
         mxAlgebra( expression=AM+CM+EM, name="VM" ),
         mxAlgebra( expression=AF+CF+EF, name="VF" ),
         mxAlgebra( expression=diag2vec(AM/VM),name="stndAM"),
         mxAlgebra( expression=diag2vec(CM/VM),name="stndCM"),
         mxAlgebra( expression=diag2vec(EM/VM),name="stndEM"),
         mxAlgebra( expression=diag2vec(AF/VF),name="stndAF"),
         mxAlgebra( expression=diag2vec(CF/VF),name="stndCF"),
         mxAlgebra( expression=diag2vec(EF/VF),name="stndEF"),
         
         # Matrices a, c, and e to store moderated a, c, and e path coefficients
         # Males 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="aIM1", lbound=-10, ubound=10 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="cIM1", lbound=-20, ubound=20 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="eIM1", lbound=-10, ubound=10 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="aIM2", lbound=-10, ubound=10 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="cIM2", lbound=-10, ubound=10 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="eIM2", lbound=-10, ubound=10 ),
         # Females
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="aIF1", lbound=-10, ubound=10 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="cIF1", lbound=-10, ubound=10 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="eIF1", lbound=-10, ubound=10 ),         
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="aIF2", lbound=-10, ubound=10 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="cIF2", lbound=-10, ubound=10 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="eIF2", lbound=-10, ubound=10 ),          
         
         # Matrix & Algebra for expected means vector
         mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values= c(.50,.20), label=c("mean_M_ed","mean_M_alc"), name="muM", lbound=-20, ubound=40  ),
         #  mxMatrix( type="Full", nrow=1, ncol=nv, free=c(F,F), values=c(0,0), label=c("bM1_ed","bM1_alc"), name="bM1", lbound=-20, ubound=20 ),
         mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(.20,-.1), label=c("bM2_ed","bM2_alc"), name="bM2", lbound=-20, ubound=20 ),
         # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM3_ed","bM3_alc"), name="bM3", lbound=-20, ubound=20 ),
         # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM4_ed","bM4_alc"), name="bM4", lbound=-20, ubound=20 ),
         # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM5_ed","bM5_alc"), name="bM5", lbound=-20, ubound=20 ),
         
         mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(.20,0), label=c("mean_F_ed","mean_F_alc"), name="muF" , lbound=-20, ubound=40 ),
         #  mxMatrix( type="Full", nrow=1, ncol=nv, free=c(F,F), values=c(0,0), label=c("bF1_ed","bF1_alc"), name="bF1", lbound=-20, ubound=20 ),
         mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,-1), label=c("bF2_ed","bF2_alc"), name="bF2", lbound=-20, ubound=20 ),
         # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bF3_ed","bF3_alc"), name="bF3", lbound=-20, ubound=20 ),
         # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bF4_ed","bF4_alc"), name="bF4", lbound=-20, ubound=20 ),
         # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bF5_ed","bF5_alc"), name="bF5", lbound=-20, ubound=20 ),
         
         # Keep at least one of these in the model even if CI=F (so comma above doesn't give error)
         #  mxCI(c("bM2","bF2","bM3","bF3","bM4","bF4","bM5","bF5",
         #         "aIM1","cIM1","eIM1","aIM2","cIM2","eIM2",
         #         "aIF1","cIF1","eIF1","aIF2","cIF2","eIF2",
         #         "aM","cM","eM","aF","cF","eF","muM","muF"))
         mxCI(c("aM","cM","eM","aF","cF","eF","muM","muF"))
         
         #  mxCI(c("aIM1","cIM1","eIM1","aIM2","cIM2","eIM2",
         #         "aIF1","cIF1","eIF1","aIF2","cIF2","eIF2"))
 ),
 
 #############
 mxModel("MZM",
         # Matrix for Moderating Variable
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[1],sep=""), name="modA1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[3],sep=""), name="modB1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[1],sep=""), name="mod2"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[2],sep=""), name="mod3"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[3],sep=""), name="mod4"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[4],sep=""), name="mod5"),
         
         # Matrices A, C, and E compute Moderated variance components (moderators for education and age only - not country)
         mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM1" ),
         mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM1" ),
         mxAlgebra((ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM1" ),
         
         mxAlgebra((ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM2" ),
         mxAlgebra((ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM2" ),
         mxAlgebra((ACE.eM + modB1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modB1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM2" ),
         
         mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM12" ),
         mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM12" ),
         
         # Algebra for expected variance/covariance matrix and expected mean vector in MZ
         mxAlgebra(rbind ( cbind(AM1+CM1+EM1 , AM12+CM12),
                           cbind(AM12+CM12  , AM2+CM2+EM2)), name="expCovMZM" ),
         
         mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanAM"),
         mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanBM"),
         mxAlgebra(cbind(meanAM,meanBM), name="expMeanM"),
         
         # Data & Objective
         mxData(observed=MZMdata, type="raw"),
         mxExpectationNormal( covariance="expCovMZM", means="expMeanM", dimnames=selvars),   
         mxFitFunctionML()
 ),
 
 
 #############    
 mxModel("DZM", 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[1],sep=""), name="modA1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[3],sep=""), name="modB1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[1],sep=""), name="mod2"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[2],sep=""), name="mod3"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[3],sep=""), name="mod4"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[4],sep=""), name="mod5"),
         
         # Matrices A, C, and E compute variance components
         mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM1" ),
         mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM1" ),
         mxAlgebra((ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM1" ),
         
         mxAlgebra((ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM2" ),
         mxAlgebra((ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM2" ),
         mxAlgebra((ACE.eM + modB1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modB1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM2" ),
         
         mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM12" ),
         mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM12" ),
         
         # Algebra for expected variance/covariance matrix in DZ
         mxAlgebra(rbind ( cbind(AM1+CM1+EM1      , 0.5%x%AM12+CM12),
                           cbind(0.5%x%AM12+CM12 , AM2+CM2+EM2)),  name="expCovDZM" ),
         
         mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanAM"),
         mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanBM"),
         mxAlgebra(cbind(meanAM,meanBM), name="expMeanM"),
         
         # Data & Objective
         mxData(observed=DZMdata, type="raw"),
         mxExpectationNormal( covariance="expCovDZM", means="expMeanM", dimnames=selvars),   
         mxFitFunctionML()
 ),
 
 #############
 mxModel("MZF",
         # Matrix for Moderating Variable
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[1],sep=""), name="modA1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[3],sep=""), name="modB1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[1],sep=""), name="mod2"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[2],sep=""), name="mod3"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[3],sep=""), name="mod4"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[4],sep=""), name="mod5"),
         
         # Matrices A, C, and E compute Moderated variance components
         mxAlgebra((ACE.aF + modA1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modA1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF1" ),
         mxAlgebra((ACE.cF + modA1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modA1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF1" ),
         mxAlgebra((ACE.eF + modA1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2) %*% t(ACE.eF + modA1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2), name="EF1" ),
         
         mxAlgebra((ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF2" ),
         mxAlgebra((ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF2" ),
         mxAlgebra((ACE.eF + modB1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2) %*% t(ACE.eF + modB1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2), name="EF2" ),
         
         mxAlgebra((ACE.aF + modA1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF12" ),
         mxAlgebra((ACE.cF + modA1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF12" ),
         
         # Algebra for expected variance/covariance matrix and expected mean vector in MZ
         mxAlgebra(rbind ( cbind(AF1+CF1+EF1 , AF12+CF12),
                           cbind(AF12+CF12  , AF2+CF2+EF2)), name="expCovMZF" ),
         
         mxAlgebra(ACE.muF + ACE.bF2 %*% mod2, name="meanAF"),
         mxAlgebra(ACE.muF + ACE.bF2 %*% mod2, name="meanBF"), 
         mxAlgebra(cbind(meanAF,meanBF), name="expMeanF"),
         
         # Data & Objective
         mxData(observed=MZFdata, type="raw"),
         mxExpectationNormal( covariance="expCovMZF", means="expMeanF", dimnames=selvars),   
         mxFitFunctionML()
 ),
 
 
 #############    
 mxModel("DZF", 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[1],sep=""), name="modA1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[3],sep=""), name="modB1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[1],sep=""), name="mod2"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[2],sep=""), name="mod3"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[3],sep=""), name="mod4"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[4],sep=""), name="mod5"),
         
         # Matrices A, C, and E compute variance components
         mxAlgebra((ACE.aF + modA1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modA1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF1" ),
         mxAlgebra((ACE.cF + modA1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modA1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF1" ),
         mxAlgebra((ACE.eF + modA1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2) %*% t(ACE.eF + modA1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2), name="EF1" ),
         
         mxAlgebra((ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF2" ),
         mxAlgebra((ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF2" ),
         mxAlgebra((ACE.eF + modB1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2) %*% t(ACE.eF + modB1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2), name="EF2" ),
         
         mxAlgebra((ACE.aF + modA1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF12" ),
         mxAlgebra((ACE.cF + modA1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF12" ),
         
         # Algebra for expected variance/covariance matrix in DZ
         mxAlgebra(rbind ( cbind(AF1+CF1+EF1      , 0.5%x%AF12+CF12),
                           cbind(0.5%x%AF12+CF12 , AF2+CF2+EF2)),  name="expCovDZF" ),
         
         mxAlgebra(ACE.muF + ACE.bF2 %*% mod2, name="meanAF"),
         mxAlgebra(ACE.muF + ACE.bF2 %*% mod2, name="meanBF"), 
         mxAlgebra(cbind(meanAF,meanBF), name="expMeanF"),
         
         # Data & Objective
         mxData(observed=DZFdata, type="raw"),
         mxExpectationNormal( covariance="expCovDZF", means="expMeanF", dimnames=selvars),   
         mxFitFunctionML()
 ),
 
 mxModel("DZO", 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[1],sep=""), name="modA1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[3],sep=""), name="modB1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[1],sep=""), name="mod2"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[2],sep=""), name="mod3"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[3],sep=""), name="mod4"), 
         # mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[4],sep=""), name="mod5"),
         
         # Matrices A, C, and E compute variance components
         mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM1" ),
         mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM1" ),
         mxAlgebra((ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM1" ),
         
         mxAlgebra((ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF2" ),
         mxAlgebra((ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF2" ),
         mxAlgebra((ACE.eF + modB1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2) %*% t(ACE.eF + modB1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2), name="EF2" ),
         
         mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AO12" ),
         mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CO12" ),
         
         # Algebra for expected variance/covariance matrix in OSDZ
         mxAlgebra(rbind ( cbind(AM1+CM1+EM1      , 0.5%x%AO12+CO12),
                           cbind(0.5%x%AO12+CO12 , AF2+CF2+EF2)),  name="expCovDZO" ),
         
         mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanAM"),
         mxAlgebra(ACE.muF + ACE.bF2 %*% mod2, name="meanBF"), 
         mxAlgebra(cbind(meanAM,meanBF), name="expMeanO"),
         
         # Data & Objective
         mxData(observed=DZOdata, type="raw"),
         mxExpectationNormal( covariance="expCovDZO", means="expMeanO", dimnames=selvars),   
         mxFitFunctionML()
 ),  
 
 mxAlgebra( expression=MZM.objective + DZM.objective + MZF.objective + DZF.objective + DZO.objective, name="neg2sumll" ),
 mxFitFunctionMultigroup(c("MZM","DZM","MZF","DZF","DZO"))
  )
  
  univModACEFit <- mxRun(univModACEModel,intervals=F)
  # univModACEFit <- mxRun(univModACEFit,intervals=F)
  # univModACEFit <- mxTryHard(univModACEFit,intervals=F)
  univModACESumm <- summary(univModACEFit)
  univModACESumm
  
  ### need to add final parts to save this output
  outfilename1 <- paste("output/modelout_univModACE_Fit_", paste(currsample, "2024", sep="_"), sep="")
  outfilename2 <- paste("output/modelout_univModACE_Summ_", paste(currsample, "2024", sep="_"), sep="")
  save(univModACEFit, file=outfilename1)
  save(univModACESumm, file=outfilename2)
  
}


loop_no_OSDZ <- c("MTSADA","OCTOTWIN","SATSA")
#loop_no_OSDZ <- c("SATSA")

for (i in loop_no_OSDZ){
  print(i)
  currsample <- i
  currdat <- filter(newtwins, SAMPLE==i)

  twinvar <- c('zISCED_trim', 'zgrams_sqrt_windsor') #variable of interest
  #modvar <- c('zage_c','Country_x1','Country_x2','Country_x3') #moderator variable
  modvar <- c('zage_c') #moderator variable
  
  # Defining Variables for OpenMx
  nv <- 2
  ntv <- nv*2
  (selvars <- paste(twinvar,c(rep("_1",nv),rep("_2",nv)),sep=""))
  (modvars <- modvar)

  # Create Separate Data Sets for MZ / DZ Twins by sex
  MZMdata <- as.data.frame(subset(currdat, zygSex==1,c(selvars,modvars))) # MZ Male Twins
  DZMdata <- as.data.frame(subset(currdat, zygSex==3,c(selvars,modvars))) # DZ Male Twins
  MZFdata <- as.data.frame(subset(currdat, zygSex==2,c(selvars,modvars))) # MZ Female Twins
  DZFdata <- as.data.frame(subset(currdat, zygSex==4,c(selvars,modvars))) # DZ Female Twins

  # Removing Twins with Missing Moderator Values   #
  ##################################################
  # only variable 1 is needed to check
  # var 2 is also based on age & everyone has country data
  MZMdata <- MZMdata[is.na(MZMdata[,modvars[1]])==0,] 
  DZMdata <- DZMdata[is.na(DZMdata[,modvars[1]])==0,]
  MZFdata <- MZFdata[is.na(MZFdata[,modvars[1]])==0,] 
  DZFdata <- DZFdata[is.na(DZFdata[,modvars[1]])==0,]

  # Remove twins with missing education variable
  MZMdata <- MZMdata[is.na(MZMdata[,selvars[3]])==0,] 
  MZMdata <- MZMdata[is.na(MZMdata[,selvars[1]])==0,] 
  DZMdata <- DZMdata[is.na(DZMdata[,selvars[3]])==0,]
  DZMdata <- DZMdata[is.na(DZMdata[,selvars[1]])==0,] 
  MZFdata <- MZFdata[is.na(MZFdata[,selvars[3]])==0,] 
  MZFdata <- MZFdata[is.na(MZFdata[,selvars[1]])==0,] 
  DZFdata <- DZFdata[is.na(DZFdata[,selvars[3]])==0,]
  DZFdata <- DZFdata[is.na(DZFdata[,selvars[1]])==0,] 	

    
  univModACEModel <- mxModel("univModACE",
   mxModel("ACE",
           # Matrices a, c, and e to store a, c, and e path coefficients
           # Males
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.7,.3,.2), name="aM", lbound=-20, ubound=20),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.1,.3,.2), name="cM", lbound=-30, ubound=30),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.7,.3,.6), name="eM", lbound=-20, ubound=20),
           # Females
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.7,.3,.6), name="aF", lbound=-20, ubound=20),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.2,.3,.2), name="cF", lbound=-20, ubound=20),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.7,.3,.6), name="eF", lbound=-20, ubound=20),
           
           # Compute heritability at mean
           mxAlgebra( expression=aM %*% t(aM), name="AM" ),
           mxAlgebra( expression=cM %*% t(cM), name="CM" ),
           mxAlgebra( expression=eM %*% t(eM), name="EM" ),
           mxAlgebra( expression=aF %*% t(aF), name="AF" ),
           mxAlgebra( expression=cF %*% t(cF), name="CF" ),
           mxAlgebra( expression=eF %*% t(eF), name="EF" ),
           # Algebra to compute total variances and standard deviations (diagonal only)
           mxAlgebra( expression=AM+CM+EM, name="VM" ),
           mxAlgebra( expression=AF+CF+EF, name="VF" ),
           mxAlgebra( expression=diag2vec(AM/VM),name="stndAM"),
           mxAlgebra( expression=diag2vec(CM/VM),name="stndCM"),
           mxAlgebra( expression=diag2vec(EM/VM),name="stndEM"),
           mxAlgebra( expression=diag2vec(AF/VF),name="stndAF"),
           mxAlgebra( expression=diag2vec(CF/VF),name="stndCF"),
           mxAlgebra( expression=diag2vec(EF/VF),name="stndEF"),
           
           # Matrices a, c, and e to store moderated a, c, and e path coefficients
           # Males 
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="aIM1", lbound=-10, ubound=10 ),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="cIM1", lbound=-20, ubound=20 ),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="eIM1", lbound=-10, ubound=10 ),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="aIM2", lbound=-10, ubound=10 ),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="cIM2", lbound=-10, ubound=10 ),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.20,0), name="eIM2", lbound=-10, ubound=10 ),
           # Females
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="aIF1", lbound=-10, ubound=10 ),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="cIF1", lbound=-10, ubound=10 ),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="eIF1", lbound=-10, ubound=10 ),         
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="aIF2", lbound=-10, ubound=10 ),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="cIF2", lbound=-10, ubound=10 ),
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="eIF2", lbound=-10, ubound=10 ),          
           
           # Matrix & Algebra for expected means vector
           mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values= c(-1,0), label=c("mean_M_ed","mean_M_alc"), name="muM", lbound=-20, ubound=40  ),
           #  mxMatrix( type="Full", nrow=1, ncol=nv, free=c(F,F), values=c(0,0), label=c("bM1_ed","bM1_alc"), name="bM1", lbound=-20, ubound=20 ),
           mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM2_ed","bM2_alc"), name="bM2", lbound=-20, ubound=20 ),
           # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM3_ed","bM3_alc"), name="bM3", lbound=-20, ubound=20 ),
           # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM4_ed","bM4_alc"), name="bM4", lbound=-20, ubound=20 ),
           # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM5_ed","bM5_alc"), name="bM5", lbound=-20, ubound=20 ),
           
           mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("mean_F_ed","mean_F_alc"), name="muF" , lbound=-20, ubound=40 ),
           #  mxMatrix( type="Full", nrow=1, ncol=nv, free=c(F,F), values=c(0,0), label=c("bF1_ed","bF1_alc"), name="bF1", lbound=-20, ubound=20 ),
           mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bF2_ed","bF2_alc"), name="bF2", lbound=-20, ubound=20 ),
           # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bF3_ed","bF3_alc"), name="bF3", lbound=-20, ubound=20 ),
           # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bF4_ed","bF4_alc"), name="bF4", lbound=-20, ubound=20 ),
           # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bF5_ed","bF5_alc"), name="bF5", lbound=-20, ubound=20 ),
           
           # Keep at least one of these in the model even if CI=F (so comma above doesn't give error)
           #  mxCI(c("bM2","bF2","bM3","bF3","bM4","bF4","bM5","bF5",
           #         "aIM1","cIM1","eIM1","aIM2","cIM2","eIM2",
           #         "aIF1","cIF1","eIF1","aIF2","cIF2","eIF2",
           #         "aM","cM","eM","aF","cF","eF","muM","muF"))
           mxCI(c("aM","cM","eM","aF","cF","eF","muM","muF"))
           
           #  mxCI(c("aIM1","cIM1","eIM1","aIM2","cIM2","eIM2",
           #         "aIF1","cIF1","eIF1","aIF2","cIF2","eIF2"))
   ),
   
   #############
   mxModel("MZM",
           # Matrix for Moderating Variable
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[1],sep=""), name="modA1"), 
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[3],sep=""), name="modB1"), 
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[1],sep=""), name="mod2"), 

           # Matrices A, C, and E compute Moderated variance components (moderators for education and age only - not country)
           mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM1" ),
           mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM1" ),
           mxAlgebra((ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM1" ),
           
           mxAlgebra((ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM2" ),
           mxAlgebra((ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM2" ),
           mxAlgebra((ACE.eM + modB1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modB1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM2" ),
           
           mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM12" ),
           mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM12" ),
           
           # Algebra for expected variance/covariance matrix and expected mean vector in MZ
           mxAlgebra(rbind ( cbind(AM1+CM1+EM1 , AM12+CM12),
                             cbind(AM12+CM12  , AM2+CM2+EM2)), name="expCovMZM" ),
           
           mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanAM"),
           mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanBM"),
           mxAlgebra(cbind(meanAM,meanBM), name="expMeanM"),
           
           # Data & Objective
           mxData(observed=MZMdata, type="raw"),
           mxExpectationNormal( covariance="expCovMZM", means="expMeanM", dimnames=selvars),   
           mxFitFunctionML()
   ),
   
   
   #############    
   mxModel("DZM", 
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[1],sep=""), name="modA1"), 
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[3],sep=""), name="modB1"), 
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[1],sep=""), name="mod2"), 

           # Matrices A, C, and E compute variance components
           mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM1" ),
           mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM1" ),
           mxAlgebra((ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM1" ),
           
           mxAlgebra((ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM2" ),
           mxAlgebra((ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM2" ),
           mxAlgebra((ACE.eM + modB1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modB1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM2" ),
           
           mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM12" ),
           mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM12" ),
           
           # Algebra for expected variance/covariance matrix in DZ
           mxAlgebra(rbind ( cbind(AM1+CM1+EM1      , 0.5%x%AM12+CM12),
                             cbind(0.5%x%AM12+CM12 , AM2+CM2+EM2)),  name="expCovDZM" ),
           
           mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanAM"),
           mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanBM"),
           mxAlgebra(cbind(meanAM,meanBM), name="expMeanM"),
           
           # Data & Objective
           mxData(observed=DZMdata, type="raw"),
           mxExpectationNormal( covariance="expCovDZM", means="expMeanM", dimnames=selvars),   
           mxFitFunctionML()
   ),
   
   #############
   mxModel("MZF",
           # Matrix for Moderating Variable
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[1],sep=""), name="modA1"), 
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[3],sep=""), name="modB1"), 
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[1],sep=""), name="mod2"), 

           # Matrices A, C, and E compute Moderated variance components
           mxAlgebra((ACE.aF + modA1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modA1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF1" ),
           mxAlgebra((ACE.cF + modA1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modA1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF1" ),
           mxAlgebra((ACE.eF + modA1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2) %*% t(ACE.eF + modA1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2), name="EF1" ),
           
           mxAlgebra((ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF2" ),
           mxAlgebra((ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF2" ),
           mxAlgebra((ACE.eF + modB1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2) %*% t(ACE.eF + modB1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2), name="EF2" ),
           
           mxAlgebra((ACE.aF + modA1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF12" ),
           mxAlgebra((ACE.cF + modA1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF12" ),
           
           # Algebra for expected variance/covariance matrix and expected mean vector in MZ
           mxAlgebra(rbind ( cbind(AF1+CF1+EF1 , AF12+CF12),
                             cbind(AF12+CF12  , AF2+CF2+EF2)), name="expCovMZF" ),
           
           mxAlgebra(ACE.muF + ACE.bF2 %*% mod2, name="meanAF"),
           mxAlgebra(ACE.muF + ACE.bF2 %*% mod2, name="meanBF"), 
           mxAlgebra(cbind(meanAF,meanBF), name="expMeanF"),
           
           # Data & Objective
           mxData(observed=MZFdata, type="raw"),
           mxExpectationNormal( covariance="expCovMZF", means="expMeanF", dimnames=selvars),   
           mxFitFunctionML()
   ),
   
   
   #############    
   mxModel("DZF", 
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[1],sep=""), name="modA1"), 
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[3],sep=""), name="modB1"), 
           mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[1],sep=""), name="mod2"), 

           # Matrices A, C, and E compute variance components
           mxAlgebra((ACE.aF + modA1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modA1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF1" ),
           mxAlgebra((ACE.cF + modA1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modA1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF1" ),
           mxAlgebra((ACE.eF + modA1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2) %*% t(ACE.eF + modA1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2), name="EF1" ),
           
           mxAlgebra((ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF2" ),
           mxAlgebra((ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF2" ),
           mxAlgebra((ACE.eF + modB1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2) %*% t(ACE.eF + modB1 %*% ACE.eIF1 + mod2 %*% ACE.eIF2), name="EF2" ),
           
           mxAlgebra((ACE.aF + modA1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2) %*% t(ACE.aF + modB1 %*% ACE.aIF1 + mod2 %*% ACE.aIF2), name="AF12" ),
           mxAlgebra((ACE.cF + modA1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2) %*% t(ACE.cF + modB1 %*% ACE.cIF1 + mod2 %*% ACE.cIF2), name="CF12" ),
           
           # Algebra for expected variance/covariance matrix in DZ
           mxAlgebra(rbind ( cbind(AF1+CF1+EF1      , 0.5%x%AF12+CF12),
                             cbind(0.5%x%AF12+CF12 , AF2+CF2+EF2)),  name="expCovDZF" ),
           
           mxAlgebra(ACE.muF + ACE.bF2 %*% mod2, name="meanAF"),
           mxAlgebra(ACE.muF + ACE.bF2 %*% mod2, name="meanBF"), 
           mxAlgebra(cbind(meanAF,meanBF), name="expMeanF"),
           
           # Data & Objective
           mxData(observed=DZFdata, type="raw"),
           mxExpectationNormal( covariance="expCovDZF", means="expMeanF", dimnames=selvars),   
           mxFitFunctionML()
   ),
   
   mxAlgebra( expression=MZM.objective + DZM.objective + MZF.objective + DZF.objective, name="neg2sumll" ),
   mxFitFunctionMultigroup(c("MZM","DZM","MZF","DZF"))
  )
  
  univModACEFit <- mxRun(univModACEModel,intervals=F)
  # univModACEFit <- mxRun(univModACEFit,intervals=T)
  #univModACEFit <- mxTryHard(univModACEFit,intervals=F)
  univModACESumm <- summary(univModACEFit)
  univModACESumm
  
  ### need to add final parts to save this output
  outfilename1 <- paste("output/modelout_univModACE_Fit_", paste(currsample, "2024", sep="_"), sep="")
  outfilename2 <- paste("output/modelout_univModACE_Summ_", paste(currsample, "2024", sep="_"), sep="")
  save(univModACEFit, file=outfilename1)
  save(univModACESumm, file=outfilename2)
  
}



loop_no_female <- c("VETSA","NASNRC")

for (i in loop_no_female){
  print(i)
  currsample <- i
  currdat <- filter(newtwins, SAMPLE==i)
  
  twinvar <- c('zISCED_trim', 'zgrams_sqrt_windsor') #variable of interest
  #modvar <- c('zage_c','Country_x1','Country_x2','Country_x3') #moderator variable
  modvar <- c('zage_c') #moderator variable
  
  # Defining Variables for OpenMx
  nv <- 2
  ntv <- nv*2
  (selvars <- paste(twinvar,c(rep("_1",nv),rep("_2",nv)),sep=""))
  (modvars <- modvar)
  
  # Create Separate Data Sets for MZ / DZ Twins by sex
  MZMdata <- as.data.frame(subset(currdat, zygSex==1,c(selvars,modvars))) # MZ Male Twins
  DZMdata <- as.data.frame(subset(currdat, zygSex==3,c(selvars,modvars))) # DZ Male Twins

  # Removing Twins with Missing Moderator Values   #
  ##################################################
  # only variable 1 is needed to check
  # var 2 is also based on age & everyone has country data
  MZMdata <- MZMdata[is.na(MZMdata[,modvars[1]])==0,] 
  DZMdata <- DZMdata[is.na(DZMdata[,modvars[1]])==0,]

  # Remove twins with missing education variable
  MZMdata <- MZMdata[is.na(MZMdata[,selvars[3]])==0,] 
  MZMdata <- MZMdata[is.na(MZMdata[,selvars[1]])==0,] 
  DZMdata <- DZMdata[is.na(DZMdata[,selvars[3]])==0,]
  DZMdata <- DZMdata[is.na(DZMdata[,selvars[1]])==0,] 

univModACEModel <- mxModel("univModACE",
 mxModel("ACE",
         # Matrices a, c, and e to store a, c, and e path coefficients
         # Males
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.7,.1,.2), name="aM", lbound=-20, ubound=20),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.1,.3,.2), name="cM", lbound=-30, ubound=30),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.7,.1,.6), name="eM", lbound=-20, ubound=20),
  
         # Matrices a, c, and e to store moderated a, c, and e path coefficients
         # Males 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="aIM1", lbound=-10, ubound=10 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="cIM1", lbound=-20, ubound=20 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=c(0,.0,0), name="eIM1", lbound=-10, ubound=10 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="aIM2", lbound=-10, ubound=10 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="cIM2", lbound=-10, ubound=10 ),
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=c(0,.0,0), name="eIM2", lbound=-10, ubound=10 ),
      
         # Compute heritability at mean
         mxAlgebra( expression=aM %*% t(aM), name="AM" ),
         mxAlgebra( expression=cM %*% t(cM), name="CM" ),
         mxAlgebra( expression=eM %*% t(eM), name="EM" ),
         # Algebra to compute total variances and standard deviations (diagonal only)
         mxAlgebra( expression=AM+CM+EM, name="VM" ),
         mxAlgebra( expression=diag2vec(AM/VM),name="stndAM"),
         mxAlgebra( expression=diag2vec(CM/VM),name="stndCM"),
         mxAlgebra( expression=diag2vec(EM/VM),name="stndEM"),
         
         # Matrix & Algebra for expected means vector
         mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values= c(0,0), label=c("mean_M_ed","mean_M_alc"), name="muM", lbound=-20, ubound=40  ),
         #  mxMatrix( type="Full", nrow=1, ncol=nv, free=c(F,F), values=c(0,0), label=c("bM1_ed","bM1_alc"), name="bM1", lbound=-20, ubound=20 ),
         mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM2_ed","bM2_alc"), name="bM2", lbound=-20, ubound=20 ),
         # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM3_ed","bM3_alc"), name="bM3", lbound=-20, ubound=20 ),
         # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM4_ed","bM4_alc"), name="bM4", lbound=-20, ubound=20 ),
         # mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM5_ed","bM5_alc"), name="bM5", lbound=-20, ubound=20 ),
    
         # Keep at least one of these in the model even if CI=F (so comma above doesn't give error)
         #  mxCI(c("bM2","bF2","bM3","bF3","bM4","bF4","bM5","bF5",
         #         "aIM1","cIM1","eIM1","aIM2","cIM2","eIM2",
         #         "aIF1","cIF1","eIF1","aIF2","cIF2","eIF2",
         #         "aM","cM","eM","aF","cF","eF","muM","muF"))
         mxCI(c("aM","cM","eM","muM"))

         #  mxCI(c("aIM1","cIM1","eIM1","aIM2","cIM2","eIM2",
         #         "aIF1","cIF1","eIF1","aIF2","cIF2","eIF2"))
 ),
 
 #############
 mxModel("MZM",
         # Matrix for Moderating Variable
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[1],sep=""), name="modA1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[3],sep=""), name="modB1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[1],sep=""), name="mod2"), 
         
         # Matrices A, C, and E compute Moderated variance components (moderators for education and age only - not country)
         mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM1" ),
         mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM1" ),
         mxAlgebra((ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM1" ),
         
         mxAlgebra((ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM2" ),
         mxAlgebra((ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM2" ),
         mxAlgebra((ACE.eM + modB1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modB1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM2" ),
         
         mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM12" ),
         mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM12" ),
         
         # Algebra for expected variance/covariance matrix and expected mean vector in MZ
         mxAlgebra(rbind ( cbind(AM1+CM1+EM1 , AM12+CM12),
                           cbind(AM12+CM12  , AM2+CM2+EM2)), name="expCovMZM" ),
         
         mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanAM"),
         mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanBM"),
         mxAlgebra(cbind(meanAM,meanBM), name="expMeanM"),
         
         # Data & Objective
         mxData(observed=MZMdata, type="raw"),
         mxExpectationNormal( covariance="expCovMZM", means="expMeanM", dimnames=selvars),   
         mxFitFunctionML()
 ),
 
 
 #############    
 mxModel("DZM", 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[1],sep=""), name="modA1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[3],sep=""), name="modB1"), 
         mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[1],sep=""), name="mod2"), 
         
         # Matrices A, C, and E compute variance components
         mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM1" ),
         mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM1" ),
         mxAlgebra((ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modA1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM1" ),
         
         mxAlgebra((ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM2" ),
         mxAlgebra((ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM2" ),
         mxAlgebra((ACE.eM + modB1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2) %*% t(ACE.eM + modB1 %*% ACE.eIM1 + mod2 %*% ACE.eIM2), name="EM2" ),
         
         mxAlgebra((ACE.aM + modA1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2) %*% t(ACE.aM + modB1 %*% ACE.aIM1 + mod2 %*% ACE.aIM2), name="AM12" ),
         mxAlgebra((ACE.cM + modA1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2) %*% t(ACE.cM + modB1 %*% ACE.cIM1 + mod2 %*% ACE.cIM2), name="CM12" ),
         
         # Algebra for expected variance/covariance matrix in DZ
         mxAlgebra(rbind ( cbind(AM1+CM1+EM1      , 0.5%x%AM12+CM12),
                           cbind(0.5%x%AM12+CM12 , AM2+CM2+EM2)),  name="expCovDZM" ),
         
         mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanAM"),
         mxAlgebra(ACE.muM + ACE.bM2 %*% mod2, name="meanBM"),
         mxAlgebra(cbind(meanAM,meanBM), name="expMeanM"),
         
         # Data & Objective
         mxData(observed=DZMdata, type="raw"),
         mxExpectationNormal( covariance="expCovDZM", means="expMeanM", dimnames=selvars),   
         mxFitFunctionML()
 ),

 mxAlgebra( expression=MZM.objective + DZM.objective, name="neg2sumll" ),
 mxFitFunctionMultigroup(c("MZM","DZM"))
)
  
  univModACEFit <- mxRun(univModACEModel,intervals=F)
  # univModACEFit <- mxRun(univModACEFit,intervals=T)
  # univModACEFit <- mxTryHard(univModACEFit,intervals=F)
  univModACESumm <- summary(univModACEFit)
  univModACESumm
  

  ### need to add final parts to save this output
  outfilename1 <- paste("output/modelout_univModACE_Fit_", paste(currsample, "2024", sep="_"), sep="")
  outfilename2 <- paste("output/modelout_univModACE_Summ_", paste(currsample, "2024", sep="_"), sep="")
  save(univModACEFit, file=outfilename1)
  save(univModACESumm, file=outfilename2)
  
  
}

# Standard errors for standardized ACEs
mxSE(ACE.stndAM,model= univModACEFit)
mxSE(ACE.stndCM,model= univModACEFit)
mxSE(ACE.stndEM,model= univModACEFit)


