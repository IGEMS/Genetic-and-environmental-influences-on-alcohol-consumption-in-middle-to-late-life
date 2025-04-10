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
source("miFunctions.R")

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
describe(newtwins$age_c) # M = 56.07 SD = 11.29
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
#load("Genetic_Data_Prepped.RData")

### #IF RUNNING ON THE CLUSTER - DO THIS ####
# sinteractive --time=8:00:00 --ntasks=8 --part=ainteractive
# ml anaconda
# conda activate /projects/lessem/software/anaconda/envs/R-latest
# R
# library(OpenMx)

########################################
#          Prepare for OpenMx          #
########################################

# Defining Variables of Interest
twinvar <- c('zISCED_trim', 'zgrams_sqrt_windsor') #variable of interest
modvar <- c('zage_c','Country_x1','Country_x2','Country_x3') #moderator variable

# Defining Variables for OpenMx
nv <- 2
ntv <- nv*2
(selvars <- paste(twinvar,c(rep("_1",nv),rep("_2",nv)),sep=""))
(modvars <- modvar)
#(modvars <- paste(modvar,c(rep("_1",4)),sep="tr"))
#modvars <- paste(modvar,c(rep("_1",4),rep("_2",4)),sep="")

# Create Separate Data Sets for MZ / DZ Twins by sex
MZMdata <- as.data.frame(subset(newtwins, zygSex==1,c(selvars,modvars))) # MZ Male Twins
DZMdata <- as.data.frame(subset(newtwins, zygSex==3,c(selvars,modvars))) # DZ Male Twins
MZFdata <- as.data.frame(subset(newtwins, zygSex==2,c(selvars,modvars))) # MZ Female Twins
DZFdata <- as.data.frame(subset(newtwins, zygSex==4,c(selvars,modvars))) # DZ Female Twins
DZOdata <- as.data.frame(subset(newtwins, zygSex==5,c(selvars,modvars))) # DZ Opp Sex Twins

############################################
### Calc full MZ pairs and full DZ pairs ###
###   (commented out when not needed)    ###
############################################
# MZMdatafull <- subset(MZMdata, is.na(MZMdata$zgrams_sqrt_windsor_1)==F & is.na(MZMdata$zgrams_sqrt_windsor_2)==F)
# MZFdatafull <- subset(MZFdata, is.na(MZFdata$zgrams_sqrt_windsor_1)==F & is.na(MZFdata$zgrams_sqrt_windsor_2)==F)
# DZMdatafull <- subset(DZMdata, is.na(DZMdata$zgrams_sqrt_windsor_1)==F & is.na(DZMdata$zgrams_sqrt_windsor_2)==F)
# DZFdatafull <- subset(DZFdata, is.na(DZFdata$zgrams_sqrt_windsor_1)==F & is.na(DZFdata$zgrams_sqrt_windsor_2)==F)
# DZOdatafull <- subset(DZOdata, is.na(DZOdata$zgrams_sqrt_windsor_1)==F & is.na(DZOdata$zgrams_sqrt_windsor_2)==F)
# describe(MZMdatafull$zgrams_sqrt_windsor_1) #5381
# describe(MZFdatafull$zgrams_sqrt_windsor_1) #6179
# describe(DZMdatafull$zgrams_sqrt_windsor_1) #3612
# describe(DZFdatafull$zgrams_sqrt_windsor_1) #4277
# describe(DZOdatafull$zgrams_sqrt_windsor_1) #7087
# 5381+6179 #11560 mz pairs
# 3612+4277+7087 #14976 dz pairs
# 11560+14976 #26536 full pairs

##################################################
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

# ########################################
# ### recalculate full MZ and DZ pairs ###
# ### (commented out when not needed)  ###
# ########################################
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
# (full_mz_pairs <- FullMZM$n +FullMZF$n) #9473 mz pairs
# (full_dz_pairs <- FullDZM$n +FullDZF$n +FullDZO$n) #14308 dz pairs
# (full_pairs <- full_mz_pairs+full_dz_pairs) #23656 full pairs


######################################
#    Print Descriptive Statistics    #
######################################
# Male Twins
colMeans(MZMdata,na.rm=TRUE)
colMeans(DZMdata,na.rm=TRUE)
cor(MZMdata,use="complete")
cor(DZMdata,use="complete")

# Female Twins
colMeans(MZFdata,na.rm=TRUE)
colMeans(DZFdata,na.rm=TRUE)
cor(MZFdata,use="complete")
cor(DZFdata,use="complete")

# OS Twins
colMeans(DZOdata,na.rm=TRUE)
cor(DZOdata,use="complete")

# Set default Mx optimizer
#mxOption(NULL, "Default optimizer", "NPSOL") 
mxOption(NULL, "Default optimizer", "SLSQP") 

# -----------------------------------------------------------------------	
# Fit ACE Model with Means and Variance Components moderation effects
# -----------------------------------------------------------------------
# NOTE: Confidence intervals took a lot of computation time so they were fit in separate models on the server

univModACEModel <- mxModel("univModACE",
mxModel("ACE",
     # Matrices a, c, and e to store a, c, and e path coefficients
     # Males
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.7,.1,.5), name="aM", lbound=-20, ubound=20),
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.1,.3,.2), name="cM", lbound=-30, ubound=30),
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.7,.1,.6), name="eM", lbound=-20, ubound=20),
     # Females
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.7,.1,.5), name="aF", lbound=-20, ubound=20),
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.1,.3,.2), name="cF", lbound=-20, ubound=20),
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=c(.7,.1,.6), name="eF", lbound=-20, ubound=20),
     
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
     mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values= c(0,0), label=c("mean_M_ed","mean_M_alc"), name="muM", lbound=-20, ubound=40  ),
   #  mxMatrix( type="Full", nrow=1, ncol=nv, free=c(F,F), values=c(0,0), label=c("bM1_ed","bM1_alc"), name="bM1", lbound=-20, ubound=20 ),
     mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM2_ed","bM2_alc"), name="bM2", lbound=-20, ubound=20 ),
     mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM3_ed","bM3_alc"), name="bM3", lbound=-20, ubound=20 ),
     mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM4_ed","bM4_alc"), name="bM4", lbound=-20, ubound=20 ),
     mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bM5_ed","bM5_alc"), name="bM5", lbound=-20, ubound=20 ),
   
     mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("mean_F_ed","mean_F_alc"), name="muF" , lbound=-20, ubound=40 ),
   #  mxMatrix( type="Full", nrow=1, ncol=nv, free=c(F,F), values=c(0,0), label=c("bF1_ed","bF1_alc"), name="bF1", lbound=-20, ubound=20 ),
     mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bF2_ed","bF2_alc"), name="bF2", lbound=-20, ubound=20 ),
     mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bF3_ed","bF3_alc"), name="bF3", lbound=-20, ubound=20 ),
     mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bF4_ed","bF4_alc"), name="bF4", lbound=-20, ubound=20 ),
     mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), label=c("bF5_ed","bF5_alc"), name="bF5", lbound=-20, ubound=20 ),
   
     # Keep at least one of these in the model even if CI=F (so comma above doesn't give error)
      #  mxCI(c("bM2","bF2","bM3","bF3","bM4","bF4","bM5","bF5",
      #         "aIM1","cIM1","eIM1","aIM2","cIM2","eIM2",
      #         "aIF1","cIF1","eIF1","aIF2","cIF2","eIF2",
      #         "aM","cM","eM","aF","cF","eF","muM","muF"))
   
        mxCI(c("aIM1","cIM1","eIM1","aIM2","cIM2","eIM2",
               "aIF1","cIF1","eIF1","aIF2","cIF2","eIF2"))
      #   mxCI(c("aM","cM","eM","aF","cF","eF","muM","muF"))
     ),

#############
mxModel("MZM",
     # Matrix for Moderating Variable
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[1],sep=""), name="modA1"), 
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",selvars[3],sep=""), name="modB1"), 
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[1],sep=""), name="mod2"), 
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[2],sep=""), name="mod3"), 
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[3],sep=""), name="mod4"), 
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[4],sep=""), name="mod5"),
     
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
     
     mxAlgebra(ACE.muM + ACE.bM2 %*% mod2 + ACE.bM3 %*% mod3+ ACE.bM4 %*% mod4 + ACE.bM5 %*% mod5, name="meanAM"),
     mxAlgebra(ACE.muM + ACE.bM2 %*% mod2 + ACE.bM3 %*% mod3+ ACE.bM4 %*% mod4 + ACE.bM5 %*% mod5, name="meanBM"),
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
      mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[2],sep=""), name="mod3"), 
      mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[3],sep=""), name="mod4"), 
      mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[4],sep=""), name="mod5"),
        
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
     
     mxAlgebra(ACE.muM + ACE.bM2 %*% mod2 + ACE.bM3 %*% mod3+ ACE.bM4 %*% mod4 + ACE.bM5 %*% mod5, name="meanAM"),
     mxAlgebra(ACE.muM + ACE.bM2 %*% mod2 + ACE.bM3 %*% mod3+ ACE.bM4 %*% mod4 + ACE.bM5 %*% mod5, name="meanBM"),
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
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[2],sep=""), name="mod3"), 
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[3],sep=""), name="mod4"), 
     mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[4],sep=""), name="mod5"),
     
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
     
     mxAlgebra(ACE.muF + ACE.bF2 %*% mod2 + ACE.bF3 %*% mod3+ ACE.bF4 %*% mod4 + ACE.bF5 %*% mod5, name="meanAF"),
     mxAlgebra(ACE.muF + ACE.bF2 %*% mod2 + ACE.bF3 %*% mod3+ ACE.bF4 %*% mod4 + ACE.bF5 %*% mod5, name="meanBF"), 
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
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[2],sep=""), name="mod3"), 
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[3],sep=""), name="mod4"), 
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[4],sep=""), name="mod5"),
        
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
     
     mxAlgebra(ACE.muF + ACE.bF2 %*% mod2 + ACE.bF3 %*% mod3+ ACE.bF4 %*% mod4 + ACE.bF5 %*% mod5, name="meanAF"),
     mxAlgebra(ACE.muF + ACE.bF2 %*% mod2 + ACE.bF3 %*% mod3+ ACE.bF4 %*% mod4 + ACE.bF5 %*% mod5, name="meanBF"), 
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
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[2],sep=""), name="mod3"), 
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[3],sep=""), name="mod4"), 
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, labels=paste("data.",modvars[4],sep=""), name="mod5"),
        
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
     
     mxAlgebra(ACE.muM + ACE.bM2 %*% mod2 + ACE.bM3 %*% mod3+ ACE.bM4 %*% mod4 + ACE.bM5 %*% mod5, name="meanAM"),
     mxAlgebra(ACE.muF + ACE.bF2 %*% mod2 + ACE.bF3 %*% mod3+ ACE.bF4 %*% mod4 + ACE.bF5 %*% mod5, name="meanBF"), 
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
univModACEFit <- mxRun(univModACEFit,intervals=T)
univModACEFit <- mxTryHard(univModACEFit,intervals=T)
univModACESumm <- summary(univModACEFit)
univModACESumm


save(univModACEFit, file="univModACEFit_4_21_2023.RData")



#############################
### Collapse 1 at a time ###
############################
## Collapse A effects across Sex
collapseMF_Aonly_sex <- mxRename(univModACEFit, "collapseMF_Aonly_sex")
collapseMF_Aonly_sex$ACE.aM <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=2, name="aM" , label=c("a11","a12","a22"), lbound=-20, ubound=20 )
collapseMF_Aonly_sex$ACE.aF <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=2, name="aF" , label=c("a11","a12","a22"), lbound=-20, ubound=20 )
collapseMF_Aonly_sex_Fit <- mxRun(collapseMF_Aonly_sex,intervals=F)
#collapseMF_Aonly_sex_Fit <- mxTryHard(collapseMF_Aonly_sex,intervals=F)
noASumm <- summary(collapseMF_Aonly_sex_Fit)
tableFitStatistics(univModACEFit,c(collapseMF_Aonly_sex_Fit))

collapseMF_A1 <- mxRename(univModACEFit, "collapseMF_A1")
collapseMF_A1$ACE.aIM1 <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=0, name="aIM1", label=c("aI1_11","aI1_12","aI1_22"), lbound=-10, ubound=10 )
collapseMF_A1$ACE.aIF1 <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=0, name="aIF1", label=c("aI1_11","aI1_12","aI1_22"), lbound=-10, ubound=10 )
collapseA1Fit <- mxRun(collapseMF_A1,intervals=F)
noASumm1 <- summary(collapseA1Fit)
tableFitStatistics(univModACEFit,c(collapseA1Fit))

collapseMF_A2y <- mxRename(univModACEFit, "collapseMF_A2y")
collapseMF_A2y$ACE.aIM2y <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=0, name="aIM2y", label=c("aI2y_11m","aI2y_12","aI2y_22"), lbound=-10, ubound=10 )
collapseMF_A2y$ACE.aIF2y <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=0, name="aIF2y", label=c("aI2y_11f","aI2y_12","aI2y_22"), lbound=-10, ubound=10 )
collapseA2yFit <- mxRun(collapseMF_A2y,intervals=F)
noASumm2y <- summary(collapseA2yFit)
tableFitStatistics(univModACEFit,c(collapseA2yFit))

collapseMF_A2o <- mxRename(univModACEFit, "collapseMF_A2o")
collapseMF_A2o$ACE.aIM2o <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=0, name="aIM2o", label=c("aI2o_11m","aI2o_12","aI2o_22"), lbound=-10, ubound=10 )
collapseMF_A2o$ACE.aIF2o <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=0, name="aIF2o", label=c("aI2o_11f","aI2o_12","aI2o_22"), lbound=-10, ubound=10 )
collapseA2oFit <- mxRun(collapseMF_A2o,intervals=F)
noASumm2o <- summary(collapseA2oFit)
tableFitStatistics(univModACEFit,c(collapseA2oFit))


## Collapse C effects across Sex
collapseMF_Conly_sex <- mxRename(univModACEFit, "collapseMF_Conly_sex")
collapseMF_Conly_sex$ACE.cM <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0, name="cM" , label=c("c11","c12","c22"), lbound=-20, ubound=20 )
collapseMF_Conly_sex$ACE.cF <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0, name="cF" , label=c("c11","c12","c22"), lbound=-20, ubound=20 )
collapseMF_Conly_sex_Fit <- mxRun(collapseMF_Conly_sex,intervals=F)
noCSumm <- summary(collapseMF_Conly_sex_Fit)
tableFitStatistics(univModACEFit,c(collapseMF_Conly_sex_Fit))

collapseMF_C1 <- mxRename(univModACEFit, "collapseMF_C1")
collapseMF_C1$ACE.cIM1 <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=0, name="cIM1", label=c("cI1_11","cI1_12","cI1_22"), lbound=-10, ubound=10 )
collapseMF_C1$ACE.cIF1 <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=0, name="cIF1", label=c("cI1_11","cI1_12","cI1_22"), lbound=-10, ubound=10 )
collapseC1Fit <- mxRun(collapseMF_C1,intervals=F)
noCSumm <- summary(collapseC1Fit)
tableFitStatistics(univModACEFit,c(collapseC1Fit))

collapseMF_C2y <- mxRename(univModACEFit, "collapseMF_C2y")
collapseMF_C2y$ACE.cIM2y <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=0, name="cIM2y", label=c("cI2y_11m","cI2y_12","cI2y_22"), lbound=-10, ubound=10 )
collapseMF_C2y$ACE.cIF2y <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=0, name="cIF2y", label=c("cI2y_11f","cI2y_12","cI2y_22"), lbound=-10, ubound=10 )
collapseC2yFit <- mxRun(collapseMF_C2y,intervals=F)
noCSumm2y <- summary(collapseC2yFit)
tableFitStatistics(univModACEFit,c(collapseC2yFit))

collapseMF_C2o <- mxRename(univModACEFit, "collapseMF_C2o")
collapseMF_C2o$ACE.cIM2o <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=0, name="cIM2", label=c("cI2o_11m","cI2o_12","cI2o_22"), lbound=-10, ubound=10 )
collapseMF_C2o$ACE.cIF2o <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=0, name="cIF2", label=c("cI2o_11f","cI2o_12","cI2o_22"), lbound=-10, ubound=10 )
collapseC2oFit <- mxRun(collapseMF_C2o,intervals=F)
noCSumm2o <- summary(collapseC2oFit)
tableFitStatistics(univModACEFit,c(collapseC2oFit))



collapseMF_Eonly_sex <- mxRename(univModACEFit, "collapseMF_Eonly_sex")
collapseMF_Eonly_sex$ACE.eM <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0, name="eM" , label=c("e11","e12","e22"), lbound=-20, ubound=20 )
collapseMF_Eonly_sex$ACE.eF <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0, name="eF" , label=c("e11","e12","e22"), lbound=-20, ubound=20 )
collapseMF_Eonly_sex_Fit <- mxRun(collapseMF_Eonly_sex,intervals=F)
#collapseMF_Eonly_sex_Fit <- mxTryHard(collapseMF_Eonly_sex,intervals=F)
noESumm <- summary(collapseMF_Eonly_sex_Fit)
tableFitStatistics(univModACEFit,c(collapseMF_Eonly_sex_Fit))

## Collapse E effects across Sex
collapseMF_E1 <- mxRename(univModACEFit, "collapseMF_E1")
collapseMF_E1$ACE.eIM1 <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=0, name="eIM1", label=c("eI1_11","eI1_12","eI1_22"), lbound=-10, ubound=10 )
collapseMF_E1$ACE.eIF1 <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,T,T), values=0, name="eIF1", label=c("eI1_11","eI1_12","eI1_22"), lbound=-10, ubound=10 )
collapseE1Fit <- mxRun(collapseMF_E1,intervals=F)
#collapseE2Fit <- mxTryHard(collapseMF_E1,intervals=F)
noESumm <- summary(collapseE1Fit)
tableFitStatistics(univModACEFit,c(collapseE1Fit))

collapseMF_E2y <- mxRename(univModACEFit, "collapseMF_E2y")
collapseMF_E2y$ACE.eIM2y <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=0, name="eIM2y", label=c("eI2y_11m","eI2y_12","eI2y_22"), lbound=-10, ubound=10 )
collapseMF_E2y$ACE.eIF2y <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=0, name="eIF2y", label=c("eI2y_11f","eI2y_12","eI2y_22"), lbound=-10, ubound=10 )
collapseE2yFit <- mxRun(collapseMF_E2y,intervals=F)
noESumm2y <- summary(collapseE2yFit)
tableFitStatistics(univModACEFit,c(collapseE2yFit))

collapseMF_E2o <- mxRename(univModACEFit, "collapseMF_E2o")
collapseMF_E2o$ACE.eIM2o <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=0, name="eIM2o", label=c("eI2o_11m","eI2o_12","eI2o_22"), lbound=-10, ubound=10 )
collapseMF_E2o$ACE.eIF2o <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,T,T), values=0, name="eIF2o", label=c("eI2o_11f","eI2o_12","eI2o_22"), lbound=-10, ubound=10 )
collapseE2oFit <- mxRun(collapseMF_E2o,intervals=F)
noESumm2o <- summary(collapseE2oFit)
tableFitStatistics(univModACEFit,c(collapseE2oFit))


## Collapse mean effects across Sex
collapseMF_mean <- mxRename(univModACEFit, "collapseMF_mean")
collapseMF_mean$ACE.muM <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), name="muM", label=c("mean_M_ed","mean_alc"), lbound=-10, ubound=10 )
collapseMF_mean$ACE.muF <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), name="muF", label=c("mean_F_ed","mean_alc"), lbound=-10, ubound=10 )
collapseMF_meanFit <- mxRun(collapseMF_mean,intervals=F)
noMeanSumm <- summary(collapseMF_meanFit)
#noASumm 
tableFitStatistics(univModACEFit,c(collapseMF_meanFit))

# collapseMF_b1 <- mxRename(univModACEFit, "collapseMF_b1")
# collapseMF_b1$ACE.bM1 <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(F,T), values=c(0,1), name="bM1", label=c("bM1_ed","b1_alc"), lbound=-10, ubound=10 )
# collapseMF_b1$ACE.bF1 <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(F,T), values=c(0,1), name="bF1", label=c("bF1_ed","b1_alc"), lbound=-10, ubound=10 )
# collapseMF_b1Fit <- mxRun(collapseMF_b1,intervals=F)
# noMeanSumm <- summary(collapseMF_b1Fit)
# #noASumm 
# tableFitStatistics(univModACEFit,c(collapseMF_b1Fit))

collapseMF_b2y <- mxRename(univModACEFit, "collapseMF_b2y")
collapseMF_b2y$ACE.bM2y <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), name="bM2y", label=c("bM2y_ed","b2y_alc"), lbound=-10, ubound=10 )
collapseMF_b2y$ACE.bF2y <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), name="bF2y", label=c("bF2y_ed","b2y_alc"), lbound=-10, ubound=10 )
collapseMF_b2yFit <- mxRun(collapseMF_b2y,intervals=F)
noMeanSumm2y <- summary(collapseMF_b2yFit)
tableFitStatistics(univModACEFit,c(collapseMF_b2yFit))

collapseMF_b2o <- mxRename(univModACEFit, "collapseMF_b2o")
collapseMF_b2o$ACE.bM2o <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), name="bM2o", label=c("bM2o_ed","b2o_alc"), lbound=-10, ubound=10 )
collapseMF_b2o$ACE.bF2o <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), name="bF2o", label=c("bF2o_ed","b2o_alc"), lbound=-10, ubound=10 )
collapseMF_b2oFit <- mxRun(collapseMF_b2o,intervals=F)
noMeanSumm2o <- summary(collapseMF_b2oFit)
tableFitStatistics(univModACEFit,c(collapseMF_b2oFit))


collapseMF_b3 <- mxRename(univModACEFit, "collapseMF_b3")
collapseMF_b3$ACE.bM3 <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), name="bM3", label=c("bM3_ed","b3_alc"), lbound=-10, ubound=10 )
collapseMF_b3$ACE.bF3 <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), name="bF3", label=c("bF3_ed","b3_alc"), lbound=-10, ubound=10 )
collapseMF_b3Fit <- mxRun(collapseMF_b3,intervals=F)
noMeanSumm3 <- summary(collapseMF_b3Fit)
tableFitStatistics(univModACEFit,c(collapseMF_b3Fit))

collapseMF_b4 <- mxRename(univModACEFit, "collapseMF_b4")
collapseMF_b4$ACE.bM4 <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), name="bM4", label=c("bM4_ed","b4_alc"), lbound=-10, ubound=10 )
collapseMF_b4$ACE.bF4 <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), name="bF4", label=c("bF4_ed","b4_alc"), lbound=-10, ubound=10 )
collapseMF_b4Fit <- mxRun(collapseMF_b4,intervals=F)
noMeanSumm4 <- summary(collapseMF_b4Fit)
tableFitStatistics(univModACEFit,c(collapseMF_b4Fit))

# collapseMF_b5 <- mxRename(univModACEFit, "collapseMF_b5")
# collapseMF_b5$ACE.bM5 <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), name="bM5", label=c("bM5_ed","b5_alc"), lbound=-10, ubound=10 )
# collapseMF_b5$ACE.bF5 <- mxMatrix( type="Full", nrow=1, ncol=nv, free=c(T,T), values=c(0,0), name="bF5", label=c("bF5_ed","b5_alc"), lbound=-10, ubound=10 )
# collapseMF_b5Fit <- mxRun(collapseMF_b5,intervals=F)
# noMeanSumm5 <- summary(collapseMF_b5Fit)
# tableFitStatistics(univModACEFit,c(collapseMF_b5Fit))


collapseMF_ACE_sex <- mxRename(univModACEFit, "collapseMF_ACE_sex")
collapseMF_ACE_sex$ACE.aM <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0, name="aM" , label=c("a11","a12","a22"), lbound=-20, ubound=20 )
collapseMF_ACE_sex$ACE.aF <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0, name="aF" , label=c("a11","a12","a22"), lbound=-20, ubound=20 )
collapseMF_ACE_sex$ACE.cM <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0, name="cM" , label=c("c11","c12","c22"), lbound=-20, ubound=20 )
collapseMF_ACE_sex$ACE.cF <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0, name="cF" , label=c("c11","c12","c22"), lbound=-20, ubound=20 )
collapseMF_ACE_sex$ACE.eM <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0, name="eM" , label=c("e11","e12","e22"), lbound=-20, ubound=20 )
collapseMF_ACE_sex$ACE.eF <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0, name="eF" , label=c("e11","e12","e22"), lbound=-20, ubound=20 )
collapseMF_ACE_sex_Fit <- mxTryHard(collapseMF_ACE_sex,intervals=F)
noACESumm <- summary(collapseMF_ACE_sex_Fit)
tableFitStatistics(univModACEFit,c(collapseMF_ACE_sex_Fit))


### Collapsing across ACE for sex
tableFitStatistics(univModACEFit,c(collapseMF_Aonly_sex_Fit, 
                                   collapseMF_Conly_sex_Fit, 
                                   collapseMF_Eonly_sex_Fit, 
                                   collapseMF_ACE_sex_Fit))
### Collapsingn effects of education across sex
tableFitStatistics(univModACEFit,c(collapseA1Fit, collapseA2Fit, collapseA3Fit, 
                                   collapseC1Fit, collapseC2Fit, collapseC3Fit, 
                                   collapseE1Fit, collapseE2Fit, collapseE3Fit,
                                   collapseACE1Fit, collapseACE2Fit, 
                                   collapseMF_meanFit, collapseMF_b2Fit, collapseMF_b3Fit, 
                                   collapseMF_b4Fit, collapseMF_b5Fit))


##########################
#### Tests within Sex ####
##########################
## Drop rA for education
noAModel_M <- mxRename(univModACEFit, "noAModel_M")
noAModel_M$ACE.aM <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,T), values=0, name="aM", lbound=-10, ubound=10 )
noAFit_M <- mxRun(noAModel_M,intervals=F)
noASumm_M <- summary(noAFit_M)
tableFitStatistics(univModACEFit,c(noAFit_M))

noAModel_F <- mxRename(univModACEFit, "noAModel_F")
noAModel_F$ACE.aF <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,T), values=0, name="aF", lbound=-10, ubound=10 )
noAFit_F <- mxRun(noAModel_F,intervals=F)
noASumm_F <- summary(noAFit_F)
tableFitStatistics(univModACEFit,c(noAFit_F))

# No A mod education
noA1Model_M <- mxRename(univModACEFit, "noA1Model_M")
noA1Model_M$ACE.aIM1 <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,F,F), values=0, name="aIM1", lbound=-10, ubound=10 )
noA1Fit_M <- mxRun(noA1Model_M,intervals=F)
noA1Summ_M <- summary(noA1Fit_M)
tableFitStatistics(univModACEFit,c(noA1Fit_M))

noA1Model_F <- mxRename(univModACEFit, "noA1Model_F")
noA1Model_F$ACE.aIF1 <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,F,F), values=0, name="aIF1", lbound=-10, ubound=10 )
noA1Fit_F <- mxRun(noA1Model_F,intervals=F)
noA1Summ_F <- summary(noA1Fit_F)
tableFitStatistics(univModACEFit,c(noA1Fit_F))


# No A mod age
noA2yModel_M <- mxRename(univModACEFit, "noA2yModel_M")
noA2yModel_M$ACE.aIM2y <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F), values=0, name="aIM2y", lbound=-10, ubound=10 )
noA2yFit_M <- mxRun(noA2yModel_M,intervals=F)
noA2ySumm_M <- summary(noA2yFit_M)
tableFitStatistics(univModACEFit,c(noA2yFit_M))

noA2yModel_F <- mxRename(univModACEFit, "noA2yModel_F")
noA2yModel_F$ACE.aIF2y <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F), values=0, name="aIF2y", lbound=-10, ubound=10 )
noA2yFit_F <- mxRun(noA2yModel_F,intervals=F)
noA2ySumm_F <- summary(noA2yFit_F)
tableFitStatistics(univModACEFit,c(noAy2Fit_F))

noA2oModel_M <- mxRename(univModACEFit, "noA2oModel_M")
noA2oModel_M$ACE.aIM2o <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F), values=0, name="aIM2o", lbound=-10, ubound=10 )
noA2oFit_M <- mxRun(noA2oModel_M,intervals=F)
noA2oSumm_M <- summary(noA2oFit_M)
tableFitStatistics(univModACEFit,c(noA2oFit_M))

noA2oModel_F <- mxRename(univModACEFit, "noA2oModel_F")
noA2oModel_F$ACE.aIF2o <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F), values=0, name="aIF2o", lbound=-10, ubound=10 )
noA2oFit_F <- mxRun(noA2oModel_F,intervals=F)
noA2oSumm_F <- summary(noA2oFit_F)
tableFitStatistics(univModACEFit,c(noA2oFit_F))


######## C ########
## Drop rC for education
noCModel_M <- mxRename(univModACEFit, "noCModel_M")
noCModel_M$ACE.cM <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,T), values=0, name="cM", lbound=-10, ubound=10 )
noCFit_M <- mxRun(noCModel_M,intervals=F)
noCSumm_M <- summary(noCFit_M)
tableFitStatistics(univModACEFit,c(noCFit_M))

noCModel_F <- mxRename(univModACEFit, "noCModel_F")
noCModel_F$ACE.cF <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,T), values=0, name="cF", lbound=-10, ubound=10 )
noCFit_F <- mxRun(noCModel_F,intervals=F)
noCSumm_F <- summary(noCFit_F)
tableFitStatistics(univModACEFit,c(noCFit_F))

# No C mod education
noC1Model_M <- mxRename(univModACEFit, "noC1Model_M")
noC1Model_M$ACE.cIM1 <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,F,F), values=0, name="cIM1", lbound=-10, ubound=10 )
noC1Fit_M <- mxRun(noC1Model_M,intervals=F)
noC1Summ_M <- summary(noC1Fit_M)
tableFitStatistics(univModACEFit,c(noC1Fit_M))

noC1Model_F <- mxRename(univModACEFit, "noC1Model_F")
noC1Model_F$ACE.cIF1 <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,F,F), values=0, name="cIF1", lbound=-10, ubound=10 )
noC1Fit_F <- mxRun(noC1Model_F,intervals=F)
noC1Summ_F <- summary(noC1Fit_F)
tableFitStatistics(univModACEFit,c(noC1Fit_F))


# No C mod ageY
noC2yModel_M <- mxRename(univModACEFit, "noC2yModel_M")
noC2yModel_M$ACE.cIM2y <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F), values=0, name="cIM2y", lbound=-10, ubound=10 )
noC2yFit_M <- mxRun(noC2yModel_M,intervals=F)
noC2ySumm_M <- summary(noC2yFit_M)
tableFitStatistics(univModACEFit,c(noC2yFit_M))

noC2yModel_F <- mxRename(univModACEFit, "noC2yModel_F")
noC2yModel_F$ACE.cIF2y <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F), values=0, name="cIF2y", lbound=-10, ubound=10 )
noC2yFit_F <- mxRun(noC2yModel_F,intervals=F)
noC2ySumm_F <- summary(noC2yFit_F)
tableFitStatistics(univModACEFit,c(noC2yFit_F))

noC2oModel_M <- mxRename(univModACEFit, "noC2oModel_M")
noC2oModel_M$ACE.cIM2o <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F), values=0, name="cIM2o", lbound=-10, ubound=10 )
noC2oFit_M <- mxRun(noC2oModel_M,intervals=F)
noC2oSumm_M <- summary(noC2oFit_M)
tableFitStatistics(univModACEFit,c(noC2oFit_M))

noC2oModel_F <- mxRename(univModACEFit, "noC2oModel_F")
noC2oModel_F$ACE.cIF2o <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F), values=0, name="cIF2o", lbound=-10, ubound=10 )
noC2oFit_F <- mxRun(noC2oModel_F,intervals=F)
noC2oSumm_F <- summary(noC2oFit_F)
tableFitStatistics(univModACEFit,c(noC2oFit_F))


######## E ########
## Drop rE for education
noEModel_M <- mxRename(univModACEFit, "noEModel_M")
noEModel_M$ACE.eM <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,T), values=0, name="eM", lbound=-10, ubound=10 )
noEFit_M <- mxRun(noEModel_M,intervals=F)
noESumm_M <- summary(noEFit_M)
tableFitStatistics(univModACEFit,c(noEFit_M))

noEModel_F <- mxRename(univModACEFit, "noEModel_F")
noEModel_F$ACE.eF <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,T), values=0, name="eF", lbound=-10, ubound=10 )
noEFit_F <- mxRun(noEModel_F,intervals=F)
noESumm_F <- summary(noEFit_F)
tableFitStatistics(univModACEFit,c(noEFit_F))

# No E mod education
noE1Model_M <- mxRename(univModACEFit, "noE1Model_M")
noE1Model_M$ACE.eIM1 <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,F,F), values=0, name="eIM1", lbound=-10, ubound=10 )
noE1Fit_M <- mxRun(noE1Model_M,intervals=F)
noE1Summ_M <- summary(noE1Fit_M)
tableFitStatistics(univModACEFit,c(noE1Fit_M))

noE1Model_F <- mxRename(univModACEFit, "noE1Model_F")
noE1Model_F$ACE.eIF1 <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(F,F,F), values=0, name="eIF1", lbound=-10, ubound=10 )
noE1Fit_F <- mxRun(noE1Model_F,intervals=F)
noE1Summ_F <- summary(noE1Fit_F)
tableFitStatistics(univModACEFit,c(noE1Fit_F))


# No E mod ageY
noE2yModel_M <- mxRename(univModACEFit, "noE2yModel_M")
noE2yModel_M$ACE.eIM2y <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F), values=0, name="eIM2y", lbound=-10, ubound=10 )
noE2yFit_M <- mxRun(noE2yModel_M,intervals=F)
noE2ySumm_M <- summary(noE2yFit_M)
tableFitStatistics(univModACEFit,c(noE2yFit_M))

noE2yModel_F <- mxRename(univModACEFit, "noE2yModel_F")
noE2yModel_F$ACE.eIF2y <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F), values=0, name="eIF2y", lbound=-10, ubound=10 )
noE2yFit_F <- mxRun(noE2yModel_F,intervals=F)
noE2ySumm_F <- summary(noE2yFit_F)
tableFitStatistics(univModACEFit,c(noE2yFit_F))

noE2oModel_M <- mxRename(univModACEFit, "noE2oModel_M")
noE2oModel_M$ACE.eIM2o <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F), values=0, name="eIM2o", lbound=-10, ubound=10 )
noE2oFit_M <- mxRun(noE2oModel_M,intervals=F)
noE2oSumm_M <- summary(noE2oFit_M)
tableFitStatistics(univModACEFit,c(noE2oFit_M))

noE2oModel_F <- mxRename(univModACEFit, "noE2oModel_F")
noE2oModel_F$ACE.eIF2o <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F), values=0, name="eIF2o", lbound=-10, ubound=10 )
noE2oFit_F <- mxRun(noE2oModel_F,intervals=F)
noE2oSumm_F <- summary(noE2oFit_F)
tableFitStatistics(univModACEFit,c(noE2oFit_F))


tableFitStatistics(univModACEFit,c(noAFit_F,  noAFit_M,  noCFit_F,  noCFit_M,  noEFit_F,  noEFit_M,
                                   noA1Fit_F, noA1Fit_M, noC1Fit_F, noC1Fit_M, noE1Fit_F, noE1Fit_M,
                                   noA2yFit_F, noA2yFit_M, noC2yFit_F, noC2yFit_M, noE2yFit_F, noE2yFit_M,
                                   noA2oFit_F, noA2oFit_M, noC2oFit_F, noC2oFit_M, noE2oFit_F, noE2oFit_M,
                                   noA3Fit_F, noA3Fit_M, noC3Fit_F, noC3Fit_M, noE3Fit_F, noE3Fit_M))



###################################################################
#####            View results ran on the server               #####
###################################################################
load(file="univModACEFit_4_26_2023.RData")
univModACESumm <- summary(univModACEFit)
univModACESumm

# Get estimates and SEs for the total ACEs
(AF <- mxEval(expression=ACE.aF %*% t(ACE.aF), model=univModACEFit))
(CF <- mxEval(expression=ACE.cF %*% t(ACE.cF), model=univModACEFit))
(EF <- mxEval(expression=ACE.eF %*% t(ACE.eF), model=univModACEFit))
(AM <- mxEval(expression=ACE.aM %*% t(ACE.aM), model=univModACEFit))
(CM <- mxEval(expression=ACE.cM %*% t(ACE.cM), model=univModACEFit))
(EM <- mxEval(expression=ACE.eM %*% t(ACE.eM), model=univModACEFit))
(VF <- AF+CF+EF)
(VM <- AM+CM+EM)
(stndAF <- mxEval(expression=diag2vec((ACE.aF %*% t(ACE.aF))/(ACE.aF %*% t(ACE.aF)+ACE.cF %*% t(ACE.cF)+ACE.eF %*% t(ACE.eF))), model=univModACEFit)[2])
(stndCF <- mxEval(expression=diag2vec((ACE.cF %*% t(ACE.cF))/(ACE.aF %*% t(ACE.aF)+ACE.cF %*% t(ACE.cF)+ACE.eF %*% t(ACE.eF))), model=univModACEFit)[2])
(stndEF <- mxEval(expression=diag2vec((ACE.eF %*% t(ACE.eF))/(ACE.aF %*% t(ACE.aF)+ACE.cF %*% t(ACE.cF)+ACE.eF %*% t(ACE.eF))), model=univModACEFit)[2])
(stndAM <- mxEval(expression=diag2vec((ACE.aM %*% t(ACE.aM))/(ACE.aM %*% t(ACE.aM)+ACE.cM %*% t(ACE.cM)+ACE.eM %*% t(ACE.eM))), model=univModACEFit)[2])
(stndCM <- mxEval(expression=diag2vec((ACE.cM %*% t(ACE.cM))/(ACE.aM %*% t(ACE.aM)+ACE.cM %*% t(ACE.cM)+ACE.eM %*% t(ACE.eM))), model=univModACEFit)[2])
(stndEM <- mxEval(expression=diag2vec((ACE.eM %*% t(ACE.eM))/(ACE.aM %*% t(ACE.aM)+ACE.cM %*% t(ACE.cM)+ACE.eM %*% t(ACE.eM))), model=univModACEFit)[2])

# Now get SEs for these computations
(se_stndAF <- mxSE(diag2vec((ACE.aF %*% t(ACE.aF))/(ACE.aF %*% t(ACE.aF)+ACE.cF %*% t(ACE.cF)+ACE.eF %*% t(ACE.eF))),model= univModACEFit)[2])
(se_stndCF <- mxSE(diag2vec((ACE.cF %*% t(ACE.cF))/(ACE.aF %*% t(ACE.aF)+ACE.cF %*% t(ACE.cF)+ACE.eF %*% t(ACE.eF))),model= univModACEFit)[2])
(se_stndEF <- mxSE(diag2vec((ACE.eF %*% t(ACE.eF))/(ACE.aF %*% t(ACE.aF)+ACE.cF %*% t(ACE.cF)+ACE.eF %*% t(ACE.eF))),model= univModACEFit)[2])
(se_stndAM <- mxSE(diag2vec((ACE.aM %*% t(ACE.aM))/(ACE.aM %*% t(ACE.aM)+ACE.cM %*% t(ACE.cM)+ACE.eM %*% t(ACE.eM))),model= univModACEFit)[2])
(se_stndCM <- mxSE(diag2vec((ACE.cM %*% t(ACE.cM))/(ACE.aM %*% t(ACE.aM)+ACE.cM %*% t(ACE.cM)+ACE.eM %*% t(ACE.eM))),model= univModACEFit)[2])
(se_stndEM <- mxSE(diag2vec((ACE.eM %*% t(ACE.eM))/(ACE.aM %*% t(ACE.aM)+ACE.cM %*% t(ACE.cM)+ACE.eM %*% t(ACE.eM))),model= univModACEFit)[2])

# Save for later
std <- c(stndAF, stndCF, stndEF, stndAM, stndCM, stndEM)
se  <- c(se_stndAF,  se_stndCF,  se_stndEF,  se_stndAM,  se_stndCM,  se_stndEM)

# Lower CI
(std - 1.95*se)
# Upper CI
(std + 1.95*se)


############################################
## Load CIs from RC to see if they worked ##
############################################
# Total ace estimates + alc/educ means
load(file="univModACE_5_2_2024.RData")
univModACESumm <- summary(univModACEFit)
univModACESumm
univModACESumm$CI

load(file="univModACE_5_4_2024.RData")
summary(univModACE_CI)$CI
load(file="univModACEsumm_5_4_2024.RData")
univModACESumm


## Standardized ACEs
load(file="univModACE_5_6_2024_stdACE.RData")
load(file="univModACEsumm_5_6_2024_stdACE.RData")
univModACESumm

## ACE moderation
load(file="univModACE_5_5_2024.RData")
load(file="univModACEsumm_5_5_2024.RData")
univModACESumm

# Betas
load(file="univModACE_5_6_2024.RData")
load(file="univModACEsumm_5_6_2024.RData")
univModACESumm
univModACESumm$CI

# Sys.setenv(OMP_NUM_THREADS=parallel::detectCores())


########################################
######        Plots       ##############
########################################

load(file="univModACE_5_6_2024_stdACE.RData")
load(file="univModACEsumm_5_6_2024_stdACE.RData")
univModACESumm

#Note you have to comment out specific parts of these scripts depending on which plot you want to output

#describe(data.frame(mergedata$zage_c_A,mergedata$age_c))
#hist(mergedata$zage_c_A)
source("Gustavson_IGEMS_Alcohol_univmodplots_Age.R") #put in working directory
univmodplots(univModACEFit, univModACESumm)

#describe(data.frame(mergedata$zeduc_A, mergedata$zeduc_B))
#table(mergedata$zeduc_A)
source("Gustavson_IGEMS_Alcohol_univmodplots_Educ.R") #put in working directory
univmodplots(univModACEFit, univModACESumm)

