# -----------------------------------------------------------------------
# Program: Analysis of Alcohol Data for IGEMS
# Author: D Gustavson
# Date: 3 30 2023
# -----------------------------------------------------------------------

## Loading Required Packages
library(psych)  	
library(dplyr)
library(lme4)
library(lmerTest)
library(nlme)
library(tidyverse)
require("splines")
library(ggExtra)
library(viridis)
library(hrbrthemes)

## Reading in Data
load(file="Harmonized_File_v3.0.RData")
load(file="Harmonized_File_byStudy_v3.0.RData")
load(file="Harmonized_File_byCountry_v3.0.RData")



################################################################################################################
###############                   Step 1:   Additional Transformations                 #########################
################################################################################################################

###############################
##    Filter Age Under 40    ##
###############################

describe(dat$AGE_1stassessed)
dat <- filter(dat, AGE_1stassessed>=40)
describe(dat$AGE_1stassessed)


##################################
##    Filter have Drink Data    ##
##################################
dat <- filter(dat, ALC_weekly_grams_1stassessed>=0)
describe(dat$AGE_1stassessed)

table(dat$SEX)

######################################
### Transform and scale weekly ALC ###
######################################
dat$grams_sqrt <- sqrt(dat$ALC_weekly_grams_1stassessed)
dat$grams_ln <- log(dat$ALC_weekly_grams_1stassessed+1)
dat$zgrams_sqrt <- scale(dat$grams_sqrt)
dat$zgrams_ln <- scale(dat$grams_ln)


####################################################
### Windsorize SQRT and LN transformed measures  ###
####################################################
(sqrt_mean <- mean(dat$grams_sqrt, na.rm=T))
(sqrt_sd   <-   sd(dat$grams_sqrt, na.rm=T))
(sqrt_high <- sqrt_mean+4*sqrt_sd)
dat$grams_sqrt_windsor <- replace(dat$grams_sqrt, dat$grams_sqrt > sqrt_high, sqrt_high) 
dat$zgrams_sqrt_windsor <- scale(dat$grams_sqrt_windsor)
#count number trimmed - N=85
sum(dat$grams_sqrt>sqrt_high, na.rm=T)

(ln_mean <- mean(dat$grams_ln, na.rm=T))
(ln_sd   <-   sd(dat$grams_ln, na.rm=T))
(ln_high <- ln_mean+4*ln_sd)
dat$grams_ln_windsor <- replace(dat$grams_ln, dat$grams_ln > ln_high, ln_high) 
dat$zgrams_ln_windsor <- scale(dat$grams_ln_windsor)

(desc <- describe(select(dat, ALC_weekly_grams_1stassessed, grams_sqrt, grams_ln, grams_sqrt_windsor, grams_ln_windsor)))
(desc2 <- select(desc, c(2:5,8:9,11:12)))
write.csv(desc2, "output/transform_windsor_descriptives.csv")

#########################
###  Plot Histograms  ###
#########################
par( mfrow= c(2,2) )
hist(dat$grams_sqrt)
hist(dat$grams_sqrt_windsor)
hist(dat$grams_ln)
hist(dat$grams_ln_windsor)
par(mfrow= c(1,1))


############################
## Centered Age Variables ##
############################
mean(dat$AGE_1stassessed, na.rm=T) # 55.57

dat$agecenter <- dat$AGE_1stassessed - mean(dat$AGE_1stassessed, na.rm=T)
dat$agecenter40 <- dat$AGE_1stassessed - 40
dat$agecenter45 <- dat$AGE_1stassessed - 45
dat$agecenter50 <- dat$AGE_1stassessed - 50
dat$agecenter55 <- dat$AGE_1stassessed - 55
dat$agecenter60 <- dat$AGE_1stassessed - 60
dat$agecenter65 <- dat$AGE_1stassessed - 65
dat$agecenter70 <- dat$AGE_1stassessed - 70
dat$agecenter75 <- dat$AGE_1stassessed - 75
dat$agecenter80 <- dat$AGE_1stassessed - 80
dat$zage <- scale(dat$AGE_1stassessed)

##########################
##    ISCED - min max   ##
##########################
summarize(dat_group, 
          N = n(),
          ISCED_mean = mean(ISCED,na.rm=T),
          ISCED_sd = sd(ISCED,na.rm=T),
          ISCED_min = min(ISCED,na.rm=T),
          ISCED_max = max(ISCED,na.rm=T))

# create alternate ISCED with minimum of 1 and maximum of 7
# also collapse across bins 4 and 5
ISCEDtemp <- replace(dat$ISCED, dat$ISCED < 1, 1) 
ISCEDtemp2 <- replace(ISCEDtemp, ISCEDtemp == 5, 4) 
ISCEDtemp3 <- replace(ISCEDtemp2, ISCEDtemp2 == 6, 5)
dat$ISCED_trim <- replace(ISCEDtemp3, ISCEDtemp >= 7, 6)

# replot ISCED by sample
dat_group <- group_by(dat,COUNTRY, SAMPLE)
ggplot(dat_group, aes(x=ISCED_trim, color=COUNTRY, fill=SAMPLE)) +
  geom_histogram(alpha=0.6, binwidth = 1) + scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) + theme_ipsum() +
  theme(legend.position="none",panel.spacing = unit(0.1, "lines"),strip.text.x = element_text(size = 8)) +
  xlim(1,7) + ylab("ISCED Trimmed") + facet_wrap(~SAMPLE, scales="free_y")

(desc <- describe(select(dat, ISCED, ISCED_trim)))

# Dichotomized variable for ISCED (2 or less vs 3 or more)
ISCEDtemp2 <- replace(dat$ISCED, dat$ISCED <= 2, 0)
dat$ISCED_dichot <- replace(ISCEDtemp2, ISCEDtemp2 >= 3, 1)
table(dat$ISCED_dichot)

dat$zISCED <- scale(dat$ISCED)
dat$zISCED_trim <- scale(dat$ISCED_trim)
dat$ISCED_center <- dat$ISCED_trim - mean(dat$ISCED_trim, na.rm=T)

table(dat$ISCED_trim)

###############################
## Factors for Sex and ISCED ##
###############################
# SEX variable: 1 = male / 2 = female
#Recode twin number for opposite sex twins (Twin 1 = Male, Twin 2 = Female)
twintemp <- replace(dat$TWIN, (dat$ZYGOS==3 & dat$SEX == 1), 1)
dat$TWIN <- replace(twintemp, (dat$ZYGOS==3 & dat$SEX == 2), 2)

# Create factor and dichotomous variables for sex and education
dat$sexfactor <- factor(dat$SEX, labels=c("male","female"))
dat$educ_factor = factor(dat$ISCED_dichot, labels=c("Low (ISCED 2 or less)","High (ISCED 3 or more)"))
dat$sexF = dat$SEX - 1 ## 0 = male, 1 = female
dat$sexM = -1*dat$SEX + 2 ## 0 = female, 1 = male


#######################
## Factor for Cohort ##
#######################
dat$DOB <- rowSums(dat[c('DOB_yy', 'DOB_YY')], na.rm=T)
describe(dat$DOB)
dat$DOB[dat$DOB == 0] <- NA

hist(dat$DOB, breaks=30, xlab = "Year of Birth", main="Cohort")
abline(v=1931)

dat$cohort <- cut(dat$DOB, breaks = c(0, 1931, Inf), labels = c("0","1"))
table(dat$cohort)

#######################################################
# Remake grouping variables with new scaled variables #
#######################################################
dat_group <- group_by(dat,COUNTRY, SAMPLE)
dat_country <- group_by(dat,COUNTRY)

##############################
## Create Demographic Table ##
##############################

(desc_all <- summarize(dat,
                           N = length(ALC_weekly_grams_1stassessed[!is.na(ALC_weekly_grams_1stassessed)]),   
                           ALC_Mean = mean(ALC_weekly_grams_1stassessed, na.rm=T),
                           ALC_SD = sd(ALC_weekly_grams_1stassessed, na.rm=T),
                           ALC_min = min(ALC_weekly_grams_1stassessed, na.rm=T),
                           ALC_max = max(ALC_weekly_grams_1stassessed, na.rm=T),
                           ALC_abs = 1-mean(ALC_currdrinker_1stassessed, na.rm=T),
                           Sex = mean(sexF, na.rm=T),
                           AGE_Mean = mean(AGE_1stassessed, na.rm=T),
                           AGE_SD = sd(AGE_1stassessed, na.rm=T),
                           AGE_min = min(AGE_1stassessed, na.rm=T),
                           AGE_max = max(AGE_1stassessed, na.rm=T),
                           EDUC_N = length(ISCED_trim[!is.na(ISCED_trim)]),
                           EDUC_Mean = mean(ISCED_trim, na.rm=T),
                           EDUC_SD = sd(ISCED_trim, na.rm=T),
                           EDUC_min = min(ISCED_trim, na.rm=T),
                           EDUC_max = max(ISCED_trim, na.rm=T),
                           cohort = mean(as.integer(cohort), na.rm=T) - 1))

(desc_country <- summarize(dat_country,
                          N = length(ALC_weekly_grams_1stassessed[!is.na(ALC_weekly_grams_1stassessed)]),   
                          ALC_Mean = mean(ALC_weekly_grams_1stassessed, na.rm=T),
                          ALC_SD = sd(ALC_weekly_grams_1stassessed, na.rm=T),
                          ALC_min = min(ALC_weekly_grams_1stassessed, na.rm=T),
                          ALC_max = max(ALC_weekly_grams_1stassessed, na.rm=T),
                          ALC_abs = 1-mean(ALC_currdrinker_1stassessed, na.rm=T),
                          Sex = mean(sexF, na.rm=T),
                          AGE_Mean = mean(AGE_1stassessed, na.rm=T),
                          AGE_SD = sd(AGE_1stassessed, na.rm=T),
                          AGE_min = min(AGE_1stassessed, na.rm=T),
                          AGE_max = max(AGE_1stassessed, na.rm=T),
                          EDUC_N = length(ISCED_trim[!is.na(ISCED_trim)]),
                          EDUC_Mean = mean(ISCED_trim, na.rm=T),
                          EDUC_SD = sd(ISCED_trim, na.rm=T),
                          EDUC_min = min(ISCED_trim, na.rm=T),
                          EDUC_max = max(ISCED_trim, na.rm=T),
                          cohort = mean(as.integer(cohort), na.rm=T) - 1))

(desc_sample <- summarize(dat_group,
                        N = length(ALC_weekly_grams_1stassessed[!is.na(ALC_weekly_grams_1stassessed)]),   
                        ALC_Mean = mean(ALC_weekly_grams_1stassessed, na.rm=T),
                        ALC_SD = sd(ALC_weekly_grams_1stassessed, na.rm=T),
                        ALC_min = min(ALC_weekly_grams_1stassessed, na.rm=T),
                        ALC_max = max(ALC_weekly_grams_1stassessed, na.rm=T),
                        ALC_abs = 1-mean(ALC_currdrinker_1stassessed, na.rm=T),
                        Sex = mean(sexF, na.rm=T),
                        AGE_Mean = mean(AGE_1stassessed, na.rm=T),
                        AGE_SD = sd(AGE_1stassessed, na.rm=T),
                        AGE_min = min(AGE_1stassessed, na.rm=T),
                        AGE_max = max(AGE_1stassessed, na.rm=T),
                        EDUC_N = length(ISCED_trim[!is.na(ISCED_trim)]),
                        EDUC_Mean = mean(ISCED_trim, na.rm=T),
                        EDUC_SD = sd(ISCED_trim, na.rm=T),
                        EDUC_min = min(ISCED_trim, na.rm=T),
                        EDUC_max = max(ISCED_trim, na.rm=T),
                        cohort = mean(as.integer(cohort), na.rm=T) - 1))

(rownames <- c(desc_country[1], desc_sample[2]))
desc1 <- rbind(desc_country[2:18], desc_sample[3:19])
desc2 <- rbind(desc_all, desc1)
rows <- c("Full Sample", "Australia", "Denmark", "Sweden", "USA", 
          "OATS", "OVER50", "LSADT", "MADT", "MIDT", "GENDER", "HARMONY", "OCTOTWIN", "SALT", "SATSA",
          "MIDUS", "MTSADA", "NAS-NRC", "VETSA")
row.names(desc2) <- rows
write.csv(desc2, "output/Descriptive_Table_7_30_2024.csv")




(desc_cohort <- summarize(dat_group,
                          N = length(ALC_weekly_grams_1stassessed[!is.na(ALC_weekly_grams_1stassessed)]),   
                          cohort = mean(as.integer(cohort, na.rm=T)) - 1
                          ))


##############################################
## Additional Descriptive Stats for Samples ##
##############################################
## First recompute year of data collection
VETSA <- filter(dat, SAMPLE=="VETSA")  ## Missing IN_yy but we can compute it
head(VETSA)
VETSA$IN_yy <- VETSA$AGE_in + VETSA$DOB

MIDUS <- filter(dat, SAMPLE=="MIDUS")  ## MIDUS used fu1
describe(MIDUS$FU1_yy)
head(MIDUS)

HARMONY <- filter(dat, SAMPLE=="HARMONY")  ## Looks like its missing IN_yy
describe(HARMONY$FU1_yy)
head(HARMONY)

NASNRC <- filter(dat, SAMPLE=="NASNRC")  ## Looks like its missing IN_yy
head(NASNRC)

# Recompute variable
(dat <- mutate(dat, yy_1st_Assessed = case_when(
  (SAMPLE=="MIDUS" | SAMPLE=="HARMONY" | SAMPLE=="NASNRC") ~ FU1_yy,
  (SAMPLE=="VETSA") ~ AGE_in + DOB,
  (SAMPLE!="MIDUS" & SAMPLE!="HARMONY" & SAMPLE!="NASNRC" & SAMPLE!="VETSA") ~ IN_yy
)))   
describe(dat$yy_1st_Assessed)
dat_group <- group_by(dat,COUNTRY, SAMPLE)
dat_country <- group_by(dat,COUNTRY)



## Additional descriptive statistics for birth year and year of data collection 

(desc_all3 <- summarize(dat,
                       N = length(ALC_weekly_grams_1stassessed[!is.na(ALC_weekly_grams_1stassessed)]),   
                       Intake_N = length(yy_1st_Assessed[!is.na(yy_1st_Assessed)]),   
                       Intake_year = mean(yy_1st_Assessed, na.rm=T),
                       Intake_min = min(yy_1st_Assessed, na.rm=T),
                       Intake_max = max(yy_1st_Assessed, na.rm=T),
                       DOB_N = length(DOB[!is.na(DOB)]),
                       DOB_mean = mean(DOB, na.rm=T),
                       DOB_min = min(DOB, na.rm=T),
                       DOB_max = max(DOB, na.rm=T),
                       EDUC_N_bin1 = sum(ISCED_trim==1, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                       EDUC_N_bin2 = sum(ISCED_trim==2, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                       EDUC_N_bin3 = sum(ISCED_trim==3, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                       EDUC_N_bin4 = sum(ISCED_trim==4, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                       EDUC_N_bin5 = sum(ISCED_trim==5, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                       EDUC_N_bin6 = sum(ISCED_trim==6, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                       cohort = mean(as.integer(cohort), na.rm=T) - 1))

(desc_country3 <- summarize(dat_country,
                        N = length(ALC_weekly_grams_1stassessed[!is.na(ALC_weekly_grams_1stassessed)]),   
                        Intake_N = length(yy_1st_Assessed[!is.na(yy_1st_Assessed)]),   
                        Intake_year = mean(yy_1st_Assessed, na.rm=T),
                        Intake_min = min(yy_1st_Assessed, na.rm=T),
                        Intake_max = max(yy_1st_Assessed, na.rm=T),
                        DOB_N = length(DOB[!is.na(DOB)]),
                        DOB_mean = mean(DOB, na.rm=T),
                        DOB_min = min(DOB, na.rm=T),
                        DOB_max = max(DOB, na.rm=T),
                        EDUC_N_bin1 = sum(ISCED_trim==1, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                        EDUC_N_bin2 = sum(ISCED_trim==2, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                        EDUC_N_bin3 = sum(ISCED_trim==3, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                        EDUC_N_bin4 = sum(ISCED_trim==4, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                        EDUC_N_bin5 = sum(ISCED_trim==5, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                        EDUC_N_bin6 = sum(ISCED_trim==6, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                        cohort = mean(as.integer(cohort), na.rm=T) - 1))

(desc_sample3 <- summarize(dat_group,
                         N = length(ALC_weekly_grams_1stassessed[!is.na(ALC_weekly_grams_1stassessed)]),   
                         Intake_N = length(yy_1st_Assessed[!is.na(yy_1st_Assessed)]),   
                         Intake_year = mean(yy_1st_Assessed, na.rm=T),
                         Intake_min = min(yy_1st_Assessed, na.rm=T),
                         Intake_max = max(yy_1st_Assessed, na.rm=T),
                         DOB_N = length(DOB[!is.na(DOB)]),
                         DOB_mean = mean(DOB, na.rm=T),
                         DOB_min = min(DOB, na.rm=T),
                         DOB_max = max(DOB, na.rm=T),
                         EDUC_N_bin1 = sum(ISCED_trim==1, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                         EDUC_N_bin2 = sum(ISCED_trim==2, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                         EDUC_N_bin3 = sum(ISCED_trim==3, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                         EDUC_N_bin4 = sum(ISCED_trim==4, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                         EDUC_N_bin5 = sum(ISCED_trim==5, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                         EDUC_N_bin6 = sum(ISCED_trim==6, na.rm=T)/length(ISCED_trim[!is.na(ISCED_trim)]),
                         cohort = mean(as.integer(cohort), na.rm=T) - 1))

(rownames <- c(desc_country3[1], desc_sample3[2]))
desc4 <- rbind(desc_country3[2:length(desc_country3)], desc_sample3[3:length(desc_sample3)])
desc5 <- rbind(desc_all3, desc4)
rows <- c("Full Sample", "Australia", "Denmark", "Sweden", "USA", 
          "OATS", "OVER50", "LSADT", "MADT", "MIDT", "GENDER", "HARMONY", "OCTOTWIN", "SALT", "SATSA",
          "MIDUS", "MTSADA", "NAS-NRC", "VETSA")
row.names(desc5) <- rows
write.csv(desc5, "output/Extra_Descriptives_Table_8_6_2024.csv")


## Save data here so we can reload if necessary
save(dat, file="Harmonized_File_v3.5.RData")
save(dat_group, file="Harmonized_File_byStudy_v3.5.RData")
save(dat_country, file="Harmonized_File_byCountry_v3.5.RData")
load(file="Harmonized_File_v3.5.RData")
load(file="Harmonized_File_byStudy_v3.5.RData")
load(file="Harmonized_File_byCountry_v3.5.RData")


################################################################################################################
###############                       Step 2:   Regression Spline Analyses             #########################
################################################################################################################

#######################################################
## Revision: 3-way interactions included from start  ##
#######################################################

dat$agesquared <- dat$agecenter*dat$agecenter
Model_full <- lmer(grams_sqrt_windsor~agecenter*sexfactor* ISCED_center + agesquared*sexfactor*ISCED_center + 
                     (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), 
                   data=dat,na.action=na.exclude, REML=F)
summary(Model_full)
anova(Model_full)
fit <- summary(Model_full)$AIC
fitquad <- c(fit[1], fit[2], fit[3])
(estsquad <- summary(Model_full)$coef)

Model_lin <- lmer(grams_sqrt_windsor~agecenter* sexfactor* ISCED_center + agecenter* ISCED_center   
                  + sexfactor * ISCED_center + (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), 
                  data=dat,na.action=na.exclude, REML=F)
summary(Model_lin)
anova(Model_lin)
fit <- summary(Model_lin)$AIC
fitlin <- c(fit[1], fit[2], fit[3])
(estslin <- summary(Model_lin)$coef)

#mean (55.58)
Model_splineM <- lmer(grams_sqrt_windsor~ (pmin(0,agecenter))*sexfactor*ISCED_center +(pmax(0,agecenter))*sexfactor*ISCED_center
                      + (pmin(0,agecenter))*ISCED_center +(pmax(0,agecenter))*ISCED_center+ sexfactor * ISCED_center +
                        (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_splineM)
est <- summary(Model_splineM)$coef
fit <- summary(Model_splineM)$AIC
fitMean <- c(fit[1], fit[2], fit[3])
estMean <- summary(Model_splineM)$coef[1:12]

#45
Model_spline45 <- lmer(grams_sqrt_windsor~ (pmin(0,agecenter45))*sexfactor*ISCED_center +(pmax(0,agecenter45))*sexfactor*ISCED_center
                       + (pmin(0,agecenter45))*ISCED_center +(pmax(0,agecenter45))*ISCED_center+ sexfactor * ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline45)
summary(Model_spline45)$coef
fit <- summary(Model_spline45)$AIC
fit45 <- c(fit[1], fit[2], fit[3])
est45 <- summary(Model_spline45)$coef[1:12]

#50
Model_spline50 <- lmer(grams_sqrt_windsor~ (pmin(0,agecenter50))*sexfactor*ISCED_center +(pmax(0,agecenter50))*sexfactor*ISCED_center
                       + (pmin(0,agecenter50))*ISCED_center +(pmax(0,agecenter50))*ISCED_center+ sexfactor * ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F,
                       control = lmerControl(optimizer ="Nelder_Mead"))  # helps with warning about convergence!
summary(Model_spline50)
summary(Model_spline50)$coef
fit <- summary(Model_spline50)$AIC
fit50 <- c(fit[1], fit[2], fit[3])
est50 <- summary(Model_spline50)$coef[1:12]

#55
Model_spline55 <- lmer(grams_sqrt_windsor~ (pmin(0,agecenter55))*sexfactor*ISCED_center +(pmax(0,agecenter55))*sexfactor*ISCED_center
                       + (pmin(0,agecenter55))*ISCED_center +(pmax(0,agecenter55))*ISCED_center+ sexfactor * ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline55)
summary(Model_spline55)$coef
fit <- summary(Model_spline55)$AIC
fit55 <- c(fit[1], fit[2], fit[3])
est55 <- summary(Model_spline55)$coef[1:12]

# 60
Model_spline60 <- lmer(grams_sqrt_windsor~ (pmin(0,agecenter60))*sexfactor*ISCED_center +(pmax(0,agecenter60))*sexfactor*ISCED_center
                       + (pmin(0,agecenter60))*ISCED_center +(pmax(0,agecenter60))*ISCED_center+ sexfactor * ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline60)
summary(Model_spline60)$coef
fit <- summary(Model_spline60)$AIC
fit60 <- c(fit[1], fit[2], fit[3])
est60 <- summary(Model_spline60)$coef[1:12]

# 65
Model_spline65 <- lmer(grams_sqrt_windsor~ (pmin(0,agecenter65))*sexfactor*ISCED_center +(pmax(0,agecenter65))*sexfactor*ISCED_center
                       + (pmin(0,agecenter65))*ISCED_center +(pmax(0,agecenter65))*ISCED_center+ sexfactor * ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline65)
summary(Model_spline65)$coef
fit <- summary(Model_spline65)$AIC
fit65 <- c(fit[1], fit[2], fit[3])
est65 <- summary(Model_spline65)$coef[1:12]

#70
Model_spline70 <- lmer(grams_sqrt_windsor~ (pmin(0,agecenter70))*sexfactor*ISCED_center +(pmax(0,agecenter70))*sexfactor*ISCED_center
                       + (pmin(0,agecenter70))*ISCED_center +(pmax(0,agecenter70))*ISCED_center+ sexfactor * ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F,
                       control = lmerControl(optimizer ="Nelder_Mead"))
summary(Model_spline70)
summary(Model_spline70)$coef
fit <- summary(Model_spline70)$AIC
fit70 <- c(fit[1], fit[2], fit[3])
est70 <- summary(Model_spline70)$coef[1:12]

#75
Model_spline75 <- lmer(grams_sqrt_windsor~ (pmin(0,agecenter75))*sexfactor*ISCED_center +(pmax(0,agecenter75))*sexfactor*ISCED_center
                       + (pmin(0,agecenter75))*ISCED_center +(pmax(0,agecenter75))*ISCED_center+ sexfactor * ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline75)
summary(Model_spline75)$coef
fit <- summary(Model_spline75)$AIC
fit75 <- c(fit[1], fit[2], fit[3])
est75 <- summary(Model_spline75)$coef[1:12]

#80
Model_spline80 <- lmer(grams_sqrt_windsor~ (pmin(0,agecenter80))*sexfactor*ISCED_center +(pmax(0,agecenter80))*sexfactor*ISCED_center
                       + (pmin(0,agecenter80))*ISCED_center +(pmax(0,agecenter80))*ISCED_center+ sexfactor * ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline80)
summary(Model_spline80)$coef
fit <- summary(Model_spline80)$AIC
fit80 <- c(fit[1], fit[2], fit[3])
est80 <- summary(Model_spline80)$coef[1:12]

fitall_3way <- rbind(fit45, fit50, fit55, fitMean,fit60, fit65, fit70, fit75, fit80, fitquad, fitlin)
write.csv(fitall_3way, "output/Regression_Spine_FitCompare_3way.csv")

estall_3way <- rbind(est45, est50, est55, estMean, est60, est65, est70, est75, est80)
colnames(estall_3way) <- c("Intercept","Age_low", "SEX", "ISCED", "Age_high",  "Age_low * SEX", "Age_low * ISCED",
                           "SEX * ISCED", "Age_high * SEX", "Age_high * ISCED",  "3-way (Age_low)","3-way (Age_high)")
write.csv(estall_3way, "output/Regression_Spine_Estimates_3way.csv")

## Add 3-way interactions to final model
Model_spline75_stnd <- lmer(scale(zgrams_sqrt_windsor)~ scale(pmin(0,agecenter75))*sexfactor*zISCED_trim +scale(pmax(0,agecenter75))*sexfactor*zISCED_trim
                            + scale(pmin(0,agecenter75))*zISCED_trim +scale(pmax(0,agecenter75))*zISCED_trim+ sexfactor * zISCED_trim +
                              (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline75_stnd)
fit <- summary(Model_spline75_stnd)$AIC
(ests<- summary(Model_spline75_stnd)$coef)
fit75_stnd <- c(fit[1], fit[2], fit[3])
est75_stnd <- summary(Model_spline75_stnd)$coef[2:10]
write.csv(ests, "output/Stand_Ests_7_30_2024_3way.csv")


#### Run some computations to help transform main effects back into reasonable metrics
# Alcohol M=81.62, SD=113.3
# Age M=55.57, SD=10.89
# Educ M=3.09, SD=1.41

# > describe(dat$grams_sqrt_windsor)
# vars     n mean   sd median trimmed  mad min   max range skew kurtosis   se
# X1    1 72371    7 5.68    6.7     6.5 6.71   0 29.84 29.84  0.6    -0.05 0.02

# Age main effect of -.05 (still based on age)
(7+-.05*5.68)^2 / 12 #expectation is 3.8 drinks per week in Europe
(7^2) / 12 # expectationof 4.1 drinks per week at mean

# Age main effect of -.09
(7+-.09*5.68)^2 / 12 #expectation is 3.5 drinks per week in Europe
(7^2) / 12 # expectationof 4.1 drinks per week at mean

# Sex main effect of -.57
(7+-.57*5.68)^2 / 12 # expectatoin is 1.2 drinks per week in Europe for females
4.1-1.2 # 2.9 fewer drinks per week for females

# Main effect for education
(7+.08*5.68)^2 / 12 # expectatoin is 4.6 drinks per week in Europe for females

# Check residuals
library(car)
shapiro.test(resid(Model_spline75_stnd)) # Doesn't work because same is too large

# Plot the residuals
plot(Model_spline75_stnd)

par(mfrow=c(2,2)) # init 4 charts in 1 panel
plot(Model_spline75_stnd)
par(mfrow=c(1,1)) # reset

# produce residual vs. fitted plot 
plot(fitted(Model_spline75_stnd), resid(Model_spline75_stnd)) 
# add a horizontal line at 0  
abline(0,0) 

# QQ plot
qqnorm(resid(Model_spline75_stnd)) 
qqline(resid(Model_spline75_stnd))

# Density
plot(density(resid(Model_spline75_stnd), na.rm=T))


# Other tests
lmtest::bptest((Model_spline75_stnd))
car::ncvTest(Model_spline75_stnd)
ncvMLM(Model_spline75_stnd)

library(influence.ME)
infl <- influence(Model_spline75_stnd, obs = TRUE)
cooks.distance(infl)

plot(infl, which = "cook")

#####################################
## Revision: test with count data  ##
#####################################

# Test with the poisson option
Model_spline75_poisson <- glmer(grams_sqrt_windsor~ AGE_1stassessed*sexfactor*ISCED_trim +
                              (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, family="poisson")
summary(Model_spline75_stnd)
fit <- summary(Model_spline75_stnd)$AIC
(ests<- summary(Model_spline75_stnd)$coef)
fit75_stnd <- c(fit[1], fit[2], fit[3])
est75_stnd <- summary(Model_spline75_stnd)$coef[2:10]
#write.csv(ests, "output/Stand_Ests_7_30_2024_3way.csv")

Model_spline75_poisson <- lmer(grams_sqrt_windsor~ scale(pmin(0,agecenter75))*sexfactor*zISCED_trim +scale(pmax(0,agecenter75))*sexfactor*zISCED_trim
                            + scale(pmin(0,agecenter75))*zISCED_trim +scale(pmax(0,agecenter75))*zISCED_trim+ sexfactor * zISCED_trim +
                              (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, family="poisson")


#################################################
## Revision: test with no interactions at all  ##
#################################################

dat$agesquared <- dat$agecenter*dat$agecenter
Model_full <- lmer(grams_sqrt_windsor~agecenter + sexfactor + ISCED_center + agesquared + (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), 
                   data=dat,na.action=na.exclude, REML=F)
summary(Model_full)
anova(Model_full)
fit <- summary(Model_full)$AIC
fitquad <- c(fit[1], fit[2], fit[3])

Model_lin <- lmer(grams_sqrt_windsor~agecenter + sexfactor + ISCED_center +  (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), 
                  data=dat,na.action=na.exclude, REML=F)
summary(Model_lin)
anova(Model_lin)
fit <- summary(Model_lin)$AIC
fitlin <- c(fit[1], fit[2], fit[3])

#mean (55.58)
Model_splineM <- lmer(grams_sqrt_windsor~  (pmin(0,agecenter)) +(pmax(0,agecenter)) + sexfactor + ISCED_center +
                        (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_splineM)
est <- summary(Model_splineM)$coef
fit <- summary(Model_splineM)$AIC
fitMean <- c(fit[1], fit[2], fit[3])
estMean <- summary(Model_splineM)$coef[1:5]

#45
Model_spline45 <- lmer(grams_sqrt_windsor~  (pmin(0,agecenter45)) +(pmax(0,agecenter45)) + sexfactor + ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline45)
summary(Model_spline45)$coef
fit <- summary(Model_spline45)$AIC
fit45 <- c(fit[1], fit[2], fit[3])
est45 <- summary(Model_spline45)$coef[1:5]

#50
Model_spline50 <- lmer(grams_sqrt_windsor~  (pmin(0,agecenter50)) +(pmax(0,agecenter50)) + sexfactor + ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F,
                       control = lmerControl(optimizer ="Nelder_Mead"))  # helps with warning about convergence!
summary(Model_spline50)
summary(Model_spline50)$coef
fit <- summary(Model_spline50)$AIC
fit50 <- c(fit[1], fit[2], fit[3])
est50 <- summary(Model_spline50)$coef[1:5]

#55
Model_spline55 <- lmer(grams_sqrt_windsor~  (pmin(0,agecenter55)) +(pmax(0,agecenter55)) + sexfactor + ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline55)
summary(Model_spline55)$coef
fit <- summary(Model_spline55)$AIC
fit55 <- c(fit[1], fit[2], fit[3])
est55 <- summary(Model_spline55)$coef[1:5]

# 60
Model_spline60 <- lmer(grams_sqrt_windsor~  (pmin(0,agecenter60)) +(pmax(0,agecenter60)) + sexfactor + ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline60)
summary(Model_spline60)$coef
fit <- summary(Model_spline60)$AIC
fit60 <- c(fit[1], fit[2], fit[3])
est60 <- summary(Model_spline60)$coef[1:5]

# 65
Model_spline65 <- lmer(grams_sqrt_windsor~  (pmin(0,agecenter65)) +(pmax(0,agecenter65)) + sexfactor + ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline65)
summary(Model_spline65)$coef
fit <- summary(Model_spline65)$AIC
fit65 <- c(fit[1], fit[2], fit[3])
est65 <- summary(Model_spline65)$coef[1:5]

#70
Model_spline70 <- lmer(grams_sqrt_windsor~  (pmin(0,agecenter70)) +(pmax(0,agecenter70)) + sexfactor + ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F,
                       control = lmerControl(optimizer ="Nelder_Mead"))
summary(Model_spline70)
summary(Model_spline70)$coef
fit <- summary(Model_spline70)$AIC
fit70 <- c(fit[1], fit[2], fit[3])
est70 <- summary(Model_spline70)$coef[1:5]

#75
Model_spline75 <- lmer(grams_sqrt_windsor~  (pmin(0,agecenter75)) +(pmax(0,agecenter75)) + sexfactor + ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline75)
summary(Model_spline75)$coef
fit <- summary(Model_spline75)$AIC
fit75 <- c(fit[1], fit[2], fit[3])
est75 <- summary(Model_spline75)$coef[1:5]

#80
Model_spline80 <- lmer(grams_sqrt_windsor~ (pmin(0,agecenter80)) +(pmax(0,agecenter80)) + sexfactor + ISCED_center +
                         (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_spline80)
summary(Model_spline80)$coef
fit <- summary(Model_spline80)$AIC
fit80 <- c(fit[1], fit[2], fit[3])
est80 <- summary(Model_spline80)$coef[1:5]

fitall_noInt <- rbind(fit45, fit50, fit55, fitMean,fit60, fit65, fit70, fit75, fit80, fitquad, fitlin)
write.csv(fitall__noInt, "output/Regression_Spine_FitCompare_noInt.csv")



#######################################
## Add cohort effects to final model ##
#######################################

Model_cohort_stnd <- lmer(scale(zgrams_sqrt_windsor)~ scale(pmin(0,agecenter75))*sexfactor*as.factor(cohort) +
                            scale(pmin(0,agecenter75))*zISCED_trim*as.factor(cohort) +
                            sexfactor*zISCED_trim*as.factor(cohort) +
                            scale(pmax(0,agecenter75))*sexfactor*zISCED_trim + 
                              (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_cohort_stnd)
fit <- summary(Model_cohort_stnd)$AIC
(ests<- summary(Model_cohort_stnd)$coef)
fit75_stnd <- c(fit[1], fit[2], fit[3])
est75_stnd <- summary(Model_cohort_stnd)$coef[1:19]
write.csv(ests, "output/Stand_Ests_cohortsensitivity_3way.csv")


Model_cohort_stnd <- lmer(scale(zgrams_sqrt_windsor)~ scale(pmin(0,agecenter75))*sexfactor*zISCED_trim*as.factor(cohort) +
                            scale(pmax(0,agecenter75))*sexfactor*zISCED_trim + 
                            (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_cohort_stnd)
fit <- summary(Model_cohort_stnd)$AIC
(ests<- summary(Model_cohort_stnd)$coef)
fit75_stnd <- c(fit[1], fit[2], fit[3])
est75_stnd <- summary(Model_cohort_stnd)$coef[1:20]
write.csv(ests, "output/Stand_Ests_cohortsensitivity_4way.csv")



Model_cohort_stnd <- lmer(scale(zgrams_sqrt_windsor)~ 
                            cohort + 
                            (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_cohort_stnd)

# count people in cohort #1 over age 75
dat75 <- dat[dat$agecenter75>0,]
table(dat75$cohort)
213/72371 #0.29% of the sample


#####################################
## Sensitivity Removing Abstainers ##
#####################################
abs_dat <- filter(dat,grams_sqrt_windsor>0)
describe(abs_dat$grams_sqrt_windsor)  # 60649
describe(abs_dat$zISCED_trim)  # 53858
abs_dat$zgrams_sqrt_windsor <- scale(abs_dat$grams_sqrt_windsor)
abs_dat$zISCED_trim <- scale(abs_dat$ISCED_trim)

abs_stnd <- lmer(scale(zgrams_sqrt_windsor)~ scale(pmin(0,agecenter75))*sexfactor*zISCED_trim +scale(pmax(0,agecenter75))*sexfactor*zISCED_trim +
                                (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=abs_dat,na.action=na.exclude, REML=F)
summary(abs_stnd)
fit <- summary(abs_stnd)$AIC
(ests<- summary(abs_stnd)$coef)
fit75_stnd <- c(fit[1], fit[2], fit[3])
est75_stnd <- summary(abs_stnd)$coef[1:12]
write.csv(ests, "output/Stand_Ests_noabstainers_3way.csv")

confint(abs_stnd)


######################################
## Sensitivity Removing Heavy Users ##
######################################

describe(dat$grams_sqrt_windsor) #72371
heavy_dat <- filter(dat,grams_sqrt_windsor<(7+5.68*1.5))
describe(heavy_dat$grams_sqrt_windsor) #60353 (1 SD)  66579 (1.5 SD)
describe(heavy_dat$zISCED_trim) #59737 (1.5 SD)

(1-66579/72371) # 8.0

heavy_dat$zgrams_sqrt_windsor <- scale(heavy_dat$grams_sqrt_windsor)
heavy_dat$zISCED_trim <- scale(heavy_dat$ISCED_trim)

heavy_stnd <- lmer(scale(zgrams_sqrt_windsor)~ scale(pmin(0,agecenter75))*sexfactor*zISCED_trim +scale(pmax(0,agecenter75))*sexfactor*zISCED_trim +
                   (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=heavy_dat,na.action=na.exclude, REML=F)
summary(heavy_stnd)
fit <- summary(heavy_stnd)$AIC
(ests<- summary(heavy_stnd)$coef)
fit75_stnd <- c(fit[1], fit[2], fit[3])
est75_stnd <- summary(heavy_stnd)$coef[1:12]
write.csv(ests, "output/Stand_Ests_noheavyusers_3way.csv")

confint(heavy_stnd)


##############################################
## Consideration of Year of Data Collection ##
##############################################

## see if we can fill in some missing info
describe(dat$yy_1st_Assessed)
dat$z_year_assessed <- scale(dat$yy_1st_Assessed)
describe(dat$DOB)

corr.test(data.frame(dat$yy_1st_Assessed, dat$AGE_1stassessed, dat$DOB))
# DOB with age - .35 / DOB with assessment .70 / AGE with DOB -.43

Model_yearcollected_stnd <- lmer(scale(zgrams_sqrt_windsor)~ scale(pmin(0,agecenter75))*sexfactor*zISCED_trim*z_year_assessed +
                            scale(pmax(0,agecenter75))*sexfactor*zISCED_trim*z_year_assessed + 
                            (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_yearcollected_stnd)
fit <- summary(Model_yearcollected_stnd)$AIC
(ests<- summary(Model_yearcollected_stnd)$coef)
fit75_stnd <- c(fit[1], fit[2], fit[3])
est75_stnd <- summary(Model_yearcollected_stnd)$coef[1:24]
write.csv(ests, "output/Stand_Ests_yearcollected_sensitivity_4way.csv")


#######################
## Main Effects Only ##
#######################

Model_ME_only <- lmer(scale(zgrams_sqrt_windsor)~ scale(pmin(0,agecenter75)) + scale(pmax(0,agecenter75)) + 
                                   sexfactor + zISCED_trim + z_year_assessed + as.factor(cohort) + 
                                   (1|COUNTRY) + (1|SAMPLE) + (1|PAIRID), data=dat,na.action=na.exclude, REML=F)
summary(Model_ME_only)
fit <- summary(Model_ME_only)$AIC
(ests<- summary(Model_ME_only)$coef)
write.csv(ests, "output/Stand_Ests_maineffect_only_sensitivity.csv")


###############
#### PLOTS ####
###############


################################
### PLOT THE SPLINE EFFECT  ####
################################
p <- ggplot(dat, aes(x = agecenter75, y = grams_sqrt_windsor)) +
  geom_point(size=1) + xlim(-30, 30) + xlab("Age (centered at age 75 years)") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
  geom_smooth(method = "lm",formula = y ~ pmin(0,x)+(pmax(0,x)), se = T)
print(p)

p <- ggplot(dat, aes(x = agecenter75, y = grams_sqrt_windsor, color=SAMPLE)) +
  geom_point(size=1) + xlim(-30, 30) + xlab("Age (centered at age 75 years)") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
  geom_smooth(method = "lm",formula = y ~ pmin(0,x)+(pmax(0,x)), se = T)
print(p)

p <- ggplot(dat, aes(x = agecenter75, y = grams_sqrt_windsor, color=COUNTRY)) +
  geom_point(size=1) + xlim(-30, 30) + xlab("Age (centered at age 75 years)") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
  geom_smooth(method = "lm",formula = y ~ pmin(0,x)+(pmax(0,x)), se = T)
print(p)

# Redo age without spline
p <- ggplot(dat, aes(x = agecenter75, y = grams_sqrt_windsor)) +
  geom_point(size=1) + xlim(-30, 30) + xlab("Age (centered at age 75 years)") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
  geom_smooth(se = T)
print(p)


### Educ
# by sample
p <- ggplot(dat, aes(x = ISCED_trim, y = grams_sqrt_windsor, colour = SAMPLE)) +
  geom_point(size=1) + xlim(1, 6) + xlab("Education (ISCED bins)") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
  geom_smooth(method = "lm", se = FALSE) 
print(p)

# by country
p <- ggplot(dat, aes(x = ISCED_trim, y = grams_sqrt_windsor, colour = COUNTRY)) +
  geom_point(size=1)+ xlim(1, 6) + xlab("Education (ISCED bins)") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
  geom_smooth(method = "lm", se = FALSE) 
 print(p)
 


###################################
### plot some 2-way interactions ##
###################################
# age x sex
p <- ggplot(dat, aes(x = agecenter75, y = grams_sqrt_windsor, color=sexfactor)) +
   geom_point(size=1) + xlim(-30, 30) + xlab("Age (centered at age 75 years)") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
   geom_smooth(method = "lm",formula = y ~ pmin(0,x)+(pmax(0,x)), se = T)
 print(p)
 
 
# education by sex
p <- ggplot(dat, aes(x = ISCED_trim, y = grams_sqrt_windsor, colour = sexfactor)) +
  geom_point(size=1)+ xlim(1, 6) + xlab("Education") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
  geom_smooth(method = "lm", se = FALSE) 
print(p)

# age x education
tempdat <- dat[!is.na(dat$ISCED_trim),]
p <- ggplot(tempdat, aes(x = agecenter75, y = grams_sqrt_windsor, color=as.factor(educ_factor))) +
  geom_point(size=1) + xlim(-30, 30) + xlab("Age (centered at age 75 years)") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
  geom_smooth(method = "lm",formula = y ~ pmin(0,x)+(pmax(0,x)), se = T)
print(p)


# Version after edits from Erik

# age x sex
ggplot(dat, aes(x = agecenter75, y = grams_sqrt_windsor, color=sexfactor)) +
  geom_point(size=1, alpha=.5, aes(fill=sexfactor), color="transparent",shape=21) + xlim(-30, 30) + xlab("Age (centered at 75 years)") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
  geom_smooth(aes(color=sexfactor), method = "lm",formula = y ~ pmin(0,x)+(pmax(0,x)), se=T)+
  scale_color_manual(values = c("red4","blue4")) + ggtitle("A. Alcohol Consumption by Age and Sex") +
  scale_fill_manual(values = c("lightcoral","steelblue1")) +
  guides(color = guide_legend(title = "Sex"),fill="none")

# education x age
levels(dat$educ_factor)
levels(dat$educ_factor) <- c("Low","High")
tempdat <- dat[!is.na(dat$ISCED_trim),]
ggplot(tempdat, aes(x = agecenter75, y = grams_sqrt_windsor, color=factor(educ_factor))) +
  geom_point(size=1, alpha=.5, aes(fill=factor(educ_factor)), color="transparent",shape=21) + xlim(-30, 30) + xlab("Age (centered at 75 years)") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
  geom_smooth(aes(color=factor(educ_factor)),method = "lm",formula = y ~ pmin(0,x)+(pmax(0,x)), se = T)+
  scale_color_manual(values = c("red4","blue4")) + ggtitle("B. Alcohol Consumption by Age and Education") +
  scale_fill_manual(values = c("lightcoral","steelblue1"))+
  guides(color = guide_legend(title = "Education"),fill="none")


########## COHORT ##############
# age x cohort
p <- ggplot(dat, aes(x = agecenter75, y = grams_sqrt_windsor, color=cohort)) +
  geom_point(size=1) + xlim(-30, 30) + xlab("Age (centered at age 75 years)") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
  geom_smooth(method = "lm",formula = y ~ pmin(0,x)+(pmax(0,x)), se = T)
print(p)

# educ x cohort
p <- ggplot(dat, aes(x = ISCED_trim, y = grams_sqrt_windsor, color=cohort)) +
  geom_point(size=1) + xlim(1, 6) + xlab("Education") + ylab("Weekly Grams of Alcohol (sqrt-transformed)") +
  geom_smooth(method = "lm", se = T)
print(p)


################################################################################################################
###############                    Step 3:   Cross-Twin Cross-Trait Corrs             #########################
################################################################################################################


##############################################
#######      Prepare Family File     #########
##############################################
names(dat)
twins <- select(dat, c(1:4,AGE_1stassessed,ALC_currdrinker_1stassessed,ALC_weekly_grams_1stassessed, 
                       COUNTRY, SAMPLE, 125:length(dat)))
names(twins)
table(twins$ZYGOS)

## Creating MZ and DZ Data Sets
# Dividing up twins
twin1 <- twins[twins$TWIN=="1",]
twin2 <- twins[twins$TWIN=="2",]

# Remerging by case to create paired data set
newtwins <- merge(twin1, twin2, by=c("PAIRID","ZYGOS","SAMPLE","COUNTRY"),all=T,suffixes=c("_1","_2"))
(length(newtwins$PAIRID)) # 49286 families

# Making data sets of Just MZ & DZ
MZdata <- newtwins[newtwins$ZYGOS==1,]
DZdata <- newtwins[newtwins$ZYGOS>=2,]
MZmale <- newtwins[newtwins$ZYGOS==1 & newtwins$sexfactor_1=="male",]
MZfem  <- newtwins[newtwins$ZYGOS==1 & newtwins$sexfactor_1=="female",]
DZmale <- newtwins[newtwins$ZYGOS==2 & newtwins$sexfactor_1=="male",]
DZfem  <- newtwins[newtwins$ZYGOS==2 & newtwins$sexfactor_1=="female",]
DZOS   <- newtwins[newtwins$ZYGOS==3,]

# check everything is accurate
table(MZmale$sexfactor_1)
table(MZmale$sexfactor_2)
table(DZmale$sexfactor_1)
table(DZmale$sexfactor_2)
table(MZfem$sexfactor_1)
table(MZfem$sexfactor_2)
table(DZfem$sexfactor_1)
table(DZfem$sexfactor_2)
table(DZOS$sexfactor_1)
table(DZOS$sexfactor_2)

##########################
#### twin correlations ###
##########################
# Full sample
MZ_all <- cor.test(MZdata$grams_sqrt_windsor_1, MZdata$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZdata$grams_sqrt_windsor_1, DZdata$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZmale$grams_sqrt_windsor_1, MZmale$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZfem$grams_sqrt_windsor_1,  MZfem$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZmale$grams_sqrt_windsor_1, DZmale$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZfem$grams_sqrt_windsor_1,  DZfem$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZOS$grams_sqrt_windsor_1,   DZOS$grams_sqrt_windsor_2)

corrs_full <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_full <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)
colnames(corrs_full) <- c("MZ All", "DZ All", "MZ Male", "MZ Female", "DZ Male", "DZ Female", "DZ Opposite")
colnames(ps_full) <- c("MZ All", "DZ All", "MZ Male", "MZ Female", "DZ Male", "DZ Female", "DZ Opposite")

#####################
# Subset by Country #
#####################
MZ_Den <- MZdata[MZdata$COUNTRY=="DEN",]
DZ_Den <- DZdata[DZdata$COUNTRY=="DEN",]
MZ_Aus <- MZdata[MZdata$COUNTRY=="AUS",]
DZ_Aus <- DZdata[DZdata$COUNTRY=="AUS",]
MZ_Swe <- MZdata[MZdata$COUNTRY=="SWE",]
DZ_Swe <- DZdata[DZdata$COUNTRY=="SWE",]
MZ_Usa <- MZdata[MZdata$COUNTRY=="USA",]
DZ_Usa <- DZdata[DZdata$COUNTRY=="USA",]

MZ_Den_M <- MZmale[MZmale$COUNTRY=="DEN",]
MZ_Den_F <- MZfem[MZfem$COUNTRY=="DEN",]
DZ_Den_M <- DZmale[DZmale$COUNTRY=="DEN",]
DZ_Den_F <- DZfem[DZfem$COUNTRY=="DEN",]
DZ_Den_O <- DZOS[DZOS$COUNTRY=="DEN",]
MZ_Aus_M <- MZmale[MZmale$COUNTRY=="AUS",]
MZ_Aus_F <- MZfem[MZfem$COUNTRY=="AUS",]
DZ_Aus_M <- DZmale[DZmale$COUNTRY=="AUS",]
DZ_Aus_F <- DZfem[DZfem$COUNTRY=="AUS",]
DZ_Aus_O <- DZOS[DZOS$COUNTRY=="AUS",]
MZ_Swe_M <- MZmale[MZmale$COUNTRY=="SWE",]
MZ_Swe_F <- MZfem[MZfem$COUNTRY=="SWE",]
DZ_Swe_M <- DZmale[DZmale$COUNTRY=="SWE",]
DZ_Swe_F <- DZfem[DZfem$COUNTRY=="SWE",]
DZ_Swe_O <- DZOS[DZOS$COUNTRY=="SWE",]
MZ_Usa_M <- MZmale[MZmale$COUNTRY=="USA",]
MZ_Usa_F <- MZfem[MZfem$COUNTRY=="USA",]
DZ_Usa_M <- DZmale[DZmale$COUNTRY=="USA",]
DZ_Usa_F <- DZfem[DZfem$COUNTRY=="USA",]
DZ_Usa_O <- DZOS[DZOS$COUNTRY=="USA",]

### Danish only
MZ_all <- cor.test(MZ_Den$grams_sqrt_windsor_1, MZ_Den$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_Den$grams_sqrt_windsor_1, DZ_Den$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_Den_M$grams_sqrt_windsor_1, MZ_Den_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_Den_F$grams_sqrt_windsor_1,  MZ_Den_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_Den_M$grams_sqrt_windsor_1, DZ_Den_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_Den_F$grams_sqrt_windsor_1,  DZ_Den_F$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_Den_O$grams_sqrt_windsor_1,   DZ_Den_O$grams_sqrt_windsor_2)
corrs_den <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_den <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

### Aus only
MZ_all <- cor.test(MZ_Aus$grams_sqrt_windsor_1, MZ_Aus$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_Aus$grams_sqrt_windsor_1, DZ_Aus$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_Aus_M$grams_sqrt_windsor_1, MZ_Aus_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_Aus_F$grams_sqrt_windsor_1,  MZ_Aus_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_Aus_M$grams_sqrt_windsor_1, DZ_Aus_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_Aus_F$grams_sqrt_windsor_1,  DZ_Aus_F$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_Aus_O$grams_sqrt_windsor_1,   DZ_Aus_O$grams_sqrt_windsor_2)
corrs_aus <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_aus <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

### Sweden only
MZ_all <- cor.test(MZ_Swe$grams_sqrt_windsor_1, MZ_Swe$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_Swe$grams_sqrt_windsor_1, DZ_Swe$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_Swe_M$grams_sqrt_windsor_1, MZ_Swe_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_Swe_F$grams_sqrt_windsor_1,  MZ_Swe_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_Swe_M$grams_sqrt_windsor_1, DZ_Swe_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_Swe_F$grams_sqrt_windsor_1,  DZ_Swe_F$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_Swe_O$grams_sqrt_windsor_1,   DZ_Swe_O$grams_sqrt_windsor_2)
corrs_swe <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_swe <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

### Usa only
MZ_all <- cor.test(MZ_Usa$grams_sqrt_windsor_1, MZ_Usa$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_Usa$grams_sqrt_windsor_1, DZ_Usa$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_Usa_M$grams_sqrt_windsor_1, MZ_Usa_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_Usa_F$grams_sqrt_windsor_1,  MZ_Usa_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_Usa_M$grams_sqrt_windsor_1, DZ_Usa_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_Usa_F$grams_sqrt_windsor_1,  DZ_Usa_F$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_Usa_O$grams_sqrt_windsor_1,   DZ_Usa_O$grams_sqrt_windsor_2)
corrs_usa <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_usa <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

## Combine across country
# (CTC_country <- rbind(corrs_full, corrs_aus, corrs_den, corrs_swe, corrs_usa))
# rownames(CTC_country) <- c("All", "Australia", "Denmark", "Sweden", "USA")
# CTC_country


table(dat$SAMPLE)
table(MZdata$SAMPLE)
table(DZdata$SAMPLE)

####################
# Subset by Sample #
####################
DZ_GENDER <- DZdata[DZdata$SAMPLE=="GENDER",] # Gender = only DZO twins
DZ_GENDER_O <- DZOS[DZOS$SAMPLE=="GENDER",]

MZ_HARMONY <- MZdata[MZdata$SAMPLE=="HARMONY",]
DZ_HARMONY <- DZdata[DZdata$SAMPLE=="HARMONY",]
MZ_HARMONY_M <- MZmale[MZmale$SAMPLE=="HARMONY",]
MZ_HARMONY_F <- MZfem[MZfem$SAMPLE=="HARMONY",]
DZ_HARMONY_M <- DZmale[DZmale$SAMPLE=="HARMONY",]
DZ_HARMONY_F <- DZfem[DZfem$SAMPLE=="HARMONY",]
DZ_HARMONY_O <- DZOS[DZOS$SAMPLE=="HARMONY",]

MZ_LSADT <- MZdata[MZdata$SAMPLE=="LSADT",]
DZ_LSADT <- DZdata[DZdata$SAMPLE=="LSADT",]
MZ_LSADT_M <- MZmale[MZmale$SAMPLE=="LSADT",]
MZ_LSADT_F <- MZfem[MZfem$SAMPLE=="LSADT",]
DZ_LSADT_M <- DZmale[DZmale$SAMPLE=="LSADT",]
DZ_LSADT_F <- DZfem[DZfem$SAMPLE=="LSADT",]
DZ_LSADT_O <- DZOS[DZOS$SAMPLE=="LSADT",]

MZ_MADT <- MZdata[MZdata$SAMPLE=="MADT",]
DZ_MADT <- DZdata[DZdata$SAMPLE=="MADT",]
MZ_MADT_M <- MZmale[MZmale$SAMPLE=="MADT",]
MZ_MADT_F <- MZfem[MZfem$SAMPLE=="MADT",]
DZ_MADT_M <- DZmale[DZmale$SAMPLE=="MADT",]
DZ_MADT_F <- DZfem[DZfem$SAMPLE=="MADT",]
DZ_MADT_O <- DZOS[DZOS$SAMPLE=="MADT",]

MZ_MIDT <- MZdata[MZdata$SAMPLE=="MIDT",]
DZ_MIDT <- DZdata[DZdata$SAMPLE=="MIDT",]
MZ_MIDT_M <- MZmale[MZmale$SAMPLE=="MIDT",]
MZ_MIDT_F <- MZfem[MZfem$SAMPLE=="MIDT",]
DZ_MIDT_M <- DZmale[DZmale$SAMPLE=="MIDT",]
DZ_MIDT_F <- DZfem[DZfem$SAMPLE=="MIDT",]
DZ_MIDT_O <- DZOS[DZOS$SAMPLE=="MIDT",]

MZ_MIDUS <- MZdata[MZdata$SAMPLE=="MIDUS",]
DZ_MIDUS <- DZdata[DZdata$SAMPLE=="MIDUS",]
MZ_MIDUS_M <- MZmale[MZmale$SAMPLE=="MIDUS",]
MZ_MIDUS_F <- MZfem[MZfem$SAMPLE=="MIDUS",]
DZ_MIDUS_M <- DZmale[DZmale$SAMPLE=="MIDUS",]
DZ_MIDUS_F <- DZfem[DZfem$SAMPLE=="MIDUS",]
DZ_MIDUS_O <- DZOS[DZOS$SAMPLE=="MIDUS",]

MZ_MTSADA <- MZdata[MZdata$SAMPLE=="MTSADA",]
DZ_MTSADA <- DZdata[DZdata$SAMPLE=="MTSADA",]
MZ_MTSADA_M <- MZmale[MZmale$SAMPLE=="MTSADA",]
MZ_MTSADA_F <- MZfem[MZfem$SAMPLE=="MTSADA",]
DZ_MTSADA_M <- DZmale[DZmale$SAMPLE=="MTSADA",]
DZ_MTSADA_F <- DZfem[DZfem$SAMPLE=="MTSADA",]
#DZ_MTSADA_O <- DZOS[DZOS$SAMPLE=="MTSADA",]

MZ_NASNRC <- MZdata[MZdata$SAMPLE=="NASNRC",]  # No females in NASNRC
DZ_NASNRC <- DZdata[DZdata$SAMPLE=="NASNRC",]
MZ_NASNRC_M <- MZmale[MZmale$SAMPLE=="NASNRC",]
DZ_NASNRC_M <- DZmale[DZmale$SAMPLE=="NASNRC",]

MZ_OATS <- MZdata[MZdata$SAMPLE=="OATS",]
DZ_OATS <- DZdata[DZdata$SAMPLE=="OATS",]
MZ_OATS_M <- MZmale[MZmale$SAMPLE=="OATS",]
MZ_OATS_F <- MZfem[MZfem$SAMPLE=="OATS",]
DZ_OATS_M <- DZmale[DZmale$SAMPLE=="OATS",]
DZ_OATS_F <- DZfem[DZfem$SAMPLE=="OATS",]
DZ_OATS_O <- DZOS[DZOS$SAMPLE=="OATS",]

MZ_OCTOTWIN <- MZdata[MZdata$SAMPLE=="OCTOTWIN",]
DZ_OCTOTWIN <- DZdata[DZdata$SAMPLE=="OCTOTWIN",]
MZ_OCTOTWIN_M <- MZmale[MZmale$SAMPLE=="OCTOTWIN",]
MZ_OCTOTWIN_F <- MZfem[MZfem$SAMPLE=="OCTOTWIN",]
DZ_OCTOTWIN_M <- DZmale[DZmale$SAMPLE=="OCTOTWIN",]
DZ_OCTOTWIN_F <- DZfem[DZfem$SAMPLE=="OCTOTWIN",] #No opposite sex twins

MZ_OVER50 <- MZdata[MZdata$SAMPLE=="OVER50",]
DZ_OVER50 <- DZdata[DZdata$SAMPLE=="OVER50",]
MZ_OVER50_M <- MZmale[MZmale$SAMPLE=="OVER50",]
MZ_OVER50_F <- MZfem[MZfem$SAMPLE=="OVER50",]
DZ_OVER50_M <- DZmale[DZmale$SAMPLE=="OVER50",]
DZ_OVER50_F <- DZfem[DZfem$SAMPLE=="OVER50",]
DZ_OVER50_O <- DZOS[DZOS$SAMPLE=="OVER50",]

MZ_SALT <- MZdata[MZdata$SAMPLE=="SALT",]
DZ_SALT <- DZdata[DZdata$SAMPLE=="SALT",]
MZ_SALT_M <- MZmale[MZmale$SAMPLE=="SALT",]
MZ_SALT_F <- MZfem[MZfem$SAMPLE=="SALT",]
DZ_SALT_M <- DZmale[DZmale$SAMPLE=="SALT",]
DZ_SALT_F <- DZfem[DZfem$SAMPLE=="SALT",]
DZ_SALT_O <- DZOS[DZOS$SAMPLE=="SALT",]

MZ_SATSA <- MZdata[MZdata$SAMPLE=="SATSA",]
DZ_SATSA <- DZdata[DZdata$SAMPLE=="SATSA",]
MZ_SATSA_M <- MZmale[MZmale$SAMPLE=="SATSA",]
MZ_SATSA_F <- MZfem[MZfem$SAMPLE=="SATSA",]
DZ_SATSA_M <- DZmale[DZmale$SAMPLE=="SATSA",]
DZ_SATSA_F <- DZfem[DZfem$SAMPLE=="SATSA",] # No opposite sex twins

MZ_VETSA <- MZdata[MZdata$SAMPLE=="VETSA",] # No females in VETSA
DZ_VETSA <- DZdata[DZdata$SAMPLE=="VETSA",]
MZ_VETSA_M <- MZmale[MZmale$SAMPLE=="VETSA",]
DZ_VETSA_M <- DZmale[DZmale$SAMPLE=="VETSA",]

### Gender only (there are only opposite sex DZ pairs in gender)
DZ_all <- cor.test(DZ_GENDER$grams_sqrt_windsor_1, DZ_GENDER$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_GENDER_O$grams_sqrt_windsor_1,   DZ_GENDER_O$grams_sqrt_windsor_2)
corrs_GENDER <- cbind("NA", DZ_all$est, "NA", "NA", "NA", "NA", DZ_o$est)
ps_GENDER <- cbind("NA", DZ_all$p.value, "NA", "NA", "NA", "NA", DZ_o$p.value)

### HARMONY only
MZ_all <- cor.test(MZ_HARMONY$grams_sqrt_windsor_1, MZ_HARMONY$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_HARMONY$grams_sqrt_windsor_1, DZ_HARMONY$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_HARMONY_M$grams_sqrt_windsor_1, MZ_HARMONY_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_HARMONY_F$grams_sqrt_windsor_1,  MZ_HARMONY_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_HARMONY_M$grams_sqrt_windsor_1, DZ_HARMONY_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_HARMONY_F$grams_sqrt_windsor_1,  DZ_HARMONY_F$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_HARMONY_O$grams_sqrt_windsor_1,   DZ_HARMONY_O$grams_sqrt_windsor_2)
corrs_HARMONY <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_HARMONY <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

### LSADT only
MZ_all <- cor.test(MZ_LSADT$grams_sqrt_windsor_1, MZ_LSADT$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_LSADT$grams_sqrt_windsor_1, DZ_LSADT$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_LSADT_M$grams_sqrt_windsor_1, MZ_LSADT_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_LSADT_F$grams_sqrt_windsor_1,  MZ_LSADT_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_LSADT_M$grams_sqrt_windsor_1, DZ_LSADT_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_LSADT_F$grams_sqrt_windsor_1,  DZ_LSADT_F$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_LSADT_O$grams_sqrt_windsor_1,   DZ_LSADT_O$grams_sqrt_windsor_2)
corrs_LSADT <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_LSADT <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

### MADT only
MZ_all <- cor.test(MZ_MADT$grams_sqrt_windsor_1, MZ_MADT$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_MADT$grams_sqrt_windsor_1, DZ_MADT$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_MADT_M$grams_sqrt_windsor_1, MZ_MADT_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_MADT_F$grams_sqrt_windsor_1,  MZ_MADT_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_MADT_M$grams_sqrt_windsor_1, DZ_MADT_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_MADT_F$grams_sqrt_windsor_1,  DZ_MADT_F$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_MADT_O$grams_sqrt_windsor_1,   DZ_MADT_O$grams_sqrt_windsor_2)
corrs_MADT <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_MADT <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

### MIDT only
MZ_all <- cor.test(MZ_MIDT$grams_sqrt_windsor_1, MZ_MIDT$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_MIDT$grams_sqrt_windsor_1, DZ_MIDT$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_MIDT_M$grams_sqrt_windsor_1, MZ_MIDT_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_MIDT_F$grams_sqrt_windsor_1,  MZ_MIDT_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_MIDT_M$grams_sqrt_windsor_1, DZ_MIDT_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_MIDT_F$grams_sqrt_windsor_1,  DZ_MIDT_F$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_MIDT_O$grams_sqrt_windsor_1,   DZ_MIDT_O$grams_sqrt_windsor_2)
corrs_MIDT <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_MIDT <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

### MIDUS only
MZ_all <- cor.test(MZ_MIDUS$grams_sqrt_windsor_1, MZ_MIDUS$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_MIDUS$grams_sqrt_windsor_1, DZ_MIDUS$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_MIDUS_M$grams_sqrt_windsor_1, MZ_MIDUS_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_MIDUS_F$grams_sqrt_windsor_1,  MZ_MIDUS_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_MIDUS_M$grams_sqrt_windsor_1, DZ_MIDUS_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_MIDUS_F$grams_sqrt_windsor_1,  DZ_MIDUS_F$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_MIDUS_O$grams_sqrt_windsor_1,   DZ_MIDUS_O$grams_sqrt_windsor_2)
corrs_MIDUS <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_MIDUS <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

### MTSADA only (no opposite sex twins)
MZ_all <- cor.test(MZ_MTSADA$grams_sqrt_windsor_1, MZ_MTSADA$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_MTSADA$grams_sqrt_windsor_1, DZ_MTSADA$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_MTSADA_M$grams_sqrt_windsor_1, MZ_MTSADA_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_MTSADA_F$grams_sqrt_windsor_1,  MZ_MTSADA_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_MTSADA_M$grams_sqrt_windsor_1, DZ_MTSADA_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_MTSADA_F$grams_sqrt_windsor_1,  DZ_MTSADA_F$grams_sqrt_windsor_2)
#DZ_o <- cor.test(DZ_MTSADA_O$grams_sqrt_windsor_1,   DZ_MTSADA_O$grams_sqrt_windsor_2)
corrs_MTSADA <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, "NA")
ps_MTSADA <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, "NA")

### NASNRC only (male only sample)
MZ_all <- cor.test(MZ_NASNRC$grams_sqrt_windsor_1, MZ_NASNRC$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_NASNRC$grams_sqrt_windsor_1, DZ_NASNRC$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_NASNRC_M$grams_sqrt_windsor_1, MZ_NASNRC_M$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_NASNRC_M$grams_sqrt_windsor_1, DZ_NASNRC_M$grams_sqrt_windsor_2)
corrs_NASNRC <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, "NA", DZ_m$est, "NA", "NA")
ps_NASNRC <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, "NA", DZ_m$p.value, "NA", "NA")

### OATS only
MZ_all <- cor.test(MZ_OATS$grams_sqrt_windsor_1, MZ_OATS$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_OATS$grams_sqrt_windsor_1, DZ_OATS$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_OATS_M$grams_sqrt_windsor_1, MZ_OATS_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_OATS_F$grams_sqrt_windsor_1,  MZ_OATS_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_OATS_M$grams_sqrt_windsor_1, DZ_OATS_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_OATS_F$grams_sqrt_windsor_1,  DZ_OATS_F$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_OATS_O$grams_sqrt_windsor_1,   DZ_OATS_O$grams_sqrt_windsor_2)
corrs_OATS <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_OATS <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

### OCTOTWIN only (no opposite sex twins)
MZ_all <- cor.test(MZ_OCTOTWIN$grams_sqrt_windsor_1, MZ_OCTOTWIN$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_OCTOTWIN$grams_sqrt_windsor_1, DZ_OCTOTWIN$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_OCTOTWIN_M$grams_sqrt_windsor_1, MZ_OCTOTWIN_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_OCTOTWIN_F$grams_sqrt_windsor_1,  MZ_OCTOTWIN_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_OCTOTWIN_M$grams_sqrt_windsor_1, DZ_OCTOTWIN_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_OCTOTWIN_F$grams_sqrt_windsor_1,  DZ_OCTOTWIN_F$grams_sqrt_windsor_2)
corrs_OCTOTWIN <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, "NA")
ps_OCTOTWIN <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

### OVER50 only
MZ_all <- cor.test(MZ_OVER50$grams_sqrt_windsor_1, MZ_OVER50$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_OVER50$grams_sqrt_windsor_1, DZ_OVER50$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_OVER50_M$grams_sqrt_windsor_1, MZ_OVER50_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_OVER50_F$grams_sqrt_windsor_1,  MZ_OVER50_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_OVER50_M$grams_sqrt_windsor_1, DZ_OVER50_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_OVER50_F$grams_sqrt_windsor_1,  DZ_OVER50_F$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_OVER50_O$grams_sqrt_windsor_1,   DZ_OVER50_O$grams_sqrt_windsor_2)
corrs_OVER50 <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_OVER50 <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

### SALT only
MZ_all <- cor.test(MZ_SALT$grams_sqrt_windsor_1, MZ_SALT$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_SALT$grams_sqrt_windsor_1, DZ_SALT$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_SALT_M$grams_sqrt_windsor_1, MZ_SALT_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_SALT_F$grams_sqrt_windsor_1,  MZ_SALT_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_SALT_M$grams_sqrt_windsor_1, DZ_SALT_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_SALT_F$grams_sqrt_windsor_1,  DZ_SALT_F$grams_sqrt_windsor_2)
DZ_o <- cor.test(DZ_SALT_O$grams_sqrt_windsor_1,   DZ_SALT_O$grams_sqrt_windsor_2)
corrs_SALT <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, DZ_o$est)
ps_SALT <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, DZ_o$p.value)

### SATSA only (no opposite sex twins)
MZ_all <- cor.test(MZ_SATSA$grams_sqrt_windsor_1, MZ_SATSA$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_SATSA$grams_sqrt_windsor_1, DZ_SATSA$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_SATSA_M$grams_sqrt_windsor_1, MZ_SATSA_M$grams_sqrt_windsor_2)
MZ_f <- cor.test(MZ_SATSA_F$grams_sqrt_windsor_1,  MZ_SATSA_F$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_SATSA_M$grams_sqrt_windsor_1, DZ_SATSA_M$grams_sqrt_windsor_2)
DZ_f <- cor.test(DZ_SATSA_F$grams_sqrt_windsor_1,  DZ_SATSA_F$grams_sqrt_windsor_2)
corrs_SATSA <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, MZ_f$est, DZ_m$est, DZ_f$est, "NA")
ps_SATSA <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, MZ_f$p.value, DZ_m$p.value, DZ_f$p.value, "NA")

### VETSA only (there are only males in VETSA)
MZ_all <- cor.test(MZ_VETSA$grams_sqrt_windsor_1, MZ_VETSA$grams_sqrt_windsor_2)
DZ_all <- cor.test(DZ_VETSA$grams_sqrt_windsor_1, DZ_VETSA$grams_sqrt_windsor_2)

MZ_m <- cor.test(MZ_VETSA_M$grams_sqrt_windsor_1, MZ_VETSA_M$grams_sqrt_windsor_2)
DZ_m <- cor.test(DZ_VETSA_M$grams_sqrt_windsor_1, DZ_VETSA_M$grams_sqrt_windsor_2)
corrs_VETSA <- cbind(MZ_all$est, DZ_all$est, MZ_m$est, "NA", DZ_m$est, "NA", "NA")
ps_VETSA <- cbind(MZ_all$p.value, DZ_all$p.value, MZ_m$p.value, "NA", DZ_m$p.value, "NA", "NA")

####################################
## Combine across country/sample ###
####################################
CTC_output <- rbind(corrs_full, corrs_aus, corrs_den, corrs_swe, corrs_usa,
                     corrs_GENDER, corrs_HARMONY, corrs_LSADT, corrs_MADT, corrs_MIDT,
                     corrs_MIDUS, corrs_MTSADA, corrs_NASNRC, corrs_OATS, corrs_OCTOTWIN,
                     corrs_OVER50, corrs_SALT, corrs_SATSA, corrs_VETSA)
rownames(CTC_output) <- c("All", "Australia", "Denmark", "Sweden", "USA",
                          "GENDER","HARMONY","LSADT","MADT","MIDT",
                          "MIDUS","MTSADA","NASNRC","OATS","OCTOTWIN",
                          "OVER50","SALT","SATSA","VETSA")
CTC_output
write.csv(CTC_output, "output/CTCT_Correlations.csv")

CTC_ps <- rbind(ps_full, ps_aus, ps_den, ps_swe, ps_usa,
                ps_GENDER, ps_HARMONY, ps_LSADT,  ps_MADT, ps_MIDT,
                ps_MIDUS,  ps_MTSADA,  ps_NASNRC, ps_OATS, ps_OCTOTWIN,
                ps_OVER50, ps_SALT,    ps_SATSA,  ps_VETSA)
rownames(CTC_ps) <- c("All", "Australia", "Denmark", "Sweden", "USA",
                          "GENDER","HARMONY","LSADT","MADT","MIDT",
                          "MIDUS","MTSADA","NASNRC","OATS","OCTOTWIN",
                          "OVER50","SALT","SATSA","VETSA")
write.csv(CTC_ps, "output/CTCT_P_Values.csv")


### Save file out for twin analyses
save(newtwins, file="Twin_File_for_OpenMX.RData")


################################################################################################################
###############                      Step 4:   Phenotypic Correlations                 #########################
################################################################################################################


###########################################################
#####    PHENO CORRS WITH AGE, SEX, AND EDUCATION     #####
###########################################################
dat$zage <- dat$AGE_1stassessed
summDat <- select(dat, c(PAIRID, zage, zISCED_trim, sexfactor, zgrams_sqrt_windsor,
                         SAMPLE, COUNTRY, ALC_currdrinker_1stassessed))

pheno_DEN <- summDat[summDat$COUNTRY=="DEN",]
pheno_USA <- summDat[summDat$COUNTRY=="USA",]
pheno_Aus <- summDat[summDat$COUNTRY=="AUS",]
pheno_Swe <- summDat[summDat$COUNTRY=="SWE",]

pheno_GENDER <- summDat[summDat$SAMPLE=="GENDER",]
pheno_HARMONY <- summDat[summDat$SAMPLE=="HARMONY",]
pheno_LSADT <- summDat[summDat$SAMPLE=="LSADT",]
pheno_MADT <- summDat[summDat$SAMPLE=="MADT",]
pheno_MIDT <- summDat[summDat$SAMPLE=="MIDT",]
pheno_MIDUS <- summDat[summDat$SAMPLE=="MIDUS",]
pheno_MTSADA <- summDat[summDat$SAMPLE=="MTSADA",]
pheno_NASNRC <- summDat[summDat$SAMPLE=="NASNRC",]
pheno_OATS <- summDat[summDat$SAMPLE=="OATS",]
pheno_OCTOTWIN <- summDat[summDat$SAMPLE=="OCTOTWIN",]
pheno_OVER50 <- summDat[summDat$SAMPLE=="OVER50",]
pheno_SALT <- summDat[summDat$SAMPLE=="SALT",]
pheno_SATSA <- summDat[summDat$SAMPLE=="SATSA",]
pheno_VETSA <- summDat[summDat$SAMPLE=="VETSA",]

## Correlations between alcohol, age, sex, and ISCED by country and sample
# Age
(regall <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID) + (1|SAMPLE) + (1|COUNTRY), data=summDat)))
(regden <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID) + (1|SAMPLE), data=pheno_DEN)))
(regusa <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID) + (1|SAMPLE), data=pheno_USA)))
(regaus <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID) + (1|SAMPLE), data=pheno_Aus)))
(regswe <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID) + (1|SAMPLE), data=pheno_Swe)))
(regGENDER  <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_GENDER)))
(regHARMONY <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_HARMONY)))
(regLSADT   <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_LSADT)))
(regMADT    <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_MADT)))
(regMIDT    <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_MIDT)))
(regMIDUS   <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_MIDUS)))
(regMTSADA   <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_MTSADA)))
(regNASNRC  <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_NASNRC)))
(regOATS    <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_OATS)))
(regOCTOTWIN<- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_OCTOTWIN)))
(regOVER50<- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_OVER50)))
(regSALT   <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_SALT)))
(regSATSA   <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_SATSA)))
(regVETSA   <- summary(lmer(zgrams_sqrt_windsor~zage + (1|PAIRID), data=pheno_VETSA)))

age_rs <- c(regall$coef[2], regaus$coef[2], regden$coef[2], regswe$coef[2], regusa$coef[2], 
            regGENDER$coef[2], regHARMONY$coef[2], regLSADT$coef[2], regMADT$coef[2], regMIDT$coef[2], regMIDUS$coef[2], regMTSADA$coef[2], 
            regNASNRC$coef[2], regOATS$coef[2], regOCTOTWIN$coef[2], regOVER50$coef[2], regSALT$coef[2], regSATSA$coef[2],regVETSA$coef[2])

age_ps <- c(regall$coef[2,5], regaus$coef[2,5], regden$coef[2,5], regswe$coef[2,5], regusa$coef[2,5], 
            regGENDER$coef[2,5], regHARMONY$coef[2,5], regLSADT$coef[2,5], regMADT$coef[2,5], regMIDT$coef[2,5], regMIDUS$coef[2,5], regMTSADA$coef[2,5], 
            regNASNRC$coef[2,5], regOATS$coef[2,5], regOCTOTWIN$coef[2,5], regOVER50$coef[2,5], regSALT$coef[2,5], regSATSA$coef[2,5],regVETSA$coef[2,5])

# education
(regall <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID) + (1|SAMPLE) + (1|COUNTRY), data=summDat, control = lmerControl(optimizer ="Nelder_Mead"))))
(regden <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID) + (1|SAMPLE), data=pheno_DEN)))
(regusa <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID) + (1|SAMPLE), data=pheno_USA)))
(regaus <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID) + (1|SAMPLE), data=pheno_Aus)))
(regswe <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID) + (1|SAMPLE), data=pheno_Swe)))
(regGENDER  <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_GENDER)))
(regHARMONY <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_HARMONY)))
(regLSADT   <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_LSADT)))
(regMADT    <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_MADT)))
(regMIDT    <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_MIDT)))
(regMIDUS   <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_MIDUS)))
(regMTSADA   <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_MTSADA)))
(regNASNRC  <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_NASNRC)))
(regOATS    <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_OATS)))
(regOCTOTWIN<- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_OCTOTWIN)))
(regOVER50  <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_OVER50)))
(regSALT   <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_SALT)))
(regSATSA   <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_SATSA)))
(regVETSA   <- summary(lmer(zgrams_sqrt_windsor~zISCED_trim + (1|PAIRID), data=pheno_VETSA)))

educ_rs <- c(regall$coef[2], regaus$coef[2], regden$coef[2], regswe$coef[2], regusa$coef[2], 
            regGENDER$coef[2], regHARMONY$coef[2], regLSADT$coef[2], regMADT$coef[2], regMIDT$coef[2], regMIDUS$coef[2], regMTSADA$coef[2], 
            regNASNRC$coef[2], regOATS$coef[2], regOCTOTWIN$coef[2], regOVER50$coef[2], regSALT$coef[2],regSATSA$coef[2], regVETSA$coef[2])

educ_ps <- c(regall$coef[2,5], regaus$coef[2,5], regden$coef[2,5], regswe$coef[2,5], regusa$coef[2,5], 
            regGENDER$coef[2,5], regHARMONY$coef[2,5], regLSADT$coef[2,5], regMADT$coef[2,5], regMIDT$coef[2,5], regMIDUS$coef[2,5], regMTSADA$coef[2,5], 
            regNASNRC$coef[2,5], regOATS$coef[2,5], regOCTOTWIN$coef[2,5], regOVER50$coef[2,5], regSALT$coef[2,5], regSATSA$coef[2,5],regVETSA$coef[2,5])



# sex
(regall <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID) + (1|SAMPLE) + (1|COUNTRY), data=summDat, control = lmerControl(optimizer ="Nelder_Mead"))))
(regden <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID) + (1|SAMPLE), data=pheno_DEN)))
(regusa <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID) + (1|SAMPLE), data=pheno_USA)))
(regaus <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID) + (1|SAMPLE), data=pheno_Aus)))
(regswe <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID) + (1|SAMPLE), data=pheno_Swe)))
(regGENDER  <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_GENDER)))
(regHARMONY <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_HARMONY)))
(regLSADT   <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_LSADT)))
(regMADT    <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_MADT)))
(regMIDT    <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_MIDT)))
(regMIDUS   <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_MIDUS)))
(regMTSADA   <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_MTSADA)))
#(regNASNRC  <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_NASNRC)))
(regOATS    <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_OATS)))
(regOCTOTWIN<- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_OCTOTWIN)))
(regOVER50  <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_OVER50)))
(regSALT    <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_SALT)))
(regSATSA   <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_SATSA)))
#(regVETSA   <- summary(lmer(zgrams_sqrt_windsor~sexfactor + (1|PAIRID), data=pheno_VETSA)))

sex_rs <- c(regall$coef[2], regaus$coef[2], regden$coef[2], regswe$coef[2], regusa$coef[2],  
             regGENDER$coef[2], regHARMONY$coef[2], regLSADT$coef[2], regMADT$coef[2], regMIDT$coef[2], regMIDUS$coef[2], regMTSADA$coef[2], 
             "NA", regOATS$coef[2], regOCTOTWIN$coef[2], regOVER50$coef[2], regSALT$coef[2], regSATSA$coef[2], "NA")

sex_ps <- c(regall$coef[2,5], regaus$coef[2,5], regden$coef[2,5], regswe$coef[2,5], regusa$coef[2,5], 
             regGENDER$coef[2,5], regHARMONY$coef[2,5], regLSADT$coef[2,5], regMADT$coef[2,5], regMIDT$coef[2,5], regMIDUS$coef[2,5],  regMTSADA$coef[2,5], 
            "NA", regOATS$coef[2,5], regOCTOTWIN$coef[2,5], regOVER50$coef[2,5], regSALT$coef[2,5], regSATSA$coef[2,5], "NA")

Ns <- c(length(summDat$zgrams_sqrt_windsor), length(pheno_Aus$zgrams_sqrt_windsor), length(pheno_DEN$zgrams_sqrt_windsor), length(pheno_Swe$zgrams_sqrt_windsor), length(pheno_USA$zgrams_sqrt_windsor),
        length(pheno_GENDER$zgrams_sqrt_windsor), length(pheno_HARMONY$zgrams_sqrt_windsor), length(pheno_LSADT$zgrams_sqrt_windsor), length(pheno_MADT$zgrams_sqrt_windsor), length(pheno_MIDT$zgrams_sqrt_windsor),
        length(pheno_MIDUS$zgrams_sqrt_windsor), length(pheno_MTSADA$zgrams_sqrt_windsor), length(pheno_NASNRC$zgrams_sqrt_windsor), length(pheno_OATS$zgrams_sqrt_windsor), length(pheno_OCTOTWIN$zgrams_sqrt_windsor), 
        length(pheno_OVER50$zgrams_sqrt_windsor),  length(pheno_SALT$zgrams_sqrt_windsor), length(pheno_SATSA$zgrams_sqrt_windsor), length(pheno_VETSA$zgrams_sqrt_windsor))


### Combine into a table
pheno_corrs_out <- cbind(Ns, age_rs, educ_rs, sex_rs, age_ps, educ_ps, sex_ps)
rownames(pheno_corrs_out) <- c("All","AUS","DEN","SWE","USA",
          "GENDER", "HARMONY", "LSADT", "MADT", "MIDT", "MIDUS", "MTSADA",
          "NASNRC", "OATS", "OCTOTWIN", "OVER50", "SALT", "SATSA", "VETSA")
write.csv(pheno_corrs_out, "output/Pheno_corrs_bySample.csv")


###############################################################
###   Some additional code that may not be needed anymore   ###
###############################################################
# 40-100 for age, 0 to 6 for education, 0-30 for sqrt-transformed alcohol
library(RColorBrewer)

# Make subsetted data
pheno_DEN <- filter(dat, COUNTRY=="DEN")
pheno_USA <- filter(dat, COUNTRY=="USA")
pheno_Swe <- filter(dat, COUNTRY=="SWE")
pheno_Aus <- filter(dat, COUNTRY=="AUS")

pheno_OATS <- filter(dat, SAMPLE=="OATS")
pheno_OVER50 <- filter(dat, SAMPLE=="OVER50")
pheno_LSADT <- filter(dat, SAMPLE=="LSADT")
pheno_MADT <- filter(dat, SAMPLE=="MADT")
pheno_MIDT <- filter(dat, SAMPLE=="MIDT")
pheno_GENDER <- filter(dat, SAMPLE=="GENDER")
pheno_HARMONY <- filter(dat, SAMPLE=="HARMONY")
pheno_OCTOTWIN <- filter(dat, SAMPLE=="OCTOTWIN")
pheno_SALT <- filter(dat, SAMPLE=="SALT")
pheno_SATSA <- filter(dat, SAMPLE=="SATSA")
pheno_MIDUS <- filter(dat, SAMPLE=="MIDUS")
pheno_MTSADA <- filter(dat, SAMPLE=="MTSADA")
pheno_NASNRC <- filter(dat, SAMPLE=="NASNRC")
pheno_VETSA <- filter(dat, SAMPLE=="VETSA")


## Age
par(mfrow=c(2,2))
hist(pheno_Aus$AGE_1stassessed,main="Australia", breaks=14, xlab="Age (Years)",xlim=c(40,100), col="mediumorchid4")
hist(pheno_DEN$AGE_1stassessed,main="Denmark", breaks=14, xlab="Age (Years)",xlim=c(40,100), col="skyblue4")
hist(pheno_Swe$AGE_1stassessed,main="Sweden", breaks=14, xlab="Age (Years)",xlim=c(40,100), col="palegreen3")
hist(pheno_USA$AGE_1stassessed,main="USA", breaks=14, xlab="Age (Years)",xlim=c(40,100), col="gold3")

par(mfrow=c(4,4))
hist(pheno_OATS$AGE_1stassessed,    main="OATS",   breaks=8, xlab=" ",xlim=c(40,100), col="mediumorchid4")
hist(pheno_OVER50$AGE_1stassessed,  main="OVER50", breaks=8, xlab=" ",xlim=c(40,100), col="mediumorchid4")
hist(pheno_LSADT$AGE_1stassessed,   main="LSADT",  breaks=8, xlab=" ",xlim=c(40,100), col="skyblue4")
hist(pheno_MADT$AGE_1stassessed,    main="MADT",   breaks=5, xlab=" ",xlim=c(40,100), col="skyblue4")
hist(pheno_MIDT$AGE_1stassessed,    main="MIDT",   breaks=8, xlab=" ",xlim=c(40,100), col="skyblue4")
hist(pheno_GENDER$AGE_1stassessed,  main="GENDER", breaks=5, xlab=" ",xlim=c(40,100), col="palegreen3")
hist(pheno_HARMONY$AGE_1stassessed, main="HARMONY",breaks=8, xlab=" ",xlim=c(40,100), col="palegreen3")
hist(pheno_OCTOTWIN$AGE_1stassessed,main="OCTOTWIN",breaks=5,xlab=" ",xlim=c(40,100), col="palegreen3")
hist(pheno_SALT$AGE_1stassessed,    main="SALT",   breaks=8, xlab=" ",xlim=c(40,100), col="palegreen3")
hist(pheno_SATSA$AGE_1stassessed,   main="SATSA",  breaks=8, xlab=" ",xlim=c(40,100), col="palegreen3")
hist(pheno_MIDUS$AGE_1stassessed,   main="MIDUS",  breaks=8, xlab=" ",xlim=c(40,100), col="gold3")
hist(pheno_MTSADA$AGE_1stassessed,  main="MTSADA", breaks=8, xlab=" ",xlim=c(40,100), col="gold3")
hist(pheno_NASNRC$AGE_1stassessed,  main="NASNRC", breaks=5, xlab=" ",xlim=c(40,100), col="gold3")
hist(pheno_VETSA$AGE_1stassessed,   main="VETSA",  breaks=3, xlab=" ",xlim=c(40,100), col="gold3")


## Education
par(mfrow=c(2,2))
hist(pheno_Aus$ISCED_trim,main="Australia", breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="mediumorchid4")
hist(pheno_DEN$ISCED_trim,main="Denmark", breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="skyblue4")
hist(pheno_Swe$ISCED_trim,main="Sweden", breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="palegreen3")
hist(pheno_USA$ISCED_trim,main="USA", breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="gold3")

par(mfrow=c(4,4))
hist(pheno_OATS$ISCED_trim,    main="OATS",   breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="mediumorchid4")
hist(pheno_OVER50$ISCED_trim,  main="OVER50", breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="mediumorchid4")
hist(pheno_LSADT$ISCED_trim,   main="LSADT",  breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="skyblue4")
hist(pheno_MADT$ISCED_trim,    main="MADT",   breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="skyblue4")
hist(pheno_MIDT$ISCED_trim,    main="MIDT",   breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="skyblue4")
hist(pheno_GENDER$ISCED_trim,  main="GENDER", breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="palegreen3")
hist(pheno_HARMONY$ISCED_trim, main="HARMONY",breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="palegreen3")
hist(pheno_OCTOTWIN$ISCED_trim,main="OCTOTWIN",breaks=seq(from=0.5, to=6.5, by=1),xlab=" ",xlim=c(0,7), col="palegreen3")
hist(pheno_SALT$ISCED_trim,    main="SALT",   breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="palegreen3")
hist(pheno_SATSA$ISCED_trim,   main="SATSA",  breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="palegreen3")
hist(pheno_MIDUS$ISCED_trim,   main="MIDUS",  breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="gold3")
hist(pheno_MTSADA$ISCED_trim,  main="MTSADA", breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="gold3")
hist(pheno_NASNRC$ISCED_trim,  main="NASNRC", breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="gold3")
hist(pheno_VETSA$ISCED_trim,   main="VETSA",  breaks=seq(from=0.5, to=6.5, by=1), xlab=" ",xlim=c(0,7), col="gold3")


## Alcohol
par(mfrow=c(2,2))
hist(pheno_Aus$grams_sqrt_windsor,main="Australia", breaks=6, xlab=" ",xlim=c(0,30), col="mediumorchid4")
hist(pheno_DEN$grams_sqrt_windsor,main="Denmark", breaks=6, xlab=" ",xlim=c(0,30), col="skyblue4")
hist(pheno_Swe$grams_sqrt_windsor,main="Sweden", breaks=6, xlab=" ",xlim=c(0,30), col="palegreen3")
hist(pheno_USA$grams_sqrt_windsor,main="USA", breaks=6, xlab=" ",xlim=c(0,30), col="gold3")

par(mfrow=c(4,4))
hist(pheno_OATS$grams_sqrt_windsor,    main="OATS",   breaks=6, xlab=" ",xlim=c(0,30), col="mediumorchid4")
hist(pheno_OVER50$grams_sqrt_windsor,  main="OVER50", breaks=6, xlab=" ",xlim=c(0,30), col="mediumorchid4")
hist(pheno_LSADT$grams_sqrt_windsor,   main="LSADT",  breaks=6, xlab=" ",xlim=c(0,30), col="skyblue4")
hist(pheno_MADT$grams_sqrt_windsor,    main="MADT",   breaks=6, xlab=" ",xlim=c(0,30), col="skyblue4")
hist(pheno_MIDT$grams_sqrt_windsor,    main="MIDT",   breaks=6, xlab=" ",xlim=c(0,30), col="skyblue4")
hist(pheno_GENDER$grams_sqrt_windsor,  main="GENDER", breaks=6, xlab=" ",xlim=c(0,30), col="palegreen3")
hist(pheno_HARMONY$grams_sqrt_windsor, main="HARMONY",breaks=6, xlab=" ",xlim=c(0,30), col="palegreen3")
hist(pheno_OCTOTWIN$grams_sqrt_windsor,main="OCTOTWIN",breaks=6,xlab=" ",xlim=c(0,30), col="palegreen3")
hist(pheno_SALT$grams_sqrt_windsor,    main="SALT",   breaks=6, xlab=" ",xlim=c(0,30), col="palegreen3")
hist(pheno_SATSA$grams_sqrt_windsor,   main="SATSA",  breaks=6, xlab=" ",xlim=c(0,30), col="palegreen3")
hist(pheno_MIDUS$grams_sqrt_windsor,   main="MIDUS",  breaks=6, xlab=" ",xlim=c(0,30), col="gold3")
hist(pheno_MTSADA$grams_sqrt_windsor,  main="MTSADA", breaks=6, xlab=" ",xlim=c(0,30), col="gold3")
hist(pheno_NASNRC$grams_sqrt_windsor,  main="NASNRC", breaks=6, xlab=" ",xlim=c(0,30), col="gold3")
hist(pheno_VETSA$grams_sqrt_windsor,   main="VETSA",  breaks=6, xlab=" ",xlim=c(0,30), col="gold3")

par(mfrow=c(1,1))

## Extra Descriptives by group


(desc_all <- summarize(dat,
                       ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                       ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

(desc_country <- summarize(dat_country,
                           ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                           ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

(desc_sample <- summarize(dat_group,
                          ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                          ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

# By Sex
datF <- filter(dat, SEX==2)
dat_countryF <- filter(dat_country, SEX==2)
(desc_all <- summarize(datF,
                       ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                       ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

(desc_country <- summarize(dat_countryF,
                           ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                           ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

datM <- filter(dat, SEX==1)
dat_countryM <- filter(dat_country, SEX==1)
(desc_all <- summarize(datM,
                       ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                       ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

(desc_country <- summarize(dat_countryM,
                           ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                           ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

# By Age
datY <- filter(dat, AGE_1stassessed<75)
dat_countryY <- filter(dat_country, AGE_1stassessed<75)
(desc_all <- summarize(datY,
                       ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                       ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

(desc_country <- summarize(dat_countryY,
                           ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                           ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

datO <- filter(dat, AGE_1stassessed>=75)
dat_countryO <- filter(dat_country, AGE_1stassessed>=75)
(desc_all <- summarize(datO,
                       ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                       ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

(desc_country <- summarize(dat_countryO,
                           ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                           ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

# By Educ
datE1 <- filter(dat, ISCED_trim==1)
dat_countryE1 <- filter(dat_country, ISCED_trim==1)
(desc_all <- summarize(datE1,
                       ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                       ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

(desc_country <- summarize(dat_countryE1,
                           ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                           ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

datE2 <- filter(dat, ISCED_trim==2)
dat_countryE2 <- filter(dat_country, ISCED_trim==2)
(desc_all <- summarize(datE2,
                       ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                       ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

(desc_country <- summarize(dat_countryE2,
                           ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                           ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))


datE3 <- filter(dat, ISCED_trim==3)
dat_countryE3 <- filter(dat_country, ISCED_trim==3)
(desc_all <- summarize(datE3,
                       ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                       ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

(desc_country <- summarize(dat_countryE3,
                           ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                           ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

datE4 <- filter(dat, ISCED_trim==4)
dat_countryE4 <- filter(dat_country, ISCED_trim==4)
(desc_all <- summarize(datE4,
                       ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                       ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

(desc_country <- summarize(dat_countryE4,
                           ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                           ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))


datE5 <- filter(dat, ISCED_trim==5)
dat_countryE5 <- filter(dat_country, ISCED_trim==5)
(desc_all <- summarize(datE5,
                       ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                       ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

(desc_country <- summarize(dat_countryE5,
                           ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                           ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

datE6 <- filter(dat, ISCED_trim==6)
dat_countryE6 <- filter(dat_country, ISCED_trim==6)
(desc_all <- summarize(datE6,
                       ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                       ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

(desc_country <- summarize(dat_countryE6,
                           ALC_Mean = mean(grams_sqrt_windsor, na.rm=T),
                           ALC_SD = sd(grams_sqrt_windsor, na.rm=T)))

