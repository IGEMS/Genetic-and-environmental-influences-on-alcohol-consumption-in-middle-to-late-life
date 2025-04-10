#univmodplots <- function(FIT,SUMM,VALUES=seq(2,6,1)){
univmodplots <- function(FIT,SUMM,VALUES=seq(-1.1,2.1,.2)){

fit <- FIT
summ <- SUMM

fit <- univModACEFit
#summ <- univModACESumm

# aM <- fit$ACE@matrices$aM@values
# cM <- fit$ACE@matrices$cM@values
# eM <- fit$ACE@matrices$eM@values
# 
# aF <- fit$ACE@matrices$aF@values
# cF <- fit$ACE@matrices$cF@values
# eF <- fit$ACE@matrices$eF@values

aF21 <- fit$ACE@matrices$aF@values[2,1]
cF21 <- fit$ACE@matrices$cF@values[2,1]
eF21 <- fit$ACE@matrices$eF@values[2,1]
aM21 <- fit$ACE@matrices$aM@values[2,1]
cM21 <- fit$ACE@matrices$cM@values[2,1]
eM21 <- fit$ACE@matrices$eM@values[2,1]

aF22 <- fit$ACE@matrices$aF@values[2,2]
cF22 <- fit$ACE@matrices$cF@values[2,2]
eF22 <- fit$ACE@matrices$eF@values[2,2]
aM22 <- fit$ACE@matrices$aM@values[2,2]
cM22 <- fit$ACE@matrices$cM@values[2,2]
eM22 <- fit$ACE@matrices$eM@values[2,2]

#should this be sqrt of squares?
aF <- sqrt(aF21*aF21+aF22*aF22)
cF <- sqrt(cF21*cF21+cF22*cF22)
eF <- sqrt(eF21*eF21+eF22*eF22)
aM <- sqrt(aM21*aM21+aM22*aM22)
cM <- sqrt(cM21*cM21+cM22*cM22)
eM <- sqrt(eM21*eM21+eM22*eM22)


aIM21 <- fit$ACE@matrices$aIM1@values[2,1]
cIM21 <- fit$ACE@matrices$cIM1@values[2,1]
eIM21 <- fit$ACE@matrices$eIM1@values[2,1]
aIF21 <- fit$ACE@matrices$aIF1@values[2,1]
cIF21 <- fit$ACE@matrices$cIF1@values[2,1]
eIF21 <- fit$ACE@matrices$eIF1@values[2,1]

aIM22 <- fit$ACE@matrices$aIM1@values[2,2]
cIM22 <- fit$ACE@matrices$cIM1@values[2,2]
eIM22 <- fit$ACE@matrices$eIM1@values[2,2]
aIF22 <- fit$ACE@matrices$aIF1@values[2,2]
cIF22 <- fit$ACE@matrices$cIF1@values[2,2]
eIF22 <- fit$ACE@matrices$eIF1@values[2,2]

amodF_21 <- rep(1,length(VALUES)) * aF21 + VALUES * aIF21 
cmodF_21 <- rep(1,length(VALUES)) * cF21 + VALUES * cIF21 
emodF_21 <- rep(1,length(VALUES)) * eF21 + VALUES * eIF21 
amodM_21 <- rep(1,length(VALUES)) * aM21 + VALUES * aIM21 
cmodM_21 <- rep(1,length(VALUES)) * cM21 + VALUES * cIM21 
emodM_21 <- rep(1,length(VALUES)) * eM21 + VALUES * eIM21

amodF_22 <- rep(1,length(VALUES)) * aF22 + VALUES * aIF22 
cmodF_22 <- rep(1,length(VALUES)) * cF22 + VALUES * cIF22 
emodF_22 <- rep(1,length(VALUES)) * eF22 + VALUES * eIF22 
amodM_22 <- rep(1,length(VALUES)) * aM22 + VALUES * aIM22 
cmodM_22 <- rep(1,length(VALUES)) * cM22 + VALUES * cIM22 
emodM_22 <- rep(1,length(VALUES)) * eM22 + VALUES * eIM22

amodF <- sqrt(amodF_21*amodF_21+amodF_22*amodF_22)
cmodF <- sqrt(cmodF_21*cmodF_21+cmodF_22*cmodF_22)
emodF <- sqrt(emodF_21*emodF_21+emodF_22*emodF_22)
amodM <- sqrt(amodM_21*amodM_21+amodM_22*amodM_22)
cmodM <- sqrt(cmodM_21*cmodM_21+cmodM_22*cmodM_22)
emodM <- sqrt(emodM_21*emodM_21+emodM_22*emodM_22)

AmodF <- amodF * amodF
CmodF <- cmodF * cmodF
EmodF <- emodF * emodF
AmodM <- amodM * amodM
CmodM <- cmodM * cmodM
EmodM <- emodM * emodM

VM <- AmodM + CmodM + EmodM
VF <- AmodF + CmodF + EmodF


ApropM <- AmodM/VM 
CpropM <- CmodM/VM 
EpropM <- EmodM/VM 

ApropF <- AmodF/VF 
CpropF <- CmodF/VF 
EpropF <- EmodF/VF 

### Add standard errors
# aM_se <- summ$parameters$Std.Error[1]
# cM_se <- summ$parameters$Std.Error[2]
# eM_se <- summ$parameters$Std.Error[3]
# aF_se <- summ$parameters$Std.Error[4]
# cF_se <- summ$parameters$Std.Error[5]
# eF_se <- summ$parameters$Std.Error[6]
# 
# aIM_se <- summ$parameters$Std.Error[7]
# cIM_se <- summ$parameters$Std.Error[8]
# eIM_se <- summ$parameters$Std.Error[9]
# aIF_se <- summ$parameters$Std.Error[16]
# cIF_se <- summ$parameters$Std.Error[17]
# eIF_se <- summ$parameters$Std.Error[18]
# 
# 
# amodF_1 <- amodF-1.96*aF_se
# amodF_2 <- amodF+1.96*aF_se
# cmodF_1 <- cmodF-1.96*cF_se
# cmodF_2 <- cmodF+1.96*cF_se
# emodF_1 <- emodF-1.96*eF_se
# emodF_2 <- emodF+1.96*eF_se
# amodM_1 <- amodM-1.96*aM_se
# amodM_2 <- amodM+1.96*aM_se
# cmodM_1 <- cmodM-1.96*cM_se
# cmodM_2 <- cmodM+1.96*cM_se
# emodM_1 <- emodM-1.96*eM_se
# emodM_2 <- emodM+1.96*eM_se
# amodF_L <- pmin(amodF_1,amodF_2)
# amodF_H <- pmax(amodF_1,amodF_2)
# cmodF_L <- pmin(cmodF_1,cmodF_2)
# cmodF_H <- pmax(cmodF_1,cmodF_2)
# emodF_L <- pmin(emodF_1,emodF_2)
# emodF_H <- pmax(emodF_1,emodF_2)
# amodM_L <- pmin(amodM_1,amodM_2)
# amodM_H <- pmax(amodM_1,amodM_2)
# cmodM_L <- pmin(cmodM_1,cmodM_2)
# cmodM_H <- pmax(cmodM_1,cmodM_2)
# emodM_L <- pmin(emodM_1,emodM_2)
# emodM_H <- pmax(emodM_1,emodM_2)
# amodF_L[amodF_L<0]<-0
# cmodF_L[cmodF_L<0]<-0
# emodF_L[emodF_L<0]<-0
# amodM_L[amodM_L<0]<-0
# cmodM_L[cmodM_L<0]<-0
# emodM_L[emodM_L<0]<-0
# 
# AmodF_L <- amodF_L * amodF_L
# AmodF_H <- amodF_H * amodF_H
# CmodF_L <- cmodF_L * cmodF_L
# CmodF_H <- cmodF_H * cmodF_H
# EmodF_L <- emodF_L * emodF_L
# EmodF_H <- emodF_H * emodF_H


### CREATE ACTUAL PLOTS ###
### Plots are created one at a time. Uncomment the plot you want, then rerun
### (There is probably a way to print all 4 at once but I didn't figure it out)
windows()
plot(VALUES, AmodF, type = "l",ylim=c(0,1.5),ylab="Variance Components",xlab="Moderating Variable (Education)",
     main="A. ACE Moderation by Education - Females  (Total Variance)", col="red2",lwd=3)
lines(VALUES, CmodF, lty=2, lwd=3, col="green4")
lines(VALUES, EmodF, lty=3, lwd=3, col="blue1")
lines(VALUES, VF, lty=4, lwd=3)
legend("topleft",c("Genetic Variance (A)","Common Environment Variance (C)",
                   "Unique Environment Variance (E)","Total Variance (V)"),
       lty=1:4, col=c("red2","green4","blue1","black"), lwd=2)

# plot(VALUES, AmodM, type = "l",ylim=c(0,1.5),ylab="Variance Components",xlab="Moderating Variable (Education)",
#      main="B. ACE Moderation by Education - Males (Total Variance)",col="red2", lwd=3)
# lines(VALUES, CmodM, lty=2, lwd=3, col="green4")
# lines(VALUES, EmodM, lty=3, lwd=3, col="blue1")
# lines(VALUES, VM, lty=4, lwd=3)
# legend("topleft",c("Genetic Variance (A)","Common Environment Variance (C)",
#                    "Unique Environment Variance (E)","Total Variance (V)"),
#        lty=1:4, col=c("red2","green4","blue1","black"),lwd=2)

# plot(VALUES, ApropF, type = "l",ylim=c(0,1),ylab="Standardized Variance Components",xlab="Moderating Variable (Education)",
#      main="C. ACE Moderation by Education - Females (% Variance)",col="red2", lwd=3)
# lines(VALUES, CpropF, lty=2, lwd=3, col="green4")
# lines(VALUES, EpropF, lty=3, lwd=3, col="blue1")
# legend("topright",c("Genetic Variance (A)","Common Environment Variance (C)","Unique Environment Variance (E)"),
#        lty=1:3, col=c("red2","green4","blue1"),lwd=2)

# plot(VALUES, ApropM, type = "l",ylim=c(0,1),ylab="Standardized Variance Components",xlab="Moderating Variable (Education)",
#      main="D. ACE Moderation by Education - Males (% Variance)",col="red2", lwd=3)
# lines(VALUES, CpropM, lty=2, lwd=3, col="green4")
# lines(VALUES, EpropM, lty=3, lwd=3, col="blue1")
# legend("topright",c("Genetic Variance (A)","Common Environment Variance (C)","Unique Environment Variance (E)"),
#        lty=1:3, col=c("red2","green4","blue1"), lwd=2)


## Prints matrix of estimated values
print(round(cbind(VALUES,AmodF,CmodF,EmodF,VF,ApropF,CpropF,EpropF),3))
print(round(cbind(VALUES,AmodM,CmodM,EmodM,VM,ApropM,CpropM,EpropM),3))

est <- round(cbind(VALUES,AmodF,CmodF,EmodF,VF,ApropF,CpropF,EpropF,
                          AmodM,CmodM,EmodM,VM,ApropM,CpropM,EpropM),3)
write.csv(est, "EducMod_estimates_5_2024.csv", row.names = F)
}