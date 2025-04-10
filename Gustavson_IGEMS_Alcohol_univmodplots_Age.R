univmodplots <- function(FIT,SUMM,VALUES=seq(-1.75,2.25,.25)){


fit <- FIT
summ <- SUMM

#fit <- univModACEFit

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


aIF_21 <- fit$ACE@matrices$aIF2@values[2,1]
cIF_21 <- fit$ACE@matrices$cIF2@values[2,1]
eIF_21 <- fit$ACE@matrices$eIF2@values[2,1]
aIM_21 <- fit$ACE@matrices$aIM2@values[2,1]
cIM_21 <- fit$ACE@matrices$cIM2@values[2,1]
eIM_21 <- fit$ACE@matrices$eIM2@values[2,1]

aIF_22 <- fit$ACE@matrices$aIF2@values[2,2]
cIF_22 <- fit$ACE@matrices$cIF2@values[2,2]
eIF_22 <- fit$ACE@matrices$eIF2@values[2,2]
aIM_22 <- fit$ACE@matrices$aIM2@values[2,2]
cIM_22 <- fit$ACE@matrices$cIM2@values[2,2]
eIM_22 <- fit$ACE@matrices$eIM2@values[2,2]


amodF_21 <- rep(1,length(VALUES)) * aF21 + VALUES * aIF_21
cmodF_21 <- rep(1,length(VALUES)) * cF21 + VALUES * cIF_21 
emodF_21 <- rep(1,length(VALUES)) * eF21 + VALUES * eIF_21 
amodM_21 <- rep(1,length(VALUES)) * aM21 + VALUES * aIM_21 
cmodM_21 <- rep(1,length(VALUES)) * cM21 + VALUES * cIM_21 
emodM_21 <- rep(1,length(VALUES)) * eM21 + VALUES * eIM_21 

amodF_22 <- rep(1,length(VALUES)) * aF22 + VALUES * aIF_22 
cmodF_22 <- rep(1,length(VALUES)) * cF22 + VALUES * cIF_22 
emodF_22 <- rep(1,length(VALUES)) * eF22 + VALUES * eIF_22 
amodM_22 <- rep(1,length(VALUES)) * aM22 + VALUES * aIM_22 
cmodM_22 <- rep(1,length(VALUES)) * cM22 + VALUES * cIM_22 
emodM_22 <- rep(1,length(VALUES)) * eM22 + VALUES * eIM_22

##alternate
# amodF_Y <- rep(1,length(VALUES1)) * aF21 + VALUES1 * aIF_Y21 + rep(1,length(VALUES1)) * aF22 + VALUES1 * aIF_Y22
# cmodF_Y <- rep(1,length(VALUES1)) * cF21 + VALUES1 * cIF_Y21 + rep(1,length(VALUES1)) * cF22 + VALUES1 * cIF_Y22
# emodF_Y <- rep(1,length(VALUES1)) * eF21 + VALUES1 * eIF_Y21 + rep(1,length(VALUES1)) * eF22 + VALUES1 * eIF_Y22
# amodF_O <- rep(1,length(VALUES2)) * aF21 + VALUES2 * aIF_O21 + rep(1,length(VALUES2)) * aF22 + VALUES2 * aIF_O22
# cmodF_O <- rep(1,length(VALUES2)) * cF21 + VALUES2 * cIF_O21 + rep(1,length(VALUES2)) * cF22 + VALUES2 * cIF_O22
# emodF_O <- rep(1,length(VALUES2)) * eF21 + VALUES2 * eIF_O21 + rep(1,length(VALUES2)) * eF22 + VALUES2 * eIF_O22
# amodM_Y <- rep(1,length(VALUES1)) * aM21 + VALUES1 * aIM_Y21 + rep(1,length(VALUES1)) * aM22 + VALUES1 * aIM_Y22
# cmodM_Y <- rep(1,length(VALUES1)) * cM21 + VALUES1 * cIM_Y21 + rep(1,length(VALUES1)) * cM22 + VALUES1 * cIM_Y22
# emodM_Y <- rep(1,length(VALUES1)) * eM21 + VALUES1 * eIM_Y21 + rep(1,length(VALUES1)) * eM22 + VALUES1 * eIM_Y22
# amodM_O <- rep(1,length(VALUES2)) * aM21 + VALUES2 * aIM_O21 + rep(1,length(VALUES2)) * aM22 + VALUES2 * aIM_O22
# cmodM_O <- rep(1,length(VALUES2)) * cM21 + VALUES2 * cIM_O21 + rep(1,length(VALUES2)) * cM22 + VALUES2 * cIM_O22
# emodM_O <- rep(1,length(VALUES2)) * eM21 + VALUES2 * eIM_O21 + rep(1,length(VALUES2)) * eM22 + VALUES2 * eIM_O22 

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



windows()

#### REVISED PLOTS FOR SPLINE METHOD
# plot(VALUES, AmodF, type = "l",ylim=c(0,max(VM)+.5),ylab="Variance Components",xlim=c(-1.75,2.25),
#      xlab="Moderating Variable (Age)",col="red2",
#      main="A. ACE Moderation by Age - Females (Total Variance)", lwd=3)
# lines(VALUES, CmodF, lty=2, lwd=3, col="green4")
# lines(VALUES, EmodF, lty=3, lwd=3, col="blue1")
# lines(VALUES, VF, lty=4, lwd=3)
# legend("topleft",c("Genetic Variance (A)","Common Environment Variance (C)",
#                    "Unique Environment Variance (E)","Total Variance (V)"),
#        col=c("red2","green4","blue1","black"),lty=1:4, lwd=2)

 # plot(VALUES, AmodM, type = "l",ylim=c(0,max(VM)+.5),ylab="Variance Components",xlim=c(-1.75,2.25),
 #      xlab="Moderating Variable (Age)",col="red2",
 #      main="B. ACE Moderation by Age - Males (Total Variance)", lwd=3)
 # lines(VALUES, CmodM, lty=2, lwd=3, col="green4")
 # lines(VALUES, EmodM, lty=3, lwd=3, col="blue1")
 # lines(VALUES, VM, lty=4, lwd=3)
 # legend("topleft",c("Genetic Variance (A)","Common Environment Variance (C)",
 #                   "Unique Environment Variance (E)","Total Variance (V)"),
 #        col=c("red2","green4","blue1","black"),lty=1:4, lwd=2)

plot(VALUES, ApropF, type = "l",ylim=c(0,1),ylab="Standardized Variance Components",xlim=c(-1.75,2.25),
     xlab="Moderating Variable (Age)",col="red2",
     main="C. ACE Moderation by Age - Females (% Variance)", lwd=3)
lines(VALUES, CpropF, lty=2, lwd=3, col="green4")
lines(VALUES, EpropF, lty=3, lwd=3, col="blue1")
legend("topleft",c("Genetic Variance (A)","Common Environment Variance (C)", "Unique Environment Variance (E)"),
       col=c("red2","green4","blue1"),lty=1:4, lwd=2)

# plot(VALUES, ApropM, type = "l",ylim=c(0,1),ylab="Standardized Variance Components",xlim=c(-1.75,2.25),
#      xlab="Moderating Variable (Age)",col="red2",
#      main="D. ACE Moderation by Age - Males (% Variance)", lwd=3)
# lines(VALUES, CpropM, lty=2, lwd=3, col="green4")
# lines(VALUES, EpropM, lty=3, lwd=3, col="blue1")
# legend("topleft",c("Genetic Variance (A)","Common Environment Variance (C)",
#                   "Unique Environment Variance (E)"),
#        col=c("red2","green4","blue1"),lty=1:4, lwd=2)


print(round(cbind(VALUES,AmodF,CmodF,EmodF,VF,
                  ApropF,CpropF,EpropF,
                  AmodM,CmodM,EmodM,VF,
                  ApropM,CpropM,EpropM),3))

est <- round(cbind(VALUES,AmodF,CmodF,EmodF,VF,
                   ApropF,CpropF,EpropF,
                   AmodM,CmodM,EmodM,VF,
                   ApropM,CpropM,EpropM),3)
write.csv(est, "AgeMod_estimates_5_2024.csv", row.names = F)

}



