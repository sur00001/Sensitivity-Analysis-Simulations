#################################################################################################
# This script estimates the data generating mechanism models that will  be used to simulate the complete data
#################################################################################################

library(tidyverse) 
library(dplyr)
library(pscl)

setwd("~/Sensitivity Analysis - Paper 1/Delta-adjusted-sensitivity-analysis")
#Open complete data 
load("C:/Users/lylae/OneDrive/Documents/Sensitivity Analysis - Paper 1/completedat4pred_6-11.rda")

#Group variables into vectors 
norms = c("norm_percent_drink","norm_num_drinks_typical","norm_most_drinks_typical","norm_percent_binge")
#intents = c("intent_freq_per_month","intent_drunk_per_month","intent_typical_drinks")
#motives = c("soc.mot","conf.mot","enh.mot","cope.mot")
dems = c("female","nonwhite","greekintent")
base.outs = c("binge0","byaacq0","max_drinks0","hs_util0")
#rcs = c("PC0","C0","A0")

#################################################################################################
# Model estimation 

# NOTE: Generating female and nonwhite from baseline distributions; predicting the rest of the variables 
#################################################################################################

#Creating log variables 

cdat$logbinge0 = log(cdat$binge0); cdat$logbinge0[cdat$binge0==0] = log(runif(sum(cdat$binge0==0),0,1))
cdat$logbinge1 = log(cdat$binge1); cdat$logbinge1[cdat$binge1==0] = log(runif(sum(cdat$binge1==0),0,1))
cdat$logbinge2 = log(cdat$binge2); cdat$logbinge2[cdat$binge2==0] = log(runif(sum(cdat$binge2==0),0,1))


#-----------------------------------------------------------------------------------------------------------
# Baseline models 
#-----------------------------------------------------------------------------------------------------------

#Greek intent
greek.form1 <- as.formula(paste0("greekintent~", 
                                 "nonwhite + female"))
greek.mod1 = glm(greek.form1,data=cdat,family="binomial"); summary(greek.mod1)

# #Norm percent drink
# norm.form1 <- as.formula(paste0("norm_percent_drink~",
#                                 "nonwhite + female + greekintent"))
# norm.mod1 = lm(norm.form1,data=cdat); summary(norm.mod1)
# 
# #Norm num drinks typical
# norm.form2 <- as.formula(paste0("norm_num_drinks_typical~", 
#                                 "norm_percent_drink+ nonwhite + female + greekintent"))
# norm.mod2 <- glm(norm.form2,data=cdat,family="quasipoisson"); summary(norm.mod2)
# 
# #Norm most drinks typical
# norm.form3 <- as.formula(paste0("norm_most_drinks_typical~", 
#                                 "norm_percent_drink + norm_num_drinks_typical + nonwhite +
#                                 female + greekintent"))
# norm.mod3 <- glm(norm.form3,data=cdat,family="quasipoisson"); summary(norm.mod3)
# 
# #Norm percent binge 
# norm.form4 <- as.formula(paste0("norm_percent_binge~", " nonwhite + female + greekintent",
#                                 "+ norm_percent_drink + norm_num_drinks_typical",
#                                 "+ norm_most_drinks_typical"))
# norm.mod4 <- lm(norm.form4,data=cdat); summary(norm.mod4)


########################### Baseline outcomes models ###################
#Models for binge0 using baseline variables 

#hs_util0
hs.form <-as.formula(paste0("hs_util0~",paste(c("nonwhite",
                                                "female","greekintent"),collapse="+")))
hs.mod <- glm(hs.form,data=cdat,family="binomial"); summary(hs.mod)


#Max_drinks0
maxd.form <- as.formula(paste0("max_drinks0~",paste(c("nonwhite",
                                                      "female","greekintent","hs_util0"),collapse="+")))
maxd.mod <- zeroinfl(maxd.form,data=cdat,dist="poisson")

#byaacq0
by.form <- as.formula(paste0("byaacq0~",paste(c("nonwhite",
                                                "female","greekintent",
                                                "max_drinks0","hs_util0"),collapse="+")))
by.mod <- zeroinfl(by.form,data=cdat,dist="poisson"); summary(by.mod)

#Binge0
binge.form <- as.formula(paste0("binge0~",paste(c("nonwhite",
                                                  "female","greekintent","hs_util0","max_drinks0","byaacq0"),collapse="+")))
#binge.modP <- glm(binge.form,data=cdat,family="poisson")
#binge.mod <- glm(binge.form,data=cdat,family="quasipoisson")
#binge.modNB <- glm.nb(binge.form,data=cdat)
binge.mod <- zeroinfl(binge.form,data=cdat,dist="poisson")

# norm.form1 <- as.formula(paste0("norm_percent_drink~",
#                                 "female + greekintent+binge0+hs_util0+max_drinks0+byaacq0"))
# norm.mod1 = lm(norm.form1,data=cdat); summary(norm.mod1)

#Predicted values of poisson, qp, and neg binom
# pred.p = predict(binge.modP, newdata=simdat,type="response")
# pred.qp = predict(binge.modQP, newdata=simdat, type="response") #same as poisson
# pred.nb = predict(binge.modNB, newdata=simdat, type="response")
# 
# 
# #Calculate dispersion parameter: 
# dp = sum(residuals(binge.mod,type ="pearson")^2)/binge.mod$df.residual
# dp
# 
# #Investigating dispersion
# plot(log(fitted(binge.mod)),log((cdat$binge0-fitted(binge.mod))^2),
#      xlab=expression(hat(mu)),ylab=expression((y-hat(mu))^2),pch=20,col="blue")
# abline(0,1) ## 'variance = mean' line
# 


################## Self monitoring measure ##############################

#sm_binge_last
smbinge.form <- as.formula(paste0("sm_binge_last~",paste(c(dems,"binge0","byaacq0","max_drinks0",'hs_util0',"A1"),
                                                         collapse="+")))

smbinge.mod <- zeroinfl(smbinge.form,data=cdat,dist="poisson"); summary(smbinge.mod)

#sm_hid_last
smhid.form <- as.formula(paste0("sm_hid_last~",paste(c("sm_binge_last",dems,"binge0",
                                                       "byaacq0","max_drinks0",'hs_util0',"A1"),
                                                     collapse="+")))

smhid.mod <- zeroinfl(smhid.form,data=cdat,dist="poisson"); summary(smhid.mod)

########################### Follow-up 1 outcomes ##########################################

#hs_util1
hs.form1 <-  as.formula(paste0("hs_util1~",paste(c('hs_util0',dems,
                                                   "A1","HD","A2","A1:A2"),collapse="+")))
hs.mod1 <- glm(hs.form1, data=cdat, family="binomial"); summary(hs.mod1)

#max_drinks1
maxd.form1 <- as.formula(paste0("max_drinks1~",paste(c("max_drinks0",dems,
                                                       "A1","HD","A2","A1:A2","hs_util1"),collapse="+")))
maxd.mod1 <- zeroinfl(maxd.form1, data=cdat,dist="poisson")

#byaacq1
by.form1 <- as.formula(paste0("byaacq1~",paste(c("byaacq0",dems,
                                                 "A1","HD","A2","A1:A2","hs_util1","max_drinks1"),collapse="+")))
by.mod1 <- zeroinfl(by.form1, data=cdat,dist="poisson")

#binge1
binge.form1 <- as.formula(paste0("binge1~",paste(c("binge0",dems,
                                                   "A1","HD","A2","A1:A2","hs_util1","max_drinks1","byaacq1"),collapse="+")))

binge.mod1 <- zeroinfl(binge.form1,data=cdat,dist="poisson")

########################## Follow-up 2 outcomes #########################################
#hs_util2
hs.form2 <-  as.formula(paste0("hs_util2~",paste(c("hs_util1","hs_util0", dems,
                                                   "A1","HD","A2","A1:A2"),collapse="+")))
hs.mod2 <- glm(hs.form2, data=cdat, family="binomial"); summary(hs.mod2)

#max_drinks2
maxd.form2 <- as.formula(paste0("max_drinks2~",paste(c("max_drinks1","max_drinks0",dems,
                                                       "A1","HD","A2","A1:A2","hs_util2"),collapse="+")))
maxd.mod2 <- zeroinfl(maxd.form2, data=cdat,dist="poisson")

#byaacq2
by.form2 <- as.formula(paste0("byaacq2~",paste(c("byaacq1","byaacq0",dems,"A1","HD","A2","A1:A2","hs_util2","max_drinks2"),collapse="+")))
by.mod2 <- zeroinfl(by.form2, data=cdat,dist="poisson")

#binge2
binge.form2 <- as.formula(paste0("binge2~",paste(c("binge1","binge0",dems,"A1","HD","A2","A1:A2","hs_util2","max_drinks2","byaacq2"),collapse="+")))
binge.mod2 <- zeroinfl(binge.form2,data=cdat,dist="poisson")


save(greek.mod1,
     smhid.mod,smbinge.mod,
     binge.mod,binge.mod1,binge.mod2,maxd.mod,maxd.mod1,maxd.mod2,
     by.mod,by.mod1,by.mod2,hs.mod,hs.mod1,hs.mod2,file="mbridgesim.ZPmodels.rda")


######################################################################################################
# Generating overdispersed data 
#https://support.bioconductor.org/p/40207/
# rqpois <- function(n, mu, theta) {
#rnbinom(n = n, mu = mu, size = mu/(theta-1))
#}
#https://stats.stackexchange.com/questions/54682/overdispersion-parameter 


