#################################################################################################
# This script estimates the data generating mechanism models that will  be used to simulate the complete data
# This script creates log normal variables and redefines the models 
# It takes out maxdrinks and byaacq to keep things simpler 
#################################################################################################

library(tidyverse) 
library(dplyr)
library(pscl)


#Open complete data 
load("C:/Users/lylae/OneDrive/Documents/Sensitivity Analysis - Paper 1/completedat4pred_lognormbinge.rda")
setwd("~/Sensitivity Analysis - Paper 1/Delta-adjusted-sensitivity-analysis")

#Group variables into vectors 
dems = c("female","nonwhite","greekintent")

#################################################################################################
# Model estimation 

# NOTE: Generating female and nonwhite from baseline distributions; predicting the rest of the variables 
#################################################################################################

#Creating log variables - now saved to completedat4pred_lognormbinge.rda file
# 
# cdat$logbinge0 = log(cdat$binge0); cdat$logbinge0[cdat$binge0==0] = log(runif(sum(cdat$binge0==0),0,1))
# cdat$logbinge1 = log(cdat$binge1); cdat$logbinge1[cdat$binge1==0] = log(runif(sum(cdat$binge1==0),0,1))
# cdat$logbinge2 = log(cdat$binge2); cdat$logbinge2[cdat$binge2==0] = log(runif(sum(cdat$binge2==0),0,1))
# cdat$logsm_binge_last = log(cdat$sm_binge_last); 
# cdat$logsm_binge_last[cdat$sm_binge_last==0 & !is.na(cdat$sm_binge_last)] = log(runif(length(which(cdat$sm_binge_last==0)),0,1))


#-----------------------------------------------------------------------------------------------------------
# Baseline models 
#-----------------------------------------------------------------------------------------------------------

#Greek intent
greek.form1 <- as.formula(paste0("greekintent~", 
                                 "nonwhite + female"))
greek.mod1 = glm(greek.form1,data=cdat,family="binomial"); summary(greek.mod1)


########################### Baseline outcomes models ###################
#Models for binge0 using baseline variables 

#hs_util0
hs.form <-as.formula(paste0("hs_util0~",paste(c("nonwhite",
                                                "female","greekintent"),collapse="+")))
hs.mod <- glm(hs.form,data=cdat,family="binomial"); summary(hs.mod)


#Binge0
logbinge.form <- as.formula(paste0("logbinge0~",paste(c("nonwhite",
                                                  "female","greekintent","hs_util0"),collapse="+")))

logbinge.mod = lm(logbinge.form,data=cdat); summary(logbinge.mod)


################## Self monitoring measure ##############################

#log sm_binge_last 
logsmbinge.form = as.formula(paste0("logsm_binge_last~",paste(c(dems,"logbinge0",'hs_util0',"A1"),collapse="+")))
logsmbinge.mod = lm(logsmbinge.form,data=cdat); summary(logsmbinge.mod)

#Null case: 
logsmbinge.mod$coefficients[7] = 0


#cdat$logHD = ifelse(cdat$logsm_binge_last<log(3) | is.na(cdat$logsm_binge_last), 0,1)


########################### Follow-up 1 outcomes ##########################################

#Add in A1:A2 interaction later 

#hs_util1
hs.form1 <-  as.formula(paste0("hs_util1~",paste(c('hs_util0',dems,
                                                   "A1","logHD","A2","logbinge0"),collapse="+")))
hs.mod1 <- glm(hs.form1, data=cdat, family="binomial"); summary(hs.mod1)
hs.mod1$coefficients[c(6,8)]=0 #null case

#logbinge1
logbinge.form1 <- as.formula(paste0("logbinge1~",paste(c("logbinge0",dems,
                                                   "A1","logHD","A2","hs_util1"),collapse="+")))
logbinge.mod1 = lm(logbinge.form1,data=cdat); summary(logbinge.mod1)
logbinge.mod1$coefficients[c(6,8)] = 0 #null case

########################## Follow-up 2 outcomes #########################################
#hs_util2
hs.form2 <-  as.formula(paste0("hs_util2~",paste(c("hs_util1","hs_util0", dems,
                                                   "A1","logHD","A2","logbinge1"),collapse="+")))
hs.mod2 <- glm(hs.form2, data=cdat, family="binomial"); summary(hs.mod2)
hs.mod2$coefficients[c(7,9)]=0


#logbinge2
logbinge.form2 <- as.formula(paste0("logbinge2~",paste(c("logbinge1","logbinge0",dems,"A1","logHD","A2","hs_util2"),collapse="+")))
logbinge.mod2 <- lm(logbinge.form2,data=cdat); summary(logbinge.mod2)
logbinge.mod2$coefficients[c(7,9)]=0

save(greek.mod1,
     logsmbinge.mod,
     logbinge.mod,logbinge.mod1,logbinge.mod2,hs.mod,hs.mod1,hs.mod2,file="mbridgeSimplesim.logbinge.Normmodels.rda")



######################################################################################################
# Generating overdispersed data 
#https://support.bioconductor.org/p/40207/
# rqpois <- function(n, mu, theta) {
#rnbinom(n = n, mu = mu, size = mu/(theta-1))
#}
#https://stats.stackexchange.com/questions/54682/overdispersion-parameter 


