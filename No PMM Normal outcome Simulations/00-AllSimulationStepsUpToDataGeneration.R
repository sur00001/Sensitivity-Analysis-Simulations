######################################################################################
# This file includes all steps of the missing data SMART simulation 
# Only non-montonone missingness, MI and CC coded in this file 
# Two analysis options: comparing all 4 regime means and first stage treatment effect 

#11/12/22 Only up to generating data is here for easier debugging of why mnar missing 
# and full data distributions are not different
######################################################################################

#Libraries 
library(dplyr); library(truncnorm)
library(MASS);library(AER); library(pscl)

setwd("~/Sensitivity Analysis - Paper 1/Sensitivity-Analysis-Simulations/No PMM Normal outcome Simulations")

#Open empirical data for estimating the simulation models (includes log outcomes)
load("C:/Users/lylae/OneDrive/Documents/Sensitivity Analysis - Paper 1/completedat4pred_lognormbinge.rda")
#Group variables into vector 
dems = c("female","nonwhite","greekintent")

# Null regime means estimated from monte carlo integration 
# using mbridge simulation models where I changed the treatment effect coefficients to be 0
#load("~/Sensitivity Analysis - Paper 1/Delta-adjusted-sensitivity-analysis/nulltruth.res.rda")

# Expit Function 
expit=function(x){
  exp(x)/(1+exp(x))
}

#-------------------------------------------------------------------------------
# SIMULATION PARAMETERS
#-------------------------------------------------------------------------------

simres = data.frame()  # initialize data frame that will store results for each rep
nrep = 1000            # Number of simulation reps
set.seed(1996) 
analysis = "timing"    # Theta of interest: "timing" or "regime means"
missTV = "MNAR"        # missingness mechanism for tailoring variable 
missY="MNAR"           # missingness mechanism for Y
missprob=.2            # probability of induced missingness in TV, Y1 and Y2\
analysis_type = "CC"   # multiple imputation or complete case analysis
n=1000                 # sample size 
m=10                   # Number of imputation datasets
truth = truemeans.null # "truemeans.null" are the "true" regime means estimated from Monte Carlo integration

#-------------------------------------------------------------------------------
# PARAMETERS IN LOGIT MODEL TO INDUCE MISSINGNESS 
#-------------------------------------------------------------------------------
alpha_m = logit(missprob-.02) #intercept to control missingness rate 
beta_y0 = 4 #effect of baseline outcome
beta_a1 = 1 #effect of first stage treatment
beta_a2 = 1 #effect of 2nd stage treatment
beta_x1 = 1 #effect of a baseline covariate (female)
beta_a1tv = 1 #effect of the interaction betw stage1 trt and intermediate tailoring var
beta_tv = 1 #effect of intermediate tailoring var
beta_y = 16 #effect of outcome at follow-up 1 and 2
beta_hd = 1 #effect of being a heavy drinker
beta_a1logbinge = 4 #effect of interaction betw stage1 trt and follow-up outcome 1 or 2
beta_a2logbinge = 4 #effect of interaction betw stage2 trt and follow-up outcome 1 or 2


#-------------------------------------------------------------------------------
# MODEL ESTIMATION

# NOTE: Generating female and nonwhite from baseline distributions; predicting the rest of the variables 
#-------------------------------------------------------------------------------

#Creating log variables - now saved to completedat4pred_lognormbinge.rda file

# cdat$logbinge0 = log(cdat$binge0); cdat$logbinge0[cdat$binge0==0] = log(runif(sum(cdat$binge0==0),0,1))
# cdat$logbinge1 = log(cdat$binge1); cdat$logbinge1[cdat$binge1==0] = log(runif(sum(cdat$binge1==0),0,1))
# cdat$logbinge2 = log(cdat$binge2); cdat$logbinge2[cdat$binge2==0] = log(runif(sum(cdat$binge2==0),0,1))
# cdat$logsm_binge_last = log(cdat$sm_binge_last); 
# cdat$logsm_binge_last[cdat$sm_binge_last==0 & !is.na(cdat$sm_binge_last)] = log(runif(length(which(cdat$sm_binge_last==0)),0,1))


########################### Baseline demographic models ###################

#Group variables into vectors 
dems = c("female","nonwhite","greekintent")

#Greek intent
greek.form1 <- as.formula(paste0("greekintent~", 
                                 "nonwhite + female"))
greek.mod1 = glm(greek.form1,data=cdat,family="binomial"); 

#Note: Female and Gender are simulated from binomial distributions 

########################### Baseline outcomes models ###################
#Models for binge0 using baseline variables 

#hs_util0
hs.form <-as.formula(paste0("hs_util0~",paste(c("nonwhite",
                                                "female","greekintent"),collapse="+")))
hs.mod <- glm(hs.form,data=cdat,family="binomial");summary(hs.mod)


#Binge0
logbinge.form <- as.formula(paste0("logbinge0~",paste(c("nonwhite",
                                                        "female","greekintent","hs_util0"),collapse="+")))
logbinge.mod = lm(logbinge.form,data=cdat);summary(logbinge.mod)

################# Tailoring variable model ##############################

#log sm_binge_last 
logsmbinge.form = as.formula(paste0("logsm_binge_last~",paste(c(dems,"logbinge0",'hs_util0',"A1"),collapse="+")))
logsmbinge.mod = lm(logsmbinge.form,data=cdat); summary(logsmbinge.mod)
#Null case: 
logsmbinge.mod$coefficients[7] = 0

#arbritary threshold chosen so there aren't over 60% nonresponders (log(2) led to 70% nonresponders)
cdat$logHD = ifelse(cdat$logsm_binge_last<log(3) | is.na(cdat$logsm_binge_last), 0,1)

########################### Follow-up 1 outcome models ##########################################

#hs_util1
hs.form1 <-  as.formula(paste0("hs_util1~",paste(c('hs_util0',dems,
                                                   "A1","logHD","A2","logbinge0"),collapse="+")))
hs.mod1 <- glm(hs.form1, data=cdat, family="binomial")
hs.mod1$coefficients[c(6,8)]=0 #null case

#logbinge1
logbinge.form1 <- as.formula(paste0("logbinge1~",paste(c("logbinge0",dems,
                                                         "A1","logHD","A2","hs_util1"),collapse="+")))
logbinge.mod1 = lm(logbinge.form1,data=cdat); summary(logbinge.mod1)
logbinge.mod1$coefficients[c(6,8)] = 0 #null case

########################## Follow-up 2 outcome models #########################################
#hs_util2
hs.form2 <-  as.formula(paste0("hs_util2~",paste(c("hs_util1", dems,
                                                   "A1","logHD","A2","logbinge1"),collapse="+")))
hs.mod2 <- glm(hs.form2, data=cdat, family="binomial"); summary(hs.mod2)
hs.mod2$coefficients[c(7,9)]=0


#logbinge2
#NOTE: #logbinge 1 and 2 are correlated in the empirical data, yet when I generate the data there is almost 0 correlation (model R2 is near 0)


#Original model 
#logbinge.form2 <- as.formula(paste0("logbinge2~",paste(c("logbinge1","logbinge0",dems,"A1","logHD","A2","hs_util2"),collapse="+"))) #original model
logbinge.form2 <- as.formula(paste0("logbinge2~",paste(c("logbinge1",dems,"A1","logHD","A2","hs_util2"),collapse="+")))
logbinge.mod2 <- lm(logbinge.form2,data=cdat); summary(logbinge.mod2) #getting an idea of how correlated logbinge1 and logbinge2 are 
logbinge.mod2$coefficients[c(7,9)]=0

#-------------------------------------------------------------------------------
# START SIMULATION
#-------------------------------------------------------------------------------

#for (rep in 1:nrep){

#-----------------------------
# Generate data  
#-----------------------------
# Simulating data using models estimated above
  
  #Initialize variables to be simulated 
  female = nonwhite = logbinge0 = logbinge1 = logbinge2 = greekintent = rep(NA,N)
  hs_util0 = hs_util1 = hs_util2 = A1 = A2 = rep(NA,N)
  #Initialize data set 
  fulldata = data.frame(id=1:N,female,nonwhite,logbinge0,logbinge1,logbinge2,greekintent,hs_util0,hs_util1,
                      hs_util2,A1,A2)
  
  ############################ Baseline data ###########################
  
  #female
  fulldata$female = rbinom(N,1,.630) #probabilities from empirical distribution 
  
  #nonwhite 
  fulldata$nonwhite = rbinom(N,1,.236) #probabilities from empirical distribution
  
  #Greek intent
  fulldata$greekintent = as.factor(ifelse(rbinom(N,1,prob = expit(predict(greek.mod1,newdata=fulldata,type="response")))==1,
                                        "Yes/Undecided","No"))
  ### Baseline outcomes ###
  fulldata$hs_util0 = rbinom(N,1,prob = expit(predict(hs.mod,newdata=fulldata,type="response")))
  
  #logBinge0
  logbinge0.mean = predict(logbinge.mod,newdata=fulldata)
  fulldata$logbinge0 = rnorm(logbinge0.mean,sigma(logbinge.mod))
  
  ### Randomization: early vs late ###
  fulldata$A1 = ifelse(rbinom(N,1,.5)==1,1,-1)
  
  ################## Intermediate variable for tailoring ##############################
  
  logsmbinge.mean = predict(logsmbinge.mod,newdata=fulldata)
  fulldata$logsm_binge_last = rnorm(logsmbinge.mean,sigma(logsmbinge.mod))
  
  
  ############################ Induce missingness in intermediate tailoring variable ###########################
  
  # Missingness rate is controlled by using an intercept that is the logit(desired rate) 
  # + centered covariates in logit missingness model
    
  obsdata = fulldata #keep fulldata and observed data with missingness separate 
   
  rand.prob = runif(N,min=0,max=1) #vector of uniform distributed values for each person 
  
   
  if (missTV=="MCAR"){yt.missprob = rep(miss_prob,N)} #probability of missingness is a constant 
   
  if(missTV=="MAR"){ #probability of missingness is a function of previously observed data 
       yt.missprob = expit(alpha_m + 
                             beta_y0*(fulldata$logbinge0-mean(fulldata$logbinge0)) +
                             beta_a1*(fulldata$A1-mean(fulldata$A1)) +
                             beta_x1*(fulldata$female-mean(fulldata$female))) 
       }
   
  if(missTV=="MNAR"){ 
       yt.missprob = expit(alpha_m + 
                             beta_tv*(fulldata$logsm_binge_last-mean(fulldata$logsm_binge_last)) +
                             beta_a1*(fulldata$A1-mean(fulldata$A1)) + 
                             #interaction of tailoring variable and first stage treatment
                             beta_a1smb*((fulldata$A1-mean(fulldata$A1))*(fulldata$logsm_binge_last-mean(fulldata$logsm_binge_last)))) 
   } 
  obsdata$logsm_binge_last[rand.prob < yt.missprob] = NA

  # Check if missingness is controlled close to desired rate 
  prop.table(table(is.na(obsdata$logsm_binge_last))) 
  
  # Check if missing data distribution of tailoring variable is significantly different than observed data distribution 
  # Should be if MNAR, should not be if MAR or MCAR
  t.test(fulldata$logsm_binge_last~is.na(obsdata$logsm_binge_last))
    
  ############################ Response status and 2nd randomization ######################
  
  #Arbitrary threshold 
  #logHD = 1 is a non-responder
  fulldata$logHD = ifelse(fulldata$logsm_binge_last< log(3),0,1) 
  obsdata$logHD = ifelse(obsdata$logsm_binge_last< log(3)|is.na(obsdata$logsm_binge_last),0,1) 
  
  #2nd randomization for full data
  fulldata$A2 = 0 
  fulldata$A2[which(fulldata$logHD==1)] = ifelse(rbinom(sum(fulldata$logHD==1),1,.5)==1, 1, -1)
  prop.table(table(fulldata$A2)) # Check distribution
  
  #2nd randomization for observed data 
  obsdata$A2 = 0 
  obsdata$A2[which(obsdata$logHD==1)] = ifelse(rbinom(sum(obsdata$logHD==1),1,.5)==1, 1, -1)
  prop.table(table(obsdata$A2)) # Check distribution
  #data with missingness should have more responders than fulldata (since missing = responder)
  
  ########################### Follow-up 1 outcomes ##########################################
  #hs util1
  fulldata$hs_util1 = rbinom(N,1,prob = expit(predict(hs.mod1,newdata=fulldata,type="response")))
  obsdata$hs_util1 = rbinom(N,1,prob = expit(predict(hs.mod1,newdata=obsdata,type="response")))
  
  #Binge1
  fulldata$logbinge1 = rnorm(predict(logbinge.mod1,newdata=fulldata),sigma(logbinge.mod1))
  obsdata$logbinge1 = rnorm(predict(logbinge.mod1,newdata=obsdata),sigma(logbinge.mod1))
  
  ########################## Follow-up 2 outcomes #########################################
  #hs util2
  fulldata$hs_util2 = rbinom(N,1,prob = expit(predict(hs.mod2,newdata=fulldata,type="response")))
  obsdata$hs_util2 = rbinom(N,1,prob = expit(predict(hs.mod2,newdata=obsdata,type="response")))
  
  #Binge2
  fulldata$logbinge2 = rnorm(predict(logbinge.mod2,newdata=fulldata),sigma(logbinge.mod2))
  obsdata$logbinge2 = rnorm(predict(logbinge.mod2,newdata=obsdata),sigma(logbinge.mod2))
 
  ############################ Induce missingness in follow-up outcomes ###########################
  rand.prob.binge1 = runif(N,min=0,max=1) #uniform distributed values for everyone's binge1 missingness 
  rand.prob.binge2 = runif(N,min=0,max=1) #uniform distributed values for everyone's  binge2 missingness 
  rand.prob.hs1 = runif(N,min=0,max=1) #uniform distributed values for everyone's hs_util1 missingness 
  rand.prob.hs2 = runif(N,min=0,max=1) #uniform distributed values for everyone's hs_util2 missingness 
  
  wide.df = obsdata #need current full data version of observed data (this full data version has tv missing)
  
  if (missY=="MCAR"){yt.missprob1 = yt.missprob2 = rep(miss_prob,N)} #probability of missingness is a constant 
  
  if(missY=="MAR"){ #probability of missingness is a function of previously observed data 
    yt.missprob1 = yt.missprob2 = expit(
                          alpha_m+
                          beta_y0*(wide.df$logbinge0-mean(wide.df$logbinge0))+
                          beta_a1*(wide.df$A1-mean(wide.df$A1))+
                          beta_a2*(wide.df$A2-mean(wide.df$A2))+ 
                          beta_hd*(wide.df$logHD-mean(wide.df$logHD)))
  }
  
  if(missY=="MNAR"){ 
    #Missingness probability at follow-up 1 relies on follow-up 1 outcome
    yt.missprob1 = expit(alpha_m+
                          beta_a1*(wide.df$A1-mean(wide.df$A1)) + beta_a2*(wide.df$A2-mean(wide.df$A2)) +
                          beta_y*(wide.df$logbinge1-mean(wide.df$logbinge1))+ 
                          beta_a1logbinge*((wide.df$A1-mean(wide.df$A1))*(wide.df$logbinge1-mean(wide.df$logbinge1)))+ 
                          beta_a2logbinge*((wide.df$A2-mean(wide.df$A2))*(wide.df$logbinge1-mean(wide.df$logbinge1))))
    
    #Missingness probability at follow-up 2 relies on follow-up 2 outcome 
    yt.missprob2 = expit(alpha_m+
                           beta_a1*(wide.df$A1-mean(wide.df$A1)) + beta_a2*(wide.df$A2-mean(wide.df$A2)) +
                           beta_y*(wide.df$logbinge2-mean(wide.df$logbinge2))+ 
                           beta_a1logbinge*((wide.df$A1-mean(wide.df$A1))*(wide.df$logbinge2-mean(wide.df$logbinge2)))+ 
                           beta_a2logbinge*((wide.df$A2-mean(wide.df$A2))*(wide.df$logbinge2-mean(wide.df$logbinge2))))
  } 
  obsdata$logbinge1[rand.prob.binge1 < yt.missprob1] = NA
  obsdata$logbinge2[rand.prob.binge2 < yt.missprob2] = NA
  obsdata$hs_util1[rand.prob.hs1 < yt.missprob1] = NA
  obsdata$hs_util2[rand.prob.hs2 < yt.missprob2] = NA
  
  #Check missing data rate
  prop.table(table(is.na(obsdata$logbinge1)))
  prop.table(table(is.na(obsdata$logbinge2)))

  # Check if missing data distribution of binge1 is significantly different than observed data distribution 
  # Should be if MNAR and MAR, should not be if MCAR
  
  par(mfrow=c(2,1))
  hist(fulldata$logbinge1); hist(obsdata$logbinge1,na.rm=TRUE)
  t.test(fulldata$logbinge1~is.na(obsdata$logbinge1)) #note that no. of non-responders, trt2 will be different for people in fulldata vs obsdata
  t.test(fulldata$logbinge2~is.na(obsdata$logbinge2))

  t.test(fulldata$hs_util1~is.na(obsdata$hs_util1)) #note that no. of non-responders, trt2 will be different for people in fulldata vs obsdata
  t.test(fulldata$hs_util2~is.na(obsdata$hs_util2))

  #MNAR missing and full distribution of outcomes are the same...



