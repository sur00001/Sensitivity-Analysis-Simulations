###################################################################################
# Data generating mechanism 
###################################################################################
library(dplyr); library(truncnorm)
library(MASS);library(AER); library(pscl)
setwd("~/Sensitivity Analysis - Paper 1/Delta-adjusted-sensitivity-analysis")
#load("~/Sensitivity Analysis - Paper 1/completedat4pred_6-11.rda")
load("~/Sensitivity Analysis - Paper 1/completedat4pred_11-6.rda")
#load("~/Sensitivity Analysis - Paper 1/Delta-adjusted-sensitivity-analysis/NULLmbridgesim.ZPmodels.rda")
#load("~/Sensitivity Analysis - Paper 1/Delta-adjusted-sensitivity-analysis/mbridgesim.ZPmodels.rda")
source("03a-SimpleLogInduceMissingness.R")
#source("01-Simulation_Models.R")
#source("01a-SimpleLogSimulations_Models.R")
load("~/Sensitivity Analysis - Paper 1/Delta-adjusted-sensitivity-analysis/mbridgeSimplesim.logbinge.Normmodels.rda")
set.seed(1996)

# Expit Function 
expit=function(x){
  exp(x)/(1+exp(x))
}

# # #Generate quasipoisson data - this is incorrect; parameterized it wrong
# rqpois <- function(n, mu, theta) {
# rnbinom(n = n, mu = mu, size = mu/(theta-1))
# }

# Simulating data function
gen.data <- function(N,missTV, missY,missprob,monotone.mis=0,rep.no){
  
  #Initialize variables
  female = nonwhite = logbinge0 = logbinge1 = logbinge2 = greekintent = rep(NA,N)
  hs_util0 = rep(NA,N)
  hs_util1 = hs_util2 = A1 = A2 = rep(NA,N)
  
  simdat = data.frame(id=1:N,female,nonwhite,logbinge0,logbinge1,logbinge2,greekintent,hs_util0,hs_util1,
                      hs_util2,A1,A2)
  # set.seed(rep.no*10)
  
  ############################################
  #Simulate baseline data 
  ############################################
  #female
  simdat$female = rbinom(N,1,.630)
  
  #nonwhite 
  simdat$nonwhite = rbinom(N,1,.236)
  
  #Greek intent
  simdat$greekintent = as.factor(ifelse(rbinom(N,1,prob = expit(predict(greek.mod1,newdata=simdat,type="response")))==1,
                                        "Yes/Undecided","No"))
  ### Baseline outcomes ###
  simdat$hs_util0 = rbinom(N,1,prob = expit(predict(hs.mod,newdata=simdat,type="response")))

  
  #logBinge0
  logbinge0.mean = predict(logbinge.mod,newdata=simdat)
  simdat$logbinge0 = rnorm(logbinge0.mean,sigma(logbinge.mod))

  ### Randomization: early vs late ###
  simdat$A1 = ifelse(rbinom(N,1,.5)==1,1,-1)
  
  ################## Self monitoring measures ##############################
  
  #sm_binge_last
  
  logsmbinge.mean = predict(logsmbinge.mod,newdata=simdat)
  simdat$logsm_binge_last = rnorm(logsmbinge.mean,sigma(logsmbinge.mod))
  
  fulldat1 = simdat 
  ############################ Induce missingness ###########################
  if (missTV!="none"){
    alpha_m = logit(missprob-.02)
    
    simdat$logsm_binge_last = induce.missYt(nperson=N,yt=simdat$logsm_binge_last,wide.df = simdat,t=.5,
                                         alpha_m= alpha_m, 
                                         beta_y= .3,
                                         beta_a1 = 1, beta_maxd = 0,beta_hs = 0,beta_by=0,beta_f = 0,
                                         beta_a1smb = .01,
                                         miss_type=missTV,miss_prob = missprob)
    
    
    
    #simdat$sm_hid_last[is.na(simdat$sm_binge_last)] = NA #if a person is missing binge_last they are also missing hid_last
    
    prop.table(table(is.na(simdat$logsm_binge_last))) 
    
  }
  
  ############################ Tailoring variable and 2nd randomization ######################
  #simdat$HD = ifelse(simdat$sm_binge_last<2 | simdat$sm_hid_last<1 | is.na(simdat$sm_binge_last),0,1) #if missing, they are a light drinker
  simdat$logHD = ifelse(simdat$logsm_binge_last< log(3)|is.na(simdat$logsm_binge_last),0,1)
  prop.table(table(simdat$logHD))
  simdat$A2 = 0 
  simdat$A2[which(simdat$logHD==1)] = ifelse(rbinom(sum(simdat$logHD==1),1,.5)==1, 1, -1)
  prop.table(table(simdat$A2))
  
  ########################### Follow-up 1 outcomes ##########################################
  #hs util1
  simdat$hs_util1 = rbinom(N,1,prob = expit(predict(hs.mod1,newdata=simdat,type="response")))

  #Binge1
  logbinge1.mean = predict(logbinge.mod1,newdata=simdat)
  simdat$logbinge1 = rnorm(logbinge1.mean,sigma(logbinge.mod1))
  

  ########################## Follow-up 2 outcomes #########################################
  
  #hs util2
  simdat$hs_util2 = rbinom(N,1,prob = expit(predict(hs.mod2,newdata=simdat,type="response")))
  
  #Binge2
  logbinge2.mean = predict(logbinge.mod2,newdata=simdat)
  simdat$logbinge2 = rnorm(logbinge2.mean,sigma(logbinge.mod2))
  
  
  fulldat= simdat
  if(missY!="none"){
    alpha_m = logit(missprob-.02)
    
    simdat$logbinge1 = induce.missYt(nperson=N,t=1,yt=simdat$logbinge1,wide.df=fulldat,miss_type=missY,
                                  alpha_m= alpha_m,  beta_y= .3, beta_a1 = 1, beta_a2 = 1,
                                  beta_maxd = 0,beta_hs = 0,beta_by=0,beta_hd = 0, beta_a1logbinge2 = .05,beta_a2logbinge2 = .05,
                                  miss_prob= missprob)  #missprob only used for mcar 
    
     prop.table(table(is.na(simdat$logbinge1)))
    
    if (monotone.mis == 1){
      fu1.misI = which(is.na(simdat$logbinge1)) #inidivudals who drop out at fu1
      #simulate fu2 missingness for remaining individuals
      simdat$logbinge2 = induce.missYt(nperson=N-length(fu1.misI),t=2,simdat$logbinge2[-c(fu1.misI)],wide.df=fulldat,miss_type=missY,miss_prob=missprob) 
      
      fu2.misI = which(is.na(simdat$logbinge2)) #individuals who drop out at fu2
      simdat[fu1.misI,names(simdat) %in% c("logbinge2","hs_util1","hs_util2","max_drinks1","max_drinks2","byaacq1","byaacq2")] =NA
      simdat[fu2.misI,names(simdat) %in% c("hs_util2","max_drinks2","byaacq2")] =NA
      
    } else{ #non-monotone missingness 
      
      simdat$logbinge2 = induce.missYt(nperson=N,t=2,yt=simdat$logbinge2,wide.df=fulldat,miss_type=missY,
                                    alpha_m= logit(missprob-.02),  beta_y= .3, beta_a1 = 1, beta_a2 = 1,
                                    beta_maxd = 0,beta_hs = 0,beta_by=0,beta_hd = 0,beta_a1logbinge2 = .01,beta_a2logbinge2 = .01, miss_prob=missprob) 
      
      
      simdat$hs_util1 = induce.missYt(nperson=N,t=1,simdat$hs_util1,wide.df=fulldat,miss_type=missY,
                                      alpha_m= logit(missprob-.02),  beta_y=.3, beta_a1 = 1, beta_a2 = 1,
                                      beta_maxd = 0,beta_hs = 0,beta_by=0,beta_hd = 0,beta_a1logbinge2 = .01,beta_a2logbinge2 = .01,miss_prob=missprob)

      simdat$hs_util2 = induce.missYt(nperson=N,t=2,simdat$hs_util2,wide.df=fulldat,miss_type=missY,
                                      alpha_m= logit(missprob-.02),  beta_y= .3, beta_a1 = 1, beta_a2 = 1,
                                      beta_maxd = 0,beta_hs = 0,beta_by=0,beta_hd = 0,beta_a1logbinge2 = .01,beta_a2logbinge2 = .01, miss_prob=missprob)
      # 
      # simdat$max_drinks1 = induce.missYt(nperson=N,t=1,simdat$max_drinks1,wide.df=fulldat,miss_type=missY,
      #                                    alpha_m= logit(missprob-.02),  beta_y= .3, beta_a1 = 1, beta_a2 = 1,
      #                                    beta_maxd = 0,beta_hs = 0,beta_by=0,beta_hd = 0,beta_a1logbinge2 = .01,beta_a2logbinge2 = .01,miss_prob=missprob) 
      # 
      # simdat$max_drinks2 = induce.missYt(nperson=N,t=2,simdat$max_drinks2,wide.df=fulldat,miss_type=missY,
      #                                    alpha_m= logit(missprob-.02),  beta_y= .3, beta_a1 = 1, beta_a2 = 1,
      #                                    beta_maxd = 0,beta_hs = 0,beta_by=0,beta_hd = 0, beta_a1logbinge2 = .01,beta_a2logbinge2 = .01, miss_prob=missprob) 
      # 
      # simdat$byaacq1 = induce.missYt(nperson=N,t=1,simdat$byaacq1,wide.df=fulldat,miss_type=missY,
      #                                alpha_m= logit(missprob-.02),  beta_y= .3, beta_a1 = 1, beta_a2 = 1,
      #                                beta_maxd = 0,beta_hs = 0,beta_by=0,beta_hd = 0, beta_a1logbinge2 = .01,beta_a2logbinge2 = .01, miss_prob=missprob) 
      # 
      # simdat$byaacq2 = induce.missYt(nperson=N,t=2,simdat$byaacq2,wide.df=fulldat,miss_type=missY,
      #                                alpha_m= logit(missprob-.02),  beta_y= .3, beta_a1 = 1, beta_a2 = 1,
      #                                beta_maxd = 0,beta_hs = 0,beta_by=0,beta_hd = 0, beta_a1logbinge2 = .01,beta_a2logbinge2 = .01, miss_prob=missprob) 
    }
    
  }
  
  return(simdat)
}

#Change model coefficients for A1 and A2 in the binge drinking models to run an alternative case 

# # Saving the missingness parameters for scenarios
# mis.betas = data.frame(miss_type="MNAR", miss_prob =.528,yt=.5,monotone=0,alpha_m=0, beta_y=0, beta_by=0,beta_hs=0,
#                                             beta_maxd=0L,beta_a1 = .1,beta_f = 0,
#                                            beta_hd=0,beta_a2=0,beta_a2binge2= 0,
#                                             beta_a1binge2=0,beta_smb=.1,beta_smh=.01,beta_a1smb=.1)
# save(mis.betas,file="missingScenarioBetas.rda")

# simdat$sm_binge_last = induce.missYt(N,t=.5,simdat$sm_binge_last,simdat,miss_type="MNAR",alpha_m=0,beta_a1 = 1,beta_smb=1,beta_smh=1,beta_a1smb=2) #~50% missingness

# simdat$sm_binge_last = induce.missYt(N,t=.5,alpha_m=alpha_ms, beta_y=beta_ys, beta_by=beta_bys,beta_hs=beta_hss,
#                                      beta_maxd=beta_maxds,beta_a1 = beta_a1s,beta_f = beta_fs,
#                                      beta_hd=beta_hds,beta_a2=beta_a2s,beta_a2binge2= beta_a2binge2s,
#                                      beta_a1binge2=beta_a1binge2s,beta_smb=beta_smbs,beta_smh=beta_smhs,beta_a1smb=beta_a1smbs,
#                                      miss_type=missTV, miss_prob =miss_prob)


# simdat$sm_binge_last = induce.missYt(N,yt=simdat$sm_binge_last,wide.df = simdat,t=.5,alpha_m=-alpha_m, beta_y=beta_y, beta_by=beta_by,beta_hs=beta_hs,
#                                      beta_maxd=beta_maxd,beta_a1 = beta_a1,beta_f = beta_f,beta_a2=beta_a2,
#                                      miss_type=missTV,miss_prob = missprob)