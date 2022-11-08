#Main function to run simulations 

setwd("~/Sensitivity Analysis - Paper 1/Delta-adjusted-sensitivity-analysis")
source("02-DataGeneration.R")
source("03- InduceMissingness.R")
source("04-MultipleImputation.R")
source("05-Analysis.R")
#load("~/Sensitivity Analysis - Paper 1/Delta-adjusted-sensitivity-analysis/truthRefMeans.rda") #true regime means
load("~/Sensitivity Analysis - Paper 1/Delta-adjusted-sensitivity-analysis/nulltruth.res.rda") #null regime means 
library(parallel)
library(doParallel)
n_cores <- detectCores(logical=FALSE) #parallel nodes

#-------------------------------------------------------------------------------
# Simulation parameters 
#-------------------------------------------------------------------------------

simres = data.frame()
nrep = 3
set.seed(1996)
analysis = "timing"
missTV = "MAR"
missY="MAR"
missprob= .2
analysis_type = "MI"
n=591
k1=0 #sensitivity parameters
k2=0
m=20 #number of imputation datasets
TV = T #tailoring variable imputation 
truth = truemeans.null
monotone.mis = 0 #indicator for monotone missingness 

#-------------------------------------------------------------------------------
# Simulation function
#-------------------------------------------------------------------------------

run.sim = function(rep,N=n,missTV=missTV, missY=missY,missprob=missprob, monotone.mis = monotone.mis,analysis=analysis,analysis_type=analysis_type,
                   k1=k1,k2=k2,m=m,TV=T,iter=35,truth=truth){
#for (rep in 1:nrep){
  
  print(paste("Nrep:",rep))
  
  #-----------------------
  # Generate data 
  #-----------------------
  obs.data = gen.data(N=n,missTV=missTV, missY=missY,missprob=missprob, monotone.mis = 0,rep.no=rep)
  
  if(analysis_type=="CC"){dat4analysis = obs.data} #obs data will be analyzed for complete case analysis
  
  #-----------------------
  #Impute data if analysis_type = "MI"
  #-----------------------
  if(analysis_type=="MI"){
    mi.dat = impute.mi(mis.dat = obs.data,k1=k1,k2=k2,m=m,TV=T,iter=35,nrep=rep)
    dat4analysis = mi.dat
  }
  
  #-----------------------
  # Run analysis 
  #-----------------------
  
    if (analysis=="regime means"){
      res = regime.means(dat=dat4analysis,analysis_type = analysis_type,truth=truth,nrep=rep, missing.dat = obs.data)
    }
  
    if(analysis=="timing"){res=timing.comp(dat=dat4analysis,analysis_type=analysis_type, missing.dat = obs.data)}

  #-----------------------
  # Combine results across reps
  #-----------------------
  #simres = rbind(simres,res)
  return(res)
}

#-------------------------------------------------------------------------------
# Start simulation 
#-------------------------------------------------------------------------------

start = Sys.time()

currentseed=.Random.seed[2]
currentseed=.Random.seed[currentseed+2]
cl <- makeCluster(n_cores)
clusterExport(cl, c("run.sim","induce.missYt","gen.data","abs","mean",
                    "impute.mi","geeglm","n","missY","missprob","missTV","monotone.mis","analysis","analysis_type","k1","k2","m","TV","iter","truth",
                    "expit","binge.mod","binge.mod1","binge.mod2","by.mod","by.mod1","by.mod2","greek.mod1","hs.mod","hs.mod1","hs.mod2",
                    "maxd.mod","maxd.mod1","maxd.mod2","smbinge.mod","smhid.mod","predict","is.na","which","logit","%>%",
                    "mutate_at","vars","starts_with","mice","complete","%nin%","timing.comp","pool","regime.means","pivot_longer",
                    "arrange","mutate","exp.est","lin.comb"))
clusterEvalQ(cl,library("pscl","tidyverse","mice"))
clusterSetRNGStream(cl=cl,currentseed)
thetas = parSapply(cl,1:nrep,FUN=function(irep){run.sim(rep=irep,N=n,missTV=missTV, missY=missY,missprob=missprob, monotone.mis = monotone.mis,analysis=analysis,analysis_type=analysis_type,
                                                        k1=k1,k2=k2,m=m,TV=T,iter=35,truth=truth)})
stopCluster(cl)

simres = data.frame(apply(data.frame(t(thetas)),2,unlist)) #takes parallel results and puts into the simres dataframe 

time = Sys.time() - start
time


#-------------------------------------------------------------------------------
# Summarize simulation results 
#-------------------------------------------------------------------------------
# 
# #-----------------------
# # Timing comparison
# #-----------------------
# 
if (analysis=="timing"){
  simres$rep = 1:nrep
  simres$covp = ifelse(0> simres$A1.lb & 0<simres$A1.ub,1,0)
  simres$rej = ifelse(simres$A1.pval<.05,1,0)
  sumsim = simres %>% summarise(theta = mean(A1.beta),theta.se = sd(A1.beta),theta.lb = theta -1.96*theta.se,
                                theta.ub = theta + 1.96*theta.se,bias = theta-0,Covp = mean(covp),rej.prob = mean(rej),
                                mean.HD.prop = mean(HD.prop),TVmissing = mean(miss.rateTV,na.rm=TRUE), Ymissing = mean(miss.rateY,na.rm=TRUE))
  #sumsim=round(sumsim,3)
}

#-----------------------
# Regime means
#-----------------------
if (analysis=="regime means"){
  
simres$pval=as.numeric(simres$pval)
simres$miss.rateTV = as.numeric(simres$miss.rateTV); simres$miss.rateY=as.numeric(simres$miss.rateY); simres$HD.prop= as.numeric(simres$HD.prop)
simres$exp.est=as.numeric(simres$exp.est); simres$se.exp.est=as.numeric(simres$se.exp.est)
simres$rej2 = simres$pval < .05 #A check to see if my rejection probability is working correctly
simres$CovP = ifelse(simres$CovP=="TRUE",1,0)
simres$reg = as.factor(simres$reg)
sumsim = simres %>% dplyr::group_by(reg) %>% dplyr::summarise(theta = mean(exp.est), theta.se = sd(exp.est),theta.lb = quantile(exp.est,.025),theta.lb2 = theta -1.96*theta.se,
                                               theta.ub = quantile(exp.est,.975), theta.ub2 = theta +1.96*theta.se, avg.se = mean(se.exp.est),covp = mean(CovP)*100,rej.prob =
                                                 mean(rej2),TVmissing = mean(miss.rateTV,na.rm=TRUE), Ymissing = mean(miss.rateY,na.rm=TRUE),
                                                 mean.HD.prop = mean(HD.prop))
sumsim$truth = truemeans.null
sumsim$bias = sumsim$theta - sumsim$truth
#rej.prob = mean(simres$rej2[simres$reg=="EE1"])
#sumsim=round(sumsim,3)
}

#-----------------------
#Add parameters that describe the missingness scenario
#-----------------------

#Adding to the individual simulation results dataframe
simres$sim.type = rep(analysis_type,dim(simres)[1]) # MI or CC
simres$missTV = rep(missTV,dim(simres)[1]) #Missingness mechanism for TV
simres$missY = rep(missY,dim(simres)[1]) #Missingness mechanism for binge outcome
simres$N = n #Sample size

#Adding same parameters to the summary of simulation results dataframe
sumsim$sim.type = rep(analysis_type,length(sumsim$theta))
sumsim$missTV = rep(missTV,length(sumsim$theta))
sumsim$missY = rep(missY,length(sumsim$theta))
sumsim$N = n
sumsim$nreps = nrep
sumsim$zstat = sumsim$bias/(sumsim$theta.se/sqrt(sumsim$nreps)) #testing if bias is significantly different from 0

#Formatting final results dataframe
sumres = data.frame(sumsim$sim.type,sumsim$missTV,sumsim$TVmissing,sumsim$Ymissing,sumsim$theta,sumsim$theta.se,sumsim$zstat,
                    sumsim$theta.lb,sumsim$theta.ub,
                    sumsim$bias,sumsim$Covp,sumsim$rej.prob,sumsim$mean.HD.prop,sumsim$N,sumsim$nreps)
#sumres[c(3,4,5,6,7,8,9,10,.11,12,13)] = round(sumres[c(3,4,5,6,7,8,9,10,.11,12,13)],3)
names(sumres) = c("Analysis","Mech",	"%TVmis","%Ymis",	"Theta-hat","SE","Zstat","Lb","Ub","Bias",
                  "CovP","Type1","%NR","N","Reps")		

head(simres) #results from each rep
sumres #summary of all reps


