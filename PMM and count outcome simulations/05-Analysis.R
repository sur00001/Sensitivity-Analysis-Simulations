# This script runs analyses for imputed data and does a complete case analysis 
# (1) Means at each time point for all four APIs, with test of any differences in proportional change (0 to 1, 0 to 2)
# (2) Pooled means at each time point by randomization (early vs late, coach vs email) with main effects

# Use multiply imputed data from mi-mbridge.R

# read in packages
library(mice)
library(dplyr)
library(geepack)
library(emmeans)
library(ggplot2)
library(readxl)
library(kableExtra)
library(tidyverse)
library(RColorBrewer)
library(grid);library(gtable)

load("~/Research with M-bridge/mi-mbridge.Rdata") #wide data - original M-bridhge imputed data 

#load("~/Sensitivity Analysis - Paper 1/Delta-adjusted-sensitivity-analysis/truthRefMeans.rda") #true regime means
load("~/Sensitivity Analysis - Paper 1/Delta-adjusted-sensitivity-analysis/nulltruth.res.rda") # null regime means

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Helper functions for inference and calculating and formatting values
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
myp.format = function(pvalue) {
  ifelse(pvalue<0.001,"<0.001",format(round(pvalue,3),nsmall=3))
}

# function to fit model
fitMod = function(imp,outcome,fam,int1=T,int2=T,allimp) {
  # subset to one imputed dataset
  mydata = allimp[allimp$.imp==imp,]
  
  
  # fit GEE model w/ exchangeable correlation
  form = ifelse(int1&int2,paste0(outcome," ~ A1*A2 + time1 + time1*A1*A2 + time2 + time2*A1*A2"),
                ifelse(int1==F,paste0(outcome," ~ A1*A2 + time1 + time2 + time2*A1*A2"),
                       ifelse(int2==F,paste0(outcome," ~ A1*A2 + time1 + time1*A1*A2 + time2"),NA)))
  mod = geeglm(as.formula(form),family=fam,data = mydata,id = id,weights = wt,
               corstr = "exchangeable")
  return(mod)
}

## functions to do inference on (exponentiated) linear combinations of parameters across imputed datasets
# function to get linear combination (default: sum) of coefficients and their variance
lin.comb = function(mymod,coef.names,contrast=1) {
  c = rep(0,length(coef(mymod))); names(c) = names(coef(mymod))
  c[coef.names] = contrast
  est = c %*% coef(mymod)
  var = t(c) %*% vcov(mymod) %*% c
  return(data.frame(est=est,var=var))
}

# function to pool linear combination across models fit to imputed datasets
pool.lin.comb = function(imp.mods,coef.names,contrast=1,adj=F,N) {
  est.vars = do.call(rbind,lapply(imp.mods,lin.comb,coef.names,contrast))
  p = ifelse(adj,16,12)
  pooled.res = pool.scalar(est.vars$est,est.vars$var,n=N,k=p)
  est = pooled.res$qbar
  var = pooled.res$t
  df = pooled.res$df
  return(c(est=est,var=var,df=df))
}
expit = function(x) {
  return(exp(x) / (1+exp(x)))
}

# function to exponentiate estimate and get SE via Delta (relies on normal approximation, but everything else uses t)
exp.est = function(est.var,level=.975) {
  est.var = unname(est.var)
  est = est.var[[1]]; var = est.var[[2]]
  exp.est = exp(est)
  se.exp.est = exp.est * sqrt(var)
  exp.lowCI = exp(est - sqrt(var)*1.96)
  exp.highCI = exp(est + sqrt(var)*1.96)
  return(c(exp.est=exp.est,se.exp.est=se.exp.est,
           exp.lowCI=exp.lowCI,exp.highCI=exp.highCI))
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#Function to get timing comparison at follow-up 1
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

timing.comp= function(dat, analysis_type, missing.dat){

if (analysis_type =="CC"){ #dat is observed data
  time1.mod = glm(binge1~A1,data=dat,family="quasipoisson")
  time1.res= summary(time1.mod)
  
  res= data.frame(A1.beta = time1.res$coefficients[2,1],A1.lb =confint(time1.mod)[2,1],
                  A1.ub = confint(time1.mod)[2,2],A1.se = time1.res$coefficients[2,2],
                  A1.pval = time1.res$coefficients[2,4])
  
  
} 
  else { #MI analysis - dat is the "imp" mice object 
   fit = lapply(1:max(dat$.imp),function(x){glm(binge1~A1,family="quasipoisson",data=dat[dat$.imp==x,])})
    time1.res = summary(pool(fit))
    res= data.frame(A1.beta = time1.res$estimate[2],A1.lb = time1.res$estimate[2] - 1.96*time1.res$std.error[2],
                    A1.ub = time1.res$estimate[2] + 1.96*time1.res$std.error[2],A1.se = time1.res$std.error[2],
                    A1.pval = time1.res$p.value[2])
  }
  res$HD.prop = prop.table(table(dat$HD))[2]
  res$miss.rateTV = prop.table(table(is.na(missing.dat$sm_binge_last)))[2] #missingness rate of TV
  res$miss.rateY = prop.table(table(is.na(missing.dat$binge1)))[2] #missingness rate of binge1 (outcome)
  
return(res)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#Function to find regime means
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

regime.means = function(dat, analysis_type,truth,nrep,missing.dat){ #"dat" will have all imputations rbinded for MI scenarios

#----------------------------------------------------
#Create long data w/ all outcomes, all randomizations 
#----------------------------------------------------
  N=max(dat$id)
  
  if (analysis_type =="CC"){ dat$.imp = rep(1, length(dat$id))}
  
  # use pivot_longer to create three rows per person in every imputed dataset
  alldat = dat %>% pivot_longer(
    cols = starts_with(c("binge","max_drinks","byaacq","hs_util")),
    names_to = c(".value", "time"),
    names_pattern = "(.*)(.)") 
  
  alldat = alldat %>% dplyr::select(.imp,id,A1,A2,HD,time,
                                    starts_with(c("binge","max_drinks","byaacq","hs_util")))
  
  alldat$time = as.numeric(alldat$time)
  
  alldat = alldat %>% arrange(.imp,id,time)
  #alldat$hs_util = as.numeric(alldat$hs_util)-1
  head(alldat)
  
  # code A1 and A2
  alldat = alldat %>%
    mutate(A2 = ifelse(HD==0,1,A2),
           wave = time+1)
  
  # make duplicates of the non-heavy drinkers; change A2
  alldat_notHD = subset(alldat,HD==0) %>%
    mutate(A2 = -1,
           wave = time+4)
  alldat = rbind.data.frame(alldat,alldat_notHD)
  
  # sort by person, wave
  alldat = alldat %>% arrange(.imp,id,wave)
  
  # make time a factor
  alldat$time1 = ifelse(alldat$time==1,1,0)
  alldat$time2 = ifelse(alldat$time==2,1,0)
  
  # code weights
  alldat = mutate(alldat,wt=ifelse(HD==1,4,2))
  
  # alldat$time1 = as.factor(alldat$time1)
  # alldat$time2 = as.factor(alldat$time2)
  # alldat$A1=as.factor(alldat$A1)
  # alldat$A2 = as.factor(alldat$A2)
#-----------------------------------------------
# Obtain regime means, SEs and p-values
#-----------------------------------------------

if (analysis_type=="CC"){
  cc.mod = geeglm(binge ~ A1*A2 + time1 + time1*A1*A2 + time2 + time2*A1*A2, family="poisson", data = alldat,
           id = id, weights = wt, corstr = "exchangeable")
  
  #cc.mod = lapply(1:max(alldat$.imp),fitMod,outcome="binge",fam="poisson") #tried to do a hack - put in list form so I can run D1() and check it with anova(). D1 could not find ubar (between variabnce so it didn't work)
 # cc.noint = lapply(1:max(alldat$.imp),fitMod,outcome="binge",fam="poisson",int1=F)
  
  cc.noint = geeglm(binge ~  A1*A2 + time1 + time2 + time2*A1*A2, family="poisson", data = alldat,
                    id = id, weights = wt, corstr = "exchangeable")
  
  EC1 = exp.est(lin.comb(cc.mod,c("time1","A1:time1","A2:time1","A1:A2:time1"),c(1,1,1,1))) # early/coach
  EE1 = exp.est(lin.comb(cc.mod,c("time1","A1:time1","A2:time1","A1:A2:time1"),c(1,1,-1,-1))) # early/email
  LC1 = exp.est(lin.comb(cc.mod,c("time1","A1:time1","A2:time1","A1:A2:time1"),c(1,-1,1,-1))) # late/coach
  LE1 = exp.est(lin.comb(cc.mod,c("time1","A1:time1","A2:time1","A1:A2:time1"),c(1,-1,-1,1))) # late/email
  
  res = data.frame(rbind(EC1,EE1,LC1,LE1))
  
  #p1 = unname(D1(cc.mod,cc.noint)$result[1,4])
 p1 = anova(cc.mod,cc.noint,test="wald")$`P(>|Chi|)`
  
} else { #analysis with imputation 
  
  mi.mod = lapply(1:max(alldat$.imp),fitMod,outcome="binge",fam="poisson",allimp=alldat)
  mi.noint = lapply(1:max(alldat$.imp),fitMod,outcome="binge",fam="poisson",int1=F,allimp=alldat)
  
  EC1 = exp.est(pool.lin.comb(mi.mod,c("time1","A1:time1","A2:time1","A1:A2:time1"),c(1,1,1,1),N=N)) # early/coach
  EE1 = exp.est(pool.lin.comb(mi.mod,c("time1","A1:time1","A2:time1","A1:A2:time1"),c(1,1,-1,-1),N=N)) # early/email
  LC1 = exp.est(pool.lin.comb(mi.mod,c("time1","A1:time1","A2:time1","A1:A2:time1"),c(1,-1,1,-1),N=N)) # late/coach
  LE1 = exp.est(pool.lin.comb(mi.mod,c("time1","A1:time1","A2:time1","A1:A2:time1"),c(1,-1,-1,1),N=N)) # late/email
  
  res = data.frame(rbind(EC1,EE1,LC1,LE1))
  # info on D1(): https://stefvanbuuren.name/fimd/sec-multiparameter.html
  # test for difference bw groups in change from baseline to fu1 
  p1 = unname(D1(mi.mod,mi.noint,dfcom=N-12)$result[1,4]) # c("A1:time1","A2:time1","A1:A2:time1")) #This function only works with mice() objects
}

#-----------------------------------------------
# Coverage prob, % reject, bias
#-----------------------------------------------
  res$truth = truth 
  res$CovP = res$truth > res$exp.lowCI & res$truth < res$exp.highCI
  res$rej = p1 < .05
  res$bias = abs(res$truth - res$exp.est)
  res$rep = nrep
  res$HD.prop = prop.table(table(dat$HD))[2]
  res$reg = c("EC1","EE1","LC1","LE1")
  res$N = max(dat$id)
  res$pval = p1
  res$miss.rateTV = prop.table(table(is.na(missing.dat$sm_binge_last)))[2] #missingness rate of TV
  res$miss.rateY = prop.table(table(is.na(missing.dat$binge1)))[2] #missingness rate of binge1 (outcome)
  res$HD.prop = prop.table(table(dat$HD))[2]
  
  
  rownames(res) = NULL
  
  return(res)

} 




















#-------------------------------------------------------------------------------
# Extra analyses and plots 
#-------------------------------------------------------------------------------

# rownames(allplotdat) = NULL
# allplotdat = mutate(allplotdat,
#                     mean = ifelse(fam=="poisson",exp.est,expit.est),
#                     lowCI = ifelse(fam=="poisson",exp.lowCI,expit.lowCI),
#                     highCI = ifelse(fam=="poisson",exp.highCI,expit.highCI))
# allplotdat = mutate(allplotdat,
#                     time = as.numeric(time),
#                     time = ifelse(outcome!="Binge drinking frequency",time+4,time),
#                     time = ifelse(outcome!="Binge drinking frequency" & time==6,7,time))
# addon = allplotdat[allplotdat$outcome=="Binge drinking frequency" & allplotdat$time==2,] %>%
#   mutate(time=3,mean=NA,lowCI=NA,highCI=NA)
# allplotdat2 = rbind.data.frame(allplotdat,addon)
# g = ggplot(allplotdat2,aes(x=time,y=mean,group=group)) + 
#   geom_point(position=position_dodge(.75),size=2.5,aes(shape=group)) +
#   geom_errorbar(aes(ymin=lowCI,ymax=highCI),width=.5,position=position_dodge(.75)) +
#   # scale_color_brewer(name="Group",palette="Dark2") +
#   # scale_linetype_manual(name="Group",values=c("solid","dashed","dotdash","dotted")) +
#   scale_shape_manual(name="Group",values=c(16,17,15,18)) +
#   ylab("Marginal mean (95% CI)") + xlab("") +
#   scale_y_continuous(n.breaks=8,limits=c(0,NA)) +
#   scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7),
#                      labels=c("Baseline","Follow-up 1","Follow-up 2\nretrospective\nmeasure","Campus\nshutdown",
#                               "Baseline","Follow-up 1","Campus\nshutdown","Follow-up 2\nconcurrent\nmeasure")) + 
#   theme(axis.text.x = element_text(size=7)) +
#   facet_wrap(~outcome,scales="free")
# g

# # functions to make plot of pooled means by timing and table w/ main effect
# combs.timing = list(e0 = data.frame(names=c("(Intercept)","A1"),cont=1),
#                     e1 = data.frame(names=c("(Intercept)","A1","time1","A1:time1"),cont=1),
#                     e2 = data.frame(names=c("(Intercept)","A1","time2","A1:time2"),cont=1),
#                     l0 = data.frame(names=c("(Intercept)","A1"),cont=c(1,-1)),
#                     l1 = data.frame(names=c("(Intercept)","A1","time1","A1:time1"),cont=c(1,-1,1,-1)),
#                     l2 = data.frame(names=c("(Intercept)","A1","time2","A1:time2"),cont=c(1,-1,1,-1)))
# 
# View(rbind(pooledTimingTable(bingefreq.itt,bingefreq.itt2,"poisson"),
#            #pooledTimingTable(anybinge.itt,anybinge.itt2,"binomial"),
#            #pooledTimingTable(maxdrink.itt,maxdrink.itt2,"poisson"),
#            pooledTimingTable(byaacq.itt,byaacq.itt2,"poisson"),
#            pooledTimingTable(hs.itt,hs.itt2,"binomial")))
# 
# # functions to make plot of pooled means by coach/email and table w/ main effect
# combs.bridge = list(c0 = data.frame(names=c("(Intercept)","A2"),cont=1),
#                     c1 = data.frame(names=c("(Intercept)","A2","time1","A2:time1"),cont=1),
#                     c2 = data.frame(names=c("(Intercept)","A2","time2","A2:time2"),cont=1),
#                     e0 = data.frame(names=c("(Intercept)","A2"),cont=c(1,-1)),
#                     e1 = data.frame(names=c("(Intercept)","A2","time1","A2:time1"),cont=c(1,-1,1,-1)),
#                     e2 = data.frame(names=c("(Intercept)","A2","time2","A2:time2"),cont=c(1,-1,1,-1)))
# 
# View(rbind(pooledBridgeTable(bingefreq.itt,bingefreq.itt2,"poisson"),
#            pooledBridgeTable(anybinge.itt,anybinge.itt2,"binomial"),
#            pooledBridgeTable(maxdrink.itt,maxdrink.itt2,"poisson"),
#            pooledBridgeTable(byaacq.itt,byaacq.itt2,"poisson"),
#            pooledBridgeTable(hs.itt,hs.itt2,"binomial")))
# 

