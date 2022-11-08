########################################################################################################################################
# This script runs a MNAR multiple imputation that incorporates sensitivity parameters for both the TV AND the outcome
########################################################################################################################################


# load libraries
library(mice);library(dplyr)
library(Hmisc)

#Read in data to impute 'misdat'
#load("~/Dissertation/Mbridge_SensAnal/dat4impute.rda") #only API group


#-------------------------------------------------------------------------------------------------------------------------------------------
# Function to do a multiple imputation with optional sensitivity parameters 

# Inputs: mis.dat (data to be imputed), k1 (late sensitivity parameter), k2 (early sensitivity parameter), m (# of imputed datasets), TV (if there 
#         should be a sensitivity parameter for the tailoring variable), nrep (the rep #; so teh random seeds don't overlap - I don't want the same imputation happening across the reps)
# Output: an complete imputed dataset 

#Note: if k1=0 and k2=0 this is the standard MAR imputation
#-------------------------------------------------------------------------------------------------------------------------------------------

impute.mi<- function(mis.dat,k1=0,k2=0,m,TV=T,iter=35,nrep){
  
  imp_long = c()
  for (i in 1:m){
    #set.seed(i*nrep+1) 
    print(paste("     Imputation:",i))
    
    #name variables 
    dems = c("female","nonwhite","greekintent")
    base.outs = c("binge0","byaacq0","max_drinks0","hs_util0")
    yr1.outs = c("binge1","byaacq1","max_drinks1","hs_util1")
    yr2.outs = c("binge2","byaacq2","max_drinks2","hs_util2")
    
    
    # code binary variables as factors
    mis.dat = mis.dat %>% 
      mutate_at(vars(starts_with("hs_util")),factor)
    
   # apply(mis.dat,2,function(x) sum(is.na(x))) #see how many missing
    
    
    ##########################################################################################################
    #Impute self-monitoring measures and tailoring outcome 
    ##########################################################################################################
    
    #Add missingness indicators for SM, FU1 and FU2 - needed for sensitivity analysis
    mis.dat$misSM <- ifelse(is.na(mis.dat$sm_binge_last),1,0)
    mis.dat$misb1 <- ifelse(is.na(mis.dat$binge1),1,0) #FU1 binge 
    mis.dat$misb2 <- ifelse(is.na(mis.dat$binge2),1,0)
    mis.dat$misby1 <- ifelse(is.na(mis.dat$byaacq1),1,0) #FU1 alcohol consequences outcome
    mis.dat$misby2 <- ifelse(is.na(mis.dat$byaacq2),1,0)
    
    
    #Creating predictor matrix
    sm.vars <- c("sm_binge_last", dems,"binge0","byaacq0","max_drinks0",'hs_util0',"A1")
    
    #Index of variables for sm imputation
    sm.i = which(names(mis.dat) %in% sm.vars)
    
    # shell imputation
    imp = mice(mis.dat[,sm.i],maxit=0)
   # imp$loggedEvents
    
    # look at method
  #  imp$method
    
    # make predictor matrix
    pred = imp$predictorMatrix
    pred[,] = 0
    
    
    pred["sm_binge_last", c(sm.vars)] = 1
    #pred["sm_hid_last", c(sm.vars[-2])] = 1
    
    
    smimp = mice(mis.dat[,c(sm.i)],m=1,method = imp$method,predictorMatrix = pred,maxit=iter,printFlag = FALSE)
    sm.imp <- complete(smimp)
    
    #New dataset with imputed self-monitoring variables
    datnosm <- mis.dat %>% dplyr::select(names(mis.dat)[names(mis.dat) %nin% names(sm.imp)])
    sm.update <- cbind(sm.imp,datnosm)
    
    
    if (TV==T){
      #Add sensitivity parameter to self-monitoring measure 
      sm.update$sm_binge_last[sm.update$A1==-1] = sm.update$sm_binge_last[sm.update$A1==-1] + sm.update$misSM[sm.update$A1==-1]*(k1/2) #late 
      sm.update$sm_binge_last[sm.update$A1==1] = sm.update$sm_binge_last[sm.update$A1==1] + sm.update$misSM[sm.update$A1==1]*(k2/2) #early 
    }
    
    #################################################################################################
    #Update HD status and 2nd randomization
    #################################################################################################
    #Note: online coach = 1, email = -1; late=-1, early=1
    
    #sm.update$HD = ifelse(sm.update$sm_binge_last>1 | sm.update$sm_hid_last>0,1,0)
    sm.update$HD = ifelse(simdat$sm_binge_last<2,0,1)
    
    #sum(sm.update$HD[sm.update$misSM==1] != sm.update$HD2[sm.update$misSM==1]) #34/74 drinkers are now heavy drinkers (46%)
    
    #Updating A2 value for those missing sm_binge and who are now heavy drinkers - 1:1 randomization
    sm.update$A2[sm.update$HD==1 & sm.update$misSM==1] = ifelse(rbinom(sum(sm.update$HD==1 & sm.update$misSM==1),1,.5)==0, -1, 1)
    
    sm.update$trajnochange = ifelse(sm.update$misSM==1 & sm.update$HD==1,0,1) #if missing (light drinker) but after imputation they are a HD, then trt trajectory changed
    #trt traj no change for people who's sm was missing
    
    ####################################################################################################
    # Follow-up 1 outcome imputation and update data with sensitivity parameter
    ####################################################################################################
    
    #Induce structural missingness: make FU1 outcomes missing if trt trajectory changed 
    sm.update[sm.update$trajnochange==0,names(sm.update) %in% yr1.outs] = NA
    sm.update[sm.update$trajnochange==0,names(sm.update) %in% yr2.outs] = NA
    
    #Using the MICE chained equations to impute all follow-up outcomes at once, so we can incorporate follow-up 2 information if available for people who didn't have their trajectory change
    
    #Creating predictor matrix
    FU.vars <- c(base.outs,yr1.outs, yr2.outs, dems,
                 "A1","A2",
                 "HD",
                 "sm_binge_last")
    
    
    #Index of variables for FU imputation
    FU.i = which(names(sm.update) %in% FU.vars)
    
    # shell imputation
    imp = mice(sm.update[,FU.i],maxit=0,printFlag = FALSE)
  #  imp$loggedEvents
    
    # look at method
    meth=imp$method
    
    # make predictor matrix
    pred = imp$predictorMatrix
    pred[,] = 0
    
    
    # Follow-up outcomes - finish putting these predictor matrices in 
    pred["binge1",c("binge0",dems,
                    "A1","HD","A2","hs_util1","max_drinks1","byaacq1")] = 1
    
    pred["binge2",c("binge1","binge0",dems,"A1","HD","A2","hs_util2","max_drinks2","byaacq2")] = 1
    
    pred["max_drinks1",c(base.outs,"max_drinks2",dems,"binge1","byaacq1","hs_util1",
                         "A1","A2","sm_binge_last","HD")] = 1
    
    pred["max_drinks2",c(base.outs,"max_drinks1",dems,"binge2","byaacq2","hs_util2",
                         "A1","A2","HD",
                         "sm_binge_last")] = 1
    
    pred["byaacq1",c(base.outs,"byaacq2",dems,"binge1","max_drinks1","hs_util1",
                     "A1","A2","HD","sm_binge_last")] = 1
    
    pred["byaacq2",c(base.outs,"byaacq1",dems,"binge2","max_drinks2","hs_util2","A1","A2",
                     "sm_binge_last","HD")] = 1
    
    
    pred["hs_util1",c(base.outs,"hs_util2",dems,"binge1","max_drinks1","byaacq1",
                      "A1","A2",
                      "sm_binge_last","HD")] = 1
    
    pred["hs_util2",c(base.outs,"hs_util1",dems,"binge2","max_drinks2","byaacq2",
                      "A1","HD","A2",
                      "sm_binge_last","HD")] = 1
    
    
    #Impute data
    fuimp = mice(sm.update[,FU.i],m=1,method = meth,predictorMatrix = pred,maxit=30,printFlag = FALSE)
  #  fuimp$loggedEvents
    
    fu.imp = complete(fuimp)
    
    #New dataset with imputed follow-up variables
    datnofu <- sm.update[,-FU.i]
    fu.update <- cbind(fu.imp,datnofu)
    
    #Add sensitivity parameter to follow-up 1
    fu.update$binge1[fu.update$A1==-1] = fu.update$binge1[fu.update$A1==-1] + fu.update$misb1[fu.update$A1==-1]*k1 #late 
    fu.update$binge1[fu.update$A1==1] = fu.update$binge1[fu.update$A1==1] + fu.update$misb1[fu.update$A1==1]*k2 #early 
    
    #Add sensitivity parameter to follow-up 2
    fu.update$binge2[fu.update$A1==-1] = fu.update$binge2[fu.update$A1==-1] + fu.update$misb2[fu.update$A1==-1]*k1 #late 
    fu.update$binge2[fu.update$A1==1] = fu.update$binge2[fu.update$A1==1] + fu.update$misb2[fu.update$A1==1]*k2 #early 
    
    
    #Add this imputed dataset to a long dataset of all m imputated datasets
    dat=fu.update
    dat$.imp = rep(i,dim(fu.update)[1]) #want length of imputed data
    imp_long = rbind(imp_long,dat)
  } 
  
  return(imp_long) #This imp_long is like Grace's mbridge_long
}





























###############################################
# Comparing my imputations to Grace's 

# library(ggplot2)
# library(esquisse)
# grace.dat = complete(grace.imp)
# df = data.frame(cbind(
#   c(grace.dat$binge1[grace.dat$API==1],
#                       dat$binge1,
#     grace.dat$binge2[grace.dat$API==1],
#     dat$binge2),
#   rep(c(rep("Grace imputation",591),rep("My imputation",591)),2),
#   c(rep("Binge1",591*2),rep("Binge2",591*2))))
# 
# names(df) = c("Binge","Imputation","Follow-up")
# esquisser(df)
# 
# ggplot(df) +
#   aes(x = Binge) +
#   geom_bar(fill = "#112446") +
#   theme_gray() +
#   facet_grid(vars(Imputation), vars(`Follow-up`))
# 
# plot(density(grace.dat$binge1[grace.dat$API==1]),
#      main="Grace's Binge1 imputation")
# 
# plot(density(dat$binge1),
#      main="My Binge1 imputation")
# 
# plot(hist(grace.dat$binge1[grace.dat$API==1]),
#      main="Grace's Binge1 imputation")
# 
# plot(hist(dat$binge1),
#      main="My Binge1 imputation")


# check for convergence with 35 more iterations
# imp2 = mice.mids(imp, maxit=35)
# plot(imp2) # pretty good mixing

# save results
#save(mbridge,imp,file="~/Desktop/MBridge/mi-mbridge.Rdata")
#save(imp_long,file="imp_longSA_k1_1_k2_.5.rda")
