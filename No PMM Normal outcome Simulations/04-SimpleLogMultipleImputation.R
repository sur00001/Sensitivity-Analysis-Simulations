########################################################################################################################################
# This script runs a MNAR multiple imputation that incorporates sensitivity parameters for both the TV AND the outcome
########################################################################################################################################

# load libraries
library(mice);library(dplyr)
library(Hmisc)


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
    base.outs = c("logbinge0","hs_util0")
    yr1.outs = c("logbinge1","hs_util1")
    yr2.outs = c("logbinge2","hs_util2")
    
    
    # code binary variables as factors
    mis.dat = mis.dat %>% 
      mutate_at(vars(starts_with("hs_util")),factor)
    
    # apply(mis.dat,2,function(x) sum(is.na(x))) #see how many missing
    
    ##########################################################################################################
    #Impute self-monitoring measures and tailoring outcome 
    ##########################################################################################################
    
    #Add missingness indicators for SM, FU1 and FU2 - needed for sensitivity analysis
    mis.dat$misSM <- ifelse(is.na(mis.dat$logsm_binge_last),1,0)
    mis.dat$misb1 <- ifelse(is.na(mis.dat$logbinge1),1,0) #FU1 logbinge 
    mis.dat$misb2 <- ifelse(is.na(mis.dat$logbinge2),1,0)
    
    #Creating predictor matrix
    sm.vars <- c("logsm_binge_last", dems,"logbinge0",'hs_util0',"A1")
    
    #Index of variables for sm imputation
    sm.i = which(names(mis.dat) %in% sm.vars)
    
    # shell imputation
    imp = mice(mis.dat[,sm.i],maxit=0)
    imp$method[7]="norm"
    # imp$loggedEvents
    
    # look at method
    #  imp$method
    
    # make predictor matrix
    pred = imp$predictorMatrix
    pred[,] = 0
    
    
    pred["logsm_binge_last", c(sm.vars[-1])] = 1
    #pred["sm_hid_last", c(sm.vars[-2])] = 1
    
    
    smimp = mice(mis.dat[,c(sm.i)],m=1,method = imp$method,predictorMatrix = pred,maxit=iter,printFlag = FALSE)
    sm.imp <- complete(smimp)
    
    #New dataset with imputed self-monitoring variables
    datnosm <- mis.dat %>% dplyr::select(names(mis.dat)[names(mis.dat) %nin% names(sm.imp)])
    sm.update <- cbind(sm.imp,datnosm)
    
    
    if (TV==T){ #CHANGE TO STANDARD DEVIATIONS
      #Add sensitivity parameter to self-monitoring measure 
      sm.update$sm_logbinge_last[sm.update$A1==-1] = sm.update$sm_logbinge_last[sm.update$A1==-1] + sm.update$misSM[sm.update$A1==-1]*(k1/2) #late 
      sm.update$sm_logbinge_last[sm.update$A1==1] = sm.update$sm_logbinge_last[sm.update$A1==1] + sm.update$misSM[sm.update$A1==1]*(k2/2) #early 
    }
    
    #################################################################################################
    #Update HD status and 2nd randomization
    #################################################################################################
    #Note: online coach = 1, email = -1; late=-1, early=1
    
    #sm.update$HD = ifelse(sm.update$sm_binge_last>1 | sm.update$sm_hid_last>0,1,0)  
    #sm.update$HD = ifelse(sm.update$logsm_binge_last>1,1,0)
    sm.update$logHD = ifelse(sm.update$logsm_binge_last< log(3),0,1)
    
    #sum(sm.update$HD[sm.update$misSM==1] != sm.update$HD2[sm.update$misSM==1]) #34/74 drinkers are now heavy drinkers (46%)
    
    #Updating A2 value for those missing sm_logbinge and who are now heavy drinkers - 1:1 randomization
    #sm.update$A2[sm.update$HD==1 & sm.update$misSM==1] = ifelse(rbinom(sum(sm.update$HD==1 & sm.update$misSM==1),1,.5)==0, -1, 1)
    
    sm.update$A2 = 0 
    #simdat$A2 = 0 
    sm.update$A2[which(sm.update$logHD==1)] = ifelse(rbinom(sum(sm.update$logHD==1),1,.5)==1, 1, -1)
    prop.table(table(sm.update$A2))
    
    
   sm.update$trajnochange = ifelse(sm.update$misSM==1 & sm.update$logHD==1,0,1) #if missing (light drinker) but after imputation they are a HD, then trt trajectory changed
    #trt traj no change for people who's sm were not missing
    
    ####################################################################################################
    # Follow-up 1 outcome imputation and update data with sensitivity parameter
    ####################################################################################################
    
    #Induce structural missingness: make FU1 outcomes missing if trt trajectory changed 
    sm.update[sm.update$trajnochange==0,names(sm.update) %in% yr1.outs] = NA
    sm.update[sm.update$trajnochange==0,names(sm.update) %in% yr2.outs] = NA
    
    #RIGHT NOW NO INTERACTION FOR SIMPLICITY - Create interaction term of A1:A2; need for congenielty with analysis model when comparing the 4 regimes 
   # sm.update$A1.A2 = sm.update$A1*sm.update$A2 
    
    
    #Using the MICE chained equations to impute all follow-up outcomes at once, so we can incorporate follow-up 2 information if available for people who didn't have their trajectory change
    
    #Creating predictor matrix
    FU.vars <- c(base.outs,yr1.outs, yr2.outs, dems,
                 "A1","A2",
                 "logHD")
    
    
    #Index of variables for FU imputation
    FU.i = which(names(sm.update) %in% FU.vars)
    
    # shell imputation
    imp = mice(sm.update[,FU.i],maxit=0,printFlag = FALSE)
    #  imp$loggedEvents
    
    # look at method
    #Change pmm to norm; can try pmm case later too
    imp$method[c(7,8)]="norm"
    meth=imp$method
    
    
    # make predictor matrix
    pred = imp$predictorMatrix
    pred[,] = 0
    
    # # code interactions
    # meth = imp$method
    # for (var in c(norms,intents)) {
    #   meth[paste0("API.",var)] = paste0("~I(API*",var,")")
    #   meth[paste0("API.A1.",var)] = paste0("~I(API.A1*",var,")")
    # }
    
    # Follow-up outcomes 
    pred["logbinge1",c("logbinge0",dems,
                    "A1","logHD","A2","hs_util1")] = 1
    
    pred["logbinge2",c("logbinge1","logbinge0",dems,"A1","logHD","A2","hs_util2")] = 1
    
    pred["hs_util1",c('hs_util0',dems,
                      "A1","logHD","A2","logbinge0")] = 1
    
    pred["hs_util2",c("hs_util1","hs_util0", dems,
                      "A1","logHD","A2","logbinge1")] = 1
    

    #Impute data
    fuimp = mice(sm.update[,FU.i],m=1,method = meth,predictorMatrix = pred,maxit=30,printFlag = FALSE)
    #  fuimp$loggedEvents
    
    fu.imp = complete(fuimp)
    
    #New dataset with imputed follow-up variables
    datnofu <- sm.update[,-FU.i]
    fu.update <- cbind(fu.imp,datnofu)
    
    #Add sensitivity parameter to follow-up 1
    fu.update$logbinge1[fu.update$A1==-1] = fu.update$logbinge1[fu.update$A1==-1] + fu.update$misb1[fu.update$A1==-1]*k1 #late 
    fu.update$logbinge1[fu.update$A1==1] = fu.update$logbinge1[fu.update$A1==1] + fu.update$misb1[fu.update$A1==1]*k2 #early 
    
    #Add sensitivity parameter to follow-up 2
    fu.update$logbinge2[fu.update$A1==-1] = fu.update$logbinge2[fu.update$A1==-1] + fu.update$misb2[fu.update$A1==-1]*k1 #late 
    fu.update$logbinge2[fu.update$A1==1] = fu.update$logbinge2[fu.update$A1==1] + fu.update$misb2[fu.update$A1==1]*k2 #early 
    
    
    #Add this imputed dataset to a long dataset of all m imputated datasets
    dat=fu.update
    dat$.imp = rep(i,dim(fu.update)[1]) #want length of imputed data
    imp_long = rbind(imp_long,dat)
  } 
  
  return(imp_long) #This imp_long is like Grace's mbridge_long
}





























###############################################
#COMPARE TO THE FULL DATA 


# Comparing my imputations to Grace's 

# library(ggplot2)
# library(esquisse)
# grace.dat = complete(grace.imp)
# df = data.frame(cbind(
#   c(grace.dat$logbinge1[grace.dat$API==1],
#                       dat$logbinge1,
#     grace.dat$logbinge2[grace.dat$API==1],
#     dat$logbinge2),
#   rep(c(rep("Grace imputation",591),rep("My imputation",591)),2),
#   c(rep("logbinge1",591*2),rep("logbinge2",591*2))))
# 
# names(df) = c("logbinge","Imputation","Follow-up")
# esquisser(df)
# 
# ggplot(df) +
#   aes(x = logbinge) +
#   geom_bar(fill = "#112446") +
#   theme_gray() +
#   facet_grid(vars(Imputation), vars(`Follow-up`))
# 
# plot(density(grace.dat$logbinge1[grace.dat$API==1]),
#      main="Grace's logbinge1 imputation")
# 
# plot(density(dat$logbinge1),
#      main="My logbinge1 imputation")
# 
# plot(hist(grace.dat$logbinge1[grace.dat$API==1]),
#      main="Grace's logbinge1 imputation")
# 
# plot(hist(dat$logbinge1),
#      main="My logbinge1 imputation")


# check for convergence with 35 more iterations
# imp2 = mice.mids(imp, maxit=35)
# plot(imp2) # pretty good mixing

# save results
#save(mbridge,imp,file="~/Desktop/MBridge/mi-mbridge.Rdata")
#save(imp_long,file="imp_longSA_k1_1_k2_.5.rda")
