#This script defines functions to induce missingness 

# Expit Function 
expit=function(x){
  exp(x)/(1+exp(x))
}

#Function to calculate missingness for a time-point based on missingness mechanism (MCAR, MAR or MNAR)
induce.missYt = function(nperson,t,yt,wide.df,alpha_m=NULL, beta_y=NULL, beta_by=NULL,beta_hs=NULL,
                         beta_maxd=NULL,beta_a1 = NULL,beta_f = NULL,
                         beta_hd=NULL,beta_a2=NULL,beta_a2binge2= NULL,
                         beta_a1binge2=NULL,beta_smb=NULL,beta_smh=NULL,beta_a1smb=NULL,
                         miss_type, miss_prob =NULL){
  
  rand.prob = runif(nperson,min=0,max=1) #vector of uniform distributed values for each person 
  
  if (miss_type=="MCAR"){yt.missprob = rep(miss_prob,nperson)} 
  
  
  if(miss_type=="MAR"){ 
    if (t<1){ #tailoring variable

      
      yt.missprob = expit(alpha_m+beta_y*(wide.df$binge0-mean(wide.df$binge0))+
                            beta_by*(wide.df$byaacq0-mean(wide.df$byaacq0))+
                            #beta_hs*(wide.df$hs_util0-mean(wide.df$hs_util0))+
                            beta_maxd*(wide.df$max_drinks0-mean(wide.df$max_drinks0)) + 
                            beta_a1*(wide.df$A1-mean(wide.df$A1))+
                            beta_f*(wide.df$female-mean(wide.df$female)))
      
    } else # follow-up times
      yt.missprob = expit(alpha_m+
                            beta_y*(wide.df$binge0-mean(wide.df$binge0))*(t==1)+
                            beta_by*(wide.df$byaacq0 - mean(wide.df$byaacq0))*(t==1)+
                            #beta_hs*(wide.df$hs_util0-mean(wide.df$hs.util0))+
                            beta_maxd*(wide.df$max_drinks0- mean(wide.df$max_drinks0)) + 
                            beta_a1*(wide.df$A1-mean(wide.df$A1))+
                            beta_a2*(wide.df$A2-mean(wide.df$A2))+ 
                            beta_hd*(wide.df$HD-mean(wide.df$HD))+
                            beta_by*(wide.df$byaacq1-mean(wide.df$byaacq1))*(t==2) + 
                            beta_y*(wide.df$binge1-mean(wide.df$binge1))*(t==2))
  }
  if(miss_type=="MNAR"){ 
    if (t<1){ #tailoring variable 
      yt.missprob = expit(alpha_m+(beta_y/2)*(wide.df$sm_binge_last-mean(wide.df$sm_binge_last))+
                            #beta_smh*(wide.df$sm_hid_last)+ 
                            beta_a1*(wide.df$A1-mean(wide.df$A1)) + beta_a1smb*((wide.df$A1-mean(wide.df$A1))*(wide.df$sm_binge_last-mean(wide.df$sm_binge_last)))) #interaction with trt and sm value
      
    } else # follow-up times
      yt.missprob = expit(alpha_m+beta_a1*(wide.df$A1-mean(wide.df$A1)) + beta_a2*(wide.df$A2-mean(wide.df$A2)) +
                            beta_y*(wide.df$binge2-mean(wide.df$binge2))+ 
                            #beta_by*wide.df$byaacq1*(t==1)+beta_by*wide.df$byaacq2*(t==2)+
                            beta_a1binge2*((wide.df$A1-mean(wide.df$A1))*(wide.df$binge2-mean(wide.df$binge2)))+ 
                            beta_a2binge2*((wide.df$A2-mean(wide.df$A2))*(wide.df$binge2-mean(wide.df$binge2))))
  }
  yt[rand.prob < yt.missprob] = NA
  return(yt)
}



