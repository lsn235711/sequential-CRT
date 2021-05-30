library(tidyverse)

setting0 = c(1,2,3,4)
x_model = "hmm" #ar or hmm

#knockoff_plus0 = 1
include_h0 = FALSE
c0 = 0.1
N = 500
NN = 5
NNN = 500

upper_quan = function(x){
  return (mean(x) + qnorm(0.95)*sd(x)/sqrt(NNN))
  #return (mean(x))
}
lower_quan = function(x){
  return (mean(x) - qnorm(0.95)*sd(x)/sqrt(NNN))
  #return (mean(x))
}

m = 18
M = 18*N*NN
results = data.frame(fdp = numeric(0), power = numeric(0),
                     include_h = logical(0), c = character(0), 
                     knockoff_plus = numeric(0), method = character(0),       
                     blackbox = character(0), time = numeric(0),
                     one_shot = logical(0), amp =numeric(0),        
                     setting_no = integer(0) 
                     )
results[,6] = as.character(results[,6])
results[,7] = as.character(results[,7])

for (arg in 1:N){
    print(arg)
    for (iii in 1:NN){
        filename = sprintf("/Users/lsn/Desktop/Can/project_CRT_seqstep/simulations/simulation_20210520_one_shot_CRT/sherlock/output_%s/output_%i_%i.Rdata",x_model, arg, iii)

        if (file.exists(filename)) {
            load(filename)
            num = dim(fdp_power)[1]
            fdp_power["amp"] = rep(amplitude,num)
            fdp_power["setting_no"] = rep(setting_no,num)
            ## take a subset of results
            fdp_power[is.na(fdp_power)] <- "NA"
            fdp_power = subset(fdp_power, (c == c0 | c == "NA")&
                     (include_h == include_h0 | c == "NA")&
                     (knockoff_plus == 1)
            )
            fdp_power[,6] = as.character(fdp_power[,6])
            fdp_power[,7] = as.character(fdp_power[,7])
            results = rbind(results, fdp_power)
        }
    }
}

results[is.na(results)] <- "NA"


if(x_model == "ar"){
  amp_multi = c(1.5,15,10,20)
} else{
  amp_multi = c(2,8,10,20)
}



k_names = c(20,20,20,20)
setting_names = c("Linear Gaussian Response", "Non-linear Gaussian Response", "Logistic Binary Response","Non-Linear Binary Response")

dat0 = results

library(dplyr)

dat = dat0 %>% group_by(amp, c, method, blackbox,include_h,knockoff_plus,setting_no, one_shot) %>%
          summarise(FDR = mean(fdp), Power = mean(power),
                    FDR_upper = upper_quan(fdp), power_upper = upper_quan(power),
                    FDR_lower = lower_quan(fdp), power_lower = lower_quan(power),
                    time = mean(time)
                    ) 
dat$c = factor(dat$c, levels = c("NA", 0.1,0.2,0.3,0.4,0.5))
dat$amp = dat$amp* amp_multi[dat$setting_no]
dat$setting = factor(setting_names[dat$setting_no], levels = setting_names[c(1,3,2,4)])

blackbox_transfer = function(a){
  if (a == "rf") {return ("Random Forest")}
  else if (a == "gb") {return ("Gradient Boosting")}
  else if (a == "lasso") {return ("Lasso/glmnet")}
}

method_transfer = function(a){
  if (a == "CRT_seqstep_exact") {return ("CRT_seqstep+ (exact)")}
  else if (a == "CRT_seqstep_inexact") {return ("CRT_seqstep+ (inexact)")}
  else if (a == "knockoffs") {return ("knockoffs")}
}

dat$Blackbox = factor(unlist(lapply(dat$blackbox, blackbox_transfer)), levels = c("Lasso/glmnet", "Random Forest","Gradient Boosting"))
dat$Method = factor(unlist(lapply(dat$method, method_transfer)), levels = c("knockoffs", "CRT_seqstep+ (inexact)","CRT_seqstep+ (exact)"))

dat_toplot = subset(dat, (c == c0 | c == "NA")&
                          (include_h == include_h0 | c == "NA")&
                          #(knockoff_plus == 1 | c == "NA")&
                          (knockoff_plus == 1)&
                          (setting_no %in% setting0[c(1,2)])&
                         (method == "CRT_seqstep_inexact")
                    )

dat_toplot1 = dat_toplot %>% group_by(method, blackbox, one_shot, setting) %>%
                         summarise(average_time = mean(time))




