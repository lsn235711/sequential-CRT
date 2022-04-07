library(tidyverse)

plot_FDR = T
plot_power = T
plot_band = F
setting0 = c(1,2,3,4)
x_model = "ar" #ar or hmm

#knockoff_plus0 = 1
include_h0 = FALSE
c0 = 0.1

upper_quan = function(x){
  return (mean(x) + qnorm(0.95)*sd(x)/sqrt(NNN))
  #return (mean(x))
}
lower_quan = function(x){
  return (mean(x) - qnorm(0.95)*sd(x)/sqrt(NNN))
  #return (mean(x))
}

save_filename = sprintf("outputs/result_%s_c%.1f.Rdata", x_model, c0)
load(save_filename)

results[is.na(results)] <- "NA"


if(x_model == "ar"){
  amp_multi = c(1.5,15,8,20)
} else{
  amp_multi = c(2,6,10,20)
}



k_names = c(20,20,20,20)
setting_names = c("Linear Gaussian Response", "Non-linear Gaussian Response", "Logistic Binary Response","Non-Linear Binary Response")

dat0 = results

library(dplyr)

dat = dat0 %>% group_by(amp, c, method, blackbox,include_h,knockoff_plus,setting_no, one_shot) %>%
          summarise(FDR = mean(fdp), Power = mean(power),
                    #FDR_upper = upper_quan(fdp), power_upper = upper_quan(power),
                    #FDR_lower = lower_quan(fdp), power_lower = lower_quan(power),
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
  if (a == "sequential_CRT_exact") {return ("sequential CRT (split)")}
  else if (a == "sequential_CRT_inexact") {return ("sequential CRT (symmetric statistics)")}
  else if (a == "knockoffs") {return ("knockoffs")}
}

dat$Blackbox = factor(unlist(lapply(dat$blackbox, blackbox_transfer)), levels = c("Lasso/glmnet", "Random Forest","Gradient Boosting"))
dat$Method = factor(unlist(lapply(dat$method, method_transfer)), levels = c("knockoffs", "sequential CRT (symmetric statistics)","sequential CRT (split)"))


#########################################
########## Power & FDR Plots ############
#########################################


dat_toplot = subset(dat, (c == c0 | c == "NA")&
                          (include_h == include_h0 | c == "NA")&
                          #(knockoff_plus == 1 | c == "NA")&
                          (knockoff_plus == 1)&
                          #(setting_no %in% setting0[c(1,2)])&
                         (method == "sequential_CRT_sym_stats")
                    )

dat_toplot1 = dat_toplot %>% group_by(method, blackbox, one_shot, setting) %>%
                         summarise(average_time = mean(time))
dat_toplot1



