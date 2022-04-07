library(tidyverse)

plot_FDR = T
plot_power = T
plot_band = F
setting0 = c(1,2,3,4)
x_model = "hmm" #ar or hmm

#knockoff_plus0 = 1
include_h0 = FALSE
c0 = 0.1

#########################################
############ Can change this ############
#########################################
N = 400
NN = 5
NNN = 100

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

track_numbers = rep(0,N)

for (arg in 1:N){
    print(arg)
    for (iii in 1:NN){
        filename = sprintf("output_%s/output_%i_%i.Rdata",x_model, arg, iii)

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
            track_numbers[arg] = iii
        }
    }
}
save_filename = sprintf("outputs/result_%s_c%.1f.Rdata", x_model, c0)
save(results, file = save_filename)
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
                          (setting_no %in% setting0)&
                         (method == "sequential_CRT_sym_stats" | 
                            (method == "knockoffs" & one_shot == TRUE)
                            )
                    )
combine_method_name = function(a,b){
  if(a == "sequential_CRT_sym_stats" & b == TRUE){
    return("sequential CRT (one−shot CRT)")
  } else if(a == "sequential_CRT_sym_stats" & b == FALSE){
    return("sequential CRT (original CRT)")
  } else{
    return("knockoffs")
  }
}

combine_method_name_vec = function(a,b){
  mapply(combine_method_name, a ,b)
}

dat_toplot = dat_toplot %>% mutate(`Variable Selection Method` = combine_method_name_vec(method, one_shot) )
dat_toplot["Variable Selection Method"] = factor(unlist(dat_toplot["Variable Selection Method"]), levels = c("sequential CRT (one−shot CRT)", "sequential CRT (original CRT)", "knockoffs"))

dat_toplot = subset(dat_toplot, (setting_no == 1 & Blackbox == "Lasso/glmnet") |
                        (setting_no == 3 & Blackbox == "Lasso/glmnet") |
                        (setting_no == 2 & Blackbox == "Random Forest") |
                        (setting_no == 4 & Blackbox == "Random Forest") |
                        (setting_no == 2 & Blackbox == "Gradient Boosting") |
                        (setting_no == 4 & Blackbox == "Gradient Boosting")
                      )



library(ggplot2)

## Colors
library(scales)
color_list = hue_pal()(4)

p_FDR = ggplot(dat_toplot, aes(x=amp, y=FDR,
                               group = interaction(`Variable Selection Method`, Blackbox),
                               linetype = Blackbox,
                               color = `Variable Selection Method`,
                               shape = Blackbox
)) +
  geom_line() + geom_point() + 
  theme_bw() + 
  theme(legend.position="none", legend.box="vertical",
        legend.margin=ggplot2::margin(c(0,0,0,0)),
        legend.key.width = grid::unit(2, "lines"))+
  geom_hline(aes(yintercept = 0.1)) +
  scale_color_manual(values = color_list[c(1,4,2)]) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted")) + 
  facet_grid(~ setting, scales = "free") + 
  xlab("") +
  ylim(c(0,0.15))

p_power = ggplot(dat_toplot, aes(x=amp, y= Power,
                                 group = interaction(`Variable Selection Method`, Blackbox),
                                 linetype = Blackbox,
                                 color = `Variable Selection Method`,
                                 shape = Blackbox
)) +
  geom_line() + geom_point() + 
  theme_bw() + 
  theme(legend.position="bottom", legend.box="vertical",
        legend.margin=ggplot2::margin(c(0,0,0,0)),
        legend.key.width = grid::unit(2, "lines"))+
  scale_color_manual(values = color_list[c(1,4,2)]) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted")) + 
  facet_grid(~ setting, scales = "free") + 
  xlab("Coefficient Amplitude") +
  ylim(c(0,1))


library(gridExtra)
pp <- grid.arrange(p_FDR, p_power, nrow=2, heights=c(1,1.5))
pp
pp %>% ggsave(file=sprintf("small_experiment_%s.pdf", x_model), width=8.2, height=4.5, units="in")


