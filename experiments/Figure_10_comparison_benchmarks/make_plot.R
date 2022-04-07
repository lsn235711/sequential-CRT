library(tidyverse)

setting0 = c(1,2,3,4)
x_model = "ar"
include_h0 = FALSE
c0 = 0.1

#########################################
###### getting data for CRT_seqstep #####
#########################################
results = data.frame(fdp = numeric(0), power = numeric(0),
                     include_h = logical(0), c = character(0), 
                     knockoff_plus = numeric(0), method = character(0),       
                     blackbox = character(0), time = numeric(0),
                     one_shot = logical(0), amp =numeric(0),        
                     setting_no = integer(0) 
                     )
results[,6] = as.character(results[,6])
results[,7] = as.character(results[,7])


for (amp1 in 1:5){
  print(amp1)
    for (setting1 in 1:4){
      print(setting1)
      for (seed1 in 1:100){
        filename = sprintf("Figure_5_FDR_power_large/output_%s/output_amp%i_setting%i_seed%i_i%i.Rdata", x_model, amp1, setting1, seed1, 1)
  
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
}
results[is.na(results)] <- "NA"

if(x_model == "ar"){
  amp_multi = c(1.5,15,6,20)
} else{
  amp_multi = c(1.6,8,6,20)
}


k_names = c(20,20,20,20)
setting_names = c("Linear Gaussian Response", "Non-linear Gaussian Response", "Logistic Binary Response","Non-Linear Binary Response")

dat0 = results

library(dplyr)

dat = dat0 %>% group_by(amp, c, method, blackbox,include_h,knockoff_plus,setting_no, one_shot) %>%
          summarise(FDR = mean(fdp), Power = mean(power),
                    time = mean(time)) 
dat$c = factor(dat$c, levels = c("NA", 0.1,0.2,0.3,0.4,0.5))
dat$amp = dat$amp* amp_multi[dat$setting_no]
dat$setting = factor(setting_names[dat$setting_no], levels = setting_names[c(1,3,2,4)])

blackbox_transfer = function(a){
  if (a == "rf") {return ("Random Forest")}
  else if (a == "gb") {return ("Gradient Boosting")}
  else if (a == "lasso") {return ("Lasso/glmnet")}
}

method_transfer = function(a){
  if (a == "sequential_CRT_split") {return ("CRT_seqstep+ (exact)")}
  else if (a == "sequential_CRT_sym_stats") {return ("sequential CRT (symmetric statistic)")}
  else if (a == "knockoffs") {return ("knockoffs")}
}

dat$Blackbox = factor(unlist(lapply(dat$blackbox, blackbox_transfer)), levels = c("Lasso/glmnet", "Random Forest","Gradient Boosting"))
dat$Method = factor(unlist(lapply(dat$method, method_transfer)), levels = c("knockoffs", "sequential CRT (symmetric statistic)","sequential CRT (split)"))

dat_toplot = subset(dat, (c == c0 | c == "NA")&
                      (include_h == include_h0 | c == "NA")&
                      #(knockoff_plus == 1 | c == "NA")&
                      (knockoff_plus == 1)&
                      (setting_no %in% setting0)
)


dat_toplot = subset(dat_toplot, (setting_no == 1 & Blackbox == "Lasso/glmnet") |
                      (setting_no == 3 & Blackbox == "Lasso/glmnet") |
                      (setting_no == 2 & Blackbox == "Random Forest") |
                      (setting_no == 4 & Blackbox == "Random Forest") |
                      (setting_no == 2 & Blackbox == "Gradient Boosting") |
                      (setting_no == 4 & Blackbox == "Gradient Boosting")
)

dat_toplot11 = subset(dat_toplot, setting_no == 1 | setting_no == 3)
dat_toplot1 = dat_toplot11[,c("Method", "FDR", "Power", "setting", "amp")]


#########################################
#### getting data for other methods #####
#########################################
results = NULL
for (amp1 in 1:5){
  print(amp1)
  for (setting1 in c(1,3)){
    for (seed1 in 1:100){
      filename = sprintf("output/output_amp%i_setting%i_seed%i.Rdata", amp1, setting1, seed1)
      
      if (file.exists(filename)) {
        load(filename)
        num = dim(fdp_power)[1]
        fdp_power["amp"] = rep(amplitude,num)
        fdp_power["setting_no"] = rep(setting_no,num)
        ## take a subset of results
        fdp_power[is.na(fdp_power)] <- "NA"
        results = rbind(results, fdp_power)
      }
    }
  }
}

dat2 = results %>% group_by(amp, method, setting_no) %>%
  summarise(FDR = mean(fdp), Power = mean(power))
dat2$amp = dat2$amp* amp_multi[dat2$setting_no]
dat2$setting = factor(setting_names[dat2$setting_no], levels = setting_names[c(1,3,2,4)])

method_transfer2 = function(a){
  if (a == "GM") {return ("Gaussian Mirror")}
  else if (a == "d0CRT") {return ("d0CRT + BH")}
  else if (a == "dICRT") {return ("dICRT + BH")}
  else if (a == "HRT") {return ("HRT")}
}

dat2$Method = factor(unlist(lapply(dat2$method, method_transfer2)), levels = c("Gaussian Mirror", "d0CRT + BH","dICRT + BH", "HRT"))

dat_toplot2 = dat2[,c("Method", "FDR", "Power", "setting", "amp")]

dat_toplot = rbind(dat_toplot1, dat_toplot2)
#########################################
########## Power & FDR Plots ############
#########################################





library(ggplot2)


## Colors
library(scales)
color_list = c(hue_pal()(4), "yellow3","#6A3D9A", "brown")

p_FDR = ggplot(dat_toplot, aes(x=amp, y=FDR,
                               group = Method,
                               color = Method
)) +
  geom_line() + geom_point() + 
  theme_bw() + 
  theme(legend.position="none")+
  geom_hline(aes(yintercept = 0.1)) +
  scale_color_manual(values = color_list[c(2,1,3,4,5,6,7)]) +
  #scale_linetype_manual(values=c("dashed", "solid", "dotted")) + 
  facet_grid(~ setting, scales = "free") + 
  xlab("") +
  ylim(c(0,1))

p_power = ggplot(dat_toplot, aes(x=amp, y= Power,
                                 group = Method,
                                 color = Method
)) +
  geom_line() + geom_point() + 
  theme_bw() + 
  theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) + 
  scale_color_manual(values = color_list[c(2,1,3,4,5,6,7)]) +
  #scale_linetype_manual(values=c("dashed", "solid", "dotted")) + 
  facet_grid(~ setting, scales = "free") + 
  xlab("Coefficient Amplitude") +
  ylim(c(0,1))


library(gridExtra)
pp <- grid.arrange(p_FDR, p_power, nrow=2, heights=c(1,1.25))
pp
pp %>% ggsave(file=sprintf("large_comparison_%s.pdf", x_model), width=7, height=7, units="in")

