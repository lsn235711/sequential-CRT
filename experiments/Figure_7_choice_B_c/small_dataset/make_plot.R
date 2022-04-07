library(tidyverse)

setting_names = c("Linear Gaussian Response", "Non-linear Gaussian Response", "Logistic Binary Response","Non-Linear Binary Response")
results = NULL
for (B in c(9,19,29,39,49)){
  for (setting1 in c(1,3)){
    for (seed1 in 1:100){
      for (amp in 1:2){
        filename = sprintf("output/output_B%i_setting%i_seed%i_amp%i.Rdata", B, setting1, seed1, amp)
        
        if (file.exists(filename)) {
          load(filename)
          num = dim(fdp_power)[1]
          fdp_power["B"] = rep(B_choice,num)
          fdp_power["setting_no"] = rep(setting_no,num)
          fdp_power["Amplitude"] = rep(amp,num)
          ## take a subset of results
          fdp_power[is.na(fdp_power)] <- "NA"
          results = rbind(results, fdp_power)
        }
      }
    }
  }
}

dat = results %>% group_by(B, method, setting_no, c, Amplitude) %>%
  summarise(FDR = mean(fdp), Power = mean(power))
dat$setting = factor(setting_names[dat$setting_no], levels = setting_names[c(1,3,2,4)])


method_transfer = function(a){
  if (a == "sequential_CRT_split") {return ("sequential CRT (split)")}
  else if (a == "sequential_CRT_sym_stats") {return ("sequential CRT (symmetric statistic)")}
  else if (a == "knockoffs") {return ("knockoffs")}
}

amp_transfer = function(a){
  if (a == 1) {return ("Low SNR")}
  else {return ("High SNR")}
}

dat$Method = factor(unlist(lapply(dat$method, method_transfer)), levels = c("knockoffs", "sequential CRT (symmetric statistic)","sequential CRT (split)"))
dat$Amplitude = factor(unlist(lapply(dat$Amplitude, amp_transfer)), levels = c("Low SNR", "High SNR"))

dat$c = factor(dat$c)
dat$B = factor(dat$B)

dat_toplot = dat[,c("Method", "FDR", "Power", "setting", "B", "c","Amplitude")]
dat_toplot = subset(dat_toplot, Method == "sequential CRT (symmetric statistic)")

#########################################
########## Power & FDR Plots ############
#########################################

library(ggplot2)


## Colors
library(scales)
color_list = c(hue_pal()(4), "yellow3","#6A3D9A", "brown")

p_FDR = ggplot(dat_toplot, aes(x=c, y=FDR,
                               group = B,
                               color = B
)) +
  geom_line() + geom_point() +
  theme_bw() +
  #theme(legend.position="none")+
  geom_hline(aes(yintercept = 0.1)) +
  #scale_color_manual(values = color_list[c(2,1,3,4,5,6,7)]) +
  #scale_linetype_manual(values=c("dashed", "solid", "dotted")) +
  facet_grid(Amplitude ~ setting, scales = "free") +
  xlab("") #+
  #ylim(c(0,1))

p_power = ggplot(dat_toplot, aes(x=c, y= Power,
                                 group = B,
                                 color = B
)) +
  geom_line() + geom_point() + 
  theme_bw() + 
  #theme(legend.position="bottom")+
  #guides(fill=guide_legend(nrow=2,byrow=TRUE)) + 
  #scale_color_manual(values = color_list[c(2,1,3,4,5,6,7)]) +
  #scale_linetype_manual(values=c("dashed", "solid", "dotted")) + 
  facet_grid(Amplitude ~ setting, scales = "free") + 
  xlab("c") #+
  #ylim(c(0,1))

p_power

# library(gridExtra)
# pp <- grid.arrange(p_FDR, p_power, nrow=2, heights=c(1,1.25))
# pp
#pp %>% ggsave(file="experiment_300.pdf", width=7, height=4, units="in")

#p_power %>% ggsave(file="experiment_300_Bc.pdf", width=7, height=4, units="in")
p_FDR %>% ggsave(file="experiment_300_Bc_FDR.pdf", width=7, height=3, units="in")

