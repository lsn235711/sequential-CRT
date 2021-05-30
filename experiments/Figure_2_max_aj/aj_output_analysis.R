library(ggplot2)
library(latex2exp)
library(dplyr)
library(tidyverse)
library(gridExtra)

########################
##  Correlation
########################
max_ajs = NULL
for (seed in 1:500){
  filename = sprintf("cor/output/output%i.Rdata",seed)
  if (file.exists(filename)) {
    load(filename)
    max_ajs = c(max_ajs, max(ajs_mean))
  }
}


q = 0.1
c = 0.3

c_plus_delta = 0.39
delta = max(max_ajs[which(max_ajs<=c_plus_delta)]) - 0.3
epsilon = mean(max_ajs>c_plus_delta)
binwidth = 0.001

## make plot
dt = data.frame(max_aj = max_ajs)

font.size = 15

p1 = ggplot(dt, aes(x=max_aj)) + 
  geom_histogram(binwidth = 0.001, fill="grey60",color="black") + 
  theme_bw() + 
  labs(x = TeX('$max_j(a_j)$')) + 
  geom_vline(aes(xintercept= c), color="blue", linetype="dashed", size=1) + 
  geom_vline(aes(xintercept= c + delta + binwidth/2), color="red", linetype="dashed", size=1) +
  theme(axis.text.x = element_text(size = font.size),
        axis.text.y = element_text(size = font.size),  
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size)) +
  geom_segment(aes(x = c + delta/3, y = 22, xend = c, yend = 22),
                 arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = c + 2*delta/3, y = 22, xend = c + delta, yend = 22),
               arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("text", x= c + delta/2, y=22, label="delta", parse=TRUE, size = 6) + 
  annotate("text", x= c + 0.008, y=27, label="c", parse=TRUE, size = 6, color = "blue") + 
  annotate("text", x= c + delta + 0.014, y=27.3, label="c + delta", parse=TRUE, size = 6, color = "red") +
  scale_x_continuous(breaks = seq(0.24, 0.44, by = 0.02), limits = c(0.24,0.44)) +
  scale_y_continuous(breaks = seq(0, 30, by = 5), limits = c(0,30)) +
  ggtitle("(a) Test statistic: absolute value of correlation")

FDR_bound_cor = q*(c+delta)/(c)*(1-c)/(1-c-delta) + epsilon

########################
##  OLS
########################

max_ajs = NULL
for (seed in 1:500){
  filename = sprintf("ols/output_ols/output%i.Rdata",seed)
  if (file.exists(filename)) {
    load(filename)
    max_ajs = c(max_ajs, max(ajs_mean))
  }
}


q = 0.1
c = 0.3

c_plus_delta2 = 0.41
delta2 = max(max_ajs[which(max_ajs<=c_plus_delta2)]) - 0.3
epsilon2 = mean(max_ajs>c_plus_delta2)
binwidth = 0.001

## make plot
dt = data.frame(max_aj = max_ajs)

font.size = 15

p2 = ggplot(dt, aes(x=max_aj)) + 
  geom_histogram(binwidth = 0.001, fill="grey60",color="black") + 
  theme_bw() + 
  labs(x = TeX('$max_j(a_j)$')) + 
  geom_vline(aes(xintercept=c), color="blue", linetype="dashed", size=1) + 
  geom_vline(aes(xintercept=c + delta2 + binwidth/2), color="red", linetype="dashed", size=1) +
  theme(axis.text.x = element_text(size = font.size),
        axis.text.y = element_text(size = font.size),  
        axis.title.x = element_text(size = font.size),
        axis.title.y = element_text(size = font.size)) +
  geom_segment(aes(x = c + delta2/3, y = 47, xend = c, yend = 47),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = c + 2*delta2/3, y = 47, xend = c + delta2, yend = 47),
               arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("text", x= c + delta2/2, y=47, label="delta", parse=TRUE, size = 6) + 
  annotate("text", x= c + 0.008, y=55, label="c", parse=TRUE, size = 6, color = "blue") + 
  annotate("text", x= c + delta2 + 0.014, y=55.3, label="c + delta", parse=TRUE, size = 6, color = "red") +
  scale_x_continuous(breaks = seq(0.24, 0.44, by = 0.02), limits = c(0.24,0.44))+
  scale_y_continuous(breaks = seq(0, 60, by = 10), limits = c(0,60))+
  ggtitle("(b) Test statistics: absolute value of regression coefficient")

FDR_bound_ols = q*(c+delta2)/(c)*(1-c)/(1-c-delta2) + epsilon2

pp <- grid.arrange(p1, p2, nrow=2)
pp
pp %>% ggsave(file="hist_max_aj.pdf", width=8, height=6, units="in")


FDR_bound_cor
FDR_bound_ols