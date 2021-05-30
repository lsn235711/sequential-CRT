set.seed(4)
m = 100
q = 0.1
c = 0.2
u1 = 0.2
u0 = 1
p_vals = c(u1*runif(m/2), runif(m/2))
type = rep(c("nonnulls", "nulls"), each = m/2)
position = 1:m

k_good = max(which((1 + cumsum(p_vals >= c))/pmax(1, cumsum(p_vals <= c)) <= (1-c)/c*q))
rejected_yes_no = (p_vals<= c)*((1:m) <= k_good)
rejected = c("rejected", "not rejected")[(p_vals<= c)*(1:m <= k_good) + 1]

type_01 = rep(c(1,2), each = m/2)
selected_type = (c("True Discoveries", "False Discoveries")[type_01])[which(as.logical(rejected_yes_no))]

##### good ordering
library(ggplot2)

library(grid)
k_hat_text <- textGrob(expression(hat(k)), gp=gpar(fontsize=13, fontface="bold"))
c_text <- textGrob("c", gp=gpar(fontsize=13, fontface="bold"))

global_size = 15
data = data.frame(p_vals = p_vals, type = type, position = 1:m)
data2 = data.frame(position = which(as.logical(rejected_yes_no)), dummy = factor(selected_type, levels = c("True Discoveries", "False Discoveries")))
p_good = ggplot(data=data, aes(x=position, y=p_vals)) +
    geom_bar(stat="identity",aes(fill = type)) +
    scale_fill_manual(values = c("lightcoral", "cornflowerblue")) + 
    geom_point(aes(x=position, y=-0.05, color = dummy,shape = dummy), data = data2, size = 2) +
    scale_colour_manual(name = "Selected Set",
                        values = c("black", "black"), 
                        labels = c("True Discoveries", "False Discoveries")) + 
    scale_shape_manual(name = "Selected Set",
                       values = c(20, 7), 
                       labels = c("True Discoveries", "False Discoveries")) + 
    geom_vline(xintercept = k_good + 0.5)+
    geom_segment(aes(x = -Inf, y = c, xend = k_good + 0.5, yend = c))+
    theme_bw() +
    ##legend
    guides(fill=guide_legend(title=""), color = guide_legend(title="Selected Set"))+
    ##labs
    labs(x = "", y = "p-values") + 
    scale_x_continuous(breaks = round(seq(0, 100, by = 20),1)) + 
    annotation_custom(c_text,xmin=-7,xmax=-7,ymin=c,ymax=c) + 
    annotation_custom(k_hat_text,xmin=k_good + 0.5,xmax=k_good + 0.5,ymin=-0.18,ymax=-0.18) + 
    coord_cartesian(clip = "off") + 
    theme(text = element_text(size=global_size))

p_good

##### bad ordering

bad_ordering = rank(rep(c(0,0.8), each = 50) + rnorm(m))
p_vals_bad = p_vals[bad_ordering]
k_bad = max(which((1 + cumsum(p_vals_bad >= c))/pmax(1, cumsum(p_vals_bad <= c)) <= (1-c)/c*q))
rejected_yes_no_bad = (p_vals_bad<= c)*((1:m) <= k_bad)

type_01 = type_01[bad_ordering]
selected_type_bad = (c("True Discoveries", "False Discoveries")[type_01])[which(as.logical(rejected_yes_no_bad))]

data = data.frame(p_vals = p_vals_bad, type = type[bad_ordering], position = (1:m))
data2 = data.frame(position = which(as.logical(rejected_yes_no_bad)), dummy = factor(selected_type_bad, levels = c("True Discoveries", "False Discoveries")))
p_bad = ggplot(data=data, aes(x=position, y=p_vals)) +
    geom_bar(stat="identity",aes(fill = type)) +
    scale_fill_manual(values = c("lightcoral", "cornflowerblue")) + 
    geom_point(aes(x=position, y=-0.05, color = dummy, shape = dummy), data = data2, size = 2) +
    scale_colour_manual(name = "Selected Set",
                        values = c("black", "black"), 
                        labels = c("True Discoveries", "False Discoveries")) + 
    scale_shape_manual(name = "Selected Set",
                        values = c(20, 7), 
                        labels = c("True Discoveries", "False Discoveries")) + 
    geom_vline(xintercept = k_bad + 0.5)+
    geom_segment(aes(x = -Inf, y = c, xend = k_bad + 0.5, yend = c))+
    theme_bw() +
    ##legend
    guides(fill=guide_legend(title=""))+
    ##labs
    labs(x = "", y = "p-values") +
    scale_x_continuous(breaks = round(seq(0, 100, by = 20),1)) +
    annotation_custom(c_text,xmin=-7,xmax=-7,ymin=c,ymax=c) +
    annotation_custom(k_hat_text,xmin=k_bad + 0.5,xmax=k_bad + 0.5,ymin=-0.18,ymax=-0.18) +
    coord_cartesian(clip = "off") + 
    theme(text = element_text(size=global_size))

p_bad

library(tidyverse)
p_good %>% ggsave(file="good_ordering.pdf", width=11, height=3.5, units="in")
p_bad %>% ggsave(file="bad_ordering.pdf", width=11, height=3.5, units="in")

