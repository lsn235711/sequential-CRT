rho = seq(0,0.5,by = 0.001)
p = 300
q = 0.1
c = 0.1
k = length(rho)
q_choices = c(0.05,0.1,0.15,0.2)


inflation_bound = NULL
qs = NULL
rhos = NULL

for (q in q_choices){
    beta = (c+(1-c)*q)/(1-c)/(1-q)
    delta = rho*(c*(1-q)+q)/c/(1-q)
    inflation_bound_one = delta/(1+beta*delta)*(c/(1-c) - (c-c*delta)/(1-(c-c*delta))*q)
    inflation_bound_one = pmin(inflation_bound_one, c*(1-q))
    
    inflation_bound = c(inflation_bound, inflation_bound_one)
    qs = c(qs, rep(q,k))
    rhos = c(rhos, rho)
}

library(ggplot2)
library(latex2exp)
data = data.frame(rho = rhos, inflation_bound = inflation_bound, q = qs)
data$q = factor(data$q)

g1 = ggplot(data, aes(x = rho, y = inflation_bound, color = q, group = q))+
    geom_line()+
    xlab(TeX("$\\rho$")) +
    ylab(TeX("FDR inflation $\\epsilon(c, q, \\rho)$")) +
    theme_bw() +
    #facet_grid( ~ bound, scales = "free", labeller = label_both) +
    scale_y_continuous(breaks = seq(0, 0.1, by = 0.02), limit = c(0,0.1)) + 
    scale_color_brewer(palette="Set1") +
    scale_linetype_manual(values=c("dotted", "solid"))
g1 
ggsave("bound.pdf",g1, width=5, height=2.7, units="in")




