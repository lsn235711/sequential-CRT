library(tidyverse)

pvals = NULL
N = 500
for (seed in 1:N){
    filename = sprintf("output_pval/output_%i.Rdata", seed)
    if (file.exists(filename)) {
        load(filename)
        pvals = rbind(pvals, p_values)
    }
}

p = 500
N = dim(pvals)[1]
corrs = NULL
corrs0 = NULL
for (which_col in 1:(p-1)){
  corrs = c(corrs, cor(pvals[,which_col], pvals[,which_col+1]))
  corrs0 = c(corrs0, cor(pvals[1:(N-1),which_col], pvals[2:N,which_col+1]))
}
correlations = corrs
correlations_control = corrs0
dat.corr = data.frame(correlation = c(correlations, correlations_control), type = rep(c("correlation", "control"), each = length(correlations)))

library(ggplot2)

pp = ggplot(data = dat.corr) + geom_histogram(aes(x=correlation, fill=type), 
                 colour="grey50", alpha=0.5, position="identity", breaks=((-20):20)/200*3) +
    theme_bw() + scale_x_continuous(breaks=((-6):6)/20)
pp
pp %>% ggsave(file="pval_corr.pdf", width=7, height=3, units="in")

