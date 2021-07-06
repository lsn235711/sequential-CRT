library(tidyverse)

correlations = NULL
correlations_control = NULL
seed = 1
filename = sprintf("output_lasso/output_%i.Rdata", seed)
load(filename)

p = 500
null_set = !(((1:(p-1)) %in% nonnulls) | (2:p %in% nonnulls))

correlations = corrs[null_set]
correlations_control = corrs0[null_set]
dat.corr = data.frame(correlation = c(correlations, correlations_control), type = rep(c("correlation", "control"), each = length(correlations)))

library(ggplot2)

p = ggplot(data = dat.corr) + geom_histogram(aes(x=correlation, fill=type), 
                 colour="grey50", alpha=0.5, position="identity", breaks=((-20):20)/100) +
    theme_bw()
p
p %>% ggsave(file="lasso_corr.pdf", width=7, height=3, units="in")

