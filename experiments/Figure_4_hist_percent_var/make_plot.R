library(tidyverse)

pvals = NULL
NseedXY = 100
NseedX2 = 4

results = NULL
for (seedXY in 1:NseedXY){
    for (seedX2 in 1:NseedX2){
        filename = sprintf("output_xvar/output_XY_%i_X2_%i.Rdata", seedXY, seedX2)
        if (file.exists(filename)) {
            load(filename)
            results = rbind(results, result)
        }
    }
}

library(tidyverse)


## within group can estimate E[Var[p1|X2, rest]]
exp_withVar = results %>% group_by(X2_group) %>%
  summarise(withVar = var(p_value)) %>%
  summarise(exp_withVar = mean(withVar))
exp_withVar = unlist(exp_withVar)

## total variance can estimate E[Var[p1, rest]]
var(results$p_value)

## percentage of variance explained by x1 itself
exp_withVar/var(results$p_value)

####################################################
pvals = NULL
NseedXY = 100
NseedX2 = 4

percentages = NULL
for (seedXY in 1:NseedXY){
  results = NULL
  for (seedX2 in 1:NseedX2){
    filename = sprintf("/Users/lsn/Desktop/Can/project_CRT_seqstep/simulations/simulation_20210621_not_too_dependent/sherlock_xvar/output_xvar/output_XY_%i_X2_%i.Rdata", seedXY, seedX2)
    if (file.exists(filename)) {
      load(filename)
      results = rbind(results, result)
    }
  }
  ## within group can estimate E[Var[p1|X2, rest]]
  exp_withVar = results %>% group_by(X2_group) %>%
    summarise(withVar = var(p_value)) %>%
    summarise(exp_withVar = mean(withVar))
  exp_withVar = unlist(exp_withVar)
  percentages = c(percentages, exp_withVar/var(results$p_value))
}

hist(percentages)

dat = data.frame(percentage = pmin(percentages,1))
pp = ggplot(data = dat) + geom_histogram(aes(x=percentage), fill = "lightblue", color = "grey50",
                                         breaks=seq(0.75,1,0.01)) +
    theme_bw()
pp
pp %>% ggsave(file="xval_percentage.pdf", width=7, height=3, units="in")


