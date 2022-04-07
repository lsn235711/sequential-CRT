args = commandArgs(TRUE)
amplitude = as.numeric(args[1])
setting_no = as.numeric(args[2])
seed = as.numeric(args[3])
#######################################################
library(GM)
library(CompQuadForm)

source("sequential_CRT_functions.R")
source("functions.R") ## distilled CRT

fdp = function(selected, nonnulls) length(setdiff(selected, nonnulls)) / max(1, length(selected))
power = function(selected, nonnulls) length(intersect(selected, nonnulls)) / length(nonnulls)


ptm <- proc.time()
n = 1000
p = 1000
rho = 0.5 ##AR model
k = 50 ## number of nonnulls
N = 1

print(amplitude)
print(setting_no)

if(setting_no == 4 | setting_no == 2){
    n = 1200
    p = 500
}

nonnulls = sort(sample(1:p,k))

set.seed(seed + 12345)
beta = rep(amplitude,k)/sqrt(n)
Sigma = toeplitz(rho^(0:(p-1)))
X = matrix(rnorm(n*p),n) %*% chol(Sigma)

if (setting_no == 1){
    beta = beta*1.5
    y_model = function(X){
        y = X %*% beta
    }
    model = "linear"
    y = y_model(X[,nonnulls]) + rnorm(n) 
    blackbox = "lasso"
}
if (setting_no == 2){
    beta = beta*15
    y_model = function(X){
        p_h = floor(dim(X)[2]/2)
        X_trans = (X[,2*(1:p_h)-1] > 0) * (X[,2*(1:p_h)] > 0)
        y = X_trans %*% beta[1:p_h]
    }
    model = "linear"
    y = y_model(X[,nonnulls]) + rnorm(n) 
    blackbox = "gb"
}
if (setting_no == 3){
    beta = beta*6
    y_model = function(X){
        y = X %*% beta
    }
    model = "logistic"
    z = y_model(X[,nonnulls])
    pr = 1/(1+exp(-z))         # pass through an inv-logit function
    y = (runif(n)<pr) + 0
    blackbox = "lasso"
}
if (setting_no == 4){
    beta = beta * 20
    Sigma = toeplitz(rho^(0:(p-1)))
    X = matrix(rnorm(n*p),n) %*% chol(Sigma)
    y_model = function(X){
        X_trans = (X > 0) - (X < 0)
        y = X_trans %*% beta
    }
    model = "logistic"
    z = y_model(X[,nonnulls])
    pr = 1/(1+exp(-z))
    y = (runif(n)<pr) + 0
    blackbox = "gb"
}


fdp_power = NULL

## Seqstep_CRT and knockoffs
# one_shot_yes = TRUE
# print(blackbox)
# for (method in c("sym_stats", "split","knockoffs")){
#     print(sprintf("Method = %s", method))
#     ##get selected sets
#     selected_sets = selected_set(X, y, Sigma = Sigma, rho = 0.5, model = model, blackbox = blackbox, one_shot = one_shot_yes,
#                                  do_sequential_CRT_sym_stats = (method == "sym_stats"),
#                                  do_sequential_CRT_split = (method == "split"),
#                                  do_knockoff = (method == "knockoffs"),
#                                  include_hs = FALSE,
#                                  c = 0.1,
#                                  knockoff_plus = 1
#                                  
#     )
#     ##find fdps and powers
#     fdps = unlist(lapply(selected_sets$selected_list, nonnulls = nonnulls, fdp))
#     powers = unlist(lapply(selected_sets$selected_list, nonnulls = nonnulls, power))
#     
#     results = data.frame(fdp = fdps, power= powers, method = selected_sets$method)
#     fdp_power = rbind(fdp_power, results)
# }

## Gaussian Mirror
if(setting_no == 1){
    fit = gm(y, X)
} else if(setting_no == 3){
    fit = gm(y, X, family = "binomial")
} 
fdp_gm = fdp(fit$gm_selected, nonnulls)
power_gm = power(fit$gm_selected, nonnulls)
fdp_power = rbind(fdp_power, data.frame(fdp = fdp_gm, power = power_gm, method = "GM"))

## d0CRT
if(setting_no == 1){
    model.dCRT = 'Gaussian_lasso'
} else if(setting_no == 3){
    model.dCRT = 'Binomial_lasso'
} else{
    model.dCRT = 'RF'
}

fit_d0CRT = dCRT(X, y, Sigma_X = Sigma, FDR = 0.1, model = model.dCRT)
fdp_d0CRT = fdp(fit_d0CRT$select_set, nonnulls)
power_d0CRT = power(fit_d0CRT$select_set, nonnulls)
fdp_power = rbind(fdp_power, data.frame(fdp = fdp_d0CRT, power = power_d0CRT, method = "d0CRT"))

## dICRT
fit_dICRT = dCRT(X, y, Sigma_X = Sigma, 
                 d.interaction = T, k = as.integer(2 * log(p)),
                 FDR = 0.1, MC_free = T, model = model.dCRT)
fdp_dICRT = fdp(fit_dICRT$select_set, nonnulls)
power_dICRT = power(fit_dICRT$select_set, nonnulls)
fdp_power = rbind(fdp_power, data.frame(fdp = fdp_dICRT, power = power_dICRT, method = "dICRT"))


## HRT
if((setting_no == 1) | (setting_no == 2)){
    model.HRT = 'gaussian'
} else{
    model.HRT = 'binomial'
} 
fit_HRT = HRT(y, X, Sigma = Sigma, FDR = 0.1, N = 50000, model = model.HRT)
fdp_HRT = fdp(fit_HRT$select_set, nonnulls)
power_HRT = power(fit_HRT$select_set, nonnulls)
fdp_power = rbind(fdp_power, data.frame(fdp = fdp_HRT, power = power_HRT, method = "HRT"))

filename = sprintf("output/output_amp%i_setting%i_seed%i.Rdata", amplitude, setting_no, seed)
save(seed,amplitude, setting_no, fdp_power, file = filename)

