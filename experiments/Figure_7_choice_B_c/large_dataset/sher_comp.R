args = commandArgs(TRUE)

B_choice = as.numeric(args[1])
setting_no = as.numeric(args[2])
seed = as.numeric(args[3])
amplitude_no = as.numeric(args[4])
#######################################################
source("sequential_CRT_functions.R")


fdp = function(selected, nonnulls) length(setdiff(selected, nonnulls)) / max(1, length(selected))
power = function(selected, nonnulls) length(intersect(selected, nonnulls)) / length(nonnulls)


ptm <- proc.time()
n = 1000
p = 1000
rho = 0.5 ##AR model
k = 50 ## number of nonnulls
N = 1

if (setting_no == 1){
    amplitude = 2 + 1.5*amplitude_no
} else{
    amplitude = 10*amplitude_no
}

nonnulls = sort(sample(1:p,k))

set.seed(seed + 12345)
beta = rep(amplitude,k)/sqrt(n)
Sigma = toeplitz(rho^(0:(p-1)))
X = matrix(rnorm(n*p),n) %*% chol(Sigma)

if (setting_no == 1){
    y_model = function(X){
        y = X %*% beta
    }
    model = "linear"
    y = y_model(X[,nonnulls]) + rnorm(n) 
    blackbox = "lasso"
}
if (setting_no == 3){
    y_model = function(X){
        y = X %*% beta
    }
    model = "logistic"
    z = y_model(X[,nonnulls])
    pr = 1/(1+exp(-z))         # pass through an inv-logit function
    y = (runif(n)<pr) + 0
    blackbox = "lasso"
}


fdp_power = NULL

## Seqstep_CRT
fdp_power = NULL
one_shot_yes = TRUE
c_choices = 0.05*(1:10)
for (method in c("sym_stats", "split")){
    print(sprintf("Method = %s", method))
    ##get selected sets
    selected_sets = selected_set(X, y, Sigma = Sigma, rho = 0.5, model = model, blackbox = blackbox, one_shot = one_shot_yes,
                                 do_sequential_CRT_sym_stats = (method == "sym_stats"),
                                 do_sequential_CRT_split = (method == "split"),
                                 do_knockoff = FALSE,
                                 include_hs = FALSE,
                                 c = c_choices,
                                 knockoff_plus = 1,
                                 B = B_choice)
    ##find fdps and powers
    fdps = unlist(lapply(selected_sets$selected_list, nonnulls = nonnulls, fdp))
    powers = unlist(lapply(selected_sets$selected_list, nonnulls = nonnulls, power))
    
    results = data.frame(fdp = fdps, power= powers, method = selected_sets$method, B = rep(B_choice, length(fdps)), c = c_choices)
    fdp_power = rbind(fdp_power, results)
}

filename = sprintf("output/output_B%i_setting%i_seed%i_amp%i.Rdata", B_choice, setting_no, seed, amplitude_no)
save(seed, B_choice, setting_no, amplitude, fdp_power, file = filename)

