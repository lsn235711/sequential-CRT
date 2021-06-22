args = commandArgs(TRUE)
amplitude = as.numeric(args[1])
setting_no = as.numeric(args[2])
seed = as.numeric(args[3])
#######################################################
source("seqstep_CRT_functions.R")
fdp = function(selected, nonnulls) length(setdiff(selected, nonnulls)) / max(1, length(selected))
power = function(selected, nonnulls) length(intersect(selected, nonnulls)) / length(nonnulls)


ptm <- proc.time()
n = 1000
p = 1000
N = 1
k = 50

print(amplitude)
print(setting_no)

if(setting_no == 4 | setting_no == 2){
    n = 1200
    p = 500
}

set.seed(1)
nonnulls = sort(sample(1:p,k))

set.seed(seed + 12345)
for (iii in 1:N){    
    
    ## HMM model
    K=5  # Number of possible states for each variable
    M=3  # Number of possible emission states for each variable
    States = 1:K
    Symbols = 1:M
    startProbs = rep(1,K)/K
    rho = 0.5
    transProbs = (1-rho)*matrix(rep(1/K,(K^2)), ncol = K) + rho*diag(rep(1,K))
    rho2 = 0.5
    Diag_g = matrix(c(rep(c(1,0.5,rep(0,(K*M-1)/(M-1)-3),0.5),M-1), 1), ncol = M)
    emissionProbs = (1-rho2)*matrix(rep(1/M,(K*M)), ncol = M) + rho2*Diag_g
    hmm1 = initHMM(States, Symbols, startProbs, transProbs, emissionProbs)
    
    
    beta = rep(amplitude,k)/sqrt(n)
    Sigma = toeplitz(rho^(0:(p-1)))
    X = matrix(rep(0,n*p), nrow = n)
    for (i in 1:n){
        X[i,] = simHMM(hmm1, length = p)$observation
    }
    
    if (setting_no == 1){
        beta = beta*1.6
        y_model = function(X){
            y = X %*% beta
        }
        model = "linear"
        y = y_model(X[,nonnulls]) + rnorm(n) 
        blackbox = "lasso"
    }
    if (setting_no == 2){
        beta = beta*8
        y_model = function(X){
            p_h = floor(dim(X)[2]/2)
            X_trans = (X[,2*(1:p_h)-1] > 1.5) * (X[,2*(1:p_h)] > 1.5)
            y = X_trans %*% beta[1:p_h]
        }
        model = "linear"
        y = y_model(X[,nonnulls]) + rnorm(n) 
        blackbox = "gb"
    }
    if (setting_no == 3){
        beta = beta*6
        y_model = function(X){
            y = (X-2) %*% beta
        }
        model = "logistic"
        z = y_model(X[,nonnulls])
        pr = 1/(1+exp(-z))         # pass through an inv-logit function
        y = (runif(n)<pr) + 0
        blackbox = "lasso"
    }
    if (setting_no == 4){
        beta = beta*20
        y_model = function(X){
            X_trans = (X > 1.5) - 2/3
            y = X_trans %*% beta
        }
        model = "logistic"
        z = y_model(X[,nonnulls])
        pr = 1/(1+exp(-z))
        y = (runif(n)<pr) + 0
        blackbox = "gb"
    }
    
    fdp_power = NULL
    for (one_shot_yes in c(TRUE)){
        print(sprintf("One_shot_CRT? %s", one_shot_yes))
        #for (blackbox in c("lasso", "gb")){
        #for (blackbox in c("lasso")){
        print(blackbox)
        for (method in c("inexact", "exact", "knockoffs")){
                print(sprintf("Method = %s", method))
                ##get selected sets
                ptm <- proc.time()
                selected_sets = selected_set(X, y, hmm1 = hmm1, X_model = "HMM", model = model, blackbox = blackbox, one_shot = one_shot_yes,
                                             do_CRT_seqstep_inexact = (method == "inexact"),
                                             do_CRT_seqstep_exact = (method == "exact"),
                                             do_knockoff = (method == "knockoffs"),
                )
                time_elapsed = (proc.time() - ptm)[3]
                ##find fdps and powers
                fdps = unlist(lapply(selected_sets$selected_list, nonnulls = nonnulls, fdp))
                powers = unlist(lapply(selected_sets$selected_list, nonnulls = nonnulls, power))
                
                results = data.frame(fdp = fdps, power= powers,include_h = selected_sets$include_h, 
                                     c = selected_sets$c, knockoff_plus = selected_sets$knockoff_plus,
                                     method = selected_sets$method, blackbox = selected_sets$blackbox,
                                     time = rep(time_elapsed, length(fdps)),
                                     one_shot = rep(one_shot_yes, length(fdps)))
                fdp_power = rbind(fdp_power, results)
            }
        
    }
    filename = sprintf("output_hmm/output_amp%i_setting%i_seed%i_i%i.Rdata", amplitude, setting_no, seed, iii)
    save(seed,amplitude, setting_no, fdp_power, iii, file = filename)
}
