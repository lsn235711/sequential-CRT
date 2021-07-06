library(glmnet)
source("seqstep_CRT_functions.R")
args = commandArgs(TRUE)
seed = as.numeric(args[1])

n = 500
p = 500
k = 30
rho = 0.5
amplitude = 4


set.seed(1)
nonnulls = sort(sample(1:p,k))

set.seed(seed)
beta = rep(amplitude,k)/sqrt(n)
Sigma = toeplitz(rho^(0:(p-1)))
y_model = function(X){
    y = X %*% beta
}
X = matrix(rnorm(n*p),n) %*% chol(Sigma)
y = y_model(X[,nonnulls]) + rnorm(n) 

sample_xj = function(X, y, j, rho = 0.5, sigma2 = 1){
    if(j == 1){
        mu_c = X[,2] * rho
        Sigmacho = sqrt(sigma2*(1 - rho^2))
    } else if(j == p){
        mu_c = X[,p-1] * rho
        Sigmacho = sqrt(sigma2*(1 - rho^2))
    } else{
        Sigma22Inv = solve(matrix(c(1,rho^2,rho^2,1),ncol = 2))
        Sigma12 = c(rho,rho)
        Sigmacho = sqrt(as.vector(sigma2*(1 - t(Sigma12) %*% Sigma22Inv %*% Sigma12)))
        SigmaProd = Sigma12 %*% Sigma22Inv
        mu_c = X[,c(j-1,j+1)] %*%  t(SigmaProd)
    }
    x_tilde = Sigmacho * rnorm(n) + mu_c
    return(x_tilde)
}

find_one_pval = function(X, y, j, rho = 0.5, sigma2 = 1, blackbox="lasso", model="linear", B = 19, one_shot = TRUE){
    betahats = abs(feature_importance(X, y, blackbox, model)[j])
    for (k in 1:B){
        Xnew = X
        Xnew[,j] = sample_xj(X,y,j)
        betahats = c(betahats, abs(feature_importance(Xnew, y, blackbox, model)[j]))
    }
    betahat0 = betahats[1]
    numgreat = sum(betahats > betahat0)
    numeq = sum(betahats == betahat0)
    numergreateq = numgreat
    if(numeq > 0){
        numergreateq = numgreat + sample(1:numeq,1)
    }
    pj = (numergreateq)/length(betahats)
    return(pj)
}

p_values = NULL
for (j in 1:p){
    print(j)
    pval_one = find_one_pval(X, y, j)
    p_values = c(p_values, pval_one)
}

filename = sprintf("output_pval/output_%i.Rdata", seed)
save(seed, p_values, file = filename)


