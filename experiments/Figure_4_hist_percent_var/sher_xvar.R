library(glmnet)
source("seqstep_CRT_functions.R")
args = commandArgs(TRUE)
seedXY = as.numeric(args[1])
seedX2 = as.numeric(args[2])

n = 500
p = 500
k = 30
rho = 0.5
amplitude = 4

N2 = 5
N1 = 20

set.seed(1)
nonnulls = sort(sample(1:p,k))

set.seed(seedXY)
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

##sample x2 conditional on x3...and Y
sample_x2 = function(X, y, rho = 0.5, sigma2 = 1){
    mu_c = X[,3] * rho
    Sigmacho = sqrt(sigma2*(1 - rho^2))
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

#pj = find_one_pval(X,y,1)
#pj

## 1 and 2 are both nulls here
set.seed(seedX2)
X2_groups = NULL
X1_groups = NULL
p_values = NULL
for (time2 in 1:N2){
    print(time2)
    ## sample X2 conditional on X3... and y
    x2new = sample_x2(X,y)
    Xnew = X
    Xnew[,2] = x2new
    for(time1 in 1:N1){
        print(time1)
        ## sample X1 conditional on the sampled X1;X3 ... and y
        x1new = sample_xj(Xnew,y,1)
        Xnew[,1] = x1new
        p2_one = find_one_pval(Xnew, y, 2)
        X2_groups = c(X2_groups, time2)
        X1_groups = c(X1_groups, time1 + time2*(N2-1))
        p_values = c(p_values, p2_one)
    }
}


result = data.frame(X2_group = X2_groups, X1_group = X1_groups, p_value = p_values)

filename = sprintf("output_xvar/output_XY_%i_X2_%i.Rdata", seedXY, seedX2)
save(seedXY, seedX2, result, file = filename)


