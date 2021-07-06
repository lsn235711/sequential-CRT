library(glmnet)
args = commandArgs(TRUE)
seed = as.numeric(args[1])

n = 500
p = 500
k = 30
rho = 0.5
amplitude = 4
N = 500

set.seed(1)
nonnulls = sort(sample(1:p,k))

set.seed(seed)
beta = rep(amplitude,k)/sqrt(n)
Sigma = toeplitz(rho^(0:(p-1)))
y_model = function(X){
    y = X %*% beta
}

beta_hat_matrix = NULL
for (time in 1:N){
    print(time)
    X = matrix(rnorm(n*p),n) %*% chol(Sigma)
    y = y_model(X[,nonnulls]) + rnorm(n) 
    fit = cv.glmnet(X,y)
    beta_hat = coef(fit, fit$lambda.min)[2:(p+1)]
    beta_hat_matrix = rbind(beta_hat_matrix, beta_hat)
}


corrs = NULL
corrs0 = NULL
for (which_col in 1:(p-1)){
    corrs = c(corrs, cor(beta_hat_matrix[,which_col], beta_hat_matrix[,which_col+1]))
    corrs0 = c(corrs0, cor(beta_hat_matrix[1:(N-1),which_col], beta_hat_matrix[2:N,which_col+1]))
}

filename = sprintf("output_lasso/output_%i.Rdata", seed)
save(seed, corrs, corrs0, nonnulls, file = filename)


