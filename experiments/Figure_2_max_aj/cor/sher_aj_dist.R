args = commandArgs(TRUE)
arg = as.numeric(args[1])
seed = arg
set.seed(seed)


##############################################
ptm <- proc.time()
# Problem parameters
n = 200 ##number of observations
p = 120 ##number of covariates
q = 0.1
amplitude = 3
B1 = 3000 ## number of sample for each one of them
##B1 = 300

rho = 0.3
c = 0.3
p_value_precision = 19
k = 30 ## number of nonnulls 30



###########################
## Create block of size 3 ####
###########################
blo_size = 3
number_blo = p/blo_size
membership = (1:p - 1)%/%blo_size + 1
neighbors_matrix = NULL
for (i in 1:p){
  neighbors_matrix = rbind(neighbors_matrix, (membership[i]-1)*blo_size + (0:(blo_size-2) + i)%%blo_size + 1)
}
Sigma = matrix(rep(0,p^2), ncol = p)
for (i in 1:number_blo){
  Sigma[(i-1)*blo_size + 1:blo_size,(i-1)*blo_size + 1:blo_size] = rho
}
diag(Sigma) = 1


###########################
## Generate Data ####
###########################
beta = rep(0,p)
set.seed(1111)
nonnulls = sample(1:p, k)
nulls = setdiff(1:p,nonnulls)
beta[nonnulls] = amplitude * (2*(rbinom(k,1,0.5))-1)/ sqrt(n)

set.seed(seed)
X = matrix(rnorm(n*p),n) %*% chol(Sigma)
Y = X %*% beta + rnorm(n)


###########################
## To compute p-values ####
###########################
small_matrix = Sigma[1:blo_size, 1:blo_size]
small_matrix22_inv = solve(small_matrix[1:(blo_size-1),1:(blo_size-1)])
small_matrix21 = small_matrix[1:(blo_size-1),blo_size]
small_matrix12 = small_matrix[blo_size, 1:(blo_size-1)]
cond_sigma2 = Sigma[blo_size,blo_size] - t(small_matrix12) %*% small_matrix22_inv %*% small_matrix21
cond_mean_12_22_inv = t(small_matrix12) %*% small_matrix22_inv

find_p_val = function(X, Y){
  ## for case where block size is 5 here
  cond_mean = cond_mean_12_22_inv[1] * (X[,neighbors_matrix[,1]] + X[,neighbors_matrix[,2]])
  cond_var = matrix(rep(cond_sigma2, each = n*p), nrow = n)
  
  corr_0 = abs(apply(X, 2, cor, Y))
  corrs = NULL
  for (l in 1:p_value_precision){
    newX = cond_mean + sqrt(cond_var)*rnorm(n*p)
    corrs = cbind(corrs, abs(apply(newX, 2, cor, Y)))
  }
  
  p_val = (apply((corrs >= corr_0), 1, sum) + 1)/ (p_value_precision)
  return(p_val)
}

find_p_val2 = function(X, Y, js){
  l = length(js)
  cond_mean = cond_mean_12_22_inv[1] * (X[,neighbors_matrix[js,1]] + X[,neighbors_matrix[js,2]])
  cond_var = matrix(rep(cond_sigma2, each = n*l), nrow = n)
  
  corr_0 = abs(apply(X[,js], 2, cor, Y))
  corrs = NULL
  for (ll in 1:p_value_precision){
    newX = cond_mean + sqrt(cond_var)*rnorm(n*l)
    corrs = cbind(corrs, abs(apply(newX, 2, cor, Y)))
  }
  
  p_val = (apply((corrs >= corr_0), 1, sum) + 1)/ (p_value_precision)
  return(p_val)
}

###################################
## To sample X condition on Y #####
###################################
sigma2_Y = as.numeric(t(beta) %*% Sigma %*% beta) + 1
Sigma_12 = Sigma %*% beta
mu_X_cond = Y %*% t(Sigma_12)/ sigma2_Y
Sigma_X_cond = Sigma - Sigma_12 %*% t(Sigma_12) / sigma2_Y
Sigma_X_cond_chol = chol(Sigma_X_cond)

#p_vals = NULL
lengths = rep(0,p)
lengths[nonnulls]  = B1
m = 0

p_val0 = find_p_val(X, Y)
p_val0_c = 2*(p_val0 <= c)-1
ajs = rep(list(NULL),p)

## to remove the NA's
for (j in nonnulls){
  ajs[[j]] = 0
}

while (min(lengths) < B1){
  m = m + 1
  if (m %% 500 == 0){
    print(m)
    print(lengths)
    print(proc.time() - ptm)
  }
  indices = which(lengths < B1)
  new_X = matrix(rnorm(p*n), nrow = n) %*% Sigma_X_cond_chol  + mu_X_cond
  p_val = rep(NA, p)
  ## we want the indices to include its neighbors
  index_ngb = unique(c(indices, as.vector(neighbors_matrix[indices,])))
  index_ngb = sort(index_ngb)
  p_val[index_ngb] = find_p_val2(new_X, Y, index_ngb)
  #p_vals = rbind(p_vals, p_val)
  
  p_val_c = 2*(p_val <= c)-1
  same_sign = p_val_c * p_val0_c
  same_sign = c(same_sign, 1, 1)
  
  for (j in indices){
    ## finding a_j here
    same_sign_subset = same_sign[neighbors_matrix[j,]]
    if(min(same_sign_subset)>0){
      a_j = (p_val_c[j] + 1)/2
      ajs[[j]] = c(ajs[[j]],a_j)
      lengths[j] = lengths[j] + 1
    }
  }
}

ajs_mean = unlist(lapply(ajs,mean))
#plot(unlist(lapply(ajs,mean)))
print(max(ajs_mean))
save(seed,ajs_mean, file = paste("output/output",seed,".Rdata", sep = ""))
