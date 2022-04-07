################################################################
## These functions run the sequential CRT
################################################################
library(knockoff)
library(glmnet)
library(xgboost)
library(randomForest)
library(HMM)
library(SNPknock)
################################################################

## Compute feature importance statistics
## Inputs: 
##    X: matrix of predictors
##    y: vector of response variable
##    blackbox: "lasso" is glmnet; "rf" is random forest; "gb" is gradient boosting
##    model: "logistic" when y is binary; "gaussian" when y is continuous
## Output:
##    A vector of feature importance (length is the same as the number of predictors in X)
feature_importance = function(X,y, blackbox, model){
    if (model == "logistic"){
        family = "binomial"
        objective = "binary:logistic"
    } else{    
        family = "gaussian"
        objective = "reg:squarederror"
    }

    if (blackbox == "lasso"){
        p = dim(X)[2]
        fit = cv.glmnet(X, y, family = family)
        lambda = fit$lambda.min
        coeffs = coef(fit,s=lambda)
        return (coeffs[2:(p+1)])
    } 
    else if(blackbox == "gb"){
        p = dim(X)[2]
        colnames(X) = paste("X.", 1:p,  sep="")
        dtrain = xgb.DMatrix(data = X, label = y) 
        params <- list(booster = "gbtree", objective = objective, eta=0.05, max_depth=2)
        xgb1 <- xgb.train(params = params, data = dtrain, nrounds = 100)
        importance_matrix = xgb.importance(feature_names = colnames(X), model = xgb1)
        importance_gb = rep(0,p)
        for (i in 1:p){
            index_of_i = which(importance_matrix$Feature == paste("X.", i,  sep=""))
            if (length(index_of_i) > 0){
                importance_gb[i] = importance_matrix$Gain[index_of_i]
            }
        }
        return (importance_gb)
    }
    else if(blackbox == "rf"){
        if (model == "logistic"){
            y = as.factor(y)
        }
        rf = randomForest(X, y, ntree=500)
        return (importance(rf))
    }
}

## Compute feature importance statistics for knockoff filter
## Inputs: 
##    X: matrix of predictors
##    X_k: matrix of knockoffs
##    y: vector of response variable
##    blackbox: "lasso" is glmnet; "rf" is random forest; "gb" is gradient boosting
##    model: "logistic" when y is binary; "gaussian" when y is continuous
## Output:
##    A vector of feature importance (length is the same as the number of predictors in X)
feature_importance_knockoffs = function(X, X_k, y, blackbox, model){
    X_2 = cbind(X, X_k)
    Z2 = feature_importance(X_2,y,blackbox, model)
    W = Z2[1:p] - Z2[(p+1):(2*p)]
    return(W)
}

## Obtain p-values from CRT and the statistics for ordering when X is from an AR(1) model
##        The statistics are computed as in the symmetric statistics method
## Inputs: 
##    X: matrix of predictors
##    y: vector of response variable
##    rho: parameter of the AR(1) model, correlation
##    sigma2: parameter of the AR(1) model, variance
##    B: number of randomizations in CRT
##    blackbox: "lasso" is glmnet; "rf" is random forest; "gb" is gradient boosting
##    model: "logistic" when y is binary; "gaussian" when y is continuous
##    one_shot: whether run CRT with "one shot" CRT
## Output:
##    A matrix, whose first column consists of the p-values, and the second column consists
##        of the statistics   
sequential_CRT_AR = function(X, y, rho, sigma2, B, blackbox, model, one_shot = TRUE){
    n = dim(X)[1]
    p = dim(X)[2]
    ps <- rep(1,p)
    betas <- rep(1,p)
    
    for (j in 1:p){
        print(j)
        ##AR model with mu = 0
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
        
        if(one_shot){
            Xnew = X
            for (k in 1:B){
                x_tilde = Sigmacho * rnorm(n) + mu_c
                Xnew = cbind(Xnew,x_tilde)
            }
            
            betahats = abs(feature_importance(Xnew, y, blackbox, model)[c(j,(p+1):(p+B))])
        } else{
            betahats = abs(feature_importance(X, y, blackbox, model)[j])
            for (k in 1:B){
                Xnew = X
                x_tilde = Sigmacho * rnorm(n) + mu_c
                Xnew[,j] = x_tilde
                betahats = c(betahats, abs(feature_importance(Xnew, y, blackbox, model)[j]))
            }
        }
        
        betahat0 = betahats[1]
        betas[j] = max(betahats)
        numgreat = sum(betahats > betahat0)
        numeq = sum(betahats == betahat0)
        numergreateq = numgreat
        if(numeq > 0){
            numergreateq = numgreat + sample(1:numeq,1)
        }
        #ps[j] = (numergreateq-runif(1))/length(betahats)
        ps[j] = (numergreateq)/length(betahats)
        betas[j] = max(abs(betahats))
    }
    return (cbind(ps,betas))
}

## Obtain p-values from CRT and the statistics for ordering when X is from a gaussian model
##        The statistics are computed as in the symmetric statistics method
## Inputs: 
##    X: matrix of predictors
##    y: vector of response variable
##    mu: parameter of the gaussian model, vector mean
##    Sigma: parameter of the gaussian model, matrix of variance
##    B: number of randomizations in CRT
##    blackbox: "lasso" is glmnet; "rf" is random forest; "gb" is gradient boosting
##    model: "logistic" when y is binary; "gaussian" when y is continuous
##    one_shot: whether run CRT with "one shot" CRT
## Output:
##    A matrix, whose first column consists of the p-values, and the second column consists
##        of the statistics  
sequential_CRT_gaussian = function(X, y, mu, Sigma, B, blackbox, model, one_shot = TRUE){
    n = dim(X)[1]
    p = dim(X)[2]
    ps <- rep(1,p)
    betas <- rep(1,p)
    
    for (j in 1:p){
        print(j)
        Sigma22 = Sigma[-j,-j]
        Sigma22Inv = solve(Sigma22)
        Sigma12 = Sigma[j,-j]
        Sigmacho = sqrt(as.vector((1 - t(Sigma12) %*% Sigma22Inv %*% Sigma12)))
        SigmaProd = Sigma12 %*% Sigma22Inv
        mu_c = (X[,-j] - mu[-j]) %*%  t(SigmaProd) + mu[j]
        
        if(one_shot){
            Xnew = X
            for (k in 1:B){
                x_tilde = Sigmacho * rnorm(n) + mu_c
                Xnew = cbind(Xnew,x_tilde)
            }
            
            betahats = abs(feature_importance(Xnew, y, blackbox, model)[c(j,(p+1):(p+B))])
        } else{
            betahats = abs(feature_importance(X, y, blackbox, model)[j])
            for (k in 1:B){
                Xnew = X
                x_tilde = Sigmacho * rnorm(n) + mu_c
                Xnew[,j] = x_tilde
                betahats = c(betahats, abs(feature_importance(Xnew, y, blackbox, model)[j]))
            }
        }
        
        betahat0 = betahats[1]
        betas[j] = max(betahats)
        numgreat = sum(betahats > betahat0)
        numeq = sum(betahats == betahat0)
        numergreateq = numgreat
        if(numeq > 0){
            numergreateq = numgreat + sample(1:numeq,1)
        }
        #ps[j] = (numergreateq-runif(1))/length(betahats)
        ps[j] = (numergreateq)/length(betahats)
        betas[j] = max(abs(betahats))
    }
    return (cbind(ps,betas))
}

## Compute the posterior probability of X_j|X_{-j} for each j
## Inputs: 
##    X: matrix of predictors
##    hmm1: an hmm object from Package HMM
## Output:
##    The posterior probability of X_j|X_{-j} for each j
compute_HMMposterior = function(X, hmm1){
    n = dim(X)[1]
    p = dim(X)[2]
    posterior = array(rep(0, n*p*M), c(n,M,p))
    for (i in 1:n){
        backward_prob = exp(backward(hmm1, X[i,]))
        forward_prob = exp(forward(hmm1, X[i,]))
        forward_prob_plus = cbind(rep(1,K), forward_prob)
        posterior_X_t = t(t(forward_prob_plus[,1:p]) %*% transProbs) * backward_prob
        posterior_X_t_scaled = scale(posterior_X_t, center = FALSE, scale = apply(posterior_X_t, 2, sum))
        posterior_O_t = t(t(posterior_X_t_scaled) %*% emissionProbs)
        posterior[i,,] = posterior_O_t
    }
    return(posterior)
}

## Sample independently from the distribution of X_j|X_{-j}
## Inputs: 
##    posterior: matrix of posterior probabilities P(X_j|X_{-j})
##    B: number of independent samples
## Output:
##    A matrix consisting of B independent samples
sample_HMM_from_posterior = function(posterior, B){
    n = dim(posterior)[1]
    sampled_HMM = matrix(rep(0, n*B), nrow = B)
    for (step in 1:B){
        u_n = runif(n)
        sampled_HMM[step,] = 1 + apply(apply(posterior,1,cumsum) <= matrix(rep(u_n, each = M), nrow = M), 2, sum)
    }    
    return (sampled_HMM)
}

## Obtain p-values from CRT and the statistics for ordering when X is from an HMM model
##        The statistics are computed as in the symmetric statistics method
## Inputs: 
##    X: matrix of predictors
##    y: vector of response variable
##    hmm1: an hmm object from the HMM package
##    B: number of randomizations in CRT
##    blackbox: "lasso" is glmnet; "rf" is random forest; "gb" is gradient boosting
##    model: "logistic" when y is binary; "gaussian" when y is continuous
##    one_shot: whether run CRT with "one shot" CRT
## Output:
##    A matrix, whose first column consists of the p-values, and the second column consists
##        of the statistics  
sequential_CRT_HMM = function(X, y, hmm1, B, blackbox, model, one_shot = TRUE){
    n = dim(X)[1]
    p = dim(X)[2]
    ps <- rep(1,p)
    betas <- rep(1,p)
    
    posterior_array = compute_HMMposterior(X, hmm1)
    
    for (j in 1:p){
        print(j)
        
        if(one_shot){
            Xnew = t(sample_HMM_from_posterior(posterior_array[,,j], B))
            Xnew = cbind(X,Xnew)
            betahats = abs(feature_importance(Xnew, y, blackbox, model)[c(j,(p+1):(p+B))])
        } else{
            betahats = abs(feature_importance(X, y, blackbox, model)[j])
            Xnewjs = t(sample_HMM_from_posterior(posterior_array[,,j], B))
            for (k in 1:B){
                Xnew = X
                x_tilde = Xnewjs[,k]
                Xnew[,j] = x_tilde
                betahats = c(betahats, abs(feature_importance(Xnew, y, blackbox, model)[j]))
            }
        }
        
        betahat0 = betahats[1]
        betas[j] = max(betahats)
        numgreat = sum(betahats > betahat0)
        numeq = sum(betahats == betahat0)
        numergreateq = numgreat
        if(numeq > 0){
            numergreateq = numgreat + sample(1:numeq,1)
        }
        #ps[j] = (numergreateq-runif(1))/length(betahats)
        ps[j] = (numergreateq)/length(betahats)
        betas[j] = max(abs(betahats))
    }
    return (cbind(ps,betas))
}



## Report discoveries from the sequential CRT
## Inputs: 
##    X: matrix of predictors
##    y: vector of response variable
##    cs: the threshold c in Selective SeqStep+
##    one_shot: whether to run the "one shot" CRT
##    include_hs: whether to use additional information in p-values, default is FALSE
##    knockoff_pluss: 1 stands for Knockoffs+/SeqStep+; 0 stands for Knockoffs/SeqStep
##    blackbox: "lasso" is glmnet; "rf" is random forest; "gb" is gradient boosting
##    model: "logistic" when y is binary; "gaussian" when y is continuous
##    B: number of randomizations in CRT
##    do_sequential_CRT_sym_stats: whether to run the symmetric statistics version of the sequential CRT
##    do_sequential_CRT_split: whether to run the split version of the sequential CRT
##    do_knockoff: whether to run knockoffs as a comparison
##    q: FDR threshold q
##    X_model: distribution of X -- "AR" stands for AR(1) model; "gaussian" stands for general
##                                  gaussian model; "hmm" stands for hmm model
##    rho: parameter of the AR(1) model, correlation
##    mu: parameter of the gaussian model, vector mean
##    Sigma: parameter of the gaussian model/AR model, matrix of variance
##    hmm1: an hmm object from the HMM package
## Output:
##    A list consisting of results from all selected settings
##              selected_list: list of selected sets
##              c: vector of thresholds c
##              include_h: vector of include_h (TRUE/FALSE)
##              knockoff_plus: vector of knockoff_plus (0/1)
##              method: vector of methods considered ("sequential_CRT_sym_stats"/"sequential_CRT_split"/"knockoffs")
##              blackbox: vector of blackboxs considered ("lasso"/"gb"/"rf")
selected_set = function(X, y, cs = 0.1*(1:5), one_shot = TRUE,
            include_hs = c(FALSE,TRUE),
            knockoff_pluss = c(0,1),
            model = "linear",
            blackbox = "lasso",
            B = 9,
            do_sequential_CRT_sym_stats = TRUE,
            do_sequential_CRT_split = TRUE,
            do_knockoff = TRUE,
            q = 0.1,
            X_model = "AR", rho = NULL, Sigma = NULL, mu = NULL,
            hmm1 = NULL
            )
            {
    
    selected_list = list()
    list_i = 1
    c_list = NULL
    include_h_list = NULL
    knockoff_plus_list = NULL
    method_list = NULL
    n = dim(X)[1]
    p = dim(X)[2]
    
    ################### the sequential CRT symmetric statistics ###########
    if (do_sequential_CRT_sym_stats){
        if(X_model == "AR"){
            pvalue_W = sequential_CRT_AR(X, y, rho, Sigma[1,1], B, blackbox, model, one_shot)
        } else if (X_model == "HMM") {
            pvalue_W = sequential_CRT_HMM(X, y, hmm1, B, blackbox, model, one_shot)
        } else if (X_model == "gaussian") {
            pvalue_W = sequential_CRT_gaussian(X, y, mu, Sigma, B, blackbox, model, one_shot)
        }
        
        for (c in cs){
            p_values = pvalue_W[,1]
            W = pvalue_W[,2]
            for (include_h in include_hs){
                if (include_h) {
                    p_values = p_values + runif(p)/(B+1)
                    h_pvalue = (p_values <= c)*(c-p_values)/c + (p_values > c)*(1-p_values)/(1-c)
                    W = W*h_pvalue
                }
                #selective Seqstep
                W_sign = (p_values <= c) - (p_values > c)
                W_mag = abs(W)
                W_eff = W_sign * W_mag
                
                ##threshold
                for (knockoff_plus in knockoff_pluss){
                    thres = knockoff.threshold(W_eff, fdr=q*(1-c)/c, offset = knockoff_plus)
                    selected = intersect(which(W_eff >= thres), which(W_eff != 0))
                    
                    selected_list[[list_i]] = selected
                    c_list = c(c_list,c)
                    include_h_list = c(include_h_list, include_h)
                    knockoff_plus_list = c(knockoff_plus_list,knockoff_plus)
                    method_list = c(method_list,"sequential_CRT_sym_stats")
                    list_i = list_i + 1
                }
            }
        }
    }
    
    ################### the sequential CRT split ###########
    if (do_sequential_CRT_split){
        pvalue_proportion = 0.5
        W_proportion = n - pvalue_proportion
        n_p = floor(n*pvalue_proportion)
        n_W = n - n_p
        p_value_indices = sample(1:n, n_p)
        X_p = X[p_value_indices,]
        X_w = X[-p_value_indices,]
        y_p = y[p_value_indices]
        y_w = y[-p_value_indices]
        
        if(X_model == "AR"){
            pvalue_W = sequential_CRT_AR(X_p, y_p, rho, Sigma[1,1], B, blackbox, model, one_shot)
        } else if (X_model == "HMM") {
            pvalue_W = sequential_CRT_HMM(X_p, y_p, hmm1, B, blackbox, model, one_shot)
        } else if (X_model == "gaussian") {
            pvalue_W = sequential_CRT_gaussian(X_p, y_p, mu, Sigma, B, blackbox, model, one_shot)
        }
        
        for (c in cs){
            p_values = pvalue_W[,1]
            W = feature_importance(X_w, y_w, blackbox,model)
            for (include_h in include_hs){
                if (include_h) {
                    p_values = p_values + runif(p)/(B+1)
                    h_pvalue = (p_values <= c)*(c-p_values)/c + (p_values > c)*(1-p_values)/(1-c)
                    W = W*h_pvalue
                }
                #selective Seqstep
                W_sign = (p_values <= c) - (p_values > c)
                W_mag = abs(W)
                W_eff = W_sign * W_mag
                
                ##threshold
                for (knockoff_plus in knockoff_pluss){
                    thres = knockoff.threshold(W_eff, fdr=q*(1-c)/c, offset = knockoff_plus)
                    selected = intersect(which(W_eff >= thres), which(W_eff != 0))
                    selected_list[[list_i]] = selected
                    c_list = c(c_list,c)
                    include_h_list = c(include_h_list, include_h)
                    knockoff_plus_list = c(knockoff_plus_list,knockoff_plus)
                    method_list = c(method_list,"sequential_CRT_split")
                    list_i = list_i + 1
                }
            }
        }
    }
    
    ###################Comparing to original knockoff ###########
    if (do_knockoff){
        if (X_model == "AR"){
            X_k = create.gaussian(X, rep(0,p), Sigma)
        }else if (X_model == "HMM"){
            K = length(hmm1$States)
            M = length(hmm1$Symbols)
            pInit = hmm1$startProbs
            Q = array(rep(0,p*K*K),c(p-1,K,K))
            for(j in 1:(p-1)) { Q[j,,] = t(hmm1$transProbs)}
            pEmit = array(rep(0,p*M*K),c(p,M,K))
            for(j in 1:p) { pEmit[j,,] = t(hmm1$emissionProbs)}
            X_k = knockoffHMM(matrix(as.integer(X-1), nrow = n), pInit, Q, pEmit) + 1
        } else if (X_model == "gaussian"){
            X_k = create.gaussian(X, mu, Sigma)
        }
        W0 = feature_importance_knockoffs(X, X_k, y, blackbox, model)
        for (knockoff_plus in knockoff_pluss){
            thres0 = knockoff.threshold(W0, fdr=q, offset = knockoff_plus)
            selected0 = intersect(which(W0 >= thres0), which(W0 != 0))
        
            selected_list[[list_i]] = selected0
            c_list = c(c_list, NA)
            include_h_list = c(include_h_list, NA)
            knockoff_plus_list = c(knockoff_plus_list, knockoff_plus)
            method_list = c(method_list,"knockoffs")
            list_i = list_i + 1
        }
    }
    
    return (list(selected_list = selected_list, c = c_list, include_h = include_h_list, 
                 knockoff_plus = knockoff_plus_list, method = method_list, blackbox = rep(blackbox,list_i-1)))
}





