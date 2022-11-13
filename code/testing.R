#' ---
#' title:  Testing code related to the project on homoegeneity pursuit in GEEs for multivariate abundance data
#' author: Originally written by FKCH
#' date: "Started September 2022"
#' ---

##---------------------
#' # Load appropriate packages and files
##---------------------
rm(list = ls())
library(tidyverse)
library(mvtnorm)
library(mvabund)
library(mgcv)
here::i_am("testing.R")
library(here)

source(here("auxfnsv0.R"))
source(here("calc_penaltyweights.R"))
source(here("create_life.R"))
source(here("familyfns.R"))
source(here("geefn.R"))
source(here("hpgeefn.R"))
source(here("updatecoeffns.R"))


##--------------------
#' # Example 1: Simulating independence response data
##---------------------
set.seed(092022)

num_units <- 400
num_resp <- 12
num_lv <- 2
num_X <- 9 # Excludes intercept
rho <- 0.5
P <- outer(1:num_X, 1:num_X, FUN = "-") %>% abs
P <- rho^P 

X <- rmvnorm(num_units, sigma = P) %>% cbind(1, .)
colnames(X) <- c("Int.", paste0("x",1:num_X))

X_coefs <- cbind(runif(num_resp, -2, 0), matrix(rnorm(num_resp*num_X, 0, sd = 0.5), nrow = num_resp))
true_loadings <- matrix(0, nrow = num_resp, ncol = num_lv)
true_disp_param <- runif(num_resp, 0, 4)


#' ## Gaussian responses
simy <- create_life(family = gaussian(), X = X, X_coefs = X_coefs, loadings = true_loadings, disp_param = true_disp_param)
X <- as.data.frame(X)
fit_goldstd <- manylm(simy ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data = X)
fit_geeind <- unpenalized_gee(y = simy, X = X, family = gaussian(), control = list(trace = TRUE))
fit_geeunstr <- unpenalized_gee(y = simy, X = X, family = gaussian(), control = list(trace = TRUE), corstr = "unstructured")
fit_geeRR <- unpenalized_gee(y = simy, X = X, family = gaussian(), control = list(trace = TRUE), corstr = "reducedrank", rank = 5)
fit_geeind_glmnet <- unpenalized_gee(y = simy, X = X, family = gaussian(), control = list(trace = TRUE), method = "glmnet")
fit_geeunstr_glmnet <- unpenalized_gee(y = simy, X = X, family = gaussian(), control = list(trace = TRUE), corstr = "unstructured", method = "glmnet")
fit_geeRR_glmnet <- unpenalized_gee(y = simy, X = X, family = gaussian(), control = list(trace = TRUE), corstr = "reducedrank", rank = 5, method = "glmnet")
# fit_geeind$coefficients_mat - t(fit_goldstd$coefficients) # Good!
# fit_geeunstr$coefficients_mat - t(fit_goldstd$coefficients) # Good!
# fit_geeRR$coefficients_mat - t(fit_goldstd$coefficients) # Good!
data.frame(true = X_coefs %>% as.vector, 
           gold = t(fit_goldstd$coefficients) %>% as.vector,
           geeind = fit_geeind$coefficients_mat %>% as.vector, 
           geeunstr = fit_geeunstr$coefficients_mat %>% as.vector, 
           geeRR = fit_geeRR$coefficients_mat %>% as.vector, 
           geeind_glmnet = fit_geeind_glmnet$coefficients_mat %>% as.vector, 
           geeunstr_glmnet = fit_geeunstr_glmnet$coefficients_mat %>% as.vector, 
           geeRR_glmnet = fit_geeRR_glmnet$coefficients_mat %>% as.vector 
           ) %>% 
   t %>% 
   dist


#' ## Binary responses
simy <- create_life(family = binomial(link = "probit"), X = X, X_coefs = X_coefs, loadings = true_loadings, disp_param = true_disp_param)
X <- as.data.frame(X)
fit_goldstd <- lapply(1:num_resp, function(x) mgcv::gam(simy[,x] ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data = X, family = binomial(link = "probit")))
fit_geeind <- unpenalized_gee(y = simy, X = X, family = binomial(link = "probit"), control = list(trace = TRUE))
fit_geeunstr <- unpenalized_gee(y = simy, X = X, family = binomial(link = "probit"), control = list(trace = TRUE), corstr = "unstructured")
fit_geeRR <- unpenalized_gee(y = simy, X = X, family = binomial(link = "probit"), control = list(trace = TRUE), corstr = "reducedrank", rank = 5)
fit_geeind_glmnet <- unpenalized_gee(y = simy, X = X, family = binomial(link = "probit"), control = list(trace = TRUE), method = "glmnet")
fit_geeunstr_glmnet <- unpenalized_gee(y = simy, X = X, family = binomial(link = "probit"), control = list(trace = TRUE), corstr = "unstructured", method = "glmnet")
fit_geeRR_glmnet <- unpenalized_gee(y = simy, X = X, family = binomial(link = "probit"), control = list(trace = TRUE), corstr = "reducedrank", rank = 5, method = "glmnet")
# fit_geeind$coefficients_mat - t(sapply(fit_goldstd, coef)) # Good!
# fit_geeunstr$coefficients_mat - t(sapply(fit_goldstd, coef)) # Good enough!
# fit_geeRR$coefficients_mat - t(sapply(fit_goldstd, coef)) # Good enough!
data.frame(true = X_coefs %>% as.vector, 
           gold = t(sapply(fit_goldstd, coef)) %>% as.vector,
           geeind = fit_geeind$coefficients_mat %>% as.vector, 
           geeind_glmnet = fit_geeind_glmnet$coefficients_mat %>% as.vector, 
           geeunstr = fit_geeunstr$coefficients_mat %>% as.vector, 
           geeunstr_glmnet = fit_geeunstr_glmnet$coefficients_mat %>% as.vector, 
           geeRR = fit_geeRR$coefficients_mat %>% as.vector, 
           geeRR_glmnet = fit_geeRR_glmnet$coefficients_mat %>% as.vector 
           ) %>% 
   t %>% 
   dist


#' ## Poisson responses 
simy <- create_life(family = poisson(), X = X, X_coefs = X_coefs, loadings = true_loadings)
X <- as.data.frame(X)
fit_goldstd <- manyglm(simy ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data = X, family = "poisson")
# fit_ind <- .manual_glm(bigX = bigX, y = c(t(simy)), num_resp = num_resp, num_units = num_units)
# fit_ind$coefficients_mat <- matrix(fit_ind$coefficients, nrow = num_resp, byrow = TRUE)
# fit_ind$coefficients_mat - fit_long$coefficients_mat 
fit_geeind <- unpenalized_gee(y = simy, X = X, family = poisson(), control = list(trace = TRUE))
fit_geeunstr <- unpenalized_gee(y = simy, X = X, family = poisson(), control = list(trace = TRUE), corstr = "unstructured")
fit_geeRR <- unpenalized_gee(y = simy, X = X, family = poisson(), control = list(trace = TRUE), corstr = "reducedrank", rank = 5)
fit_geeind_glmnet <- unpenalized_gee(y = simy, X = X, family = poisson(), control = list(trace = TRUE), method = "glmnet")
fit_geeunstr_glmnet <- unpenalized_gee(y = simy, X = X, family = poisson(), control = list(trace = TRUE), corstr = "unstructured", method = "glmnet")
fit_geeRR_glmnet <- unpenalized_gee(y = simy, X = X, family = poisson(), control = list(trace = TRUE), corstr = "reducedrank", rank = 5, method = "glmnet")
# fit_geeind$coefficients_mat - t(fit_goldstd$coefficients) # Good!
# fit_geeunstr$coefficients_mat - t(fit_goldstd$coefficients) # Good enough!
# fit_geeRR$coefficients_mat - t(fit_goldstd$coefficients) # Good enough!
data.frame(true = X_coefs %>% as.vector, 
           gold = t(fit_goldstd$coefficients) %>% as.vector,
           geeind = fit_geeind$coefficients_mat %>% as.vector, 
           geeind_glmnet = fit_geeind_glmnet$coefficients_mat %>% as.vector, 
           geeunstr = fit_geeunstr$coefficients_mat %>% as.vector, 
           geeunstr_glmnet = fit_geeunstr_glmnet$coefficients_mat %>% as.vector, 
           geeRR = fit_geeRR$coefficients_mat %>% as.vector, 
           geeRR_glmnet = fit_geeRR_glmnet$coefficients_mat %>% as.vector 
           ) %>% 
   t %>% 
   dist


#' ## Negative binomial responses 
simy <- create_life(family = nb2(), X = X, X_coefs = X_coefs, loadings = true_loadings, disp_param = true_disp_param)
X <- as.data.frame(X)
fit_goldstd <- manyglm(simy ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data = X, family = "negative.binomial")
fit_geeind <- unpenalized_gee(y = simy, X = X, family = nb2(), control = list(trace = TRUE, step_halving = TRUE))
fit_geeunstr <- unpenalized_gee(y = simy, X = X, family = nb2(), control = list(trace = TRUE, step_halving = TRUE), corstr = "unstructured")
fit_geeRR <- unpenalized_gee(y = simy, X = X, family = nb2(), control = list(trace = TRUE, step_halving = TRUE), corstr = "reducedrank", rank = 5)
fit_geeind_glmnet <- unpenalized_gee(y = simy, X = X, family = nb2(), control = list(trace = TRUE), method = "glmnet")
fit_geeunstr_glmnet <- unpenalized_gee(y = simy, X = X, family = nb2(), control = list(trace = TRUE), corstr = "unstructured", method = "glmnet")
fit_geeRR_glmnet <- unpenalized_gee(y = simy, X = X, family = nb2(), control = list(trace = TRUE), corstr = "reducedrank", rank = 5, method = "glmnet")
# fit_geeind$coefficients_mat - t(fit_goldstd$coefficients) # Good enough? I think not estimating the dispersion parameter well is hurting
# fit_geeunstr$coefficients_mat - t(fit_goldstd$coefficients) # Good enough! I think not estimating the dispersion parameter well is hurting
# fit_geeRR$coefficients_mat - t(fit_goldstd$coefficients) # Good enough! I think not estimating the dispersion parameter well is hurting
data.frame(true = X_coefs %>% as.vector, 
           gold = t(fit_goldstd$coefficients) %>% as.vector,
           geeind = fit_geeind$coefficients_mat %>% as.vector, 
           geeind_glmnet = fit_geeind_glmnet$coefficients_mat %>% as.vector, 
           geeunstr = fit_geeunstr$coefficients_mat %>% as.vector, 
           geeunstr_glmnet = fit_geeunstr_glmnet$coefficients_mat %>% as.vector, 
           geeRR = fit_geeRR$coefficients_mat %>% as.vector, 
           geeRR_glmnet = fit_geeRR_glmnet$coefficients_mat %>% as.vector 
           ) %>% 
   t %>% 
   dist


#' ## Tweedie responses 
simy <- create_life(family = tweedielogfam(), X = X, X_coefs = X_coefs, loadings = true_loadings, 
                    disp_param = true_disp_param, power_param = runif(num_resp, 1.05, 1.95))
X <- as.data.frame(X)
fit_goldstd <- lapply(1:num_resp, function(x) gam(simy[,x] ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data = X, family = tw()))
fit_geeind <- unpenalized_gee(y = simy, X = X, family = tweedielogfam(), control = list(trace = TRUE))
fit_geeunstr <- unpenalized_gee(y = simy, X = X, family = tweedielogfam(), control = list(trace = TRUE), corstr = "unstructured")
fit_geeRR <- unpenalized_gee(y = simy, X = X, family = tweedielogfam(), control = list(trace = TRUE), corstr = "reducedrank", rank = 5)
fit_geeind_glmnet <- unpenalized_gee(y = simy, X = X, family = tweedielogfam(), control = list(trace = TRUE), method = "glmnet")
fit_geeunstr_glmnet <- unpenalized_gee(y = simy, X = X, family = tweedielogfam(), control = list(trace = TRUE), corstr = "unstructured", method = "glmnet")
fit_geeRR_glmnet <- unpenalized_gee(y = simy, X = X, family = tweedielogfam(), control = list(trace = TRUE), corstr = "reducedrank", rank = 5, method = "glmnet")
# fit_geeind$coefficients_mat - t(sapply(fit_goldstd, coef)) # Good enough? I think not estimating the dispersion parameter well is hurting
# fit_geeunstr$coefficients_mat - t(sapply(fit_goldstd, coef)) # Performance starts to hurt a bit more here
# fit_geeRR$coefficients_mat - t(sapply(fit_goldstd, coef)) # Performance starts to hurt a bit more here
data.frame(true = X_coefs %>% as.vector, 
           gold = t(sapply(fit_goldstd, coef)) %>% as.vector,
           geeind = fit_geeind$coefficients_mat %>% as.vector, 
           geeind_glmnet = fit_geeind_glmnet$coefficients_mat %>% as.vector, 
           geeunstr = fit_geeunstr$coefficients_mat %>% as.vector, 
           geeunstr_glmnet = fit_geeunstr_glmnet$coefficients_mat %>% as.vector, 
           geeRR = fit_geeRR$coefficients_mat %>% as.vector, 
           geeRR_glmnet = fit_geeRR_glmnet$coefficients_mat %>% as.vector 
           ) %>% 
   t %>% 
   dist
# fit_gee$power_param
# sapply(fit_goldstd, function(x) x$family$getTheta(TRUE))
# fit_gee$disp_param
# sapply(fit_goldstd, function(x) x$scale)


#' ## Beta responses 
simy <- create_life(family = betalogitfam(), X = X, X_coefs = X_coefs, loadings = true_loadings, disp_param = true_disp_param)
X <- as.data.frame(X)
fit_goldstd <- lapply(1:num_resp, function(x) mgcv::gam(simy[,x] ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data = X, family = betar()))
fit_geeind <- unpenalized_gee(y = simy, X = X, family = betalogitfam(), control = list(trace = TRUE))
fit_geeunstr <- unpenalized_gee(y = simy, X = X, family = betalogitfam(), control = list(trace = TRUE), corstr = "unstructured")
fit_geeRR <- unpenalized_gee(y = simy, X = X, family = betalogitfam(), control = list(trace = TRUE), corstr = "reducedrank", rank = 5)
fit_geeind_glmnet <- unpenalized_gee(y = simy, X = X, family = betalogitfam(), control = list(trace = TRUE), method = "glmnet")
fit_geeunstr_glmnet <- unpenalized_gee(y = simy, X = X, family = betalogitfam(), control = list(trace = TRUE), corstr = "unstructured", method = "glmnet")
fit_geeRR_glmnet <- unpenalized_gee(y = simy, X = X, family = betalogitfam(), control = list(trace = TRUE), corstr = "reducedrank", rank = 5, method = "glmnet")
data.frame(true = X_coefs %>% as.vector, 
           gold = t(sapply(fit_goldstd, coef)) %>% as.vector,
           geeind = fit_geeind$coefficients_mat %>% as.vector, 
           geeind_glmnet = fit_geeind_glmnet$coefficients_mat %>% as.vector, 
           geeunstr = fit_geeunstr$coefficients_mat %>% as.vector, 
           geeunstr_glmnet = fit_geeunstr_glmnet$coefficients_mat %>% as.vector, 
           geeRR = fit_geeRR$coefficients_mat %>% as.vector, 
           geeRR_glmnet = fit_geeRR_glmnet$coefficients_mat %>% as.vector 
           ) %>% 
   t %>% 
   dist
# fit_gee$disp_param
# sapply(fit_goldstd, function(x) x$family$getTheta(TRUE))


##--------------------
#' # Example 2a: Simulating independence Gaussian response data, but now add some homogeneity and sparsity in the regression coefficients for each covariate. Also start to test homogeneity pursuit in GEEs, given a specific value of lambda
#' 
#' 24 September 2022: SCAD/MCP no longer pursued as in additional preliminary work, it was found to be both slower and less stable (typically leading worst results) using ncvreg under the hood. By contrast adaptive lasso using glmnet continues to work well and be relatively efficient 
##---------------------
set.seed(092022)

num_units <- 400
num_resp <- 12
num_lv <- 2
num_X <- 9 # Excludes intercept
rho <- 0.5
P <- outer(1:num_X, 1:num_X, FUN = "-") %>% abs
P <- rho^P 

X <- rmvnorm(num_units, sigma = P) %>% cbind(1, .)
colnames(X) <- c("Int.", paste0("x",1:num_X))

X_coefs <- cbind(runif(num_resp, -2, 0), 
                 rnorm(num_resp,0,sd=0.25), # Unique values for each response
                 replicate(3, sample(c(-0.5,0,0.5),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 3 unique values, including zero
                 replicate(3, sample(c(-1,-0.5,0.5,1),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 4 unique non-zero values
                 replicate(2, numeric(num_resp))) 
# Note not covariate is such that it has the same value for all responses -- this is not really believable in ecology =P
true_loadings <- matrix(0, nrow = num_resp, ncol = num_lv)
true_disp_param <- runif(num_resp, 0, 4)
rm(P, rho, num_lv)


#' ## Simulate data and fit some unpenalized GEEs
simy <- create_life(family = gaussian(), X = X, X_coefs = X_coefs, loadings = true_loadings, disp_param = rep(1, num_resp))
X <- as.data.frame(X)
fit_goldstd <- manylm(simy ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data = X)
fit_geeind_glmnet <- unpenalized_gee(y = simy, X = X, family = gaussian(), control = list(trace = TRUE), method = "glmnet")
fit_geeunstr_glmnet <- unpenalized_gee(y = simy, X = X, family = gaussian(), control = list(trace = TRUE), corstr = "unstructured", method = "glmnet")
fit_geeRR_glmnet <- unpenalized_gee(y = simy, X = X, family = gaussian(), control = list(trace = TRUE), corstr = "reducedrank", rank = 5, method = "glmnet")


# ## Calculate and interpret the penalty weight matrix
# calc_penaltyweights(X_coefs = X_coefs, kappa = rep(c(0,1), c(1,num_X)))
# penweights_ADL <- calc_penaltyweights(X_coefs = fit_geeind_glmnet$coefficients_mat, kappa = rep(c(0,1), c(1,num_X)))
# penweights_ADL

# ## Calculate and check the reparametrized coefficients
RSmatrix <- .reparametrizationmatrix(X_coefs = X_coefs)
reparamcoefs <- (RSmatrix %*% as.vector(t(X_coefs))) %>% 
    matrix(., nrow = num_resp, byrow = TRUE) 
# reparamcoefs - penweights_ADL$reparametrized_coefficients # Good!
# rm(reparamcoefs) 
# 
originalcoefs <- solve(RSmatrix) %*% as.vector(t(reparamcoefs)) %>% 
    matrix(., nrow = num_resp, byrow = TRUE)
# originalcoefs - fit_geeind_glmnet$coefficients_mat # Good!
# rm(originalcoefs)


#' ## Start testing hpgee
fit_hpgee_path <- hpgee_regpath(y = simy, X = X, family = gaussian(), corstr = "independence", multicore_path = TRUE)

fit_hpgee_BIC <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$BIC)], 
                           family = gaussian(), corstr = "independence")
fit_hpgee_BIClogNloglogp <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$BIC_logNloglogp)], 
                           family = gaussian(), corstr = "independence")
fit_hpgee_EBIC1 <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$EBIC1)], 
                           family = gaussian(), corstr = "independence")
fit_hpgee_EBIC2 <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$EBIC2)], 
                           family = gaussian(), corstr = "independence")
fit_hpgee_ERIC <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$ERIC1)], 
                           family = gaussian(), corstr = "independence")


data.frame(true = X_coefs %>% as.vector, 
           gold = t(fit_goldstd$coefficients) %>% as.vector,
           geeind = fit_geeind_glmnet$coefficients_mat %>% as.vector, 
           geeunstr = fit_geeunstr_glmnet$coefficients_mat %>% as.vector, 
           geeRR = fit_geeRR_glmnet$coefficients_mat %>% as.vector,
           hpgee_BIC = fit_hpgee_BIC$coefficients_mat %>% as.vector,
           hpgee_BIClogNloglogp = fit_hpgee_BIClogNloglogp$coefficients_mat %>% as.vector,
           hpgee_EBIC1 = fit_hpgee_EBIC1$coefficients_mat %>% as.vector,
           hpgee_EBIC2 = fit_hpgee_EBIC2$coefficients_mat %>% as.vector,
           hpgee_ERIC = fit_hpgee_ERIC$coefficients_mat %>% as.vector
           ) %>% 
   t %>% 
   dist


##--------------------
#' # Example 2b: Simulating independence binary response data, but now add some homogeneity and sparsity in the regression coefficients for each covariate. Also start to test homogeneity pursuit in GEEs, given a specific value of lambda
##---------------------
set.seed(092022)

num_units <- 400
num_resp <- 12
num_lv <- 2
num_X <- 9 # Excludes intercept
rho <- 0.5
P <- outer(1:num_X, 1:num_X, FUN = "-") %>% abs
P <- rho^P 

X <- rmvnorm(num_units, sigma = P) %>% cbind(1, .)
colnames(X) <- c("Int.", paste0("x",1:num_X))

X_coefs <- cbind(runif(num_resp, -2, 0), 
                 rnorm(num_resp,0,sd=0.25), # Unique values for each response
                 replicate(3, sample(c(-0.5,0,0.5),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 3 unique values, including zero
                 replicate(3, sample(c(-1,-0.5,0.5,1),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 4 unique non-zero values
                 replicate(2, numeric(num_resp))) 
# Note not covariate is such that it has the same value for all responses -- this is not really believable in ecology =P
true_loadings <- matrix(0, nrow = num_resp, ncol = num_lv)
true_disp_param <- runif(num_resp, 0, 4)
rm(P, rho, num_lv)


#' ## Simulate data and fit some unpenalized GEEs
simy <- create_life(family = binomial(link = "probit"), X = X, X_coefs = X_coefs, loadings = true_loadings)
X <- as.data.frame(X)
fit_goldstd <- lapply(1:num_resp, function(x) mgcv::gam(simy[,x] ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data = X, family = binomial(link = "probit")))
fit_geeind_glmnet <- unpenalized_gee(y = simy, X = X, family = binomial(link = "probit"), control = list(trace = TRUE), method = "glmnet")
fit_geeunstr_glmnet <- unpenalized_gee(y = simy, X = X, family = binomial(link = "probit"), control = list(trace = TRUE), corstr = "unstructured", method = "glmnet")
fit_geeRR_glmnet <- unpenalized_gee(y = simy, X = X, family = binomial(link = "probit"), control = list(trace = TRUE), corstr = "reducedrank", rank = 5, method = "glmnet")


# ## Calculate and interpret the penalty weight matrix
# calc_penaltyweights(X_coefs = X_coefs, kappa = rep(c(0,1), c(1,num_X)))$penalty
# penweights_ADL <- calc_penaltyweights(X_coefs = fit_geeind_glmnet$coefficients_mat, kappa = rep(c(0,1), c(1,num_X)))
# penweights_ADL

# ## Calculate and check the reparametrized coefficients
# RSmatrix <- .reparametrizationmatrix(X_coefs = fit_geeind_glmnet$coefficients_mat)
# reparamcoefs <- (RSmatrix %*% as.vector(t(fit_geeind_glmnet$coefficients_mat))) %>% 
#    matrix(., nrow = num_resp, byrow = TRUE) 
# reparamcoefs - penweights_ADL$reparametrized_coefficients # Good!
# rm(reparamcoefs) 
# 
# originalcoefs <- solve(RSmatrix) %*% as.vector(t(penweights_ADL$reparametrized_coefficients)) %>% 
#    matrix(., nrow = num_resp, byrow = TRUE)
# originalcoefs - fit_geeind_glmnet$coefficients_mat # Good!
# rm(originalcoefs)


#' ## Start testing hpgee
fit_hpgee_path <- hpgee_regpath(y = simy, X = X, family = binomial(link = "probit"), corstr = "unstructured", multicore_path = TRUE)

fit_hpgee_BIC <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$BIC)], 
                           family = binomial(link = "probit"), corstr = "unstructured")
fit_hpgee_BIClogNloglogp <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$BIC_logNloglogp)], 
                           family = binomial(link = "probit"), corstr = "unstructured")
fit_hpgee_EBIC1 <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$EBIC1)], 
                           family = binomial(link = "probit"), corstr = "unstructured")
fit_hpgee_EBIC2 <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$EBIC2)], 
                           family = binomial(link = "probit"), corstr = "unstructured")
fit_hpgee_ERIC <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$ERIC1)], 
                           family = binomial(link = "probit"), corstr = "unstructured")


data.frame(true = X_coefs %>% as.vector, 
           gold = t(sapply(fit_goldstd, coef)) %>% as.vector,
           geeind = fit_geeind_glmnet$coefficients_mat %>% as.vector, 
           geeunstr = fit_geeunstr_glmnet$coefficients_mat %>% as.vector, 
           geeRR = fit_geeRR_glmnet$coefficients_mat %>% as.vector,
           hpgee_BIC = fit_hpgee_BIC$coefficients_mat %>% as.vector,
           hpgee_BIClogNloglogp = fit_hpgee_BIClogNloglogp$coefficients_mat %>% as.vector,
           hpgee_EBIC1 = fit_hpgee_EBIC1$coefficients_mat %>% as.vector,
           hpgee_EBIC2 = fit_hpgee_EBIC2$coefficients_mat %>% as.vector,
           hpgee_ERIC = fit_hpgee_ERIC$coefficients_mat %>% as.vector
           ) %>% 
   t %>% 
   dist

fit_hpgee_BIC$coefficients_matrix


##--------------------
#' # Example 2c: Simulating independence Poisson response data, but now add some homogeneity and sparsity in the regression coefficients for each covariate. Also start to test homogeneity pursuit in GEEs, given a specific value of lambda
##---------------------
set.seed(092022)

num_units <- 400
num_resp <- 12
num_lv <- 2
num_X <- 9 # Excludes intercept
rho <- 0.5
P <- outer(1:num_X, 1:num_X, FUN = "-") %>% abs
P <- rho^P 

X <- rmvnorm(num_units, sigma = P) %>% cbind(1, .)
colnames(X) <- c("Int.", paste0("x",1:num_X))

X_coefs <- cbind(runif(num_resp, -2, 0), 
                 rnorm(num_resp,0,sd=0.25), # Unique values for each response
                 replicate(3, sample(c(-0.5,0,0.5),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 3 unique values, including zero
                 replicate(3, sample(c(-1,-0.5,0.5,1),size=num_resp,replace=TRUE)), # 3 covariates where in each covariate there are at most 4 unique non-zero values
                 replicate(2, numeric(num_resp))) 
# Note not covariate is such that it has the same value for all responses -- this is not really believable in ecology =P
true_loadings <- matrix(0, nrow = num_resp, ncol = num_lv)
true_disp_param <- runif(num_resp, 0, 4)
rm(P, rho, num_lv)


#' ## Simulate data and fit some unpenalized GEEs
simy <- create_life(family = poisson(), X = X, X_coefs = X_coefs, loadings = true_loadings)
X <- as.data.frame(X)
fit_goldstd <- lapply(1:num_resp, function(x) mgcv::gam(simy[,x] ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data = X, family = poisson()))
fit_geeind_glmnet <- unpenalized_gee(y = simy, X = X, family = poisson(), control = list(trace = TRUE), method = "glmnet")
fit_geeunstr_glmnet <- unpenalized_gee(y = simy, X = X, family = poisson(), control = list(trace = TRUE), corstr = "unstructured", method = "glmnet")
fit_geeRR_glmnet <- unpenalized_gee(y = simy, X = X, family = poisson(), control = list(trace = TRUE), corstr = "reducedrank", rank = 5, method = "glmnet")


# ## Calculate and interpret the penalty weight matrix
# calc_penaltyweights(X_coefs = X_coefs, kappa = rep(c(0,1), c(1,num_X)))$penalty
# penweights_ADL <- calc_penaltyweights(X_coefs = fit_geeind_glmnet$coefficients_mat, kappa = rep(c(0,1), c(1,num_X)))
# penweights_ADL

# ## Calculate and check the reparametrized coefficients
# RSmatrix <- .reparametrizationmatrix(X_coefs = fit_geeind_glmnet$coefficients_mat)
# reparamcoefs <- (RSmatrix %*% as.vector(t(fit_geeind_glmnet$coefficients_mat))) %>% 
#    matrix(., nrow = num_resp, byrow = TRUE) 
# reparamcoefs - penweights_ADL$reparametrized_coefficients # Good!
# rm(reparamcoefs) 
# 
# originalcoefs <- solve(RSmatrix) %*% as.vector(t(penweights_ADL$reparametrized_coefficients)) %>% 
#    matrix(., nrow = num_resp, byrow = TRUE)
# originalcoefs - fit_geeind_glmnet$coefficients_mat # Good!
# rm(originalcoefs)


#' ## Start testing hpgee
fit_hpgee_path <- hpgee_regpath(y = simy, X = X, family = poisson(), corstr = "reducedrank", rank = 5, multicore_path = TRUE)

fit_hpgee_BIC <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$BIC)], 
                           family = poisson(), corstr = "reducedrank", rank = 5)
fit_hpgee_BIClogNloglogp <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$BIC_logNloglogp)], 
                           family = poisson(), corstr = "reducedrank", rank = 5)
fit_hpgee_EBIC1 <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$EBIC1)], 
                           family = poisson(), corstr = "reducedrank", rank = 5)
fit_hpgee_EBIC2 <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$EBIC2)], 
                           family = poisson(), corstr = "reducedrank", rank = 5)
fit_hpgee_ERIC <- hpgee_fit(hpgee_regpath = fit_hpgee_path, lambda = fit_hpgee_path$lambda[which.min(fit_hpgee_path$ICs$ERIC1)], 
                           family = poisson(), corstr = "reducedrank", rank = 5)


data.frame(true = X_coefs %>% as.vector, 
           gold = t(sapply(fit_goldstd, coef)) %>% as.vector,
           geeind = fit_geeind_glmnet$coefficients_mat %>% as.vector, 
           geeunstr = fit_geeunstr_glmnet$coefficients_mat %>% as.vector, 
           geeRR = fit_geeRR_glmnet$coefficients_mat %>% as.vector,
           hpgee_BIC = fit_hpgee_BIC$coefficients_mat %>% as.vector,
           hpgee_BIClogNloglogp = fit_hpgee_BIClogNloglogp$coefficients_mat %>% as.vector,
           hpgee_EBIC1 = fit_hpgee_EBIC1$coefficients_mat %>% as.vector,
           hpgee_EBIC2 = fit_hpgee_EBIC2$coefficients_mat %>% as.vector,
           hpgee_ERIC = fit_hpgee_ERIC$coefficients_mat %>% as.vector
           ) %>% 
   t %>% 
   dist

