#' ---
#' title: Applying HPGEEs to Great Barrier Reef (GBR) dataset; see Section 5 of the associated manuscript "Simultaneous homogeneity pursuit and variable selection in regression models for multivariate abundance data" for more details
#' abstract: The original GBR data is sourced from Pichler et al. 2007 <http://www.frdc.com.au/Archived-Reports/FRDC%20Projects/2003-021-DLD.pdf>. We provide an example dataset which possess the same structure as the actual data used for analysis in the manuscript, but the values are altered to mask their original values.
#' author: Originally written by FKCH
#' date: "Code started Oct 2022"
#' ---

##---------------------
#' # Load appropriate packages and files
##---------------------
rm(list = ls())
library(tidyverse)
library(mvtnorm)
library(mvabund)
library(mgcv)
library(blockCV)
library(patchwork)
library(cluster)
library(sf)
library(ROCR)

library(colorspace)
library(plotly)
library(withr)
#library(htmlwidgets)

here::i_am("application_GBR/application.R")
library(here)

source(here("code","auxfnsv0.R"))
source(here("code","calc_penaltyweights.R"))
source(here("code","create_life.R"))
source(here("code","familyfns.R"))
source(here("code","geefn.R"))
source(here("code","hpgeefn.R"))
source(here("code","updatecoeffns.R"))
source(here("code","ADLgeefn.R"))


##--------------------
#' # Load in data
##---------------------
load(file = here("application_GBR", "gbr_exampledata.RData"))
dim(resp_dat)
colSums(resp_dat)

dim(environ_dat)
summary(environ_dat)

# Scale covariates
environ_dat[,-c(1:2)] <- scale(environ_dat[,-c(1:2)]) %>% 
   as.data.frame

summary(environ_dat)



##-------------------
#' # Form spatial blocks in preparation for cross-validation
##-------------------
num_resp <- ncol(resp_dat)
num_units <- nrow(resp_dat)

longlat_dat <- environ_dat %>% 
   dplyr::select(GRD_LON:GRD_LAT)

pa_data <- st_as_sf(longlat_dat, coords = c("GRD_LON", "GRD_LAT"), crs = 4326)

sb <- spatialBlock(speciesData = pa_data,
                   theRange = 150000, # size of the blocks
                   k = 5,
                   selection = "random",
                   seed = 102022,
                   progress = TRUE)

sb$plots + geom_sf(data = pa_data, alpha = 0.5, color = "grey")


##-------------------
#' # Perform five-fold cross-validation
##---------------------
X <- model.matrix(~ GBR_BATHY + GBR_SLOPE + GBR_ASPECT + GBR_TS_BSTRESS + GA_GRAVEL + GA_MUD + + GA_CRBNT + CRS_O2_AV + CRS_T_AV + SW_CHLA_AV, data = environ_dat) %>% 
     as.data.frame() 
summary(X)


all_predictions <- vector("list", 5)

for(k0 in 1:5) {
   message("Starting fold ", k0)
   train_index <- sb$folds[[k0]][[1]]
   test_index <- sb$folds[[k0]][[2]]
   train_resp_dat <- resp_dat[train_index,]
   train_X <- X[train_index,]
   test_X <- X[test_index,] %>% as.matrix
   
   # HPGEE with independence working correlation 
   fit_hpgee <- hpgee_regpath(y = train_resp_dat, X = train_X, family = binomial(link = "probit"), corstr = "independence", rank = 3, multicore_path = TRUE, num_cores = 5)
   fit_hpgee_IC <- hpgee_fit(hpgee_regpath = fit_hpgee, lambda = fit_hpgee$lambda[which.min(fit_hpgee$ICs$ERIC)], family = binomial(link = "probit"), corstr = "independence", rank = 3)
   all_predictions[[k0]]$HPGEE_ERIC <- pnorm(tcrossprod(test_X, fit_hpgee_IC$coefficients_matrix))

   
   # HPGEE with reduced rank working correlation 
   fit_hpgee <- hpgee_regpath(y = train_resp_dat, X = train_X, family = binomial(link = "probit"), corstr = "reducedrank", rank = 3, multicore_path = TRUE, num_cores = 5)
   fit_hpgee_IC <- hpgee_fit(hpgee_regpath = fit_hpgee, lambda = fit_hpgee$lambda[which.min(fit_hpgee$ICs$ERIC)], family = binomial(link = "probit"), corstr = "reducedrank", rank = 3)
   all_predictions[[k0]]$HPGEE_RR_ERIC <- pnorm(tcrossprod(test_X, fit_hpgee_IC$coefficients_matrix))

      
   ## Adaptive lasso GEEs with independence working correlation 
   fit_adlgee <- adlgee_regpath(y = train_resp_dat, X = train_X, family = binomial(link = "probit"), corstr = "independence", multicore_path = TRUE)
   fit_adlgee_IC <- adlgee_fit(adlgee_regpath = fit_adlgee, lambda = fit_adlgee$lambda[which.min(fit_adlgee$ICs$ERIC)], family = binomial(link = "probit"), corstr = "independence")
   all_other_predictions[[k0]]$ADLGEE_RR_ERIC <- pnorm(tcrossprod(test_X, fit_adlgee_IC$coefficients_matrix))
   
   
   ## Adaptive lasso GEEs with reduced-rank working correlation 
   fit_adlgee <- adlgee_regpath(y = train_resp_dat, X = train_X, family = binomial(link = "probit"), corstr = "reducedrank", rank = 3, multicore_path = TRUE)
   fit_adlgee_IC <- adlgee_fit(adlgee_regpath = fit_adlgee, lambda = fit_adlgee$lambda[which.min(fit_adlgee$ICs$ERIC)], family = binomial(link = "probit"), corstr = "reducedrank", rank = 3)
   all_other_predictions[[k0]]$ADLGEE_ERIC <- pnorm(tcrossprod(test_X, fit_adlgee_IC$coefficients_matrix))
   
   
   # Stacked penalized glms using glmnet, using 10-fold cross-validation with deviance as the loss function and selecting lambda.1se
   fit_glmnet <- list(coefficient_matrix = matrix(0, nrow = num_resp, ncol = ncol(train_X)))
   for(k1 in 1:num_resp) {
      find_cv <- cv.glmnet(y = train_resp_dat[,k1], x = as.matrix(train_X[,-1,drop=FALSE]), parallel = TRUE, type.measure = "default", family = binomial(link = "probit"))
      fit_glmnet$coefficient_matrix[k1,] <- as.vector(coef(find_cv, s = "lambda.1se", exact = TRUE))
      }
   all_other_predictions[[k0]]$glmnet <- pnorm(tcrossprod(test_X, fit_glmnet$coefficient_matrix)) 
   
   
   # Unpenalized GEEs 
   fit_geeind <- unpenalized_gee(y = train_resp_dat, X = train_X, family = binomial(link = "probit"), corstr = "independence")
   fit_geeRR <- unpenalized_gee(y = train_resp_dat, X = train_X, family = binomial(link = "probit"), corstr = "reducedrank", rank = 3)
   all_predictions[[k0]]$GEE_ind <- pnorm(tcrossprod(test_X, fit_geeind$coefficients_matrix))
   all_predictions[[k0]]$GEE_RR <- pnorm(tcrossprod(test_X, fit_geeRR$coefficients_matrix))

   
   # Two-stage unpenalized GEEs followed by K-means clustering, assuming either independence or reduced-rank working correlation
   fit_geeind <- unpenalized_gee(y = train_resp_dat, X = train_X, family = binomial(link = "probit"), control = list(trace = TRUE))
   fit_kmeansind <- list(coefficient_matrix = fit_geeind$coefficients_matrix)
   registerDoParallel(cores = detectCores() - 1)
   clusteringcoefficient_fn <- function(l) {
      cw_coefficients <- matrix(fit_geeind$coefficients_matrix[,l], ncol = 1) 
      gap_stat <- clusGap(cw_coefficients, FUN = kmeans, nstart = 25, K.max = num_resp-1, verbose = FALSE)
      cw_coefficients_cluster <- kmeans(cw_coefficients, centers = which.max(gap_stat$Tab[,3]), nstart = 25, iter.max = 50)
      out <- cw_coefficients_cluster$centers[cw_coefficients_cluster$cluster]
      return(out)
      }
   fit_kmeansind$coefficient_matrix[,-1] <- foreach(s0 = 2:ncol(X), .combine = "cbind") %dopar% clusteringcoefficient_fn(l = s0)
   all_other_predictions[[k0]]$kmeans_ind <- pnorm(tcrossprod(test_X, fit_kmeansind$coefficient_matrix))
   
   fit_geeRR <- unpenalized_gee(y = train_resp_dat, X = train_X, family = binomial(link = "probit"), control = list(trace = TRUE), corstr = "reducedrank", rank = 3)
   fit_kmeansRR <- list(coefficient_matrix = fit_geeRR$coefficients_matrix)
   registerDoParallel(cores = detectCores() - 1)
   clusteringcoefficient_fn <- function(l) {
      cw_coefficients <- matrix(fit_geeRR$coefficients_matrix[,l], ncol = 1) 
      gap_stat <- clusGap(cw_coefficients, FUN = kmeans, nstart = 25, K.max = num_resp-1, verbose = FALSE)
      cw_coefficients_cluster <- kmeans(cw_coefficients, centers = which.max(gap_stat$Tab[,3]), nstart = 25, iter.max = 50)
      out <- cw_coefficients_cluster$centers[cw_coefficients_cluster$cluster]
      return(out)
      }
   fit_kmeansRR$coefficient_matrix[,-1] <- foreach(s0 = 2:ncol(X), .combine = "cbind") %dopar% clusteringcoefficient_fn(l = s0)
   all_other_predictions[[k0]]$kmeans_RR <- pnorm(tcrossprod(test_X, fit_kmeansRR$coefficient_matrix))
   }

save.image(file = "GBR_spatialCVpredictions.RData")


##-------------------
#' # Assess cross-validation performance 
##-------------------
load(file = "GBR_spatialCVpredictions.RData")

all_predictions %>% 
   str

all_RMSE <- all_AUC <- all_mce <- array(NA, dim = c(num_resp, length(all_predictions[[1]]), 5), 
                                dimnames = list(spp = 1:num_resp, methods = names(all_predictions[[1]]), fold = 1:5) )

for(k0 in 1:5) {
   test_index <- sb$folds[[k0]][[2]]
   test_resp_dat <- resp_dat[test_index,]
   
   all_RMSE[,1:length(all_predictions[[1]]),k0] <- sapply(1:length(all_predictions[[k0]]), function(k1) {
      sqrt(colMeans((test_resp_dat - all_predictions[[k0]][[k1]])^2))
      })
   
   all_AUC[,1:length(all_predictions[[1]]),k0] <- sapply(1:length(all_predictions[[k0]]), function(k1) {
      sapply(1:num_resp, function(j) { pred <- prediction(all_predictions[[k0]][[k1]][,j], test_resp_dat[,j]); performance(pred, "auc")@y.values[[1]] } ) 
      })

   all_mce[,1:length(all_predictions[[1]]),k0] <- sapply(1:length(all_predictions[[k0]]), function(k1) {
      sapply(1:num_resp, function(j) { pred <- prediction(all_predictions[[k0]][[k1]][,j], test_resp_dat[,j]); performance(pred, "mxe")@y.values[[1]] } ) 
      })
   }



sel_methods <- c("HPGEE_RR_ERIC", "GEE_RR", "glmnet", "kmeans_RR", "ADLGEE_RR_ERIC")

RMSE_results <- abind::abind(lapply(1:5, function(k) all_RMSE[,,k]-all_RMSE[,"HPGEE_RR_ERIC",k]), along = 3) %>% 
   apply(., c(1,2), mean) %>% 
   as.data.frame %>% 
   rownames_to_column(var = "species") %>% 
   dplyr::select(-GEE_ind) %>% 
   pivot_longer(-species, names_to = "method") %>% 
   dplyr::filter(method %in% sel_methods) %>% 
   mutate(method = fct_recode(method, HPGEE = "HPGEE_RR_ERIC", GEE = "GEE_RR", "K-means" = "kmeans_RR", ADLGEE = "ADLGEE_RR_ERIC")) %>% 
   mutate(method = fct_relevel(method, "HPGEE", "ADLGEE", "glmnet", "GEE", "K-means")) 

ggplot(RMSE_results %>% dplyr::filter(method != "HPGEE"), aes(x = method, y = value)) +
   geom_boxplot() +
   geom_hline(yintercept = 0, linetype = 2) +
   labs(y = "Difference in RMSE (relative to HPGEE)", x = "Method", title = "RMSE") +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))


crossentropy_results <- abind::abind(lapply(1:5, function(k) all_mce[,,k]-all_mce[,"HPGEE_RR_ERIC",k]), along = 3) %>% 
   apply(., c(1,2), mean) %>% 
   as.data.frame %>% 
   rownames_to_column(var = "species") %>% 
   dplyr::select(-GEE_ind) %>% 
   pivot_longer(-species, names_to = "method") %>% 
   dplyr::filter(method %in% sel_methods) %>% 
   mutate(method = fct_recode(method, HPGEE = "HPGEE_RR_ERIC", GEE = "GEE_RR", "K-means" = "kmeans_RR", ADLGEE = "ADLGEE_RR_ERIC")) %>% 
   mutate(method = fct_relevel(method, "HPGEE", "ADLGEE", "glmnet", "GEE", "K-means")) 

ggplot(crossentropy_results %>% dplyr::filter(method != "HPGEE"), aes(x = method, y = value)) +
   geom_boxplot() +
   geom_hline(yintercept = 0, linetype = 2) +
   labs(y = "Difference in MCE (relative to HPGEE)", x = "Method", title = "Mean cross-entropy") +
   scale_y_continuous(limits = c(-0.025, 0.025)) +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))


AUC_results <- abind::abind(lapply(1:5, function(k) all_AUC[,,k]-all_AUC[,"HPGEE_RR_ERIC",k]), along = 3) %>% 
   apply(., c(1,2), mean) %>% 
   as.data.frame %>% 
   rownames_to_column(var = "species") %>% 
   dplyr::select(-GEE_ind) %>% 
   pivot_longer(-species, names_to = "method") %>% 
   dplyr::filter(method %in% sel_methods) %>% 
   mutate(method = fct_recode(method, HPGEE = "HPGEE_RR_ERIC", GEE = "GEE_RR", "K-means" = "kmeans_RR", ADLGEE = "ADLGEE_RR_ERIC")) %>% 
   mutate(method = fct_relevel(method, "HPGEE", "ADLGEE", "glmnet", "GEE", "K-means")) 

ggplot(AUC_results %>% dplyr::filter(method != "HPGEE"), aes(x = method, y = value)) +
   geom_boxplot() +
   geom_hline(yintercept = 0, linetype = 2) +
   labs(y = "Difference in AUC (relative to HPGEE)", x = "Method", title = "AUC") +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))


##-------------------
#' # Fit HPGEEs to full dataset, and examine results
##-------------------
X <- model.matrix(~ GBR_BATHY + GBR_SLOPE + GBR_ASPECT + GBR_TS_BSTRESS + GA_GRAVEL + GA_MUD + + GA_CRBNT + CRS_O2_AV + CRS_T_AV + SW_CHLA_AV, data = environ_dat) %>% 
     as.data.frame()
cor(X[,-1]) %>% round(3)
summary(X)


# HPGEE with reduced rank working correlation
fit_hpgee_RR <- hpgee_regpath(y = resp_dat, X = X, family = binomial(link = "probit"), corstr = "reducedrank", rank = 3, multicore_path = TRUE)
fit_hpgee_RR_ERIC <- hpgee_fit(hpgee_regpath = fit_hpgee_RR, lambda = fit_hpgee_RR$lambda[which.min(fit_hpgee_RR$ICs$ERIC)], family = binomial(link = "probit"), corstr = "reducedrank", rank = 3)


getcoef_mat <- fit_hpgee_RR_ERIC$coefficients_matrix
get_num_unique_values <- apply(getcoef_mat, 2, function(x) length(unique(x)))
covariate_names <- paste0(c("Intercept", "Bathymetry", "Slope", "Aspect", "Bottom stress", "Gravel", "Mud", "Carbonate", "Oxygen", "Temperature", "Chlorophyll-a"), " (", get_num_unique_values, ")")
colnames(getcoef_mat) <- covariate_names
abbreviated_spp_names <- colnames(resp_dat)
spp_scientific_names <- c(
   "Scyllarus_martensii", 
   "Parthenope_longispinus", 
   "Lophiotoma_acuta", 
   "Leucosia_ocellata", 
   "Lovenia_elongata", 
   "Xenophora_solarioides", 
   "Lomopsis_sp",
   "Luidia_hardwicki", 
   "Leionucula_superba", 
   "Murex_tenuirostrum", 
   "Myrine_kesslerii", 
   "Sorsogona_tuberculata", 
   "Portunus_gracilimanus", 
   "Amusium_pleuronectes", 
   "Ophiacantha_indica", 
   "Temnopleuridae_sp", 
   "Hippospongia_elastica", 
   "Brissopsis_luzonica", 
   "Myra_tumidospina", 
   "Charybdis_truncata")

spp_common_names <- c( # As best as possible via Google searching!
   "Striated locust lobster", 
   "White long-armed crab", 
   "Marbled turrid", 
   "Pebble crab sp A", 
   "Heart urchin", 
   "Minute carrier shell", 
   "Bivalve sp A",
   "Seastar", 
   "Saltwater clam sp A", 
   "Sea snail sp A", 
   "Pebble crab sp B", 
   "Tuberculate flathead", 
   "Swimming crab", 
   "Bivalve sp B", 
   "Brittle stars", 
   "Sea urchin sp A", 
   "Horny sponges", 
   "Sea urchin sp B", 
   "Crab sp A", 
   "Blunt-toothed crab")

rownames(getcoef_mat) <- spp_common_names

getcoef_mat


ggplot(getcoef_mat[,-1] %>% 
          as.data.frame() %>% 
          rownames_to_column(var = "species") %>% 
          pivot_longer(-species) %>% 
          mutate(species = fct_inorder(species)) %>% 
          mutate(species = fct_relabel(species, ~ gsub("\\.", " ", .x))) %>% 
          mutate(iszero = 1 + 1*(value == 0)) %>% 
          mutate(name = fct_inorder(name)), 
       aes(x = name, y = species, fill = factor(value))) + 
   geom_tile(color = "gray", show.legend = FALSE) +
   geom_text(aes(label = round(value,2), fontface = iszero), color = "black", show.legend = FALSE) + #color = factor(value)), 
   labs(x = "Covariate (number of unique estimates)", y = "Species", fill = "Estimate") +
   scale_fill_discrete_diverging() +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank())


apply(getcoef_mat, 2, unique)


getcoef_long <- getcoef_mat[,-1] %>% 
   as.data.frame() %>% 
   rownames_to_column(var = "Species") %>% 
   pivot_longer(-Species, names_to = "Covariate", values_to = "Estimate") 

for(k1 in 1:length(unique(getcoef_long$Covariate))) {
   aggregate(Species ~ Estimate, data = getcoef_long %>% dplyr::filter(Covariate == unique(getcoef_long$Covariate)[k1]), FUN = print) %>% 
      print
   }



#' # Construct a regularization path for each covariate, which effectively acts as tree showing the clustering taking place across species
covariate_names <- c("Intercept", "Bathymetry", "Slope", "Aspect", "Bottom stress", "Gravel", "Mud", "Carbonate", "Oxygen", "Temperature", "Chlorophyll-a")

all_coefficients_array <- array(NA, 
                                dim = c(num_resp, ncol(X), length(fit_hpgee_RR$lambda)),
                                dimnames = list(species = spp_common_names, covariates = covariate_names, lambda = 1:length(fit_hpgee_RR$lambda)))

for(k0 in 1:length(fit_hpgee$lambda)) {
   all_coefficients_array[,,k0] <- matrix(fit_hpgee_RR$coefficients_path[,k0], nrow = num_resp, byrow = TRUE)
   }


for(k0 in 2:length(covariate_names)) {
   cw_coefficients <- all_coefficients_array[,k0,] %>% 
      t %>% 
      as.data.frame %>% 
      mutate(lambda = fit_hpgee_RR$lambda) %>% 
      relocate(lambda) %>% 
      pivot_longer(-lambda, names_to = "species")
   
   p <- ggplot(cw_coefficients, aes(x = lambda, y = value, color = species)) +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_line() +
      geom_vline(xintercept = fit_hpgee_RR$lambda[which.min(fit_hpgee_RR$ICs$ERIC1)], col = "darkred", linetype = 2) +
      scale_x_log10() +
      scale_color_discrete_qualitative() +
      theme_bw() +
      labs(x = "Tuning parameter (lambda)", y = "Coefficient value", title = covariate_names[k0], color = "Species") +
         coord_flip() 
   
   with_options(
      list(digits = 4), ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.2))
      )
   }








##-------------------
sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Linux Mint 19
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#  [1] LC_CTYPE=en_AU.UTF-8       LC_NUMERIC=C               LC_TIME=en_AU.UTF-8        LC_COLLATE=en_AU.UTF-8     LC_MONETARY=en_AU.UTF-8   
#  [6] LC_MESSAGES=en_AU.UTF-8    LC_PAPER=en_AU.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_AU.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] tweedie_2.3.3        MASS_7.3-58.1        glmnet_4.1-4         Matrix_1.2-18        doParallel_1.0.17    iterators_1.0.14     foreach_1.5.2       
#  [8] here_1.0.1           withr_2.5.0          plotly_4.9.2.1       colorspace_2.0-3     ROCR_1.0-7           gplots_3.0.1.2       sf_0.9-5            
# [15] cluster_2.1.0        patchwork_1.0.1.9000 blockCV_2.1.4        mgcv_1.8-40          nlme_3.1-149         mvabund_4.2.1        mvtnorm_1.1-3       
# [22] forcats_0.5.2        stringr_1.4.0        dplyr_1.0.9          purrr_0.3.4          readr_2.1.2          tidyr_1.2.0          tibble_3.1.8        
# [29] ggplot2_3.4.0        tidyverse_1.3.2     
# 
# loaded via a namespace (and not attached):
#  [1] bitops_1.0-6        fs_1.5.0            lubridate_1.8.0     progress_1.2.2      httr_1.4.2          rprojroot_2.0.2     tools_3.6.3        
#  [8] backports_1.2.1     utf8_1.2.2          R6_2.5.1            KernSmooth_2.23-17  lazyeval_0.2.2      DBI_1.1.0           raster_3.4-5       
# [15] sp_1.5-0            tidyselect_1.1.2    prettyunits_1.1.1   compiler_3.6.3      cli_3.4.1           rvest_1.0.3         xml2_1.3.3         
# [22] caTools_1.17.1.4    scales_1.2.1        classInt_0.4-3      digest_0.6.29       pkgconfig_2.0.3     htmltools_0.5.1.1   dbplyr_2.2.1       
# [29] htmlwidgets_1.5.1   rlang_1.0.6         readxl_1.3.1        rstudioapi_0.13     shape_1.4.4         generics_0.1.3      jsonlite_1.7.2     
# [36] gtools_3.8.1        googlesheets4_1.0.1 magrittr_2.0.3      Rcpp_1.0.8.3        munsell_0.5.0       fansi_1.0.3         abind_1.4-5        
# [43] lifecycle_1.0.3     stringi_1.7.6       grid_3.6.3          gdata_2.18.0        crayon_1.5.1        lattice_0.20-41     haven_2.5.1        
# [50] splines_3.6.3       hms_1.1.2           pillar_1.8.1        codetools_0.2-16    reprex_2.0.2        glue_1.6.2          data.table_1.14.2  
# [57] modelr_0.1.9        vctrs_0.5.1         tzdb_0.3.0          cellranger_1.1.0    gtable_0.3.1        assertthat_0.2.1    broom_1.0.1        
# [64] e1071_1.7-3         survival_3.2-7      viridisLite_0.4.1   class_7.3-17        googledrive_2.0.0   gargle_1.2.1        units_0.6-7        
# [71] statmod_1.4.36      ellipsis_0.3.2     
