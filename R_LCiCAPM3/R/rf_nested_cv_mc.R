# Load libraries
library(caret)
library(SummarizedExperiment)
library(parallel)
library(doParallel)
library(tidyverse) #includes dplyr and purr
library(randomForest)
library(glmnet)

#set up parallel cluster unless already defined
start_cluster <- function(){
  cl <<- makePSOCKcluster(15)
  registerDoParallel(cl)
}
if(!exists("cl"))
  start_cluster()

#set seed for reproducibility
set.seed(1)

#set number of folds for cross validation
num_folds = 5


#' rf_nestedcv_stable
#'
#' Performs stability-based feature selection with nested cross-validation using a random forest classifier.
#' @param .ex SummarizedExperiment object
#' @param form Formula for classification
#' @param partition_var Column in colData indicating training partition
#' @param num_folds Number of outer folds
#' @param use_matrix Name of assay (default: log2_normalized_counts)
#' @param num_subsamples Number of subsamples for stability selection
#' @param stability_cutoff Proportion threshold for stable feature inclusion
#' @param seed Random seed
#' @param stability_sample_fraction Fraction of data to sample in each iteration
#' @param ntree Number of trees for Random Forest
#' @param preprocess_method List of preprocessing steps (center, scale, etc.)
#' @return A list containing accuracy, AUC, stable features, final model, etc.


rf_nestedcv_stable <- function(
    .ex, 
    form=status~., 
    partition_var='partition_new', 
    num_folds=5,
    use_matrix="log2_normalized_counts", 
    num_subsamples = 50, 
    stability_cutoff=0.75, 
    seed=1,
    stability_sample_fraction = 0.5,
    ntree = 100,
    preprocess_method = c("center", "scale")
  ){
  
  set.seed(seed)
  seeds <- sample(1:1000, num_subsamples)
  
  alpha_values <- seq(0.1, 1, by=0.1)
  
  class_var = as.character(form[[2]])
  
  ex <- .ex[, colData(.ex)[,partition_var] == 'train']
  
  x <- data.frame(t(assay(ex, use_matrix)))
  
  # filter zero variance columns
  x <- x[,apply(x, MARGIN=2, var)!=0]
  
  # remove highly correlated variables
  #z <- caret::findCorrelation(cor(x), cutoff = .9)
  #if(length(z) > 0)
  #  x <- x[,-z]
  
  # add status to train
  x[,class_var] <- factor(colData(ex)[,class_var], levels=c('benign', 'malignant'))
  x$smoking_status <- colData(ex)$smoking_status
  #tr$gender = ex_tr$sex
  #tr$age <- ex_tr$age_at_collection
  #tr$nodule_size <- ex_tr$mayo_nodule_size_mm
  data <- x
  # create model matrix
  #x <- model.matrix(form, x, NULL)
  #x <- data.frame(x[,-grepl("(Intercept)", colnames(x))])
  #x[,class_var] <- tr[,class_var]
  
  # create folds
  outer_folds <- createFolds(x[,class_var], k = num_folds, list = TRUE, returnTrain = TRUE)
  
  all_stable_features <- c()
  
  # outer loop
  # for debug
  #outer_train_indices = outer_folds[[1]]
  #outer_res <- lapply(outer_folds, function(outer_train_indices){
  
  outer_res <- list()
  
  for(jj in 1:length(outer_folds)){
    
    outer_train_indices = outer_folds[[jj]]
      
    outer_train <- x[outer_train_indices,]
    outer_test  <- x[-outer_train_indices,]
    
    n_cores = parallel::detectCores()-1
    
    feature_selection_results <- mclapply(mc.cores=n_cores, 1:num_subsamples, function(i){
      
      set.seed(seeds[i])
      # this is imbalanced
      #subsample <- sample_n(outer_train, size=nrow(outer_train)/2)
      
      # this corrects the imbalance
      #strata = c(class_var)
      strata = c(class_var, 'smoking_status')
      subsample <- outer_train |> 
        dplyr::mutate(n = dplyr::row_number()) |>
        dplyr::select(n, dplyr::everything()) |>
        dplyr::group_by(across(all_of(strata))) |>
        dplyr::sample_frac(stability_sample_fraction) |>
        dplyr::select(-n)
        
      #x_subsample <- model.matrix(~., subsample[,-which(names(subsample) == class_var)])
      x_subsample <- model.matrix(form, subsample)
      x_subsample <- x_subsample[,-1] # remove intercept 
      
      # filter zero var genes
      x_subsample <- x_subsample[,apply(x_subsample, MARGIN=2, var)!=0]

      # scale and center values
      # ** todo try including/excluding dummy variables (ie smoking status) then add back in
      preProcValues <- preProcess(x_subsample, method = c("center", "scale"))
      x_subsample <- predict(preProcValues, x_subsample)
    
      y_subsample <- as.factor(as.numeric(subsample[[class_var]]))
      
      ## inner loop
      cv_results <- list()
      for (alpha in alpha_values){
        cvfit <- cv.glmnet(x_subsample, y_subsample, family="binomial", alpha=alpha, nfolds=5, standardize=F)
        cv_results[[paste(alpha)]] <- cvfit
      }
      cv_error = sapply(cv_results, function(x) min(x$cvm))
      
      # select the best model based on lowest cv error
      best_alpha <- alpha_values[which.min(cv_error)]
      best_model <- cv_results[which.min(cv_error)][[1]]
      
      fit <- glmnet(x_subsample, y=y_subsample, alpha=best_alpha, lambda = best_model$lambda.min, family='binomial', standardize=F)
      
      # get coefficients
      selected_features = c()
      coefficients <- coef(fit)
      ind <- which(coefficients!=0)
      ind <- ind[ind > 1]
      if (length(ind) > 0)
        selected_features <- rownames(coefficients)[ind]
        
      #feature_selection_results[[i]] <- selected_features
      
      return( selected_features )
      
    })
    
    all_features <- unlist(feature_selection_results)
    stability_scores <- table(all_features) / num_subsamples
    
    stable_features <- names(stability_scores[stability_scores > quantile(stability_scores, stability_cutoff)])
    # Fixed proportion rule (most common interpretation)
    #stable_features <- names(stability_scores[stability_scores >= stability_cutoff])

    
    if(length(stable_features) == 0){
      message("no stable features found, try increasing num_subsamples or decreasing stability_cutoff")
      return(
        list(
          accuracy = NA,
          auc = NA,
          stable_features = NA,
          outer_res = NA,
          final_model = NA,
          preprocess_values = NA,
          error_code = 1
        )
      )
    }
    
    # pre process
    outer_train_p <- outer_train
    outer_test_p <- outer_test
    if ( ! is.null(preprocess_method) ){
      preProcValues <- preProcess(outer_train, method = preprocess_method)
    
      outer_train_p <- predict(preProcValues, outer_train)
      outer_test_p <- predict(preProcValues, outer_test)
    }
      
    # create model matrix
    outer_train_m <- model.matrix(form, outer_train_p)
    outer_train_m <- outer_train_m[,-1] # remove intercept 
    outer_train_m <- cbind(outer_train_m, outer_train_p[,class_var])
    colnames(outer_train_m)[ncol(outer_train_m)] = class_var

    
    outer_test_m <- model.matrix(form, outer_test_p)
    outer_test_m <- outer_test_m[,-1] # remove intercept 
    outer_test_m <- cbind(outer_test_m, outer_test_p[,class_var])
    colnames(outer_test_m)[ncol(outer_test_m)] = class_var
    
    # select just stable features
    outer_train_f <- outer_train_m[, c(stable_features, class_var)] |> as.data.frame()
    outer_test_f <- outer_test_m[, c(stable_features, class_var)] |> as.data.frame()
    
    # reassign factors for train/test
    #   sort levels for consistency
    for (name in colnames(outer_train_f)[colnames(outer_train_f) %in% colnames(colData(ex))])
      outer_train_f[,name] <- factor(outer_train_f[,name], levels = sort(unique(outer_train_f[,name])))
      outer_test_f[,name] <- factor(outer_test_f[,name], levels = sort(unique(outer_test_f[,name])))
    
    # rename columns for randomForest compatability
    outer_train_r <- outer_train_f
    outer_test_r <- outer_test_f
    colnames(outer_train_r)  <- colnames(outer_train_r) |> str_replace(":", ".")
    colnames(outer_test_r)  <- colnames(outer_test_r) |> str_replace(":", ".")
    
    model <- randomForest(status ~., data=outer_train_r, importance=T, ntree=ntree)
    pred <- predict(model, newdata = outer_test_r,type="prob")
    
    outer_res[[jj]] <- list(
      accuracy = confusionMatrix(
        factor(ifelse(pred[,2] > .5, 2, 1), levels=c(1,2)), 
        outer_test_f[,class_var]
      )$overall["Accuracy"],
      auc = pROC::roc(outer_test_f[,class_var], pred[,2], quiet=T)$auc[[1]],
      stable_features = stable_features,
      model = model
    )
    
    # return (
    #   list(
    #     accuracy = confusionMatrix(pred, outer_test_f[,class_var])$overall["Accuracy"],
    #     stable_features = stable_features,
    #     model = model
    #   )
    # )
    
  }
  # })
  
  # features from stability analysis
  stable_features <- unique(unlist( lapply(outer_res, FUN=function(x) x$stable_features)))
  
  # features for final model
  primary_terms <- stable_features |> str_split("\\:") |> unlist() |> unique()
  primary_terms <- primary_terms[primary_terms %in% rownames(ex)]
  clinical_features <- colnames(data)[! colnames(data) %in% rownames(ex)]
  
  tr <- data[,c(primary_terms, clinical_features)]
  
  # pre process
  preProcValues = NULL
  if ( ! is.null(preprocess_method) ){
    preProcValues <- preProcess(tr, method = preprocess_method)
    tr <- predict(preProcValues, tr)
  }
    
  # create model matrix
  tr <- model.matrix(form, tr)
  tr <- tr[,-1] |> as.data.frame()# remove intercept 
  tr[,class_var] <- colData(ex[,rownames(tr)])[,class_var]
  assertthat::assert_that(all(tr[,class_var] == data[,class_var]))
  
  tr <- tr[,c(stable_features, class_var)]
 
  colnames(tr) <- colnames(tr) |> str_replace("\\:", "\\.")
  
  final_model  <- randomForest(status ~., data=tr, importance=T, ntree=ntree)

  return(
    list(
      accuracy = mean(unlist(lapply(outer_res, FUN=function(x) x$accuracy))),
      auc = mean(unlist(lapply(outer_res, FUN=function(x) x$auc))),
      stable_features = stable_features,
      outer_res = outer_res,
      final_model = final_model,
      preprocess_values = preProcValues,
      # inside return(...)
      seed = seed,
      seeds_used = seeds, #jjs added
      error_code = 0 #jjs added
    )
  )
}



