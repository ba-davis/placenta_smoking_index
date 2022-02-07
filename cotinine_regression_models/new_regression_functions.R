



# INPUT: colname: name of column in pd to be added to beta matrix (usually the response variable)

prep_input <- function(beta_matrix, pd, colname) {
  # remove rows from metadata table if the sample is not found as a column in
  #   the beta matrix
  coldata <- pd[rownames(pd) %in% colnames(beta_matrix), ]

  # take transpose to put CpGs as columns and samples as rows
  beta_t <- data.frame(t(beta_matrix))

  # add status column to represent smoking status
  # first, add Sample_Name to beta_t for merging
  beta_t$Sample_Name <- rownames(beta_t)
  beta_t <- cbind(beta_t,
    coldata[ , colname][match(rownames(beta_t), rownames(coldata))]
  )
  colnames(beta_t)[ncol(beta_t)] <- colname

  return(beta_t)
}


#---------------------------------------------------------#
# FUNCTION                                                #
# "gen_recipe"                                            #
# recipe of pre-processing steps on dataset for modeling  #
#---------------------------------------------------------#

# for now, using this for all ML methods

# recipe function to perform on datasets
# include the prep step last
gen_recipe <- function(dataset) {
  recipe(Late_Cotinine_ng_ml ~ ., data=dataset) %>%
    update_role(Sample_Name, new_role="ID") %>%
    step_zv(all_numeric(), -all_outcomes()) %>%
    step_normalize(all_numeric(), -all_outcomes())
}


#-----------------------------------------------------------------------#
# FUNCTION                                                              #
# "penalized_logreg_var_imp"                                            #
# function to obtain nonzero cpgs for LASSO or Elastic Net fitted model #
#-----------------------------------------------------------------------#
# INPUTS:        mod_fit: a fitted model object
#      hyperparam_string: string with hyperparameters and values for naming exported file
#                     id: my_folds$id, passed in via execute_modeling function argument
penalized_logreg_var_imp <- function(mod_fit, hyperparam_string, id) {
  # collect variable importance
  imp_vars <- tidy(mod_fit)
  # subset to features with nonzero coefficients
  nonzero_imp_vars <- imp_vars[imp_vars$estimate != 0, ]

  # grab intercept
  my_int <- nonzero_imp_vars[1,]
  # order by absolute value of coefficients
  nonzero_imp_vars <- nonzero_imp_vars[-1,]
  nonzero_imp_vars2 <- nonzero_imp_vars[order(abs(nonzero_imp_vars$estimate), decreasing=T), ]
  # add intercept
  nonzero_imp_vars <- rbind(my_int, nonzero_imp_vars2)

  # export variable importance
  write.table(as.data.frame(nonzero_imp_vars),
    paste0(hyperparam_string, "_", id, "_nonzero_predictors.txt"),
    sep="\t",
    col.names=T,
    row.names=F,
    quote=F
  )
}


# New execute modeling function without subsetting nCpGs (just use input CpGs)

#-------------------------------------------------------------------#
# FUNCTION                                                          #
# "execute_modeling"                                                #
# general function to perform modeling and assessment on a CV fold  #
#-------------------------------------------------------------------#

# This function calls several functions in order to perform modeling and
#   assessment
# INPUTS: splits: CV fold (my_folds$splits[[1]])
#                id: my_folds$id, fold of the CV
#                pd: metadata table object
#          mod_spec: the model specification object
#            method: ML method used, for proper variable importance collection
#                      and export
#                    one of ("penalized_logreg", "random_forest", "gbm")
#                    penalized_logreg is used for lasso and elastic net
# hyperparam_string: the string of hyperparam values to be used in file name
#                      export
#      pred_results: a tibble with the hyperparameters as columns

execute_modeling <- function(splits, id, pd, mod_spec, method, hyperparam_string, pred_results) {
  # prepare the analysis cv split set to be used with limma
  #limma_input_mat <- prep_for_limma(splits)

  # get analysis set of CV split
  df <- analysis(splits)
  
  # pre-process with recipe
  analysis_prepped <- gen_recipe(df) %>%
    prep(strings_as_factors = FALSE)

  # bake the data
  analysis_baked <- analysis_prepped %>%
    bake(new_data=df)

  # fit the input model on pre-processed analysis set
  # note that we need to exclude our ID variable (Sample_Name)
  #   because we are not using a workflow()
  # the fit function doesn't know that the ID column of the tibble has a
  #   special role
  mod_fit <- mod_spec %>%
    fit(Late_Cotinine_ng_ml ~ .,
    data=analysis_baked[ ,!(colnames(analysis_baked) %in% "Sample_Name")])

  # get assessment set of this split
  assessment_set <- assessment(splits)
  # pre-process the assessment set using the statistics computed from the
  #   analysis set
  # subset assessment set to include the same CpGs as used in the cleaned
  #   analysis set
  assessment_set_sub <- assessment_set[,colnames(assessment_set) %in%
    colnames(df)]
  # preprocess with same prepped data
  assessment_baked <- analysis_prepped %>%
    bake(new_data=assessment_set_sub)

  #----- Variable Importance -----#

  print("collecting feature importance for penalized logistic regression model.")
  penalized_logreg_var_imp(mod_fit=mod_fit,
    hyperparam_string=hyperparam_string,
    id=id
  )
  #----- Store Prediction Results -----#

  # First repeat hyperparameter columns to have correct number of rows for
  #   merging
  rep_num <- length(assessment_baked$Sample_Name)
  newtib <- pred_results %>%
    dplyr::slice(rep(1:n(), each = rep_num))

  # Add the prediction results to the input pred_results tibble
  mod_results <- newtib %>%
    add_column("id"=id) %>%
    add_column("sample_name"=assessment_baked$Sample_Name) %>%
    add_column("truth"=assessment_baked$Late_Cotinine_ng_ml) %>%
    add_column("prediction"=mod_fit %>%
      predict(new_data=assessment_baked) %>%
        unlist())

  #----- Store Performance Metrics -----#
  my_mets <- tibble("rmse"=mod_results %>%
                      metrics(truth, prediction) %>%
                        filter(.metric == "rmse") %>%
                          pull(.estimate),
                    "rsq"=mod_results %>%
                      metrics(truth, prediction) %>%
                        filter(.metric == "rsq") %>%
                          pull(.estimate),
                    "mae"=mod_results %>%
                      metrics(truth, prediction) %>%
                        filter(.metric == "mae") %>%
                          pull(.estimate)
  )
 
  #----- Combine Prediction Results and Metrics -----#

  # add metrics to results tibble
  my_res <- mod_results %>%
    add_column(rmse=my_mets$rmse) %>%
    add_column(rsq=my_mets$rsq) %>%
    add_column(mae=my_mets$mae)

  return(my_res)
}
