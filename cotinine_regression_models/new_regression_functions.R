
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
#              padj: for "find_assoc_features", padj cutoff for top CpGs to
#                      return
#                 n: for "find_assoc_features", number of top CpGs to return
#            rm_cor: for "find_assoc_features", whether or not to remove
#                      correlated top CpGs
#        cor_cutoff: for "find_assoc_features", cutoff value to remove
#                      correlated CpGs
#          mod_spec: the model specification object
#            method: ML method used, for proper variable importance collection
#                      and export
#                    one of ("penalized_logreg", "random_forest", "gbm")
#                    penalized_logreg is used for lasso and elastic net
# hyperparam_string: the string of hyperparam values to be used in file name
#                      export
#      pred_results: a tibble with the hyperparameters as columns

execute_modeling <- function(splits, response_var, id, pd, padj=0.05, n=100, rm_cor=TRUE,
  cor_cutoff=0.75, mod_spec, method, hyperparam_string, pred_results) {
  # prepare the analysis cv split set to be used with limma
  #limma_input_mat <- prep_for_limma(splits)

  # get analysis set of CV split
  df <- analysis(splits)
  
  # pre-process with recipe
  analysis_prepped <- gen_recipe(df) %>%
    prep(strings_as_factors = FALSE)

  # bake the data
  analysis_baked <- analysis_prepped %>%
    bake(new_data=df

  # fit the input model on pre-processed analysis set
  # note that we need to exclude our ID variable (Sample_Name)
  #   because we are not using a workflow()
  # the fit function doesn't know that the ID column of the tibble has a
  #   special role
  mod_fit <- mod_spec %>%
    fit(Late_Cotinine_ng_ml ~ ., data=analysis_baked[ ,!(colnames(analysis_baked) %in% "Sample_Name")])

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





)


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



# Function to prepare a cv split object (from rsample's vfold_cv) by prepping
#   analysis set for Limma
#   - remove columns: Sample_Name and status
#   - transponse to have samples as columns and CpGs as rows
prep_for_limma2 <- function(splits, colname) {
  # obtain the data for a split from the analysis set
  df <- analysis(splits)

  # remove Sample_Name and status
  df2 <- subset(df, select=-c(Sample_Name, Late_Cotinine_ng_ml))

  # transpose
  df_t <- as.data.frame(t(df2))

  # return(df_t)
}

