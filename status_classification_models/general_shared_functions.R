
# GENERAL SHARED FUNCTIONS
# functions useful for various ML methods

##############
# FUNCTIONS: #
##############

# 1. prepare_mat: clean up beta matrix for modeling input (tranpose to put CpGs #
#      as columns, add Sample_Name and Status
# 2. prep_for_limma: prepare cv analysis split object for Limma DMC analysis
#      (basically undo what above function did)
# 3. find_assoc_features: Limma DMC Analysis to find smoking associated CpGs
# 4. gen_recipe: recipe of pre-processing steps on dataset for modeling
#      (set Sample_Name as ID, center and scale, then prep)
# 5. plot_metric_heatmap: create a heatmap of two hyperparameters (used in Lasso)
# 5. plot_three_hyperparams_heatmap: create heatmaps of two hyperparameters,
#      faceted across a third hyperparameter (elastic net, others?)
# 6. get_best_hps: get the best hyperparameter combo according to performance
#      metric
# 7. penalized_logreg_var_imp: collect and export variable importance for
#      lasso and elastic net regression models
# 8. rf_var_imp: collect and export variable importance for rf model
# 9. gbm_var_imp: collect and export variable importance for gbm model
# 10. execute_modeling: general function to perform modeling and assessment on a
#       CV fold
# 11. important_cpg_barplot: plot a barplot of the top 20 CpGs from a var importance txt file
# 12. extract_mean_sd_table: after fitting model and obtaining cpgs with nonzero coefs,
#       extract the mean and sd of the nonzero coef cpgs from the training data.
#       This table can then be used to standardize new data when predicting. 

#-----------------------------------------------------------------------------#
#-------------------------------------------#
# FUNCTION                                  #
# "prepare_mat"                             #
# set up the beta matrix for modeling input #
#-------------------------------------------#

# Function to prepare new beta matrix for ML modeling
#  - transpose to add CpGs as columns and samples as rows
#  - add "Sample_Name" and "status" columns from pd
prepare_mat <- function(beta_matrix, pd) {
  # remove rows from metadata table if the sample is not found as a column in
  #   the beta matrix
  coldata <- pd[rownames(pd) %in% colnames(beta_matrix), ]

  # take transpose to put CpGs as columns and samples as rows
  beta_t <- data.frame(t(beta_matrix))

  # add status column to represent smoking status
  # first, add Sample_Name to beta_t for merging
  beta_t$Sample_Name <- rownames(beta_t)
  beta_t <- cbind(beta_t,
    coldata[ ,"status"][match(rownames(beta_t), rownames(coldata))]
  )
  colnames(beta_t)[ncol(beta_t)] <- "status"

  return(beta_t)
}

#--------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------#
#---------------------------------------------------------#
# FUNCTION                                                #
# "prep_for_limma"                                        #
# prepare cv analysis split object for Limma DMC analysis #
#---------------------------------------------------------#

# Function to prepare a cv split object (from rsample's vfold_cv) by prepping
#   analysis set for Limma
#   - remove columns: Sample_Name and status
#   - transponse to have samples as columns and CpGs as rows
prep_for_limma <- function(splits) {
  # obtain the data for a split from the analysis set
  df <- analysis(splits)

  # remove Sample_Name and status
  df2 <- subset(df, select=-c(Sample_Name,status))

  # transpose
  df_t <- as.data.frame(t(df2))

  # return(df_t)
}

#---------------------------------------------------------------------------------#

#----------------------------------------------------#
# FUNCTION                                           #
# "find_assoc_features"                              #
# Limma DMC Analysis to find smoking associated CpGs #
#----------------------------------------------------#

# Function to perform diff analysis with Limma to find DMCs
#  perform "Smoker vs Nonsmoker" comparison
#  returns the input beta matrix now with only sig CpGs from comparison
find_assoc_features <- function(beta_matrix, pd, padj=0.05, n=100,
  rm_cor=TRUE, cor_cutoff=0.75) {
  # remove rows from metadata table if the sample is not found as a column in the
  #   beta matrix
  coldata <- pd[rownames(pd) %in% colnames(beta_matrix), ]

  # unadjusted model
  design <- model.matrix(~0+status, data=coldata)

  # create contrast matrix of smoker vs nonsmoker
  contMatrix <- makeContrasts(statussmoker - statusnonsmoker, levels=design)

  # fit linear models for each CpG
  fit <- lmFit(beta_matrix, design)

  # perform contrast
  fit2 <- contrasts.fit(fit, contMatrix)

  # compute:
  #  t-statistics,
  #  moderated F-statistic,
  #  and log-odds by empirical Bayes moderation of standard errors towards a
  #    common value
  fit2 <- eBayes(fit2)

  # extract top-ranked CpGs (all CpGs, ordered by rank of abs(t-statistic)
  DMPs <- topTable(fit2, num=Inf, coef=1, sort.by="P")
  DMPs$probe <- row.names(DMPs)

  # select CpGs with padj < assigned value (default 0.05)
  #DMPs_sig <- DMPs[DMPs$adj.P.Val < padj, ]
  # select the top n CpGs based on lowest unadjusted p-value
  DMPs_sig <- DMPs[c(1:n), ]

  # subset beta matrix to sig CpGs
  include <- DMPs_sig$probe
  beta_keep <- beta_matrix[include,]
  print(paste0("Found ", nrow(beta_keep), " sig smoking associated CpGs."))

  if (rm_cor==FALSE) {
    return(beta_keep)
  }
  else if (rm_cor==TRUE) {
    # get a pairwise pearson correlation matrix for each pair of sig CpGs
    cor <- cor(t(beta_keep))

    # find the highly correlated sig diff CpGs
    # use pairwise absolute value correlation cutoff of assigned value
    #   (default 0.75)
    # return column names
    highlyCorrelated <- findCorrelation(cor, cutoff=cor_cutoff, names=TRUE)

    # remove highly correlated CpGs from sig diff beta df
    beta_keep2 <- beta_keep[!row.names(beta_keep) %in% highlyCorrelated, ]
    print(paste0("Removed ",
      nrow(beta_keep) - nrow(beta_keep2),
      " correlated CpGs.")
    )
    print(paste0("Returning ", nrow(beta_keep2), " CpGs."))

    return(beta_keep2)
  }
}

#-------------------------------------------------------------------------------#

#---------------------------------------------------------#
# FUNCTION                                                #
# "gen_recipe"                                            #
# recipe of pre-processing steps on dataset for modeling  #
#---------------------------------------------------------#

# for now, using this for all ML methods

# recipe function to perform on datasets
# include the prep step last
gen_recipe <- function(dataset) {
  recipe(status ~ ., data=dataset) %>%
    update_role(Sample_Name, new_role="ID") %>%
    step_zv(all_numeric(), -all_outcomes()) %>%
    step_normalize(all_numeric(), -all_outcomes())
}


#-------------------------------------------------------------------------------#

#-------------------------------#
# FUNCTION                      #
# "plot_metric_heatmap"         #
# plot heatmap of 2 hyperparams #
#-------------------------------#

# Function to create heatmaps of two hyperparameters
# useful for Lasso regression

# df: dataframe output (mydf) of the loop running CV on each hyperparameter
#       combo
# param1: first hyperparameter, must be colname of df
#         the values in this column should be just the value of that
#           hyperparameter (for plotting)
# param2: second hyperparameter, must be colname of df
#         the values in this column should be just the value of that
#           hyperparameter (for plotting)
# num_metrics: number of performance metrics, usually 3 with all my hardcode
plot_metric_heatmap <- function(df, param1, param2, num_metrics=3,
  plotfile="grid_metric_heatmap.png", output_df=TRUE,
  outfile="grid_metric_summary.txt") {
  # create a tibble of the results and add unqID column (unique combo of
  #   lambda and fold)
  mydf.tib <- as_tibble(df)
  # create unique id (combo of the two hyperparameters and fold)
  mydf.tib <- mydf.tib %>%
    add_column(unqID=paste0(mydf.tib[[param1]],
      "_", mydf.tib[[param2]], "_", mydf.tib$id)
    )

  # create tibble of metrics per lambda/alpha/fold
  myres <- mydf.tib %>%
    group_by(unqID) %>%
    summarize(kap = unique(kap),
      accuracy = unique(accuracy),
      roc_auc = unique(roc_auc))

  # create unique key for just the two hyperparameters
  myres$param_key <- gsub("_Fold[0-9]*", "", myres$unqID)

  # create standard error function
  std <- function(x) sd(x)/sqrt(length(x))

  # Create metric plot df for plotting
  foo <- myres %>%
    group_by(param_key) %>%
    summarize(mean_kap = mean(kap),
      std_kap = std(kap),
      mean_accuracy = mean(accuracy),
      std_accuracy = std(accuracy),
      mean_roc_auc = mean(roc_auc),
      std_roc_auc = std(roc_auc))

  # create a numeric columns containing the value of each hyperparameter
  foo <- foo %>%
    separate(param_key, c(param1, param2), sep="_")

  penalty_vals <- num_metrics * nrow(foo)

  # Create plot df
  plotdf <- data.frame(p1=as.factor(as.numeric(c(rep(foo[[param1]], num_metrics)))),
    p2=as.factor(as.numeric(c(rep(foo[[param2]], num_metrics)))),
    metric=c(rep("kap", penalty_vals), rep("accuracy", penalty_vals),
      rep("roc_auc", penalty_vals)),
    mean=c(foo$mean_kap, foo$mean_accuracy, foo$mean_roc_auc),
    std=c(foo$std_kap, foo$std_accuracy, foo$std_roc_auc)
  )

  #plotdf <- data.frame(p1=as.factor(as.numeric(formatC(as.numeric(c(rep(
  #  foo[[param1]], num_metrics))),format="e", digits=2))),
  #  p2=as.factor(as.numeric(formatC(as.numeric(c(rep(
  #    foo[[param2]], num_metrics))), format="e", digits=2))),
  #  metric=c(rep("kap", penalty_vals), rep("accuracy", penalty_vals),
  #    rep("roc_auc", penalty_vals)),
  #  mean=c(foo$mean_kap, foo$mean_accuracy, foo$mean_roc_auc),
  #  std=c(foo$std_kap, foo$std_accuracy, foo$std_roc_auc)
  #)

  # Plot heatmap for each metric combined in one plot
  ggplot(plotdf, aes(p1, p2, fill=mean)) +
    geom_tile() +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~metric, scales = "free", nrow = 3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab(param1) +
    ylab(param2)
  ggsave(plotfile)

  # if specified, write plotdf to txt file
  if (output_df==TRUE) {
    colnames(plotdf)[c(1,2)] <- c(param1, param2)
    write.table(plotdf, outfile, sep="\t", col.names=T, row.names=F, quote=F)

    return(plotdf)
  }

}

#------------------------------------------------------------------------------#

#-------------------------------------------------------------#
# FUNCTION                                                    #
# "plot_three_hyperparams_heatmap"                            #
# plot heatmap of 2 hyperparams faceted across 3rd hyperparam #
#-------------------------------------------------------------#

# Function to create heatmaps of two hyperparameters, faceted across a third
#   hyperparameter
# useful for elastic net regression
plot_three_hyperparams_heatmap <- function(mytib, p1, p2, p3, metric="kap",
  num_metrics=3, plotfile="grid_metric_heatmap.png", output_df=TRUE,
  outfile="grid_metric_summary.txt") {

  # create unique id (combo of the three hyperparameters and fold)
  mytib <- mytib %>%
    add_column(unqID=paste0(mytib[[p1]],
      "_", mytib[[p2]], "_", mytib[[p3]], "_", mytib$id)
    )

  # create tibble of metrics per each unique hyperparam/fold combo
  myres <- mytib %>%
    group_by(unqID) %>%
    summarize(kap = unique(kap),
      accuracy = unique(accuracy),
      roc_auc = unique(roc_auc))


  # create unique key for just the hyperparameters
  myres$param_key <- gsub("_Fold[0-9]*", "", myres$unqID)

  # create standard error function
  std <- function(x) sd(x)/sqrt(length(x))

  # Create metric plot df for plotting
  foo <- myres %>%
    group_by(param_key) %>%
    summarize(mean_kap = mean(kap),
      std_kap = std(kap),
      mean_accuracy = mean(accuracy),
      std_accuracy = std(accuracy),
      mean_roc_auc = mean(roc_auc),
      std_roc_auc = std(roc_auc))

  # create numeric columns containing the value of each hyperparameter
  foo <- foo %>%
    separate(param_key, c(p1, p2, p3), sep="_")

  # Create a new df in "long" format for plotting
  #repeat_num <- num_metrics*nrow(foo)
  # Create plot df
  plotdf <- data.frame(p1=as.factor(as.numeric(formatC(as.numeric(c(rep(
    foo[[p1]], num_metrics))), format="e", digits=2))),
    p2=as.factor(as.numeric(formatC(as.numeric(c(rep(
      foo[[p2]], num_metrics))), format="e", digits=2))),
    p3=as.factor(as.numeric(c(rep(foo[[p3]], num_metrics)))),
    metric=c(rep("kap", nrow(foo)), rep("accuracy", nrow(foo)),
      rep("roc_auc", nrow(foo))),
    mean=c(foo$mean_kap, foo$mean_accuracy, foo$mean_roc_auc),
    std=c(foo$std_kap, foo$std_accuracy, foo$std_roc_auc)
  )

  # subset to metric of choice
  plotdf2 <- plotdf[plotdf$metric==metric, ]
  # change colname of mean to include metric
  #colnames(plotdf2)[5] <- paste0("mean_", metric)

  # Plot heatmap for each metric combined in one plot
  ggplot(plotdf2, aes(p1, p2, fill=mean)) +
    geom_tile() +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~p3, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab(p1) +
    ylab(p2)
    #guides(fill=guide_legend(paste0("mean ", metric)))
  ggsave(plotfile)

  # if specified, write plotdf to txt file
  if (output_df==TRUE) {
    colnames(plotdf2)[c(1,2,3)] <- c(p1, p2, p3)
    write.table(plotdf2, outfile, sep="\t", col.names=T, row.names=F, quote=F)

    return(plotdf2)
  }

}

#-----------------------------------------------------------------------------#

# get the best hyperparameter combo according to performance metric
# takes the output of the above function as input (a "plotdf")
get_best_hps <- function(df, metric, n) {
  # subset to only specified metric
  df.sub <- df[df$metric==metric, ]

  # sort by descending mean value
  df.sub2 <- df.sub[order(df.sub$mean, decreasing=T), ]

  print(df.sub2[c(1:n),])
  return(df.sub2)
}

#------------------------------------------------------------------------------#

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

  # export variable importance
  write.table(as.data.frame(nonzero_imp_vars),
    paste0(hyperparam_string, "_", id, "_nonzero_predictors.txt"),
    sep="\t",
    col.names=T,
    row.names=F,
    quote=F
  )
}

#--------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------#
# FUNCTION                                                           #
# "rf_var_imp"                                                       #
# function to obtain var importance for random forest fitted model   #
# currently from the engine "ranger" using "impurity" for importance #
#--------------------------------------------------------------------#
# INPUTS:        mod_fit: a fitted model object (random forest)
#      hyperparam_string: string with hyperparameters and values for naming exported file
#                     id: my_folds$id, passed in via execute_modeling function argument
rf_var_imp <- function(mod_fit, hyperparam_string, id) {
  # gather the ranger object
  ranger_obj <- mod_fit$fit

  # get variable importance
  important_vars <- data.frame(CpG=names(ranger_obj$variable.importance),
    imp=ranger_obj$variable.importance)

  # export the importance for all CpGs for this fold
  write.table(important_vars,
    paste0(hyperparam_string, "_", id, "_var_importance.txt"),
    sep="\t",
  	col.names=T,
  	row.names=F,
  	quote=F
  )
}

#--------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# FUNCTION                                                                  #
# "gbm_var_imp"                                                             #
# function to obtain var importance for gradien boost machine fitted model  #
# currently from the engine "xgboost"                                       #    
#---------------------------------------------------------------------------#
# INPUTS:        mod_fit: a fitted model object (random forest)
#      hyperparam_string: string with hyperparameters and values for naming exported file
#                     id: my_folds$id, passed in via execute_modeling function argument
gbm_var_imp <- function(mod_fit, hyperparam_string, id) {
  # get variable importance from xgboost fitted model
  # default for binary logistic regression is logloss (not "error")
  imp_vars <- xgb.importance(model=mod_fit$fit)

  # export
  write.table(imp_vars,
    paste0(hyperparam_string, "_", id, "_var_importance.txt"),
    sep="\t",
  	col.names=T,
  	row.names=F,
  	quote=F
  )
}

#--------------------------------------------------------------------------------------------#

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

execute_modeling <- function(splits, id, pd, padj=0.05, n=100, rm_cor=TRUE,
  cor_cutoff=0.75, mod_spec, method, hyperparam_string, pred_results) {
  # prepare the analysis cv split set to be used with limma
  limma_input_mat <- prep_for_limma(splits)

  # run diff analysis to select top smoking associated CpGs
  df <- find_assoc_features(beta_matrix=limma_input_mat,
    pd=pd,
    padj=padj,
    n=n,
    rm_cor=rm_cor,
    cor_cutoff=cor_cutoff
  )

  # prepare the format of the subsetted beta matrix for modeling
  df_clean <- prepare_mat(beta_matrix=df,
    pd=pd
  )

  # pre-process with recipe
  analysis_prepped <- gen_recipe(df_clean) %>%
    prep(strings_as_factors = FALSE)

  # bake the data
  analysis_baked <- analysis_prepped %>%
    bake(new_data=df_clean)

  # fit the input model on pre-processed analysis set
  # note that we need to exclude our ID variable (Sample_Name)
  #   because we are not using a workflow()
  # the fit function doesn't know that the ID column of the tibble has a
  #   special role
  mod_fit <- mod_spec %>%
    fit(status ~ ., data=analysis_baked[ ,!(colnames(analysis_baked) %in% "Sample_Name")])

  # get assessment set of this split
  assessment_set <- assessment(splits)
  # pre-process the assessment set using the statistics computed from the
  #   analysis set
  # subset assessment set to include the same CpGs as used in the cleaned
  #   analysis set
  assessment_set_sub <- assessment_set[,colnames(assessment_set) %in%
    colnames(df_clean)]
  # preprocess with same prepped data
  assessment_baked <- analysis_prepped %>%
    bake(new_data=assessment_set_sub)

  #----- Variable Importance -----#

  if (method=="penalized_logreg") {
    print("collecting feature importance for penalized logistic regression model.")
    penalized_logreg_var_imp(mod_fit=mod_fit,
      hyperparam_string=hyperparam_string,
      id=id
    )
  }
  else if (method=="random_forest") {
    print("collecting feature importance for random forest model.")
    rf_var_imp(mod_fit=mod_fit,
      hyperparam_string=hyperparam_string,
      id=id
    )
  }
  else if (method=="gbm") {
    print("collecting feature importance for gradient boost machine model.")
    gbm_var_imp(mod_fit=mod_fit,
      hyperparam_string=hyperparam_string,
      id=id
    )  
  }

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
    add_column("truth"=assessment_baked$status) %>%
    add_column("prediction"=mod_fit %>%
      predict(new_data=assessment_baked) %>%
        unlist()) %>%
    add_column("prob_nonsmoker"=predict(mod_fit,
      new_data=assessment_baked, type="prob")[,1] %>%
        unlist()) %>%
    add_column("prob_smoker"=predict(mod_fit,
      new_data=assessment_baked, type="prob")[,2] %>%
        unlist())
      

  #----- Store Performance Metrics -----#
  my_mets <- tibble("accuracy"=mod_results %>%
    metrics(truth, prediction) %>%
    filter(.metric == "accuracy") %>%
    pull(.estimate),
    "kap"=mod_results %>%
    metrics(truth, prediction) %>%
    filter(.metric == "kap") %>%
    pull(.estimate),
    "roc_auc"=mod_results %>%
    roc_auc(truth, prob_nonsmoker) %>%
    pull(.estimate)
  )

  #----- Combine Prediction Results and Metrics -----#

  # add metrics to results tibble
  my_res <- mod_results %>%
    add_column(accuracy=my_mets$accuracy) %>%
    add_column(kap=my_mets$kap) %>%
    add_column(roc_auc=my_mets$roc_auc)

  return(my_res)
}


#--------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------#
# FUNCTION                                                          #
# "important_cpg_barplot"                                           #
# general function to plot barplot of var importance from txt file  #
#-------------------------------------------------------------------#

# INPUTS:
# var_imp_file: a "nonzero_predictors.txt" or "variable_importance.txt" file
#               must be ordered by descending abs(estimate)
#               can have Intercept row (will be removed if "Intercept" is present)
#               assumes 2nd column is desired estimate/importance column
#     filename: desired name of plot file
#   plot_title: desired plot title
#    bar_color: desired bar color

important_cpg_barplot <- function(var_imp_file, filename="important_cpgss_barplot.png", plot_title="Important CpGs", bar_color="gray") {
  # read in important cpgs file ordered by abs(estimate), and
  # remove intercept row if present, and
  # select the top 20 rows (in terms of abs(estimate))
  df <- read_delim(var_imp_file, delim="\t") %>%
    filter(!grepl('Intercept', term))

  # if there are more than 20 rows, filter to top 20
  # assuming data is ordered by descending abs(estimate)
  if (nrow(df) > 20) {
    print("Found more than 20 rows, subsetting to first 20.")
    df <- df %>%
      slice_head(n = 20)
  }

  # add an absolute value of estimate column for plotting
  # (because lasso and elastic net can have negative values)
  # assumes estimate column is 2nd column
  df <- df %>%
    mutate(abs_est = abs(.[[2]]))

  # clean up colnames for plotting
  colnames(df)[1] <- "term"

  myplot <- ggplot(df, aes(x = reorder(term, -abs_est), y = abs_est)) +
    labs(x="CpG ID", y="Importance", title=plot_title) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_bar(stat="identity", fill=bar_color) +
    theme_minimal() +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename=filename, plot=myplot)
}

#--------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------#
# FUNCTION                                                          #
# "extract_mean_sd_table"                                           #
# obtain mean and sd of cpgs with nonzero coefs from training data  #
#-------------------------------------------------------------------#

# INPUTS:
#         df_clean: training data, rows are samples and columns are cpgs
#                   rownames and colnames should be present
# nonzero_imp_vars: dataframe, column "term" is cpgs with nonzero coefs
#           export: if TRUE (default), export the mean sd cpg table to a txt file
#          outfile: name of desired output txt file, if export is TRUE

extract_mean_sd_table <- function(df_clean=df_clean, nonzero_imp_vars=nonzero_imp_vars, export=TRUE,
  outfile="mean_sd_table.txt") {

  # subset df_clean to contain only cpgs with nonzero coefs
  mydf <- as.data.frame(df_clean[ ,colnames(df_clean) %in% nonzero_imp_vars$term])

  # to clean a df if it only has 1 cpg column
  if (ncol(mydf)==1) {
    rownames(mydf) <- rownames(df_clean)
    colnames(mydf) <- colnames(df_clean)[colnames(df_clean) %in% nonzero_imp_vars$term]
  }

  # For each nonzero coef cpg in the training data, calculate and store mean and sd
  mean_vec <- vector()
  sd_vec <- vector()
  for (i in 1:ncol(mydf)) {
    mean_val <- mean(mydf[ ,i])
    sd_val <- sd(mydf[ ,i])
    mean_vec <- c(mean_vec, mean_val)
    sd_vec <- c(sd_vec, sd_val)
  }

  # create a new df with columns for cpg, mean, and sd
  newdf <- data.frame(cpg=colnames(mydf),
    mean=mean_vec,
    sd=sd_vec)

  # if export is true, write the mean sd cpg table to a txt file
  if (export==TRUE) {
    write.table(newdf,
    outfile,
    sep="\t",
    col.names=T,
    row.names=F,
    quote=F)
  }

  return(newdf)
}
