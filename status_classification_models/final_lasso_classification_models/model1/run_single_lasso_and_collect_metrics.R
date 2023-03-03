#!/home/groups/hoolock2/u0/bd/miniconda3/envs/tidymodels/bin/Rscript

########################
# For LASSO Regression #
########################

# LIBRARIES
library(tidymodels) # framework for machine learning models
library(optparse)   # for command line arguments
library(limma)      # for performing differential analysis on CpG beta values
library(caret)      # to remove highly correlated sig diff CpGs
library(tidyverse)  # for manipulations

#------------------------------------------------------------------------------#

# Source functions
source(paste0("/home/groups/hoolock2/u0/bd/Projects/lyndsey_shorey_project/",
  "smoking_index_github/placenta_smoking_index/status_classification_models/",
  "general_shared_functions.R"))

#------------------------------------------------------------------------------#

# Read in data and metadata

# full beta matrix: 714666 CpGs x 96 samples
beta <- read.delim(gzfile(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/prelim_data/cleaned_data/new_beta.txt.gz")))
# pd table for 96 samples as rows
pd <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/prelim_data/cleaned_data/pd.txt"))
pd$status <- as.factor(pd$status)

# Prepare beta matrix for analysis
#   - transpose to place CpGs as columns and samples as rows
#   - add Sample_Name column to use as ID
#   - add status column to use as response variable
mydata <- prepare_mat(beta, pd)

#-------------------------------------------------------------------------------------------------------------------#

# Store hyperparameter values
penalty_val <- 3.79269019073225e-09
nCpG_val <- 120
# Store hyperparameter value string for file names
hyper_string <- paste0("nCpG_", nCpG_val, "_lambda_", penalty_val)

# set random seed
set.seed(42)

#-------------------------------------------------------------------------------------------------------------------#

# Prepare input data for limma

# remove Sample_Name and status
df2 <- subset(mydata, select=-c(Sample_Name,status))
# transpose
df_t <- as.data.frame(t(df2))

# run diff analysis to select top smoking associated CpGs
df <- find_assoc_features(beta_matrix=df_t,
                          pd=pd,
                          padj=0.05,
                          n=nCpG_val,
                          rm_cor=TRUE,
                          cor_cutoff=0.75
)

# prepare the subsetted beta matrix for LASSO
df_clean <- prepare_mat(beta_matrix=df,
                        pd=pd
)

# export this final "clean" input cpg table for future reference
write.table(df_clean,
  "model1.clean_input_betas.txt",
  sep="\t",
  col.names=T,
  row.names=T,
  quote=F)
  
# pre-process with recipe
analysis_prepped <- gen_recipe(df_clean) %>%
  prep(strings_as_factors = FALSE)

# bake the data
analysis_baked <- analysis_prepped %>%
  bake(new_data=df_clean)
    
# specify LASSO logistic regression
lasso_spec <- logistic_reg(mode="classification", penalty=penalty_val, mixture=1) %>%
  set_engine("glmnet")

# fit lasso logsitic reg model on pre-processed analysis set
lasso_fit <- lasso_spec %>%
  fit(status ~ ., data=analysis_baked[,!colnames(analysis_baked) %in% c("Sample_Name")])

# export the fitted model
saveRDS(lasso_fit, "model1_lasso_fit.RDS")
# can be re-loaded with: my_fitted_model <- readRDS("model1_lasso_fit.RDS")

# obtain predictors with nonzero coefficients from fitted model
imp_vars <- tidy(lasso_fit)
nonzero_imp_vars <- imp_vars[imp_vars$estimate != 0, ]

#----------------------#
# extract the mean and sd for each cpg with a nonzero coef

# subset df_clean df to contain only cpgs with nonzero coefs
mydf <- df_clean[ ,colnames(df_clean) %in% nonzero_imp_vars$term]
# do
mean_vec <- vector()
sd_vec <- vector()
for (i in 1:ncol(mydf)) {
  mean_val <- mean(mydf[ ,i])
  sd_val <- sd(mydf[ ,i])
  mean_vec <- c(mean_vec, mean_val)
  sd_vec <- c(sd_vec, sd_val)
}
newdf <- data.frame(cpg=colnames(mydf),
  mean=mean_vec,
  sd=sd_vec)
write.table(newdf,
  "model1.mean_sd_table.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

#------------------------#

# export all imp vars (102 cpgs + intercept)
write.table(imp_vars,
  "imp_vars.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

# remove intercept row if present,
# and order rows by descending absolute value of estimate
# then add intercept row to the bottom
nonzero_imp_vars2 <- nonzero_imp_vars %>%
  filter(!grepl('Intercept', term)) %>%
  arrange(desc(abs(estimate))) %>%
  rbind(nonzero_imp_vars[1,])

# export important predictor vars with nonzero coefficients for this fold
write.table(as.data.frame(nonzero_imp_vars2),
            paste0(hyper_string, "_nonzero_predictors.txt"),
            sep="\t",
            col.names=T,
            row.names=F,
            quote=F
)

# Plot variable importance of top 20 CpGs
important_cpg_barplot(var_imp_file = paste0(hyper_string, "_nonzero_predictors.txt"),
    filename = "lasso_important_cpgs_barplot.png",
    plot_title = "LASSO Model 1 EPIC Importance CpGs",
    bar_color = "cornflowerblue"
)

#------------------------------------#

# subset a vcsip beta matrix of nonzero cpgs
mydf <- df_clean[ ,colnames(df_clean) %in% nonzero_imp_vars$term]
write.table(mydf,
  "model1.vcsip.raw_beta_matrix.txt",
  sep="\t",
  col.names=T,
  row.names=T,
  quote=F)

# subset the normalized vcsip beta matrix of nonzero cpgs
mydfnorm <- analysis_baked[ ,colnames(analysis_baked) %in% nonzero_imp_vars$term]
write.table(mydfnorm,
  "model1.vcsip.norm_beta_matrix.txt",
  sep="\t",
  col.names=T,
  row.names=T,
  quote=F)

#------------------------------------#

# predict and store predictions in a new dataframe
pred_res <- data.frame(Sample_Name=analysis_baked$Sample_Name,
  status=analysis_baked$status) %>%
  add_column("prediction" = predict(lasso_fit, analysis_baked) %>%
  unlist()) %>%
  add_column("prob_nonsmoker" = predict(lasso_fit, analysis_baked, type = "prob")[,1] %>%
  unlist()) %>%
  add_column("prob_smoker" = predict(lasso_fit, analysis_baked, type = "prob")[,2] %>%
  unlist())

# export for future reference
write.table(pred_res,
  "model1.prediction_table.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

# convert the "truth" (status) to factor
pred_res$status <- as.factor(pred_res$status)

#-------#

# Get accuracy, sensitivity, specificity, cohen's kappa, AUC_ROC, and partial AUC from 0.9-1
mod_results <- pred_res
colnames(mod_results)[2] <- "truth"
my_mets <- tibble("accuracy"=mod_results %>%
                    metrics(truth, prediction) %>%
                    filter(.metric == "accuracy") %>%
                    pull(.estimate),
		  "sensitivity"=mod_results %>%
                    sens(truth, prediction) %>%
                    pull(.estimate),
                  "specificity"=mod_results %>%
                    yardstick::spec(truth, prediction) %>%
                    pull(.estimate),
                  "kap"=mod_results %>%
                    metrics(truth, prediction) %>%
                    filter(.metric == "kap") %>%
                    pull(.estimate),
                  "roc_auc"=mod_results %>%
                    roc_auc(truth, prob_nonsmoker) %>%
                    pull(.estimate),
		  "p_auc"=pROC::auc(as.numeric(mod_results$truth),
		    as.numeric(mod_results$prediction),
		    partial.auc=c(0.9, 1))[1]
)

# export metrics table
write.table(my_mets,
  "model1_vcsip_performance_metrics.txt",
  sep="\t",
  col.names=TRUE,
  row.names=FALSE,
  quote=FALSE)

# calculate the roc curve and plot the roc curve

p <- roc_curve(pred_res, status, prob_nonsmoker) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  ggtitle("Model 1 EPIC ROC")
ggsave(filename="model1_vcsip_roc.pdf", p)
