#!/home/groups/hoolock2/u0/bd/miniconda3/envs/tidymodels/bin/Rscript

######################################
# Model 5: LASSO Logistic Regression #
######################################

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

#------------------------------------------------------------------------------#

# Store hyperparameter values
penalty_val <- 0.0885866790410083
nCpG_val <- 120
# Store hyperparameter value string for file names
hyper_string <- paste0("nCpG_", nCpG_val, "_lambda_", penalty_val)

# set random seed
set.seed(42)

#------------------------------------------------------------------------------#

# Prepare input data for limma

# remove Sample_Name and status
df2 <- subset(mydata, select = -c(Sample_Name, status))
# transpose
df_t <- as.data.frame(t(df2))

# run diff analysis to select top smoking associated CpGs
df <- find_assoc_features(beta_matrix = df_t,
  pd = pd,
  padj = 0.05,
  n = nCpG_val,
  rm_cor = TRUE,
  cor_cutoff = 0.75
)

# prepare the subsetted beta matrix for LASSO
df_clean <- prepare_mat(beta_matrix = df,
  pd = pd
)

# export this final "clean" input cpg table for future reference
write.table(df_clean,
  "model5.clean_input_betas.txt",
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE)

# pre-process with recipe
analysis_prepped <- gen_recipe(df_clean) %>%
  prep(strings_as_factors = FALSE)

# bake the data
analysis_baked <- analysis_prepped %>%
  bake(new_data = df_clean)

# specify LASSO logistic regression
lasso_spec <- logistic_reg(mode = "classification",
  penalty = penalty_val, mixture = 1) %>%
  set_engine("glmnet")

# fit lasso logsitic reg model on pre-processed analysis set
lasso_fit <- lasso_spec %>%
  fit(status ~ .,
  data = analysis_baked[, !colnames(analysis_baked) %in% c("Sample_Name")])

# export the fitted model
saveRDS(lasso_fit, "model5_lasso_fit.RDS")

# obtain predictors with nonzero coefficients from fitted model
imp_vars <- tidy(lasso_fit)
nonzero_imp_vars <- imp_vars[imp_vars$estimate != 0, ]

#----------------------#

# extract mean sd table for cpgs with nonzero coefs from training data
mean_sd_tab <- extract_mean_sd_table(df_clean,
  nonzero_imp_vars,
  export = TRUE,
  outfile = "model5.mean_sd_table.txt")

#---------------------#

# export all imp vars
write.table(imp_vars,
  "imp_vars.txt",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE)

# remove intercept row if present,
# and order rows by descending absolute value of estimate
# then add intercept row to the bottom
nonzero_imp_vars2 <- nonzero_imp_vars %>%
  filter(!grepl("Intercept", term)) %>%
  arrange(desc(abs(estimate))) %>%
  rbind(nonzero_imp_vars[1,])

# export important predictor vars with nonzero coefficients for this fold
write.table(as.data.frame(nonzero_imp_vars2),
  paste0(hyper_string, "_nonzero_predictors.txt"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# Plot variable importance of top 20 CpGs
important_cpg_barplot(var_imp_file = paste0(hyper_string,
  "_nonzero_predictors.txt"),
  filename = "lasso_important_cpgs_barplot.png",
  plot_title = "LASSO Importance CpGs",
  bar_color = "cornflowerblue"
)

# END






#------------------------------------#

# rough explore:

# subset the beta matrix to cpgs in the model (zero or nonzero coef)
beta.sub <- df_t[rownames(df_t) %in% imp_vars$term, ]

# take transpose to put cpgs as columns
newdf <- as.data.frame(t(beta.sub))

# add Sample_Name and status columns
newdf$Sample_Name <- pd$Sample_Name
newdf$status <- pd$status

# Pretend this is a new dataset
# pre-process it
foo <- gen_recipe(newdf) %>%
  prep(strings_as_factors = FALSE) %>%
  bake(new_data=newdf)

# Use the fitted model to predict on new data (all 96 vcsip samples)
res <- predict(lasso_fit, new_data=foo)
