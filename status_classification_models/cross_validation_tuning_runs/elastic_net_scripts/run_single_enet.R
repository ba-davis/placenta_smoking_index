#!/home/groups/hoolock2/u0/bd/miniconda3/envs/tidymodels/bin/Rscript

##############################
# For Elastic Net Regression #
##############################

# LIBRARIES
library(tidymodels) # framework for machine learning models
library(optparse)   # for command line arguments
library(limma)      # for performing differential analysis on CpG beta values
library(caret)      # to remove highly correlated sig diff CpGs
library(tidyverse)  # for manipulations

#-------------------------------------------------------------------------------------------------------------------#

# Source functions
source("/home/groups/hoolock2/u0/bd/Projects/lyndsey_shorey_project/smoking_index_github/placenta_smoking_index/status_classification_models/general_shared_functions.R")

#-------------------------------------------------------------------------------------------------------------------#

# Read in data and metadata
beta <- read.delim(gzfile("/home/groups/hoolock2/u0/bd/Projects/lyndsey_shorey_project/prelim_data/cleaned_data/beta.txt.gz"))
pd <- read.delim("/home/groups/hoolock2/u0/bd/Projects/lyndsey_shorey_project/prelim_data/cleaned_data/pd.txt")
pd$status <- as.factor(pd$status)

# Prepare beta matrix for analysis
#   - transpose to place CpGs as columns and samples as rows
#   - add Sample_Name column to use as ID
#   - add status column to use as response variable
mydata <- prepare_mat(beta, pd)

#-------------------------------------------------------------------------------------------------------------------#

# Store hyperparameter values
penalty_val <- 0.0000000166810053720006
mixture_val <-0.666666666666667
nCpG_val <- 287
# Store hyperparameter value string for file names
hyper_string <- paste0("nCpG_", nCpG_val, "_lambda_", penalty_val, "_alpha_", mixture_val)

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

# pre-process with recipe
analysis_prepped <- gen_recipe(df_clean) %>%
  prep(strings_as_factors = FALSE)

# bake the data
analysis_baked <- analysis_prepped %>%
  bake(new_data=df_clean)
    
# specify Elastic Net logistic regression
mod_spec <- logistic_reg(mode="classification", penalty=penalty_val, mixture=mixture_val) %>%
  set_engine("glmnet")

# fit lasso logsitic reg model on pre-processed analysis set
mod_fit <- mod_spec %>%
  fit(status ~ ., data=analysis_baked[,!colnames(analysis_baked) %in% c("Sample_Name")])

# obtain predictors with nonzero coefficients from fitted model
imp_vars <- tidy(mod_fit)
nonzero_imp_vars <- imp_vars[imp_vars$estimate != 0, ]

# remove intercept row if present,
# and order rows by descending absolute value of estimate
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
df.new <- "nCpG_287_lambda_1.66810053720006e-08_alpha_0.666666666666667_nonzero_predictors.txt"

important_cpg_barplot(var_imp_file = df.new,
    filename = "enet_important_cpgs_barplot.png",
    plot_title = "Elastic Net Importance CpGs",
    bar_color = "coral1"
)
