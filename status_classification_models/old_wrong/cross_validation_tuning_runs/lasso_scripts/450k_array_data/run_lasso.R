#!/home/groups/hoolock2/u0/bd/miniconda3/envs/tidymodels/bin/Rscript

####################
# LASSO Regression #
####################
# Script to perform lasso regression to predict smoking status
# Uses specified hyperparameter values for the model
# Performs 10-fold CV
# Input beta matrix is EPIC microarray data
# Model is defined with hyperparameters, then 10-fold CV is performed
#   to estimate model performance with the given hyperparameter values
# This script will be launched in parallel via SLURM with a different job
#   for each unique hyperparameter value combination
# Output will be variable importance and model predictions/performance metrics

#------------------------------------------------------------------------------#

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
  "general_shared_functions.R")
)

#------------------------------------------------------------------------------#

# Read in the hyperparameter values from command line
option_list = list(
  make_option(c("-a", "--hyperparameter1"), type  = "numeric", default = NULL, 
              help = "value for hyperparameter1", metavar = "character"),
  make_option(c("-b", "--hyperparameter2"), type = "numeric", default = NULL,
              help = "value for hyperparameter2", metavar = "character"),
  make_option(c("-s", "--seed"), type="numeric", default = NULL,
              help = "value for random seed", metavar = "character")	      
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
if (is.null(opt$hyperparameter1)){
  print_help(opt_parser)
  stop("Missing hyperparameter1 value.n", call. = FALSE)
}
if (is.null(opt$hyperparameter2)){
  print_help(opt_parser)
  stop("Missing hyperparameter2 value.n", call. = FALSE)
}
if (is.null(opt$seed)){
  print_help(opt_parser)
  stop("Missing random seed value.n", call. = FALSE)
}

# Store hyperparameter values
penalty_val <- opt$hyperparameter1
nCpG_val <- opt$hyperparameter2
# Store hyperparameter value string for file names
hyper_string <- paste0("nCpG_", nCpG_val, "_lambda_", penalty_val)

#------------------------------------------------------------------------------#

# Read in data and metadata
beta <- read.delim(gzfile("/home/groups/hoolock2/u0/bd/Projects/lyndsey_shorey_project/models_450k_only/450k_probes/beta_450k.txt.gz"))

pd <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/prelim_data/cleaned_data/pd.txt")
)
pd$status <- as.factor(pd$status)

# Prepare beta matrix for analysis
#   - transpose to place CpGs as columns and samples as rows
#   - add Sample_Name column to use as ID
#   - add status column to use as response variable
mydata <- prepare_mat(beta, pd)

#------------------------------------------------------------------------------#

# Prepare data for general modeling function

# set random seed
set.seed(opt$seed)

# perform 10-Fold stratified CV
my_folds <- vfold_cv(mydata, strata = status, v = 10)

# specify the model with the currently defined hyperparameter(s)
# here, we get the hyperparameters from the command line options
mod_spec <- logistic_reg(mode = "classification",
  penalty = penalty_val,
  mixture = 1) %>%
    set_engine("glmnet")

# Create pred_resultsv1 tibble to hold hyperparameter values
pred_resultsv1 <- tibble(
  "n_CpG" = nCpG_val,
  "lambda" = penalty_val
)    

#---------------------------------------------------------------------------------------------------------------------#

# Perform modeling function using these hyperparameter values

# Map the formula to a dataframe with map2_df
# perform the modeling and assessment function on each CV split and fold
# use the 
res <- map2_df(
  .x = my_folds$splits,
  .y = my_folds$id,
  ~ execute_modeling(split = .x,
    id = .y,
		pd = pd,
		n = nCpG_val,
		rm_cor = TRUE,
		cor_cutoff = 0.75,
		mod_spec = mod_spec,
		method = "penalized_logreg",
		hyperparam_string = hyper_string,
		pred_results = pred_resultsv1)
)

# export results df
write.table(as.data.frame(res),
  paste0("results_", hyper_string, ".txt"),
	sep = "\t",
	col.names = T,
	row.names = F,
	quote = F
)

print("SESSION_INFO")
print(sessionInfo())