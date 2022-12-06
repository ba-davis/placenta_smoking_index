
# conda activate tidymodels

#-------------------------------------------------------------------------#

# Load Libraries
library(tidymodels)
library(vroom)
library(writexl)

#-------------------------------------------------------------------------#

# source functions
source(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/smoking_index_github/placenta_smoking_index/",
  "status_classification_models/general_shared_functions.R")
)

#-------------------------------------------------------------------------#

#########
# LASSO #
#########

# path to folder containing performance metrics txt files
dir <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/using_pmid_34429407_cpgs/",
  "lasso/try1/performance_metrics")

# read in all perormance metric txt files
#mydf <- vroom(list.files(dir, pattern="^results_", full=TRUE))

# read in all perormance metric txt files
# Use "col_types" parameter for vroom to set the lambda column to "character"
#   this is to avoid truncating the long decimal when reading into R
mydf <- vroom(list.files(dir, pattern="^results_", full=TRUE), col_types="dcccccddddd")


# Plot performance metric average for each hyperparameter combo

outplot_name <- "lasso_hyperparam_heatmap.pdf"
outfile_name <- "lasso_hyperparam_heatmap.txt"

metric_df <- plot_metric_heatmap(df=mydf,
                                 param1="n_CpG",
		                 param2="lambda",
		                 num_metrics=3,
		                 plotfile=outplot_name,
		                 output_df=TRUE,
		                 outfile=outfile_name
)

# get best performing hyperparams
lasso_kap_metrics <- get_best_hps(metric_df, "kap", 5)
# Ensure no duplicates
lasso_kap_metrics <- lasso_kap_metrics[!duplicated(lasso_kap_metrics[c(1,2)]),]

# export lasso best kap table
#lasso_kap_metrics$lambda <- as.character(lasso_kap_metrics$lambda)
write.table(lasso_kap_metrics, "lasso_kap_metrics.txt", sep="\t", col.names=T, row.names=F, quote=F)
