

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
  "lyndsey_shorey_project/classification_models/EPIC_array_data/lasso/try1")

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

# export the cv hyper tuning performance metrics
write.table(metric_df,
  "epic_lasso_hyper_tuning_metrics.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

# get best performing hyperparams
lasso_kap_metrics <- get_best_hps(metric_df, "kap", 5)
# Ensure no duplicates
lasso_kap_metrics <- lasso_kap_metrics[!duplicated(lasso_kap_metrics[c(1,2)]),]

# export lasso best kap table
#lasso_kap_metrics$lambda <- as.character(lasso_kap_metrics$lambda)
write.table(lasso_kap_metrics, "lasso_kap_metrics.txt", sep="\t", col.names=T, row.names=F, quote=F)

#-------------------------------------------------------------------------#

###############
# Elastic Net #
###############

# path to folder containing performance metrics txt files
dir <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/EPIC_array_data/",
  "elastic_net/try1")

# read in all perormance metric txt files
#mydf <- vroom(list.files(dir, pattern="^results_", full=TRUE))
# read in all perormance metric txt files
# note vroom was getting "too many open files error"
myfiles <- list.files(dir, pattern="^results_", full=TRUE)
for (i in 1:length(myfiles)) {
  if (i==1) {
    mydf <- as_tibble(read.delim(myfiles[i], header=T))
  }
  else if (i > 1) {
    newdf <- as_tibble(read.delim(myfiles[i], header=T))
    mydf <- rbind(mydf, newdf)
  }
}

outplot_name <- "elastic_net_hyperparam_heatmap.pdf"
outfile_name <- "elastic_net_hyperparam_heatmap.txt"

# Crete the plot of 3 hyperparams heatmap
metric_df <- plot_three_hyperparams_heatmap(mytib=mydf,
                                            p1="lambda",
					    p2="alpha",
					    p3="n_CpG",
					    metric="kap",
					    num_metrics=3,
					    plotfile=outplot_name,
					    output_df=TRUE,
					    outfile=outfile_name
)

# export the cv hyper tuning performance metrics
write.table(metric_df,
  "epic_elasticnet_hyper_tuning_metrics.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

# get best performing hyperparams
enet_kap_metrics <- get_best_hps(metric_df, "kap", 5)
# Ensure no duplicates
enet_kap_metrics <- enet_kap_metrics[!duplicated(enet_kap_metrics[c(1,2,3)]),]

# get accuracy metrics

outplot_name <- "elastic_net_hyperparam_accuracy_heatmap.pdf"
outfile_name <- "elastic_net_hyperparam_accuracy_heatmap.txt"

# Crete the plot of 3 hyperparams heatmap
metric_df <- plot_three_hyperparams_heatmap(mytib=mydf,
                                            p1="lambda",
                                            p2="alpha",
                                            p3="n_CpG",
                                            metric="accuracy",
                                            num_metrics=3,
                                            plotfile=outplot_name,
                                            output_df=TRUE,
                                            outfile=outfile_name
)

# export the cv hyper tuning performance metrics
write.table(metric_df,
  "epic_elasticnet_hyper_tuning_accuracy_metrics.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

# get roc_auc metrics

outplot_name <- "elastic_net_hyperparam_rocauc_heatmap.pdf"
outfile_name <- "elastic_net_hyperparam_rocauc_heatmap.txt"

# Crete the plot of 3 hyperparams heatmap
metric_df <- plot_three_hyperparams_heatmap(mytib=mydf,
                                            p1="lambda",
                                            p2="alpha",
                                            p3="n_CpG",
                                            metric="roc_auc",
                                            num_metrics=3,
                                            plotfile=outplot_name,
                                            output_df=TRUE,
                                            outfile=outfile_name
)

# export the cv hyper tuning performance metrics
write.table(metric_df,
  "epic_elasticnet_hyper_tuning_rocauc_metrics.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

#-------------------------------------------------------------------------#

#################
# Random Forest #
#################

# First round of hyperparameters

# path to folder containing performance metrics txt files
dir <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/EPIC_array_data/",
  "random_forest/trya")

# read in all perormance metric txt files
mydf <- vroom(list.files(dir, pattern="^results_", full=TRUE))

# Plot performance metric average for each hyperparameter combo

outplot_name <- "rf_a_hyperparam_heatmap.pdf"
outfile_name <- "rf_a_hyperparam_heatmap.txt"

metric_df <- plot_metric_heatmap(df=mydf,
                                 param1="nCpG",
                                 param2="ntree",
                                 num_metrics=3,
                                 plotfile=outplot_name,
                                 output_df=TRUE,
                                 outfile=outfile_name
)

# export the cv hyper tuning performance metrics
write.table(metric_df,
  "epic_rf1_hyper_tuning_metrics.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

# get best performing hyperparams
rf_a_kap_metrics <- get_best_hps(metric_df, "kap", 5)
# Ensure no duplicates
rf_a_kap_metrics <- rf_a_kap_metrics[!duplicated(rf_a_kap_metrics[c(1,2)]),]

# Second round of hyperparameters

# path to folder containing performance metrics txt files
dir <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/EPIC_array_data/",
  "random_forest/tryb")

# read in all perormance metric txt files
mydf <- vroom(list.files(dir, pattern="^results_", full=TRUE))

# Plot performance metric average for each hyperparameter combo

outplot_name <- "rf_b_hyperparam_heatmap.pdf"
outfile_name <- "rf_b_hyperparam_heatmap.txt"

metric_df <- plot_metric_heatmap(df=mydf,
                                 param1="mtry",
                                 param2="min_n",
                                 num_metrics=3,
                                 plotfile=outplot_name,
                                 output_df=TRUE,
                                 outfile=outfile_name
)

# export the cv hyper tuning performance metrics
write.table(metric_df,
  "epic_rf2_hyper_tuning_metrics.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

# get best performing hyperparams
rf_b_kap_metrics <- get_best_hps(metric_df, "kap", 5)
# Ensure no duplicates
rf_b_kap_metrics <- rf_b_kap_metrics[!duplicated(rf_b_kap_metrics[c(1,2)]),]

#-------------------------------------------------------------------------#

#############################
# Gradient Boosting Machine #
#############################

# path to folder containing performance metrics txt files
dir <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/EPIC_array_data/",
  "gradient_boost/try1")

# read in all perormance metric txt files
# note vroom was getting "too many open files error"
myfiles <- list.files(dir, pattern="^results_", full=TRUE)
for (i in 1:length(myfiles)) {
  if (i==1) {
    mydf <- as_tibble(read.delim(myfiles[i], header=T))
  }
  else if (i > 1) {
    newdf <- as_tibble(read.delim(myfiles[i], header=T))
    mydf <- rbind(mydf, newdf)
  }
}

# Plot performance metric average for each hyperparameter combo

outplot_name <- "gbm_hyperparam_heatmap.pdf"
outfile_name <- "gbm_hyperparam_heatmap.txt"

# Crete the plot of 3 hyperparams heatmap
metric_df <- plot_three_hyperparams_heatmap(mytib=mydf,
                                            p1="depth",
                                            p2="ntree",
                                            p3="nCpG",
                                            metric="kap",
                                            num_metrics=3,
                                            plotfile=outplot_name,
                                            output_df=TRUE,
                                            outfile=outfile_name
)

# export the cv hyper tuning performance metrics
write.table(metric_df,
  "epic_gbm_hyper_tuning_metrics.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

# get best performing hyperparams
gbm_kap_metrics <- get_best_hps(metric_df, "kap", 5)
# Ensure no duplicates
gbm_kap_metrics <- gbm_kap_metrics[!duplicated(gbm_kap_metrics[c(1,2,3)]),]

# get accuracy metrics

outplot_name <- "gbm_hyperparam_accuracy_heatmap.pdf"
outfile_name <- "gbm_hyperparam_accuracy_heatmap.txt"

# Crete the plot of 3 hyperparams heatmap
metric_df <- plot_three_hyperparams_heatmap(mytib=mydf,
                                            p1="depth",
                                            p2="ntree",
                                            p3="nCpG",
                                            metric="accuracy",
                                            num_metrics=3,
                                            plotfile=outplot_name,
                                            output_df=TRUE,
                                            outfile=outfile_name
)

# export the cv hyper tuning performance metrics
write.table(metric_df,
  "epic_gbm_hyper_tuning_accuracy_metrics.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

# get roc_auc metrics

outplot_name <- "gbm_hyperparam_rocauc_heatmap.pdf"
outfile_name <- "gbm_hyperparam_rocauc_heatmap.txt"

# Crete the plot of 3 hyperparams heatmap
metric_df <- plot_three_hyperparams_heatmap(mytib=mydf,
                                            p1="depth",
                                            p2="ntree",
                                            p3="nCpG",
                                            metric="roc_auc",
                                            num_metrics=3,
                                            plotfile=outplot_name,
                                            output_df=TRUE,
                                            outfile=outfile_name
)

# export the cv hyper tuning performance metrics
write.table(metric_df,
  "epic_gbm_hyper_tuning_rocauc_metrics.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

#-------------------------------------------------------------------------#

# output the performance metrics

# create list of dfs
perfs <- list(lasso_kap_metrics, enet_kap_metrics, rf_a_kap_metrics, rf_b_kap_metrics, gbm_kap_metrics)
names(perfs) <- c("lasso", "enet", "rf_a", "rf_b", "gbm")

write_xlsx(perfs, "EPIC_hyperparam_res.xlsx")
