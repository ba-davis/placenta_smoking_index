#!/home/groups/hoolock2/u0/bd/miniconda3/envs/tidymodels/bin/Rscript

library(tidyverse)
library(vroom)



#-------------------------------------------------------------------------------------------------------------------#

# Source functions
source("/home/groups/hoolock2/u0/bd/Projects/lyndsey_shorey_project/smoking_index_github/placenta_smoking_index/status_classification_models/general_shared_functions.R")

#-------------------------------------------------------------------------------------------------------------------#

# path to folder containing performance metrics txt files
#dir="/home/groups/hoolock2/u0/bd/Projects/lyndsey_shorey_project/models_450k_only/lasso/new4/performance_metrics"
dir = getwd()

# read in all perormance metric df's
mydf <- vroom(list.files(dir, pattern="^results_", full=TRUE))


# Plot performance metric average for each hyperparameter combo

outplot_name <- "lasso_hyperparam_heatmap.png"
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
kap_metrics <- get_best_hps(metric_df, "kap", 5)

write.table(kap_metrics,
            "best_kap_hyperparams_table.txt",
	    sep="\t",
	    col.names=T,
	    row.names=F,
	    quote=F
)
