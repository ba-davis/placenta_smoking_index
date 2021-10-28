#!/home/groups/hoolock2/u0/bd/miniconda3/envs/tidymodels/bin/Rscript

library(tidyverse)
library(vroom)



#-------------------------------------------------------------------------------------------------------------------#

# Source functions
source("/home/groups/hoolock2/u0/bd/Projects/lyndsey_shorey_project/smoking_index_github/placenta_smoking_index/status_classification_models/general_shared_functions.R")

#-------------------------------------------------------------------------------------------------------------------#

# path to folder containing performance metrics txt files
dir="/home/groups/hoolock2/u0/bd/Projects/lyndsey_shorey_project/classification_models/EPIC_array_data/elastic_net/try1"

# read in all perormance metric df's
mytib <- vroom(list.files(dir, pattern="^results_", full=TRUE))

outplot_name <- "elastic_net_hyperparam_heatmap.png"
outfile_name <- "elastic_net_hyperparam_heatmap.txt"

# Crete the plot of 3 hyperparams heatmap
metric_df <- plot_three_hyperparams_heatmap(mytib=mytib,
                                            p1="lambda",
					    p2="alpha",
					    p3="n_CpG",
					    metric="kap",
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
