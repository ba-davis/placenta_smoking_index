
# conda activate tidymodels

library(tidyverse)
library(ggrepel)
source("psi_v_cot_plot_fxns.R")

# read in the pd file
pd <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/prelim_data/cleaned_data/pd.txt"))

# read in the cotinine metadata
cotinine_dat <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
    "lyndsey_shorey_project/metadata/Cotinine_ng_ml_pd.txt"), header = TRUE)
colnames(cotinine_dat)[1] <- "Sample_Name"

# path to the d2 file showing psi score per sample
#mypath <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
# "lyndsey_shorey_project/classification_models/fall_2022_final_models/",
 # "summarize_models/psi_scores")

#--------------------------------------------------------------------#

# model 1
plot_psi_v_cot(infile = paste0("../../psi_scores/",
    "model1_vcsip_psi_scores_and_norm_psi_scores.txt"),
    pd, cotinine_dat,
    mod_num = "model1_EPIC",
    y_rho = -0.8,
    y_pval = -1)

# model 2
plot_psi_v_cot(infile = paste0("../../psi_scores/",
    "model2_vcsip_psi_scores_and_norm_psi_scores.txt"),
    pd, cotinine_dat,
    mod_num = "model2_450k",
    y_rho = -1,
    y_pval = -1.2)

# model 3
plot_psi_v_cot(infile = paste0("../../psi_scores/",
    "model3_vcsip_psi_scores_and_norm_psi_scores.txt"),
    pd, cotinine_dat,
    mod_num = "model3_PACE",
    y_rho = -2.1,
    y_pval = -2.4)
