# VCSIP data
# Plot ROC AUC and partial AUC curves for
#   - prob_nonsmoker
#   - PSI
#   - norm PSI
# Use summary tables for each model that includes
#   prob_nonsmoker, psi, norm_psi for each sample

# conda activate tidymodels

library(tidymodels)
source("roc_plot_fxns.R")

#-------------------#
###########
# MODEL 1 #
###########

# read in model 1 vcsip summary table
df <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/",
  "fall_2022_final_models/summarize_models/summary_tables/",
  "model1_vcsip_summary_table.txt"), header = TRUE)

# prob_nonsmoker
plot_clean_roc(filename_prefix = "vcsip_model1_auc_prob_nonsmoker_plot",
    actual = df$status, predicted = df$prob_nonsmoker)
plot_clean_pauc(filename_prefix = "vcsip_model1_pauc_prob_nonsmoker_plot",
    actual = df$status, predicted = df$prob_nonsmoker)

# psi
plot_clean_roc(filename_prefix = "vcsip_model1_auc_psi_plot",
    actual = df$status, predicted = df$psi)
plot_clean_pauc(filename_prefix = "vcsip_model1_pauc_psi_plot",
    actual = df$status, predicted = df$psi)

# norm psi
plot_clean_roc(filename_prefix = "vcsip_model1_auc_norm_psi_plot",
    actual = df$status, predicted = df$norm_psi)
plot_clean_pauc(filename_prefix = "vcsip_model1_pauc_norm_psi_plot",
    actual = df$status, predicted = df$norm_psi)

#-------------------#
###########
# MODEL 2 #
###########

# read in model 2 vcsip summary table
df <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/",
  "fall_2022_final_models/summarize_models/summary_tables/",
  "model2_vcsip_summary_table.txt"), header = TRUE)

# prob_nonsmoker
plot_clean_roc(filename_prefix = "vcsip_model2_auc_prob_nonsmoker_plot",
    actual = df$status, predicted = df$prob_nonsmoker)
plot_clean_pauc(filename_prefix = "vcsip_model2_pauc_prob_nonsmoker_plot",
    actual = df$status, predicted = df$prob_nonsmoker)

# psi
plot_clean_roc(filename_prefix = "vcsip_model2_auc_psi_plot",
    actual = df$status, predicted = df$psi)
plot_clean_pauc(filename_prefix = "vcsip_model2_pauc_psi_plot",
    actual = df$status, predicted = df$psi)

# norm psi
plot_clean_roc(filename_prefix = "vcsip_model2_auc_norm_psi_plot",
    actual = df$status, predicted = df$norm_psi)
plot_clean_pauc(filename_prefix = "vcsip_model2_pauc_norm_psi_plot",
    actual = df$status, predicted = df$norm_psi)

#-------------------#
###########
# MODEL 2 #
###########

# read in model 3 vcsip summary table
df <- read.delim(paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/",
  "fall_2022_final_models/summarize_models/summary_tables/",
  "model3_vcsip_summary_table.txt"), header = TRUE)

# prob_nonsmoker
plot_clean_roc(filename_prefix = "vcsip_model3_auc_prob_nonsmoker_plot",
    actual = df$status, predicted = df$prob_nonsmoker)
plot_clean_pauc(filename_prefix = "vcsip_model3_pauc_prob_nonsmoker_plot",
    actual = df$status, predicted = df$prob_nonsmoker)

# psi
plot_clean_roc(filename_prefix = "vcsip_model3_auc_psi_plot",
    actual = df$status, predicted = df$psi)
plot_clean_pauc(filename_prefix = "vcsip_model3_pauc_psi_plot",
    actual = df$status, predicted = df$psi)

# norm psi
plot_clean_roc(filename_prefix = "vcsip_model3_auc_norm_psi_plot",
    actual = df$status, predicted = df$norm_psi)
plot_clean_pauc(filename_prefix = "vcsip_model3_pauc_norm_psi_plot",
    actual = df$status, predicted = df$norm_psi)
