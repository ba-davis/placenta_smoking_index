
# Get ROC Curve and AUC ROC
# Using Tidymodels (yardstick)
# "roc_curve" and "roc_auc" functions


library(tidymodels)
library(tidyverse)

# Define function to make this quicker
plot_yardstick_roc <- function(infile, outfile="roc_curve.png", title=NULL) {
  # read in file with truth and predicted probability columns
  mod_res <- read.delim(infile, header=T)

  # make sure truth column is a factor
  mod_res$truth <- as.factor(mod_res$truth)

  # calculate the roc curve and plot the roc curve
  p <- roc_curve(mod_res, truth, prob_nonsmoker) %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity)) +
    geom_path() +
    geom_abline(lty = 3) +
    coord_equal() +
    theme_bw() +
    ggtitle(title)
  ggsave(filename=outfile, p)

  # calculate the global auc_roc
  global_auc_roc <- roc_auc(mod_res, truth, "prob_nonsmoker")$.estimate

  # calculate the mean auc_roc across folds
  mean_auc_roc <- mean(mod_res[!(duplicated(mod_res$id)) ,
    c(3,11)]$roc_auc)

  #summary_res <- data.frame(global_auc_roc = global_auc_roc,
  #  mean_auc_roc = mean_auc_roc)

  res <- c(global_auc_roc, mean_auc_roc)
  return(res)
}

#----------------------------------------------#

# Read in the prediction results

# model 1
infile1 <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/EPIC_array_data/lasso/",
  "try1/results_nCpG_450_lambda_0.000206913808111479.txt")

# model 2
infile2 <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/EPIC_array_data/lasso/",
  "try1/results_nCpG_1000_lambda_0.297635144163131.txt")

# model 3
infile3 <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/450k_array_data/lasso/",
  "try1/results_nCpG_120_lambda_1.12883789168469e-09.txt")

# model 4
infile4 <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/using_pmid_34429407_cpgs/",
  "lasso/try1/performance_metrics/",
  "results_nCpG_98_lambda_3.35981828628379e-10.txt")

# model 5
infile5 <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/EPIC_array_data/lasso/",
  "try1/results_nCpG_120_lambda_0.0885866790410083.txt")

#---------------------------#

# Plot with function
mod1_auc <- plot_yardstick_roc(infile1, "model1_vcsip_roc_curve.pdf", "Model 1 ROC")
mod2_auc <- plot_yardstick_roc(infile2, "model2_vcsip_roc_curve.pdf", "Model 2 ROC")
mod3_auc <- plot_yardstick_roc(infile3, "model3_vcsip_roc_curve.pdf", "Model 3 ROC")
mod4_auc <- plot_yardstick_roc(infile4, "model4_vcsip_roc_curve.pdf", "Model 4 ROC")
mod5_auc <- plot_yardstick_roc(infile5, "model5_vcsip_roc_curve.pdf", "Model 5 ROC")

#----------------------------#

# combine the auc dfs for each model
newdf <- as.data.frame(t(data.frame(model1=mod1_auc,
  model2=mod2_auc,
  model3=mod3_auc,
  model4=mod4_auc,
  model5=mod5_auc)))

newdf <- newdf[ ,c(3,1,2)]
colnames(newdf) <- c("model", "global_auc_roc", "mean_auc_roc")

write.table(newdf,
  "final_models_cv_auc_roc_table.txt",
  sep="\t",
  col.names=T,
  row.names=T)
