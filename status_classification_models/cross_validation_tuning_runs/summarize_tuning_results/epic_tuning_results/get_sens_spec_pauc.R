

# conda activate tidymodels

#-------------------------------------------------------------------------#

# Load Libraries
library(tidymodels)
library(vroom)
library(writexl)



# path to folder containing performance metrics txt files
dir <- paste0("/home/groups/hoolock2/u0/bd/Projects/",
  "lyndsey_shorey_project/classification_models/EPIC_array_data/lasso/",
  "try1_again_collect_more_metrics")

# read in all perormance metric txt files
# Use "col_types" parameter for vroom to set the lambda column to "character"
#   this is to avoid truncating the long decimal when reading into R
mydf <- vroom(list.files(dir, pattern="^results_", full=TRUE), col_types="dcccccdddddddd")
mydf.tib <- as_tibble(mydf)

param1 <- "n_CpG"
param2 <- "lambda"
# create unique id (combo of the two hyperparameters and fold)
mydf.tib <- mydf.tib %>%
  add_column(unqID=paste0(mydf.tib[[param1]],
  "_", mydf.tib[[param2]], "_", mydf.tib$id))

# create tibble of metrics per lambda/alpha/fold
myres <- mydf.tib %>%
  group_by(unqID) %>%
    summarize(kap = unique(kap),
    accuracy = unique(accuracy),
    roc_auc = unique(roc_auc),
    sensitivity = unique(sensitivity),
    specificity = unique(specificity),
    pauc = unique(pauc))

# create unique key for just the two hyperparameters
myres$param_key <- gsub("_Fold[0-9]*", "", myres$unqID)

# create standard error function
std <- function(x) sd(x)/sqrt(length(x))

# Create metric plot df for plotting
foo <- myres %>%
  group_by(param_key) %>%
  summarize(mean_kap = mean(kap),
    std_kap = std(kap),
    mean_accuracy = mean(accuracy),
    std_accuracy = std(accuracy),
    mean_roc_auc = mean(roc_auc),
    std_roc_auc = std(roc_auc),
    mean_sens = mean(sensitivity),
    std_sens = std(sensitivity),
    mean_spec = mean(specificity),
    std_spec = std(specificity),
    mean_pauc = mean(pauc),
    std_pauc = std(pauc))

# create a numeric columns containing the value of each hyperparameter
foo <- foo %>%
  separate(param_key, c(param1, param2), sep="_")

num_metrics <- 6

# Create plot df
plotdf <- data.frame(p1=as.factor(as.numeric(c(rep(foo[[param1]], num_metrics)))),
  p2=as.factor(as.numeric(c(rep(foo[[param2]], num_metrics)))),
  metric=c(rep("kap", nrow(foo)), rep("accuracy", nrow(foo)),
    rep("roc_auc", nrow(foo)), rep("sensitivity", nrow(foo)),
    rep("specificity", nrow(foo)), rep("pauc", nrow(foo))),
  mean=c(foo$mean_kap, foo$mean_accuracy, foo$mean_roc_auc,
    foo$mean_sens, foo$mean_spec, foo$mean_pauc),
  std=c(foo$std_kap, foo$std_accuracy, foo$std_roc_auc,
    foo$std_sens, foo$std_spec, foo$std_pauc))

# export
write.table(plotdf,
  "epic_cv_lasso_6_metrics.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)
