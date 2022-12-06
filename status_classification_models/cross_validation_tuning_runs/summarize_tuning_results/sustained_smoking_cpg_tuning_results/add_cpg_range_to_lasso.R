

# conda activate tidymodels

#-------------------------------------------------------------------------#

# Load Libraries
library(tidymodels)
library(vroom)
library(writexl)

#-------------------------------------------------------------------------#

# path to feature/predictor files
dir <- "/home/groups/hoolock2/u0/bd/Projects/lyndsey_shorey_project/classification_models/using_pmid_34429407_cpgs/lasso/try1/var_importance"

#mydf <- vroom(list.files(dir, pattern=".txt$", full=TRUE))
# note vroom was getting "too many open files error"
myfiles <- list.files(dir, pattern=".txt$", full=TRUE)
for (i in 1:length(myfiles)) {
  name <- gsub(paste0(dir, "/"), "", myfiles[i])
  name <- gsub("_nonzero_predictors.txt", "", name)
  if (i==1) {
    mydf <- as_tibble(read.delim(myfiles[i], header=T))
    mydf$name <- name
  }
  else if (i > 1) {
    newdf <- as_tibble(read.delim(myfiles[i], header=T))
    newdf$name <- name
    mydf <- rbind(mydf, newdf)
  }
}

# remove the intercept rows
#mydf2 <- mydf[!(mydf$term=="(Intercept)"), ]

#----------------------------------------------------------#

test <- as_tibble(as.data.frame(table(mydf$name)))
test$combo <- gsub("_Fold[0-9][0-9]", "", test$Var1)

# get the min "Freq" value for each combo
test2 <- test %>%
  group_by(combo) %>%
  slice(which.min(Freq))

# get the max "Freq" value for each combo
test3 <- test %>%
  group_by(combo) %>%
  slice(which.max(Freq))

# newdf
newdf <- data.frame(combo=test2$combo,
  min.cpg=test2$Freq,
  max.cpg=test3$Freq
)

# subtract 1 from each min and max column to remove the "intercept" count of that combo
newdf$min.cpg <- newdf$min.cpg - 1
newdf$max.cpg <- newdf$max.cpg - 1

#----------------------------------------------------------#

# read in existing kap metrics file
metrics <- read.delim("lasso_kap_metrics.txt", header=T)
metrics$combo <- paste0("nCpG_", metrics$n_CpG, "_lambda_", metrics$lambda)

# Merge
foo <- merge(newdf, metrics, by="combo")
foo2 <- foo[ ,c(4,5,6,7,8,2,3)]
# sort by mean kap value
foo3 <- foo2[order(foo2$mean, decreasing=T), ]

write.table(foo3, "lasso_kaps_and_nonzero_cpg_range_table.txt", sep="\t", col.names=T, row.names=F, quote=F)
