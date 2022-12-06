
# conda activate rvenns

library(tidyverse)
library(UpSetR)
library(readxl)


# read in excel file of shared cpgs across 5 models
mydat <- read_xlsx("../cpg_tables/all_model_cpgs_annot.xlsx")

colnames(mydat)[1] <- "term"

# Get vectors of CpGs with nonzero coef in each model
mod1 <- mydat %>%
  drop_na(model1_EPIC_coef) %>%
  filter(model1_EPIC_coef != 0) %>%
  select(term) %>%
  unlist(., use.names=FALSE)

mod2 <- mydat %>%
  drop_na(model2_450k_coef) %>%
  filter(model2_450k_coef != 0) %>%
  select(term) %>%
  unlist(., use.names=FALSE)

mod3 <- mydat %>%
  drop_na(model3_PACE_coef) %>%
  filter(model3_PACE_coef != 0) %>%
  select(term) %>%
  unlist(., use.names=FALSE)

# define input list of vectors
input_list <- list(mod1, mod2, mod3)

# define category names
cat_names <- c("Model1_EPIC", "Model2_450k", "Model3_PACE")
title <- "Shared CpGs among 3 models"
names(input_list) <- cat_names

png("shared_cpgs_upset_plot.png", width=5, height=5, units='in', res=600)
upset(fromList(input_list), order.by = "freq")
dev.off()

