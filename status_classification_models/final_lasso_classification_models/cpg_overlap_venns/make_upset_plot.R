
# conda activate rvenns

library(tidyverse)
library(UpSetR)
library(readxl)

# read in excel file of shared cpgs across 5 models
mydat <- read_xlsx("../all_model_cpgs.xlsx")

colnames(mydat)[1] <- "term"

# Get vectors of CpGs with nonzero coef in each model
mod1 <- mydat %>%
  drop_na(model1_coef) %>%
  filter(model1_coef != 0) %>%
  select(term) %>%
  unlist(., use.names=FALSE)

mod2 <- mydat %>%
  drop_na(model2_coef) %>%
  filter(model2_coef != 0) %>%
  select(term) %>%
  unlist(., use.names=FALSE)

mod3 <- mydat %>%
  drop_na(model3_coef) %>%
  filter(model3_coef != 0) %>%
  select(term) %>%
  unlist(., use.names=FALSE)

mod4 <- mydat %>%
  drop_na(model4_coef) %>%
  filter(model4_coef != 0) %>%
  select(term) %>%
  unlist(., use.names=FALSE)

mod5 <- mydat %>%
  drop_na(model5_coef) %>%
  filter(model5_coef != 0) %>%
  select(term) %>%
  unlist(., use.names=FALSE)

# define input list of vectors
input_list <- list(mod1, mod2, mod3, mod4, mod5)

# define category names
cat_names <- c("Model_1", "Model_2", "Model_3", "Model_4", "Model_5")
title <- "Shared CpGs among 5 models"
names(input_list) <- cat_names

png("shared_cpgs_upset_plot.png", width=5, height=5, units='in', res=600)
upset(fromList(input_list), order.by = "freq")
dev.off()
