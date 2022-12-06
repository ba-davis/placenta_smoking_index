
# conda activate rvenns

library(tidyverse)
library(VennDiagram)
library(UpSetR)
library(readxl)
#library(writexl)

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

# calculate overlaps
overlap <- calculate.overlap(input_list)

# venn
venn.plot <- draw.triple.venn(area1 = length(mod1),
  area2 = length(mod2),
  area3 = length(mod3),
  n12 = length(overlap[["a2"]]) + length(overlap[["a5"]]),
  n23 = length(overlap[["a6"]]) + length(overlap[["a5"]]),
  n13 = length(overlap[["a4"]]) + length(overlap[["a5"]]),
  n123 = length(overlap[["a5"]]),
  category=cat_names,
  euler.d=FALSE,
  scaled=FALSE,
  fill=c("red", "yellow", "blue"),
  cex=3,
  cat.fontface=2,
  lty = "blank",
  main=title)

pdf("shared_cpgs_5models_venn.pdf")
grid.draw(venn.plot)
dev.off()




# actual numbers
#venn.plot <- draw.triple.venn(area1 = 18,
  area2 = 21,
  area3 = 18,
  n12 =6,
  n23 = 8,
  n13 = 5,
  n123 = 5,
  category=cat_names,
  euler.d=FALSE,
  scaled=FALSE,
  fill=c("red", "yellow", "blue"),
  cex=3,
  cat.fontface=2,
  #cat.pos=cat_pos,
  #ext.percent=ext_percent,
  lty = "blank",
  main=title)
