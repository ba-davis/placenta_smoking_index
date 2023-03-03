
# conda activate tidymodels

library(tidymodels)

# read in summary tables
mod1 <- read.delim("../summary_tables/model1_vcsip_summary_table.txt")
mod2 <- read.delim("../summary_tables/model2_vcsip_summary_table.txt")
mod3 <- read.delim("../summary_tables/model3_vcsip_summary_table.txt")

# get confusion matrix
mod1_conf <- conf_mat(mod1, status, prediction)
mod2_conf <- conf_mat(mod2, status, prediction)
mod3_conf <- conf_mat(mod3, status, prediction)
