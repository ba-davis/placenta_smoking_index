
# conda activate tidymodels

library(tidyverse)
library(readxl)
library(ggrepel)

myfiles <- list.files("..", pattern="psi_scores.txt$", full=T)

dfs <- lapply(myfiles, read.delim, header=T)
names(dfs) <- c("model1_EPIC", "model2_450k", "model3_PACE")

#--------------------------------------#

# normalized psi

# add model name to norm_psi colname for each df
for (i in 1:length(dfs)) {
  colnames(dfs[[i]])[3] <- names(dfs)[i]
}

# combine dfs into one df
my_table <- Reduce(function(df1, df2) merge(df1, df2, by = "Sample_Name", all = TRUE), dfs)
my_table2 <- my_table[ ,c(1,10,3,6,9)]

# melt to long format
plotdf <- pivot_longer(my_table2, model1_EPIC:model3_PACE, names_to="model", values_to="norm_psi")
colnames(plotdf)[2] <- "status"

#-------------#

# Plot
p <- ggplot(plotdf, aes(x=model, y=norm_psi, fill=status)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Model", y="PSI", title="VCSIP 3 Models PSI") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size=9)) +
  theme(axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("combined_PSI_by_status_vcsip_boxplot.pdf", p, width=3.5, height=3.5, dpi=300)
