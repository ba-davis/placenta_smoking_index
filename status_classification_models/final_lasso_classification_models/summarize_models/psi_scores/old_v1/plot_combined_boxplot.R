
library(tidyverse)
library(readxl)
library(ggrepel)

myfiles <- list.files(getwd(), pattern="psi_vcsip.txt$", full=T)

dfs <- lapply(myfiles, read.delim, header=T)
names(dfs) <- c("model1_EPIC", "model2_450k", "model3_PACE")

# add model name to psi colname for each df
for (i in 1:length(dfs)) {
  colnames(dfs[[i]])[3] <- names(dfs)[i]
}

# combine dfs into one df
my_table <- Reduce(function(df1, df2) merge(df1, df2, by = "pd.Sample_Name", all = TRUE), dfs)
my_table2 <- my_table[ ,c(1,2,3,5,7)]

# melt to long format
plotdf <- pivot_longer(my_table2, model1_EPIC:model3_PACE, names_to="model", values_to="psi")

# label outliers
# outlier function to label outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

plotdf <- plotdf %>%
  group_by(model) %>%
  group_by(pd.status.x) %>%
  mutate(outlier = ifelse(is_outlier(psi), pd.Sample_Name, NA))

p <- ggplot(plotdf, aes(x=model, y=psi, fill=pd.status.x)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("combined_PSIscore_by_status_vcsip_boxplot_outliers.png", p, dpi=600)

p <- ggplot(plotdf, aes(x=model, y=psi, fill=pd.status.x)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI Score", title="VCSIP Data PSI Scores") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("combined_PSIscore_by_status_vcsip_boxplot.png", p, dpi=600)

#ggsave("combined_PSIscore_by_status_vcsip_boxplot_outliers.pdf", p)
