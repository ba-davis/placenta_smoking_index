
# conda activate tidymodels

library(tidyverse)
library(readxl)
library(ggrepel)

myfiles <- list.files("..", pattern="psi_scores.txt$", full=T)

dfs <- lapply(myfiles, read.delim, header=T)
names(dfs) <- c("model1_EPIC", "model2_450k", "model3_PACE")

#---------------------------------------#

# label outliers

# outlier function to label outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

#-------------------------------#
# Model 1

# add outlier column to df
d2 <- dfs[[1]]
colnames(d2)[1] <- "Sample_Name"

d2 <- d2 %>%
  group_by(status) %>%
  mutate(outlier = ifelse(is_outlier(norm_psi), Sample_Name, NA))

p <- ggplot(d2, aes(x=status, y=norm_psi, fill=status)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI", title="VCSIP Model1 EPIC PSI") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size=11)) +
  theme(axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=17))
ggsave("model1_PSI_by_status_vcsip_boxplot_outliers.pdf", p, width=3.5, height=3.5, dpi=300)

p <- ggplot(d2, aes(x=status, y=norm_psi, fill=status)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI", title="VCSIP Model1 EPIC PSI") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size=11)) +
  theme(axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=17))
ggsave("model1_PSI_by_status_vcsip_boxplot.pdf", p, width=3.5, height=3.5, dpi=300)

#-----------------------------#
# Model 2

# add outlier column
d2 <- dfs[[2]]
colnames(d2)[1] <- "Sample_Name"

d2 <- d2 %>%
  group_by(status) %>%
  mutate(outlier = ifelse(is_outlier(norm_psi), Sample_Name, NA))

p <- ggplot(d2, aes(x=status, y=norm_psi, fill=status)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI", title="VCSIP Model2 450k PSI") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size=11)) +
  theme(axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=17))
ggsave("model2_PSI_by_status_vcsip_boxplot_outliers.pdf", p, width=3.5, height=3.5, dpi=300)

p <- ggplot(d2, aes(x=status, y=norm_psi, fill=status)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI", title="VCSIP Model2 450k PSI") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size=11)) +
  theme(axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=17))
ggsave("model2_PSI_by_status_vcsip_boxplot.pdf", p, width=3.5, height=3.5, dpi=300)

#-----------------------------#
# Model 3

# add outlier
d2 <- dfs[[3]]
colnames(d2)[1] <- "Sample_Name"

d2 <- d2 %>%
  group_by(status) %>%
  mutate(outlier = ifelse(is_outlier(norm_psi), Sample_Name, NA))

p <- ggplot(d2, aes(x=status, y=norm_psi, fill=status)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  theme_classic() +
  labs(x="Smoking Status", y="PSI", title="VCSIP Model3 PACE PSI") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size=11)) +
  theme(axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=17))
ggsave("model3_PSI_by_status_vcsip_boxplot_outliers.pdf", p, width=3.5, height=3.5, dpi=300)

p <- ggplot(d2, aes(x=status, y=norm_psi, fill=status)) +
  geom_boxplot() +
  theme_classic() +
  labs(x="Smoking Status", y="PSI", title="VCSIP Model3 PACE PSI") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size=11)) +
  theme(axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=17))
ggsave("model3_PSI_by_status_vcsip_boxplot.pdf", p, width=3.5, height=3.5, dpi=300)
