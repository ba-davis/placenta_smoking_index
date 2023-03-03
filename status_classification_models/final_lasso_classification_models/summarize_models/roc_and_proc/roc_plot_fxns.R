
library(tidymodels)

# Function to plot ROC plot
plot_clean_roc <- function(filename_prefix, width = 3.5, height = 3.5,
    actual, predicted) {

    # plot ROC plot
    pdf(paste0(filename_prefix, ".pdf"), width = width, height = height)
    pROC::plot.roc(actual, predicted,
        percent = TRUE,
        print.auc = TRUE,
        print.auc.col = "white",
        print.auc.cex = 0.9,
        print.auc.x = 50,
        print.auc.y = 10,
        auc.polygon = TRUE,
        auc.polygon.col = "#1c61b6",
        max.auc.polygon = TRUE,
        max.auc.polygon.col = "#1c61b622",
        main = "AUC")
    dev.off()
}


# Function to plot partial ROC plot
plot_clean_pauc <- function(filename_prefix, width = 3.5, height = 3.5,
    actual, predicted) {

    # plot pROC plot
    pdf(paste0(filename_prefix, ".pdf"), width = width, height = height)
    pROC::plot.roc(actual, predicted,       # data
        percent = TRUE,                     # show all values in percent
        partial.auc = c(100, 90),
        partial.auc.correct = TRUE,         # define a partial AUC (pAUC)
        print.auc = TRUE,
        #display pAUC value on the plot with following options:
        print.auc.pattern = "Corrected pAUC\n(100-90%% SP):\n%.1f%%",
        print.auc.col = "#1c61b6",
        print.auc.cex = 0.9,
        print.auc.x = 55,
        print.auc.y = 25,
        auc.polygon = TRUE,
        auc.polygon.col = "#1c61b6",       # show pAUC as a polygon
        max.auc.polygon = TRUE,
        max.auc.polygon.col = "#1c61b622", # also show the 100% polygon
        main = "Partial AUC (pAUC)")
    dev.off()
}
