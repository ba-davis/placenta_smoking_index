

#library(ggplot2)

# plot psi v cotinine function
plot_psi_v_cot <- function(infile, pd, cotinine_dat, mod_num,
    y_rho, y_pval, width = 3.5, height = 3.5) {

    d2 <- read.delim(infile, header = TRUE)
    colnames(d2)[1] <- "Sample_Name"
    colnames(d2)[4] <- "smoking_status"
    colnames(d2)[2] <- "PSI"

    # put order of samples in d2 in same order as samples in pd
    d2 <- d2[match(pd$Sample_Name, d2$Sample_Name), ]

    # add cotinine data to d2
    d2$Rand_Cotinine_ng_ml <- cotinine_dat$Rand_Cotinine_ng_ml
    d2$Mid_Cotinine_ng_ml <- cotinine_dat$Mid_Cotinine_ng_ml
    d2$Late_Cotinine_ng_ml <- cotinine_dat$Late_Cotinine_ng_ml

    #-------------------------------------------------------------------#
    # Rand Cotinine

    # get pearson correlation coefficient and pval
    rho <- paste0("r = ", round(cor.test(d2$PSI, d2$Rand_Cotinine_ng_ml,
        method = "pearson")$estimate, 4))
    pval <- paste0("p = ", round(cor.test(d2$PSI, d2$Rand_Cotinine_ng_ml,
        method = "pearson")$p.value, 7))

    # plot PSI vs Rand cotinine
    p <- ggplot(d2, aes(x = Rand_Cotinine_ng_ml, y = PSI,
        color = smoking_status)) +
        geom_point(size = 3) +
        geom_smooth(method = lm) +
        annotate("text", x = 50, y = y_rho, label = rho, size = 5) +
        annotate("text", x = 70, y = y_pval, label = pval, size = 5) +
        theme_classic() +
        labs(x = "Rand Cotinine ng/ml",
            y = "PSI Score", title = paste0(mod_num, " VCSIP")) +
        theme(plot.title = element_text(hjust = 0.5))+
        theme(axis.text.x = element_text(size = 13)) +
        theme(axis.text.y = element_text(size = 13)) +
        theme(axis.title = element_text(size = 13)) +
        guides(color = FALSE)
    ggsave(paste0(mod_num, "_PSIscore_by_randCot_vcsip_scatterplot.pdf"), p,
        width = width, height = height, dpi = 300)

    #-----------#
    # Mid Cotinine

    rho <- paste0("r = ", round(cor.test(d2$PSI, d2$Mid_Cotinine_ng_ml,
        method = "pearson")$estimate, 4))
    pval <- paste0("p = ", round(cor.test(d2$PSI, d2$Mid_Cotinine_ng_ml,
        method = "pearson")$p.value, 7))

    # plot PSI vs Rand cotinine
    p <- ggplot(d2, aes(x = Mid_Cotinine_ng_ml, y = PSI,
        color = smoking_status)) +
        geom_point(size = 3) +
        geom_smooth(method = lm) +
        annotate("text", x = 50, y = y_rho, label = rho, size = 5) +
        annotate("text", x = 70, y = y_pval, label = pval, size = 5) +
        theme_classic() +
        labs(x = "Mid Cotinine ng/ml",
            y = "PSI Score", title = paste0(mod_num, " VCSIP")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(size = 13)) +
        theme(axis.text.y = element_text(size = 13)) +
        theme(axis.title = element_text(size = 13)) +
        guides(color = FALSE)
    ggsave(paste0(mod_num, "_PSIscore_by_midCot_vcsip_scatterplot.pdf"), p,
        width = width, height = height, dpi = 300)

    #-----------#
    # Late Cotinine

    rho <- paste0("r = ", round(cor.test(d2$PSI, d2$Late_Cotinine_ng_ml,
        method = "pearson")$estimate, 4))
    pval <- paste0("p = ", round(cor.test(d2$PSI, d2$Late_Cotinine_ng_ml,
        method = "pearson")$p.value, 7))

    # plot PSI vs Rand cotinine
    p <- ggplot(d2, aes(x = Late_Cotinine_ng_ml, y = PSI,
        color = smoking_status)) +
        geom_point(size = 3) +
        geom_smooth(method = lm) +
        annotate("text", x = 50, y = y_rho, label = rho, size = 5) +
        annotate("text", x = 70, y = y_pval, label = pval, size = 5) +
        theme_classic() +
        labs(x = "Late Cotinine ng/ml",
            y = "PSI Score", title = paste0(mod_num, " VCSIP")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(size = 13)) +
        theme(axis.text.y = element_text(size = 13)) +
        theme(axis.title = element_text(size = 13)) +
        guides(color = FALSE)
    ggsave(paste0(mod_num, "_PSIscore_by_lateCot_vcsip_scatterplot.pdf"), p,
        width = width, height = height, dpi = 300)
}