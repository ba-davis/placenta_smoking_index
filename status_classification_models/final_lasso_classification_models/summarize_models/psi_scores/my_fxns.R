
# normalize betas function
#   beta_sub: samples as rows, cpgs as columns
#             samples are rownames(beta_sub), no sample column
#   mytab: model reference table with cols <cpg> <estimate> <mean> <sd> <y_int>
norm_betas <- function(beta_sub, mytab) {
    # centering (subtract the mean from each cpg) and
    # scaling (divide by standard deviation)

    # for each column in the input beta matrix
    #   - subtract the mean for that cpg in mytab from
    #     each row in the cpg column in user input beta matrix.
    #   - Then, divide this new value by the sd value for the cpg
    for (i in seq_len(ncol(beta_sub))) {
        # subtract the mean
        beta_sub[[colnames(beta_sub)[i]]] <- beta_sub[[colnames(beta_sub)[i]]] -
        mytab[mytab$cpg == colnames(beta_sub)[i], "mean"]
        # divide by the standard deviation
        beta_sub[[colnames(beta_sub)[i]]] <- beta_sub[[colnames(beta_sub)[i]]] /
        mytab[mytab$cpg == colnames(beta_sub)[i], "sd"]
    }

    return(beta_sub)
}

# multiply beta value by coefficient and store in matrix
#   beta_sub: samples as rows, cpgs as columns
#             samples are rownames(beta_sub), no sample column
#   mytab: model reference table with cols <cpg> <estimate> <mean> <sd> <y_int>
calculate_psi <- function(beta_sub, mytab) {
    # put the cpg rows of transposed beta_sub in the same order as
    # the cpg rows in ref table mytab$cpg
    B <- t(beta_sub)
    B <- B[match(mytab$cpg, rownames(B)), ]
    # sanity check
    if (!(identical(rownames(B), mytab$cpg))) {
    print(paste0("WARNING: cpgs are not in proper order.",
    " Problem with PSI calculation."))
    }

    y <- NULL
    i <- NULL
    A <- mytab$estimate
    y <- matrix(nrow = nrow(B), ncol = ncol(B))
    # for each row (cpg) of B, multiply the entire row by
    #   the coefficient for that cpg
    for (i in 1:nrow(B)){
        y[i, ] <- rbind(as.numeric(B[i, ]) * A[i])
    }

    s <- colSums(y)
    # store normalized psi score per sample in dataframe
    d2 <- data.frame(colnames(B), s)

    return(d2)
}
