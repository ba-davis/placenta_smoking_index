


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

# convert logit probability value to probability value
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# take input beta matrix and predict to get output summary table
# input_beta: samples as rows, cpgs as columns
#             samples are rownames(beta_sub), no sample column
# mytab: model reference table with cols <cpg> <estimate> <mean> <sd> <y_int>
#
# returns output df of results with columns:
#   <sample> <psi> <norm_psi> <logit_prob>
#     <probability_smoker> <probability_nonsmoker>
predict_on_new_data <- function(input_beta, mytab, prob_cutoff = 0.5) {
    # check how many of the required cpgs for the model are present as columns
    # in user input beta matrix
    print("Checking input CpGs for those required by the model.")
    required_cpg_num <- length(mytab$cpg)
    present_cpg_num <- length(mytab$cpg[mytab$cpg %in% colnames(input_beta)])
    print(paste0(present_cpg_num, " out of ", required_cpg_num,
        " required CpGs are present."))
    if (present_cpg_num != required_cpg_num) {
        print(paste0(required_cpg_num - present_cpg_num, " CpGs are missing."))
    }

    # subset input beta matrix to contain only present required cpgs
    beta_sub <- input_beta[, colnames(input_beta) %in% mytab$cpg]

    # Calculate PSI from unnormalized beta matrix
    d2 <- calculate_psi(beta_sub, mytab)
    colnames(d2) <- c("sample", "psi")

    # normalize the beta matrix according to training data
    beta_norm <- norm_betas(beta_sub, mytab)

    # Calculate PSI score from the normalized betas
    d2$norm_psi <- calculate_psi(beta_norm, mytab)$s

    # add y-intercept to the normalized psi
    y_int <- unique(mytab$y_int)
    d2$logit_prob <- d2$norm_psi + y_int

    # add predicted probability columns
    d2$probability_smoker <- sapply(d2$logit_prob, logit2prob)
    d2$probability_nonsmoker <- 1 - d2$probability_smoker

    # add predicted class column
    d2$predicted_class <- ifelse(d2$probability_smoker >= prob_cutoff,
        "smoker", "nonsmoker")

    return(d2)
}