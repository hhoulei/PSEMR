#' Causal Mediation Analysis with Multiple Causally Ordered Mediators
#' @name Ordered_Egger
#' @param betaYG -- a vector, genetic summary statistics for outcome (coefficient)
#' @param sebetaYG -- a vector, genetic summary statistics for outcome (standard error)
#' @param betaXG -- a vector, genetic summary statistics for exposure (coefficient)
#' @param sebetaXG -- a vector,genetic summary statistics for exposure (standard error)
#' @param betaMG -- matrix or dataframe with n mediators, genetic summary statistics for mediators (coefficient)
#' @param sebetaMG -- matrix or dataframe with n mediators, genetic summary statistics formediators (standard error)
#' @return list
#' @rdname Ordered_Egger
#' @export

Ordered_Egger <- function(betaYG, sebetaYG, betaXG, sebetaXG, betaMG, sebetaMG) {
    
    betaMG <- as.matrix(betaMG)
    sebetaMG <- as.matrix(sebetaMG)
    
    # total effect
    fit1 <- lm(betaYG ~ betaXG, weights = sebetaYG^(-2))
    total_effect <- fit1$coef[2]
    se_total_effect <- summary(fit1)$coef[2, 2]/min(summary(fit1)$sigma, 1)
    ci_upper_total <- total_effect + qt(0.975, df = length(betaXG) - 2) * se_total_effect
    ci_lower_total <- total_effect - qt(0.975, df = length(betaXG) - 2) * se_total_effect
    
    data1 <- data.frame(betaYG, betaMG, betaXG)
    se_data1 <- data.frame(sebetaYG, sebetaMG, sebetaXG)
    result <- list()
    h = 1
    data2 <- data1
    p_test <- NULL
    while (ncol(data1) > 2) {
        
        data2 <- data1[, -1]
        reg1 <- lm(data1[, 1] ~ ., weights = se_data1[, 1]^-2, data = data2)
        reg1_result <- data.frame(coef = reg1$coef, se = summary(reg1)$coef[, 2]/min(summary(reg1)$sigma, 
            1), pvalue = summary(reg1)$coef[, 4])
        reg1_result$low <- reg1_result[, 1] - qt(0.975, df = length(betaXG) - ncol(data1)) * reg1_result[, 
            2]
        reg1_result$up <- reg1_result[, 1] + qt(0.975, df = length(betaXG) - ncol(data1)) * reg1_result[, 
            2]
        
        p_test <- rbind(p_test, reg1_result[1, ])
        
        result[[h]] <- as.matrix(reg1_result[-1, ])
        h = h + 1
        data1 <- data1[, -1]
        se_data1 <- se_data1[, -1]
    }
    reg2 <- lm(data1[, 1] ~ data1[, 2], weights = se_data1[, 1]^-2)
    reg2_result <- data.frame(coef = reg2$coef, se = summary(reg2)$coef[, 2]/min(summary(reg2)$sigma, 1), 
        pvalue = summary(reg2)$coef[, 4])
    reg2_result$low <- reg2_result[, 1] - qt(0.975, df = length(betaXG) - ncol(data1)) * reg2_result[, 2]
    reg2_result$up <- reg2_result[, 1] + qt(0.975, df = length(betaXG) - ncol(data1)) * reg2_result[, 2]
    p_test <- rbind(p_test, reg2_result[1, ])
    result[[h]] <- as.matrix(reg2_result[-1, ])
    
    zero <- matrix(0, nrow = (ncol(betaMG) + 2), ncol = (ncol(betaMG) + 2))
    nna <- c("X", paste0(rep("M", ncol(betaMG)), 1:ncol(betaMG)), "Y")
    colnames(zero) <- nna
    rownames(zero) <- nna
    M_result <- list(M_effect = zero, M_se = zero, M_pvalue = zero, M_low = zero, M_up = zero)
    
    for (k in 1:ncol(result[[1]])) {
        M_result[[k]][1, ncol(M_result[[k]])] <- result[[1]][nrow(result[[1]]), k]
        M_result[[k]][2:(ncol(M_result[[k]]) - 1), ncol(M_result[[k]])] <- result[[1]][1:(nrow(result[[1]]) - 
            1), k]
    }
    
    if (length(result) == 2) {
        
        M_result$M_effect[1, 2] <- result[[2]][1, 1]
        M_result$M_se[1, 2] <- result[[2]][1, 2]
        M_result$M_pvalue[1, 2] <- result[[2]][1, 3]
        M_result$M_low[1, 2] <- result[[2]][1, 4]
        M_result$M_up[1, 2] <- result[[2]][1, 5]
        
    } else {
        for (m in 2:(length(result) - 1)) {
            for (n in 1:ncol(result[[m]])) {
                M_result[[n]][1, m] <- result[[m]][nrow(result[[m]]), n]
                M_result[[n]][m, (m + 1):(ncol(M_result[[n]]) - 1)] <- result[[m]][1:(nrow(result[[m]]) - 
                  1), n]
            }
        }
        
        for (l in 1:ncol(result[[length(result)]])) {
            M_result[[l]][1, (ncol(M_result[[l]]) - 1)] <- result[[length(result)]][1, l]
        }
    }
    
    to <- c(total_effect, se_total_effect, ci_lower_total, ci_upper_total)
    names(to) <- c("total_effect", "se_total_effect", "ci_lower_total", "ci_upper_total")
    M_result[["Total"]] <- to
    
    
    rownames(p_test) <- c("G-Y", paste0("G-M", 1:ncol(betaMG)))
    M_result[["Pleiotry"]] <- p_test
    
    return(M_result)
}
