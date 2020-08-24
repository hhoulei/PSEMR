#' Causal Mediation Analysis with Multiple Causally Ordered Mediators
#' @name Path_effect
#' @param coef_m -- pathway martrix
#' @param effect_m -- estimated effect in Multi_IVW or Multi_MR_Egger
#' @param se_m -- estimated se in Multi_IVW or Multi_MR_Egger
#' @param straps -- repeated times for MCMC
#' @rdname Path_effect
#' @export

Path_effect <- function(coef_m, effect_m, se_m, straps) {
    
    bb <- c(coef_m * se_m)
    effect <- c(coef_m * effect_m)[bb != 0]
    se <- bb[bb != 0]
    
    indirect_effect <- prod(effect)
    
    indirect_boot = NULL
    for (k in 1:straps) {
        once <- NULL
        for (v in 1:length(effect)) {
            once <- c(once, rnorm(1, effect[v], se[v]))
        }
        indirect_boot[k] = prod(once)
    }
    
    se = sd(indirect_boot)
    ci_lower = sort(indirect_boot)[0.025 * straps + 1]
    ci_upper = sort(indirect_boot)[0.975 * straps]
    
    result <- c(indirect_effect, se, ci_lower, ci_upper)
    names(result) <- c("indirect_effect", "se", "ci_lower", "ci_upper")
    return(result)
    
}
