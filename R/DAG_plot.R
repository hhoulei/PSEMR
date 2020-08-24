
#' Causal Mediation Analysis with Multiple Causally Ordered Mediators
#' @name DAG_plot
#' @param va -- variable's name
#' @param pse_result -- results of PSE-MR
#' @param ylow -- lowest value of ylab
#' @param yup -- uppest value of ylab
#' @param tt -- label name
#' @importFrom dagitty dagitty
#' @importFrom ggdag ggdag
#' @importFrom ggplot2 ggplot
#' @return list
#' @rdname DAG_plot
#' @export

library(ggdag)
library(dagitty)
# library(lavaan)
library(png)

DAG_plot <- function(va, pse_result, ylow, yup, tt) {
    
    # ylow=0.6,yup=3.2
    
    nM <- length(va) - 2
    px <- paste0(va[1], " [exposure,pos=", "\"", 0, ",", (1 + nM)/2, "\"", "] ")
    py <- paste0(va[2], " [outcome,pos=", "\"", 2, ",", (1 + nM)/2, "\"", "] ")
    
    vamp <- NULL
    oo1 <- which(pse_result$M_pvalue[1, ] < 0.05 & pse_result$M_pvalue[1, ] > 0) - 1
    for (i1 in oo1) {
        vamp <- c(vamp, paste0(" ", va[1], " -> ", va[2 + i1], " "))
    }
    oo2 <- which(pse_result$M_pvalue[, nM + 2] < 0.05 & pse_result$M_pvalue[, nM + 2] > 0) - 1
    for (i2 in oo2) {
        vamp <- c(vamp, paste0(va[2 + i2], " -> ", va[2], " "))
    }
    if (pse_result$M_pvalue[1, nM + 2] < 0.05 & pse_result$M_pvalue[1, nM + 2] > 0) 
        vamp <- c(vamp, paste0(va[1], " -> ", va[2], " "))
    
    vamp_arc <- NULL
    name_vamp_arc <- NULL
    ty <- combn(nM, 2)
    for (j in 1:ncol(ty)) {
        ord <- diff(ty[, j])
        if (pse_result$M_pvalue[ty[1, j] + 1, ty[2, j] + 1] < 0.05 & pse_result$M_pvalue[ty[1, j] + 1, ty[2, 
            j] + 1] > 0) {
            io <- va[2 + ty[, j]]
            if (ord == 1) {
                vamp <- c(vamp, paste0(io[1], " <- ", io[2], " "))
            } else {
                vamp_arc <- c(vamp_arc, paste0(io[1], " <- ", io[2], " "))
                name_vamp_arc <- rbind(name_vamp_arc, c(io[1], io[2]))
            }
        }
    }
    
    varb_ALLp <- ""
    varb_ALLp1 <- ""
    varb_ALLp <- paste0(varb_ALLp, px, py)
    for (b in 1:nM) {
        varb_ALLp <- paste0(varb_ALLp, va[2 + b], " [pos=", "\"", 1, ",", b, "\"", "] ")
        if (va[2 + b] %in% name_vamp_arc) {
            varb_ALLp1 <- paste0(varb_ALLp1, va[2 + b], " [pos=", "\"", 1, ",", b, "\"", "] ")
        }
    }
    
    for (k in 1:length(vamp)) {
        varb_ALLp <- paste0(varb_ALLp, vamp[k])
    }
    
    for (c in 1:length(vamp_arc)) {
        varb_ALLp1 <- paste0(varb_ALLp1, vamp_arc[c])
    }
    
    
    # g1 <- dagitty(paste0('dag { ',paste0(varb_ALLp,varb_ALLp1),'}'))
    g01 <- dagitty(paste0("dag { ", varb_ALLp, "}"))
    
    # tidy_dagitty(g1) rr <- tidy_dagitty(g01)$data[which(is.na(tidy_dagitty(g01)$data[,8])),]
    # tidy_dagitty(g1)$data[which(unlist(tidy_dagitty(g1)$data[,1]) %in% rr[,1]),8] <- TRUE
    
    g2 <- ggdag(g01, edge_type = "link_arc", node_size = 18, text_size = 6)
    g3 <- g2 + xlab(NULL) + ylab(NULL) + theme_bw() + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        panel.border = element_blank()) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        panel.border = element_blank()) + scale_y_continuous(limits = c(ylow, yup), breaks = NULL) + scale_x_continuous(limits = c(-0.2, 
        2.2), breaks = NULL) + geom_dag_edges(edge_width = 1.3, arrow_directed = grid::arrow(length = grid::unit(8, 
        "pt"), type = "closed")) + annotate("text", x = -0.15, y = yup, label = tt, size = 10)
    
    if (is.null(name_vamp_arc)) {
        g31 <- NULL
    } else {
        g02 <- dagitty(paste0("dag { ", varb_ALLp1, "}"))
        g21 <- ggdag(g02, edge_type = "arc", node_size = 18, text_size = 6)
        g31 <- g21 + xlab(NULL) + ylab(NULL) + theme_bw() + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
            panel.border = element_blank()) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
            panel.border = element_blank()) + scale_y_continuous(limits = c(ylow, yup), breaks = NULL) + scale_x_continuous(limits = c(-0.2, 
            2.2), breaks = NULL) + geom_dag_edges_arc(edge_width = 1.3, arrow = grid::arrow(length = grid::unit(8, 
            "pt"), type = "closed")) + theme(panel.background = element_rect(fill = "transparent", colour = NA), 
            panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.background = element_rect(fill = "transparent", 
                colour = NA)) + annotate("text", x = -0.15, y = yup, label = tt, size = 10)
    }
    
    return(list(g3, g31))
}

# va <- c('BMI','CVD','HDL','TG') DAG_plot(va,pse_result=ivw_result) DAG_plot(va,pse_result=egger_result)
# g1 <- dagitty(paste0('dag { ',paste0(varb_ALL,varb_ALLp),'}')) g1 <- dagitty(paste0('dag {
# ',varb_ALLp,'}')) plot(g1)







