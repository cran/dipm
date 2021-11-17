#' @title Capturing Effects for Survival Outcomes 
#' with Two Treatment Groups
#' @description This function calculates the effect of a splitting
#' node in a Cox proportional hazards model with two treatment groups
#'
#' @param X A list containing: 
#'          the vector of outcomes Y, 
#'          the vector of censoring indicators C, 
#'          the vector of treatment assignments treatment,
#'          and the vector of split indicators split
#' @return A squared z-statistic
#'         of the split by treatment interaction term
#'         in a Cox proportional hazards model
#'         
#' @references Chen, V., Li, C., and Zhang, H. (2021). The dipm R 
#'             package: implementing the depth importance in 
#'             precision medicine (DIPM) tree and forest based method.
#'             \emph{Manuscript}.
#' 
#'             Chen, V. and Zhang, H. (2020). Depth importance in 
#'             precision medicine (DIPM): a tree and forest based method. 
#'             In \emph{Contemporary Experimental Design, 
#'             Multivariate Analysis and Data Mining}, 243-259.
#' 
#'             Chen, V. and Zhang, H. (2020). Depth importance in 
#'             precision medicine (DIPM): A tree-and forest-based 
#'             method for right-censored survival outcomes. 
#'             \emph{Biostatistics}.
#'             
#' @importFrom survival coxph
#' @importFrom utils capture.output
#' @noRd

coxph_R_to_C = function(X){

    old = options()         
    on.exit(options(old))  
    
    Y = X[, 1]
    C = X[, 2]
    treatment = X[, 3]
    split = X[, 4]

#    if a treatment group only has censored observations, return 0
    if(length(table(C[which(treatment==0)])) == 1 ||
         length(table(C[which(treatment==1)])) == 1){
        return(0)
    }

#    save "coxph" warnings in x
    options(warn = 1) 

    x = capture.output({
        fit0 = coxph(Surv(Y, C) ~ treatment + split + treatment * split)
    }, type = "message")

#    if there is a warning, return 0
    if(length(x) > 0){  # if there is at least 1 warning
        return(0)
    }

#    get z-statistic of split by treatment interaction term
    fit = summary(fit0)
    val = fit$"coefficients"["treatment:split", "z"]

    return(val^2)
}