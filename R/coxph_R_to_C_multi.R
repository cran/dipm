#' @title Capturing Effects for Survival Outcomes 
#' with Multiple Treatment Groups
#' @description This function calculates the effect of a splitting
#' node in a Cox proportional hazards model with multiple treatment
#' groups
#'
#' @param X A list containing: 
#'          the vector of outcomes Y, 
#'          the vector of censoring indicators C, 
#'          the vector of treatment assignments treatment,
#'          and the vector of split indicators split; 
#'          Note that here, the treatment variable is 
#'          not binary, i.e., there may be more than 2 
#'          treatment groups
#' @return A squared z-statistic
#'         with the largest magnitude among all of the 
#'         split by treatment interaction terms in a 
#'         Cox proportional hazards model
#'         
#' @references Chen, V., Li, C., and Zhang, H. (2022). dipm: an 
#'             R package implementing the Depth Importance in 
#'             Precision Medicine (DIPM) tree and Forest-based method.
#'             \emph{Bioinformatics Advances}, \strong{2}(1), vbac041.
#' 
#'             Chen, V. and Zhang, H. (2020). Depth importance in 
#'             precision medicine (DIPM): a tree and forest based method. 
#'             In \emph{Contemporary Experimental Design, 
#'             Multivariate Analysis and Data Mining}, 243-259.
#' 
#'             Chen, V. and Zhang, H. (2022). Depth importance in 
#'             precision medicine (DIPM): A tree-and forest-based 
#'             method for right-censored survival outcomes. 
#'             \emph{Biostatistics} \strong{23}(1), 157-172.
#'             
#' @importFrom survival coxph
#' @importFrom utils capture.output
#' @noRd

coxph_R_to_C_multi = function(X){
    
    old = options()         
    on.exit(options(old))  

    Y = X[, 1]
    C = X[, 2]
    treatment = as.factor(X[, 3])
    split = X[, 4]

#    if a treatment group only has censored observations, return 0
    trts = unique(X[, 3])
    ntrt = length(trts)

    for(i in 1:ntrt){
        if(length(table(C[which(treatment==trts[i])])) == 1){
            return(0)
        }
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

#    get the z-statistic with the largest magnitude of all 
#    split by treatment interaction terms
    fit = summary(fit0)
    i_int = grep(":split", rownames(fit$"coefficients"), value = FALSE)
    val = max(abs(fit$"coefficients"[i_int, "z"]))

    return(val^2)
}
