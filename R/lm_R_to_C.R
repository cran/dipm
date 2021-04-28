#' @title Capturing Effects for Continuous Outcomes 
#' with Multiple Treatment Groups
#' @description This function calculates the effect of a splitting
#' node in a linear regression model with multiple treatment
#' groups
#'
#' @param X A list containing: 
#'          the vector of outcomes Y, 
#'          the vector of treatment assignments treatment,
#'          and the vector of split indicators split;
#'          Note that here, the treatment variable is 
#'          not binary, i.e., there may be more than 
#'          2 treatment groups
#' @return A squared z-statistic 
#'         with the largest magnitude among all of the 
#'         split by treatment interaction terms in a 
#'         multiple linear regression model
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
#' @importFrom stats lm
#' @importFrom utils capture.output
#' @noRd

lm_R_to_C <- function(X) {
    
    old <- options()         
    on.exit(options(old))  

    Y=X[,1]
    treatment=as.factor(X[,2])
    split=X[,3]

#    save warnings in x
    options(warn=1) 

    x=capture.output({
        fit0=lm(Y~treatment+split+treatment*split)
    },type="message")

#    if there is a warning, return 0
    if ( length(x) > 0 ) {  # if there is at least 1 warning
        return(0)
    }

#    get the z-statistic with the largest magnitude of all 
#    split by treatment interaction terms
    fit=summary(fit0)
    i_int=grep(":split",rownames(fit$"coefficients"),value=FALSE)
    val=max(abs(fit$"coefficients"[i_int,"t value"]))

    return(val^2)
}
