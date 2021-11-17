
#' @title Depth Importance in Precision Medicine (DIPM)
#' @description This function creates a classification tree
#'              designed to identify subgroups in which subjects
#'              perform especially well or especially poorly in a
#'              given treatment group.
#'
#' @param formula A description of the model to be fit with format
#'                \code{Y ~ treatment | X1 + X2} for data with a
#'                continuous outcome variable Y and 
#'                \code{Surv(Y, C) ~ treatment | X1 + X2} for data with
#'                a right-censored survival outcome variable Y and
#'                a censoring indicator C
#' @param data A matrix or data frame of the data
#' @param types A vector, data frame, or matrix of the types
#'              of each variable in the data; if left blank, the
#'              default is to assume all of the candidate split
#'              variables are ordinal; otherwise, all variables in 
#'              the data must be specified, and the possible variable 
#'              types are: "response", "treatment", "C", "binary", 
#'              "ordinal", and "nominal" for outcome variable Y, the 
#'              treatment variable, the censoring indicator (if 
#'              applicable), binary candidate split variables, ordinal
#'              candidate split variables, and nominal candidate split
#'              variables respectively
#' @param nmin An integer specifying the minimum node size of
#'             the overall classification tree
#' @param nmin2 An integer specifying the minimum node size of
#'              embedded trees
#' @param ntree An integer specifying the number of embedded trees
#'              to construct at each node of the overall 
#'              classification tree; if left blank, the default value
#'              of \code{ceiling(min(max(sqrt(n), sqrt(nc)), 1000))}
#'              will be used if \code{mtry = Inf} below and 
#'              \code{ceiling(min(max(n, nc), 1000))} otherwise; 
#'              \code{n} is the total sample size of the data, and
#'              \code{nc} is the total number of candidate split
#'              variables
#' @param mtry An integer specifying the number of candidate split
#'             variables to randomly select at each node of 
#'             embedded trees; if \code{mtry} is set equal to the
#'             default value of Inf, then all possible splits of all
#'             candidate split variables are considered at the nodes
#'             of the embedded trees; otherwise, a recommended value
#'             of \code{mtry} to use is the square root of the total
#'             number of candidate split variables rounded up to the
#'             nearest integer
#' @param maxdepth An integer specifying the maximum depth of the
#'                 overall classification tree; this argument is 
#'                 optional but useful for shortening computation 
#'                 time; if left blank, the default is to grow the 
#'                 full tree until the minimum node size \code{nmin} 
#'                 is reached
#' @param maxdepth2 An integer specifying the maximum depth of 
#'                  embedded trees; this argument is optional but 
#'                  useful for shortening computation time; if left 
#'                  blank, the default is to grow each full, 
#'                  embedded tree until the minimum node size
#'                  \code{nmin2} is reached
#' @param print A boolean (TRUE/FALSE) value, where TRUE prints
#'              a more readable version of the final tree to the
#'              screen
#' @param dataframe A boolean (TRUE/FALSE) value, where TRUE returns
#'                  the final tree as a dataframe
#' @param prune A boolean (TRUE/FALSE) value, where TRUE prunes
#'              the final tree using \code{pmprune} function
#'
#' @details This function creates a classification tree to identify
#'          subgroups relevant to the precision medicine setting.
#'          At each node of the classification tree, a random forest
#'          of so-called embedded trees are fit and used to
#'          calculate a depth variable importance score for each
#'          candidate split variable in the data. The candidate
#'          split variable with the largest variable importance score
#'          is identified as the best split variable of the node.
#'          Then, all possible splits of the selected split variable
#'          are considered, and the split with the greatest split
#'          criteria value is finally selected as the best split of
#'          the best variable.
#' 
#'          The depth variable importance score was originally
#'          proposed by Chen et al. (2007), and the score has been
#'          adapted to the precision medicine setting here.
#'          The depth variable importance score is a relatively
#'          simple measure that takes into account two components:
#'          the depth of a split in a tree and the strength of the
#'          split. The strength of the split is captured with a G
#'          test statistic that may be modified depending on the
#'          type of analysis at hand. When the outcome variable is
#'          continuous, G is the test statistic that tests the
#'          significance of the split by treatment interaction term
#'          in a linear regression model. When the outcome variable
#'          is a right-censored survival time, G is the test statistic
#'          that tests the significance of the split by interaction
#'          term in a Cox proportional hazards model.
#'
#'          When using \code{dipm}, note the following
#'          requirements for the supplied data. First, the dataset
#'          must contain an outcome variable Y and a treatment
#'          variable. If Y is a right-censored survival time
#'          outcome, then there must also be a censoring indicator
#'          C, where values of 1 denote the occurrence of the 
#'          (harmful) event of interest, and values of 0 denote
#'          censoring. If there are only two treatment groups, then
#'          the two possible values must be 0 or 1. If there are
#'          more than two treatment groups, then the possible values
#'          must be integers starting from 1 to the total number of
#'          treatment assignments. In regard to the candidate split
#'          variables, if a variable is binary, then the variable
#'          must take values of 0 or 1. If a variable is nominal,
#'          then the values must be integers starting from 1 to the
#'          total number of categories. There cannot be any missing
#'          values in the dataset. For candidate split variables
#'          with missing values, the missings together (MT) method
#'          proposed by Zhang et al. (1996) is helpful.
#'
#' @return \code{dipm} returns the final classification tree as a 
#'         \code{party} object by default or a data frame. See
#'         Hothorn and Zeileis (2015) for details. The data 
#'         frame contains the following columns of information:
#'         \item{node}{Unique integer values that identify each node
#'                     in the tree, where all of the nodes are
#'                     indexed starting from 1}
#'         \item{splitvar}{Integers that represent the candidate split
#'                         variable used to split each node, where
#'                         all of the variables are indexed starting
#'                         from 1; for terminal nodes, i.e., nodes
#'                         without child nodes, the value is set 
#'                         equal to NA}
#'         \item{splitvar_name}{The names of the candidate split 
#'                              variables used to split each node
#'                              obtained from the column names of the
#'                              supplied data; for terminal nodes,
#'                              the value is set equal to NA}
#'         \item{type}{Characters that denote the type of each 
#'                     candidate split variable; "bin" is for binary
#'                     variables, "ord" for ordinal, and "nom" for
#'                     nominal; for terminal nodes, the value is set
#'                     equal to NA}
#'         \item{splitval}{Values of the left child node of the 
#'                         current split/node; for binary variables,
#'                         a value of 0 is printed, and subjects with
#'                         values of 0 for the current \code{splitvar}
#'                         are in the left child node, while subjects
#'                         with values of 1 are in the right child
#'                         node; for ordinal variables,
#'                         \code{splitval} is numeric and implies
#'                         that subjects with values of the current
#'                         \code{splitvar} less than or equal to
#'                         \code{splitval} are in the left child 
#'                         node, while the remaining subjects with 
#'                         values greater than \code{splitval} are in 
#'                         the right child node; for nominal
#'                         variables, the \code{splitval} is a set of
#'                         integers separated by commas, and subjects
#'                         in that set of categories are in the left
#'                         child node, while the remaining subjects
#'                         are in the right child node; for terminal
#'                         nodes, the value is set equal to NA}
#'         \item{lchild}{Integers that represent the index (i.e.,
#'                       \code{node} value) of each node's left
#'                       child node; for terminal nodes, the value is
#'                       set equal to NA}
#'         \item{rchild}{Integers that represent the index (i.e.,
#'                       \code{node} value) of each node's right
#'                       child node; for terminal nodes, the value is
#'                       set equal to NA}
#'         \item{depth}{Integers that specify the depth of each
#'                      node; the root node has depth 1, its 
#'                      children have depth 2, etc.}
#'         \item{nsubj}{Integers that count the total number of
#'                      subjects within each node}
#'         \item{besttrt}{Integers that denote the identified best 
#'                        treatment assignment of each node}
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
#'             Chen, X., Liu, C.-T., Zhang, M., and Zhang, H. (2007).
#'             A forest-based approach to identifying gene and
#'             gene-gene interactions. \emph{Proceedings of the 
#'             National Academy of Sciences of the United States of 
#'             America} \strong{204}, 19199-19203.
#'
#'             Zhang, H., Holford, T., and Bracken, M.B. (1996).
#'             A tree-based method of analysis for prospective
#'             studies. \emph{Statistics in Medicine} \strong{15},
#'             37-49.
#'             
#'             Hothorn, T. and Zeileis, A. (2015). partykit: 
#'             a modular toolkit for recursive partytioning in R. 
#'             \emph{The Journal of Machine Learning Research} 
#'             \strong{16}(1), 3905-3909.
#'
#' @seealso \code{\link{spmtree}}
#'
#' @examples
#' 
#' #
#' # ... an example with a continuous outcome variable
#' #     and two treatment groups
#' #
#'
#' N = 100
#' set.seed(123)
#'
#' # generate binary treatments
#' treatment = rbinom(N, 1, 0.5)
#'
#' # generate candidate split variables
#' X1 = rnorm(n = N, mean = 0, sd = 1)
#' X2 = rnorm(n = N, mean = 0, sd = 1)
#' X3 = rnorm(n = N, mean = 0, sd = 1)
#' X4 = rnorm(n = N, mean = 0, sd = 1)
#' X5 = rnorm(n = N, mean = 0, sd = 1)
#' X = cbind(X1, X2, X3, X4, X5)
#' colnames(X) = paste0("X", 1:5)
#'
#' # generate continuous outcome variable
#' calculateLink = function(X, treatment){
#'
#'     ((X[, 1] <= 0) & (X[, 2] <= 0)) *
#'         (25 * (1 - treatment) + 8 * treatment) + 
#'
#'     ((X[, 1] <= 0) & (X[, 2] > 0)) *
#'         (18 * (1 - treatment) + 20 * treatment) +
#'
#'     ((X[, 1] > 0) & (X[, 3] <= 0)) *
#'         (20 * (1 - treatment) + 18 * treatment) + 
#'
#'     ((X[, 1] > 0) & (X[, 3] > 0)) *
#'         (8 * (1 - treatment) + 25 * treatment)
#' }
#'
#' Link = calculateLink(X, treatment)
#' Y = rnorm(N, mean = Link, sd = 1)
#'
#' # combine variables in a data frame
#' data = data.frame(X, Y, treatment)
#' 
#' # fit a dipm classification tree 
#' tree1 = dipm(Y ~ treatment | ., data, mtry = 1, maxdepth = 3)
#' # predict optimal treatment for new subjects
#' predict(tree1, newdata = head(data), 
#' FUN = function(n)  as.numeric(n$info$opt_trt))
#'
#'\donttest{
#' #
#' # ... an example with a continuous outcome variable
#' #     and three treatment groups
#' #
#' 
#' N = 600
#' set.seed(123)
#' 
#' # generate treatments
#' treatment = sample(1:3, N, replace = TRUE)
#' 
#' # generate candidate split variables
#' X1 = round(rnorm(n = N, mean = 0, sd = 1), 4)
#' X2 = round(rnorm(n = N, mean = 0, sd = 1), 4)
#' X3 = sample(1:4, N, replace = TRUE)
#' X4 = sample(1:5, N, replace = TRUE)
#' X5 = rbinom(N, 1, 0.5)
#' X6 = rbinom(N, 1, 0.5)
#' X7 = rbinom(N, 1, 0.5)
#' X = cbind(X1, X2, X3, X4, X5, X6, X7)
#' colnames(X) = paste0("X", 1:7)
#' 
#' # generate continuous outcome variable
#' calculateLink = function(X,treatment){
#' 
#'     10.2 - 0.3 * (treatment == 1) - 0.1 * X[, 1] + 
#'     2.1 * (treatment == 1) * X[, 1] +
#'     1.2 * X[, 2]
#' }
#' 
#' Link = calculateLink(X, treatment)
#' Y = rnorm(N, mean = Link, sd = 1)
#' 
#' # combine variables in a data frame
#' data = data.frame(X, Y, treatment)
#' 
#' # create vector of variable types
#' types = c(rep("ordinal", 2), rep("nominal", 2), rep("binary", 3),
#'         "response", "treatment")
#' 
#' # fit a dipm classification tree 
#' tree2 = dipm(Y ~ treatment | ., data, types = types, maxdepth = 2)
#'
#' #
#' # ... an example with a survival outcome variable
#' #     and two treatment groups
#' #
#'
#' N = 300
#' set.seed(321)
#'
#' # generate binary treatments
#' treatment = rbinom(N, 1, 0.5)
#'
#' # generate candidate split variables
#' X1 = rnorm(n = N, mean = 0, sd = 1)
#' X2 = rnorm(n = N, mean = 0, sd = 1)
#' X3 = rnorm(n = N, mean = 0, sd = 1)
#' X4 = rnorm(n = N, mean = 0, sd = 1)
#' X5 = rnorm(n = N, mean = 0, sd = 1)
#' X = cbind(X1, X2, X3, X4, X5)
#' colnames(X) = paste0("X", 1:5)
#'
#' # generate survival outcome variable
#' calculateLink = function(X, treatment){
#'
#'     X[, 1] + 0.5 * X[, 3] + (3 * treatment - 1.5) * (abs(X[, 5]) - 0.67)
#' }
#'
#' Link = calculateLink(X, treatment)
#' T = rexp(N, exp(-Link))
#' C0 = rexp(N, 0.1 * exp(X[, 5] + X[, 2]))
#' Y = pmin(T, C0)
#' C = (T <= C0)
#'
#' # combine variables in a data frame
#' data = data.frame(X, Y, C, treatment)
#' 
#' # fit a dipm classification tree 
#' tree3 = dipm(Surv(Y, C) ~ treatment | .,data, ntree = 1, maxdepth = 2,
#'            maxdepth2 = 6)
#'
#' #
#' # ... an example with a survival outcome variable
#' #     and four treatment groups
#' #
#' 
#' N = 800
#' set.seed(321)
#' 
#' # generate treatments
#' treatment = sample(1:4, N, replace = TRUE)
#' 
#' # generate candidate split variables
#' X1 = round(rnorm(n = N, mean = 0, sd = 1), 4)
#' X2 = round(rnorm(n = N, mean = 0, sd = 1), 4)
#' X3 = sample(1:4, N, replace = TRUE)
#' X4 = sample(1:5, N, replace = TRUE)
#' X5 = rbinom(N, 1, 0.5)
#' X6 = rbinom(N, 1, 0.5)
#' X7 = rbinom(N, 1, 0.5)
#' X = cbind(X1, X2, X3, X4, X5, X6, X7)
#' colnames(X) = paste0("X", 1:7)
#' 
#' # generate survival outcome variable
#' calculateLink = function(X, treatment, noise){
#' 
#'     -0.2 * (treatment == 1) +
#'     -1.1 * X[, 1] + 
#'     1.2 * (treatment == 1) * X[, 1] +
#'     1.2 * X[, 2]
#' }
#' 
#' Link = calculateLink(X, treatment)
#' T = rweibull(N, shape=2, scale = exp(Link))
#' Cnoise = runif(n = N) + runif(n = N)
#' C0 = rexp(N, exp(0.3 * -Cnoise))
#' Y = pmin(T, C0)
#' C = (T <= C0)
#' 
#' # combine variables in a data frame
#' data = data.frame(X, Y, C, treatment)
#' 
#' # create vector of variable types
#' types = c(rep("ordinal", 2), rep("nominal", 2), rep("binary", 3), 
#'         "response", "C", "treatment")
#' 
#' # fit two dipm classification trees 
#' tree4 = dipm(Surv(Y, C) ~ treatment | ., data, types = types, ntree = 1, 
#'            maxdepth = 2, maxdepth2 = 6)
#' tree5 = dipm(Surv(Y, C) ~ treatment | X3 + X4, data, types = types, ntree = 1, 
#'            maxdepth = 2, maxdepth2 = 6)
#' }
#' 
#' @export
#' @import partykit
#' @import survival
#' @import stats

dipm = function(formula,
                 data,
                 types = NULL,
                 nmin = 5,
                 nmin2 = 5,
                 ntree = NULL,
                 mtry = Inf,
                 maxdepth = Inf,
                 maxdepth2 = Inf,
                 print = TRUE,
                 dataframe = FALSE,
                 prune = FALSE){

#    check inputs
    if(missing(formula)){
        stop("The formula input is missing.")
    }

    if(missing(data)){
        stop("The data input is missing.")
    }

#    coerce data input to R "data.frame" object
    data = as.data.frame(data)

#    coerce formula input to R "formula" object
    form = as.formula(formula)

#    if not missing, coerce types input to R "data.frame" object
    if(missing(types) == FALSE){

        types = as.data.frame(types)

        if(nrow(types) != 1){

            types = t(types)
            types = as.data.frame(types)
        }

        if(ncol(types) != ncol(data)){
            stop("The number of variables in types does not equal the number of variables in the data.")
        }

        colnames(types) = colnames(data)
    }

#    get names of variables in the formula
    form_vars = all.vars(form)      # all variables
    form_lhs = all.vars(form[[2]])  # variables to the left of ~
    form_rhs = all.vars(form[[3]])  # variables to the right of ~

#    response variable should always be (first) in lhs
    Y = data[, form_lhs[1]]

#    get censoring variable if applicable
    if(length(form_lhs) == 1){

        C = rep(0, nrow(data))
        surv = 0
    }

    if(length(form_lhs) == 2){

        C = data[, form_lhs[2]]
        surv = 1
    }

#    treatment variable should always be first in rhs
    treatment = data[, form_rhs[1]]

#    determine appropriate method from data and value of "mtry"
    ntrts = nlevels(as.factor(treatment))
    
    if(mtry == Inf){
        mtry = 0
    }
    
    if(maxdepth == Inf){
        maxdepth = -7
    }
    
    if(maxdepth2 == Inf){
        maxdepth2 = -7
    }

    if(ntrts <= 1){
        stop("At least 2 treatment groups are required.")
    }

    if(ntrts == 2){

        if(surv == 0){

            if(mtry == 0) method = 6
            else method = 7

        }else if(surv == 1){

            if(mtry == 0) method = 12
            else method = 13
        }

    }else if(ntrts > 2){

        if(surv == 0){

            if(mtry == 0) method = 22
            else method = 23

        }else if(surv == 1){

            if( mtry == 0 ) method=20
            else method=21
        }
    }

#    get matrix of candidate split variables X
    if( form_rhs[2] == "." ){  # account for Y ~ treatment | . 
                                 # formula

        exclude = c(which(colnames(data) == form_lhs[1]), # Y variable
                  which(colnames(data) == form_rhs[1])) # treatment

        if(length(form_lhs) == 2){ # exclude censoring indicator
            exclude = c(exclude,which(colnames(data) == form_lhs[2]))
        }

        X = data.frame(data[, -exclude])
        types = types[, -exclude]

    }else{
        include = which(colnames(data) %in% form_rhs[-1])
        X = data.frame(data[, include])
        types = types[, include]
    }

#    calculate number of observations n and variables nc
    n = nrow(X)
    nc = ncol(X)
    if(nc == 1){
        names(X) = names(data)[-exclude]
    }

#    use recommended value of total number of embedded trees
#    if ntree argument is blank
    if( missing(ntree)){
        if(method %in% c(6, 12, 20, 22)){
            ntree = ceiling(min(max(sqrt(n), sqrt(nc)), 1000)) 
        }else if(method %in% c(7, 13, 21, 23)){
            ntree = ceiling(min(max(n, nc), 1000)) 
        }
    }

#    prepare types
    if(is.null(types)){  

        types = rep(2, nc) # default is to assume all candidate
                        # split variables are ordinal

        message("Note that all candidate split variables are assumed to be ordinal.")

    }else{
        lll = ncol(types)
        for(i in 1:lll){
            if(types[i] == "binary") types[i] = 1
            if(types[i] == "ordinal") types[i] = 2
            if(types[i] == "nominal") types[i] = 3
        }
    }

#    create array of number of categories for nominal variables
    ifnominal = any(types == 3)
    if(ifnominal == TRUE){

        inom = which(types == 3)
        for(i in 1:length(inom)){
            X[, inom[i]] = factor(X[, inom[i]])
            data[, colnames(X)[inom[i]]] = X[, inom[i]]
        }

        ncat = sapply(X, function(x) 
                    if(is.null(levels(x))) -7 
                    else max(as.numeric(levels(x)[x])))

    } else {

        ncat = rep(-7, nc)
    }

#    prepare covariate data
    XC = t(X)

#    set types of R arguments to C
    storage.mode(ntree) = "integer"
    storage.mode(n) = "integer"
    storage.mode(nc) = "integer"
    storage.mode(Y) = "double"
    storage.mode(XC) = "double"
    storage.mode(types) = "integer"
    storage.mode(ncat) = "integer"
    storage.mode(treatment) = "integer"
    storage.mode(C) = "integer"
    storage.mode(nmin) = "integer"
    storage.mode(nmin2) = "integer"
    storage.mode(mtry) = "integer"
    storage.mode(maxdepth) = "integer"
    storage.mode(maxdepth2) = "integer"
    storage.mode(method) = "integer"

    tree = .Call("maketree",
               ntree = ntree,
               n = n,
               nc = nc,
               Y = Y,
               X = XC,
               types = types,
               ncat = ncat,
               treat = treatment,
               censor = C,
               nmin = nmin,
               nmin2 = nmin2,
               mtry = mtry,
               maxdepth = maxdepth,
               maxdepth2 = maxdepth2,
               method = method,
               environment(lm_R_to_C))

    rm(XC)

#    reformat tree
    tree_txt = data.frame(as.vector(tree[[1]]),
                        as.vector(tree[[2]]),
                        as.vector(tree[[3]]),
                        as.vector(tree[[4]]),
                        as.vector(tree[[5]]),
                        as.vector(tree[[7]]),
                        as.vector(tree[[8]]),
                        as.vector(tree[[9]]),
                        as.vector(tree[[6]]),
                        as.vector(tree[[10]]),
                        as.vector(tree[[11]]),
                        as.vector(tree[[12]]),
                        as.vector(tree[[13]]),
                        as.vector(tree[[14]]),
                        as.vector(tree[[15]]),
                        as.vector(tree[[16]]),
                        as.vector(tree[[17]]))

    colnames(tree_txt) = c("node",
                         "splitvar",
                         "type",
                         "sign",
                         "splitval",
                         "parent",
                         "lchild",
                         "rchild",
                         "depth",
                         "nsubj",
                         "ntrt0",
                         "ntrt1",
                         "r0",
                         "r1",
                         "p0",
                         "p1",
                         "besttrt")

#    process tree output and/or print tree to screen
    if(form_rhs[2] != "."){
        splitvar_include = t(data.frame(include))
        colnames(splitvar_include) = colnames(X)
    }else{
        splitvar_include = NULL
    }
    tree_txt = print.dipm(tree_txt, X, Y, C, treatment,
                        types, ncat, method, ntree, print,
                        splitvar_include)
    if(prune){
        tree_txt = pmprune(tree_txt)
    }
    
    if(dataframe){
        return(tree_txt)
    }else{
        tree_pn = ini_node(1, tree_txt, data, form_rhs[1], surv)
        if(surv){
            tree_py = party(tree_pn, data,
                            fitted = data.frame(
                                 "(fitted)" = fitted_node(tree_pn, data = data),
                                 "(response)" = Surv(Y, C), check.names = F),
                             terms = terms(form))
        }else{
            tree_py = party(tree_pn, data,
                            fitted = data.frame(
                                 "(fitted)" = fitted_node(tree_pn, data = data),
                                 "(response)" = Y, check.names = F),
                             terms = terms(form))
        }
        
        return(tree_py) 
    }

}
