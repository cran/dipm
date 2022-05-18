
#' @title Pruning A Precision Medicine Tree
#' @description This function prunes classification trees designed
#'              for the precision medicine setting.
#'
#' @param tree A data frame object returned from either the 
#'             \code{dipm()} or \code{spmtree()} function
#'
#' @details This function implements the simple pruning strategy
#'          proposed and used in Tsai et al. (2016). Terminal
#'          sister nodes, i.e., nodes with no child nodes that
#'          share the same parent node, are removed if they have
#'          the same identified optimal treatment assignment.
#'
#' @return \code{pmprune} returns the pruned classification tree as a 
#'         data frame. The data frame contains the following columns 
#'         of information:
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
#' @references Tsai, W.-M., Zhang, H., Buta, E., O'Malley, S., 
#'             Gueorguieva, R. (2016). A modified classification
#'             tree method for personalized medicine decisions.
#'             \emph{Statistics and its Interface} \strong{9}, 
#'             239-253.
#'
#' @seealso \code{\link{dipm}}, \code{\link{spmtree}}
#'
#' @examples
#' 
#' #
#' # ... an example with a continuous outcome variable
#' #     and three treatment groups
#' #
#' 
#' \donttest{
#' N = 100
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
#' calculateLink = function(X, treatment){
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
#'             "response", "treatment")
#' 
#' # fit a classification tree
#' tree = spmtree(Y ~ treatment | ., data, types = types, dataframe = TRUE)
#'
#' # prune the tree
#' ptree = pmprune(tree)
#' }
#' 
#' @export

pmprune = function(tree){
#
#  This function accepts as an argument a tree data
#  frame object and returns a pruned tree. Terminal
#  node pairs are pruned when they have the same best
#  treatment class.
#
    
    if(!inherits(tree, "data.frame")){
        stop("The input must be a data frame.")
    }
    
    if(nrow(tree) == 1){
        return(tree)
    }

    lll = seq(2, nrow(tree), 2)
    lll = lll[order(-lll)]

    for(i in lll){

#        do not consider non-terminal nodes
        if(is.na(tree$lchild[i]) == FALSE) next()

#        the pair of nodes must be two terminal nodes
        if(is.na(tree$lchild[i + 1]) == FALSE) next()

#        the pair of nodes must have the same best treat to be pruned
        if(tree$besttrt[i] != tree$besttrt[i + 1]) next()

        j = which(tree$lchild == tree$node[i])
        tree$lchild[j] = NA
        tree$rchild[j] = NA
        tree$splitvar[j] = NA

        # to be compatible with party
        # tree = tree[-c(i, i + 1), ] 
    }

    return(tree)
}
