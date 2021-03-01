#' @import partykit
ini_node <- function(j, tree, data, trt, surv){
  if(!is.na(tree$splitvar[j])){
    if(tree$type[j]=="ord"||tree$type[j]=="bin"){
      split <- partysplit(as.integer(tree$splitvar[j]),
                          breaks = as.numeric(tree$splitval[j]))
    }else if(tree$type[j]=="nom"){
      temp = unlist(strsplit(tree$splitval[j], "{|,|}", perl = T))
      nom_splitvar <- temp[which(temp!="")]
      
      nom_levels <- sort(unique(data[,tree$splitvar_name[j]]))
      nom_index <- rep(2, length(nom_levels))
      nom_index[nom_levels %in% nom_splitvar] <- 1
      
      split <- partysplit(as.integer(tree$splitvar[j]),
                          index = as.integer(nom_index))
    }
    
    node <- partynode(split$varid, split = split,
                      # info = paste("trt: ", 
                      # tree$besttrt[j], sep=""),
    kids = list(ini_node(tree$lchild[j], tree, data, trt, surv),
            ini_node(tree$rchild[j], tree, data, trt, surv)))
  }else{
    if(class(data[,trt])=="integer"||class(data[,trt])=="numeric"){
      node <- partynode(as.integer(tree$node[j]),
                        info = list(
    opt_trt = as.character(unique(data[,trt])[
      which(unique(data[,trt])==tree$besttrt[j])]), 
    nobs = tree$nsubj[j]))
    }else{
      if(!surv){
        node <- partynode(as.integer(tree$node[j]),
                          info = list(
        opt_trt = as.character(levels(data[,trt])[tree$besttrt[j]]),
                          nobs = tree$nsubj[j]))
      }else{
        node <- partynode(as.integer(tree$node[j]),
                          info = list(
                          opt_trt = as.character(tree$besttrt[j]),
                          nobs = tree$nsubj[j]))
      }

    }

  }
  return(node)
}

#' @title Prediction of Optimal Treatment 
#' @description This function predicts the optimal treatment 
#' for a new subject based on either \code{dipm} or \code{spmtree}.
#' 
#' @param tree A \code{party} tree object returned from either the 
#' \code{dipm()} or \code{spmtree()} function
#' @param newdata A matrix or data frame of the new data
#' 
#' @details This function implements the optimal treatment 
#' prediction proposed in Chen and Zhang (2020a, b).
#' 
#' @return \code{predict.dipm} returns the predicted optimal 
#' treatment for a new subject 
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
#' @seealso \code{\link{dipm}}, \code{\link{spmtree}}
#' 
#' @examples
#' 
#' #
#' # ... an example with a continuous outcome variable
#' #     and two treatment groups
#' #
#'
#' N=100
#' set.seed(123)
#'
#' # generate binary treatments
#' treatment=rbinom(N,1,0.5)
#'
#' # generate candidate split variables
#' X1=rnorm(n=N,mean=0,sd=1)
#' X2=rnorm(n=N,mean=0,sd=1)
#' X3=rnorm(n=N,mean=0,sd=1)
#' X4=rnorm(n=N,mean=0,sd=1)
#' X5=rnorm(n=N,mean=0,sd=1)
#' X=cbind(X1,X2,X3,X4,X5)
#' colnames(X)=paste0("X",1:5)
#'
#' # generate continuous outcome variable
#' calculateLink <- function(X,treatment) {
#'
#'     ( (X[,1] <= 0) & (X[,2] <= 0) )*
#'         ( 25*(1-treatment) + 8*treatment) + 
#'
#'     ( (X[,1] <= 0) & (X[,2] > 0) )*
#'         ( 18*(1-treatment) + 20*treatment ) +
#'
#'     ( (X[,1] > 0) & (X[,3] <= 0) )*
#'         ( 20*(1-treatment) + 18*treatment ) + 
#'
#'     ( (X[,1] > 0) & (X[,3] > 0) )*
#'         ( 8*(1-treatment) + 25*treatment )
#' }
#'
#' Link=calculateLink(X,treatment)
#' Y=rnorm(N,mean=Link,sd=1)
#'
#' # combine variables in a data frame
#' data=data.frame(X,Y,treatment)
#' 
#' # fit a dipm classification tree
#' tree=dipm(Y~treatment | .,data,mtry=1,maxdepth=3) 
#' predict_dipm(tree, newdata = head(data))
#'           
#'                                     
#' @export
#' @import partykit
#' @importFrom utils head
#' @importFrom stats predict

predict_dipm <- function(tree, newdata){
  predict(tree, newdata=head(newdata),
          FUN = function(n) paste("Optimal treatment is:", 
                                  n$info$opt_trt))
}
