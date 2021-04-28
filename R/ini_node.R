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

