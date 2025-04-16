
#' @title Panel-Generator for Visualization of A Precision Medicine Tree 
#' @description This function provides a new plot method for \code{dipm}
#'  and \code{spmtree}. It visualizes stratified treatment groups through 
#'  boxplots for a continuous outcome and survival plots for a survival outcome, 
#'  respectively.
#' 
#' @param obj A \code{party} tree object returned from either the 
#' \code{dipm()} or \code{spmtree()} function
#' @param ... Arguments passed on to plotfun
#' 
#' @details This function visualizes the precision medicine trees
#'  proposed in Chen and Zhang (2020a, b).
#' 
#' @return No return value, called for plot
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
#'             Seibold, H., Zeileis, A., and Hothorn, T. (2019). 
#'             model4you: An R package for personalised treatment 
#'             effect estimation. \emph{Journal of Open Research 
#'             Software} \strong{7}(1).
#'             
#'             Hothorn, T. and Zeileis, A. (2015). partykit: 
#'             a modular toolkit for recursive partytioning in R. 
#'             \emph{The Journal of Machine Learning Research} 
#'             \strong{16}(1), 3905-3909.
#'             
#' @seealso \code{\link{dipm}}, \code{\link{spmtree}}
#' 
#' @examples
#' 
#' #' #
#' # ... an example with a continuous outcome variable
#' #     and two treatment groups
#' #
#'
#' N = 100
#' set.seed(123)
#' 
#' if (!identical(tolower(Sys.getenv("NOT_CRAN")), "true")){
#' Sys.setenv(OMP_THREAD_LIMIT = "2")
#' }
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
#'     ((X[,1 ] > 0) & (X[, 3] <= 0)) *
#'         (20 * (1 - treatment) + 18 * treatment) + 
#'
#'     ((X[,1] > 0) & (X[,3] > 0)) *
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
#' tree = dipm(Y ~ treatment | ., data, mtry = 1, maxdepth = 3) 
#' plot(tree, terminal_panel = node_dipm)
#'             
#' 
#'                                     
#' @export
#' @import partykit
#' @import ggplot2
#' @import grid
#' @importFrom survival survfit
#' @importFrom stats as.formula 

node_dipm = function(obj, ...)
{
  ## old <- options()
  ## extract response
  y = obj$fitted[["(response)"]]
  is_surv = inherits(y, "Surv") || !is.null(obj$dat$C)
  is_num = is.numeric(y) || is.null(obj$dat$C)
  stopifnot(is_surv || is_num)
  
  ## panel function for nodes
  rval = function(node, 
         # .nid = function(node) paste0(nam[id_node(node)], 
         # ", n = ", node$info$nobs, ", treatment = ", 
         # node$info$trt)
         .nid = function(node) paste0( 
         "n = ", node$info$nobs, ", optimal trt =  ", 
         node$info$opt_trt)){

    ## extract data
    nam = names(obj)
    nid = id_node(node)
    dat = data_party(obj, nid)
    nid_text = .nid(node)
    form_lhs = all.vars(obj$terms[[2]])
    form_rhs = all.vars(obj$terms[[3]])
    
    top_vp = viewport(layout = grid.layout(nrow = 2,
                      heights = unit(c(0.1, 1), "null")))
    pushViewport(top_vp)
    
    node_vp = viewport(
      layout.pos.row = 2,
      width = unit(0.95, "npc"),
      height = unit(0.95, "npc")
    )
    pushViewport(node_vp)
    
    grid.rect(gp = gpar(fill = "white", cex = 0.5))
    
    ## viewport plot
    plotvp = viewport(
      width = unit(0.95, "npc"),
      height = unit(0.95, "npc")
    )
    pushViewport(plotvp)
    
    if(is_surv){
      
      ymax = max(dat[, form_lhs[1]])
      yrange = c(0, ymax)
      
      form = as.formula(paste0("Surv(",
                               form_lhs[1],
                               ",", form_lhs[2],
                               ")~", form_rhs[1]))
      dat[, form_rhs[1]] = as.factor(dat[, form_rhs[1]])
      mod = survfit(form, data = dat)
      mod_strata = mod$strata
      trt_ind = unlist(sapply(1:length(mod_strata), 
                 function(x){rep(levels(dat[, form_rhs[1]])[x],
                 mod_strata[x])}, simplify = T))
      pr = data.frame("treatment" = trt_ind,
                       "Y" = mod$time, 
                       "Probability" = mod$surv)
      colnames(pr) = c(form_rhs[1], form_lhs[1], "Probability")
      pr[[1]] = as.factor(pr[[1]])
      
      ## plot
      pl = ggplot(data = pr, 
                  aes_string(x = form_lhs[1], y = "Probability", 
                              group = form_rhs[1],
                              color = form_rhs[1])) +
        geom_step() + 
        coord_cartesian(xlim = yrange, ylim = 0:1) +
        theme_classic() + 
        theme(legend.position="top")
    }else{
      pr = data.frame("treatment" = dat[, form_rhs[1]],
                      "Y" = dat[, form_lhs[1]])
      colnames(pr) = c(form_rhs[1], form_lhs[1])
      pr[[1]] = as.factor(pr[[1]])
      
      pl = ggplot(data = pr, 
                  aes_string(x = form_rhs[1], y = form_lhs[1], 
                             color = form_rhs[1])) +
        geom_boxplot() + 
        theme_classic() + 
        theme(legend.position = "none")
    }
   
    print(pl, vp = plotvp)
    popViewport()
    
    nodeIDvp = viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
             width = max(unit(1, "lines"), 
                        unit(1.3, "strwidth", nid_text)),
             height = max(unit(1, "lines"), 
                          unit(1.3, "strheight", nid_text)))
    pushViewport(nodeIDvp)
    grid.rect(gp = gpar(fill = "white", cex = 0.5))
    grid.text(nid_text)
    popViewport()
    
    upViewport(n = 2)
  }

  ## new <- options()
  ## print(all.equal(old, new))
  return(rval)
}
class(node_dipm) = "grapcon_generator"

