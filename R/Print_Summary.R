#' Print for a fitted model represented by an \code{eefAnalytics} object.
#' @rdname print
#' @param x  Object of class \code{eefAnalytics}
#' @param ... Additional arguments of \code{\link[base]{print}}
#' @return  Print conditional and unconditional effect sizes.
#' @export
print.eefAnalytics <- function(x,...) {
  Checks <- sum(x$Function %in% c("srtBayes","crtBayes","mstBayes") )
  if(Checks==0){Approach="Frequentist"}else{Approach="Bayesian"}
  cat("\nModel Info:")
  cat("\n Method:    ", x$Method)
  cat("\n Design:    ", toupper(substr(x$Function,1,3)))
  cat("\n Approach:  ", Approach )
  cat("\n Function:  ", x$Function)
  if(is.null(x$Bootstrap)==FALSE) {
    if(x$Function!="srtFREQ") {cat("\n Bootstrap: ", x$Type)}
    cat("\n CI:        ", x$CI)
  }
  cat("\n---------\n")
  cat("\n")
  ES0=x$ES
  ES1= x$Unconditional$ES
  cat("Result for: Conditional effect size")
  cat("\n")
  print(ES0)
  cat("\n")
  cat("Result for: Unconditional effect size")
  cat("\n")
  print(ES1)
  cat("\n")
  if(sum(x$Function %in% c("srtBayes","crtBayes","mstBayes") )==0){
    cat("Please use summary to get more results \n")
  }else{
    cat("Please use summary to get more results")
    cat("\nAnd use the model object to check for convergence")
  }
}

####################################################################################################



#' Summary for a fitted model represented by an \code{eefAnalytics} object.
#' @rdname summary
#' @param object  Object of class \code{eefAnalytics}
#' @param ... Additional arguments of \code{\link[base]{summary}}
#' @return Returns relevant summary including Rhat and effective sample sizes.
#' @export
summary.eefAnalytics <- function(object,...){
  Checks <- sum(object$Function %in% c("srtBayes","crtBayes","mstBayes") )
  cat("\n method:       ", object$Method)
  cat("\n Design:       ", object$Function)
  if(Checks>0){cat("\n observations: ", length(object$Model$y))}

  res <- object
  if(Checks>0){
    Beta1 <- data.frame( summary(object$Model,pars=c("alpha","beta")))
    res$Beta <- cbind(object$Beta,round(Beta1[,c("sd","n_eff","Rhat")],2))
  }
  cat("\n")
  Beta <- res$Beta
  print(Beta)
  cat("\n")
  ES0=object$ES
  ES1= object$Unconditional$ES
  cat("Result for: Conditional effect size")
  cat("\n")
  print(ES0)
  cat("\n")
  cat("Result for: Unconditional effect size")
  cat("\n")
  print(ES1)
  cat("\n")

  class(res) <- "eefAnalyticssummary"
  invisible(res)
}

