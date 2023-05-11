#############################################################################
############# SRT Bayesian main functions ################################################


#' Analysis of Simple Randomised Education Trials using Bayesian Linear Regression Model with Vague Priors.
#'
#' \code{srtBayes} performs analysis of educational trials under the assumption of independent errors among pupils using Bayesian framework with Stan.
#' This can also be used with schools as fixed effects.
#'
#' @export
#' @param formula The model to be analysed is of the form y~x1+x2+.... Where y is the outcome variable and Xs are the independent variables.
#' @param intervention A string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.
#' @param baseln A string variable allowing the user to specify the reference category for intervention variable. When not specified, the first level will be used as a reference.
#' @param threshold a scalar or vector of pre-specified threshold(s) for estimating Bayesian posterior probability such that the observed effect size is greater than or equal to the threshold(s).
#' @param adaptD As this function uses rstanarm, this term provides the target average proposal acceptance probability during Stan’s adaptation period. Default is NULL.
#' @param nsim number of MCMC iterations per chain. Default is 2000.
#' @param condopt additional arguments of \code{\link[rstanarm]{stan_glm}} to be passed only to the conditional model specification (for example, defining priors only for the conditional model, etc.).
#' @param uncopt additional arguments of \code{\link[rstanarm]{stan_glm}} to be passed only to the unconditional model specification (for example, defining priors only for the unconditional model, etc.).
#' @param ... Additional arguments of \code{\link[rstanarm]{stan_glm}} to be passed both to the conditional and unconditional model specifications.
#' @param data Data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and credible intervals for the variables specified in the model. Use \code{summary.eefAnalytics} to get Rhat and effective sample size for each estimate.
#' \item \code{ES}: Conditional Hedges' g effect size and its 95% credible intervals.
#' \item \code{sigma2}: Residual variance.
#' \item \code{ProbES}: A matrix of Bayesian posterior probabilities such that the observed effect size is greater than or equal to a pre-specified threshold(s).
#' \item \code{Model}: A stan_glm object used in ES computation, this object can be used for convergence diagnostic.
#' \item \code{Unconditional}: A list of unconditional effect sizes, sigma2 and ProbES obtained based on residual variance from the unconditional model (model with only the intercept as a fixed effect).
#'  }
#' @example inst/examples/srtBExample.R
srtBayes <- function(formula,intervention,baseln,adaptD,nsim=2000,condopt,uncopt,data,threshold=1:10/10,...)UseMethod("srtBayes")

#' @export
srtBayes.default <- function(formula,intervention,baseln,adaptD,nsim=2000,condopt,uncopt,data,threshold=1:10/10,...){stop("No correct formula input given.")}

#' @export
srtBayes.formula <- function(formula,intervention,baseln,adaptD,nsim=2000,condopt,uncopt,data,threshold=1:10/10,...){

  #if(nsim<2000){stop("nsim must be at least 2000")}

  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),intervention))])
  tmp3 <- which(colnames(data)==intervention)
  if( missing(baseln)){data[,tmp3] <- as.factor(data[,tmp3])}
  if(!missing(baseln)){data[,tmp3] <- relevel(as.factor(data[,tmp3]),baseln)}
  mf <- model.frame(formula=formula, data=data)
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  posttest <- model.response(mf)
  intervention <- intervention
  trt <- data[,which(colnames(data)==intervention)]
  btp  <- which(tmp %in% paste0(intervention, unique(data[, intervention]) ))
  if(missing(condopt)){condopt=NULL}
  if(missing(uncopt)){uncopt=NULL}
  LHS.formula <- "post"
  RHS.formula <-paste(tmp[-1], collapse = "+")
  new.formula0 <- as.formula(post~1)
  new.formula  <- as.formula(paste0(LHS.formula,"~",RHS.formula ))
  new.data <- data.frame(post=posttest, fixedDesignMatrix[,-1])

  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}
  if(missing(adaptD)){adaptD=NULL}
  output <- srtB(formula=new.formula, unc.formula=new.formula0,intervention=intervention,data1=new.data,threshold=threshold,tmp=tmp,btp=btp,adaptD=adaptD,nsim=nsim,condopt=condopt,uncopt=uncopt,...)
  output$Method <- "LM"
  output$Function <- "srtBayes"
  class(output) <- "eefAnalytics"
  return(output)
}



################################################################################################################################
################################################################################################################################


## - internal
srtB <- function(formula,unc.formula,intervention,tmp,btp,data1,threshold,adaptD=adaptD,nsim=nsim,condopt,uncopt,...){

  freqFit <- do.call(stan_glm, c(list(formula = formula, adapt_delta=adaptD,iter=nsim,seed=1234, data = data1, ...), condopt))
  freqFit0 <-do.call(stan_glm, c(list(formula = unc.formula, adapt_delta=adaptD,iter=nsim, seed=1234, data = data1, ...), uncopt))

  Sim_Beta <- data.frame(as.matrix(freqFit))
  #btp <- which(substring(names(Sim_Beta),1,nchar(intervention))==intervention &
                 #nchar(names(Sim_Beta))==(nchar(intervention)+1))

  output1 <- NULL
  output2 <- NULL
  ProbES  <- list()
  for( i in 1:length(btp )){

    Output0=Bsrt.Summary(freqFit=freqFit,freqFit0=freqFit0,Sim_Beta= Sim_Beta, btp= btp[i])
    Beta <-  Output0$Beta
    Sim_Beta_t <-  Output0$Sim_Beta_t
    sigma2.1 <-  Output0$sigma2.1
    sigma2.2 <-  Output0$sigma2.2
    Sim_sigma2.1 <-  Output0$Sim_sigma2.1
    Sim_sigma2.2 <-  Output0$Sim_sigma2.2
    sim_ES1 <-  Output0$sim_ES1
    sim_ES2 <-  Output0$sim_ES2

    ProbES0 <- data.frame(cbind(Cond=sapply(threshold, function(j) mean(sim_ES1 > j)),
                                Uncond=sapply(threshold, function(j) mean(sim_ES2 > j))))

    row.names(ProbES0) <- paste0("P(ES>",threshold,")")
    ProbES[[i]] <- round(ProbES0,3)


    gtmp1 <- round(c("ES"=mean(sim_ES1), quantile(sim_ES1,probs=c(0.025,0.975))),2)
    gtmp2 <- round(c("ES"=mean(sim_ES2), quantile(sim_ES2,probs=c(0.025,0.975))),2)
    output1 <- rbind(output1,gtmp1)
    output2 <- rbind(output2,gtmp2)
  }

  row.names(output1) <- names(Sim_Beta)[btp]
  row.names(output2) <- names(Sim_Beta)[btp]
  names(ProbES) <- names(Sim_Beta)[btp]
  output1<- data.frame(output1)
  output2<- data.frame(output2)
  colnames(output1)<- c("Estimate","95% LB","95% UB")
  colnames(output2)<- c("Estimate","95% LB","95% UB")
  output <- list(Beta=round(Beta,2),ES=round(output1,2),sigma2=sigma2.1, ProbES=lapply(ProbES,function(x) x[, "Cond",drop=FALSE]), Model=freqFit,
                 Unconditional= list(ES=round(output2,2),sigma2=sigma2.2,ProbES=lapply(ProbES,function(x) x[, "Uncond",drop=FALSE])))
  return(output)
}

## summarise beta parameters - internal
Bsrt.Summary <- function(freqFit,freqFit0, Sim_Beta, btp){

  Sim_Beta_t <- as.matrix(freqFit,pars=names(Sim_Beta[btp])) #treatment effect
  Beta <- data.frame( summary(freqFit,pars=c("alpha","beta"), probs = c(0.025, 0.975), digits=2 ))[c("mean","X2.5.","X97.5.")]
  names(Beta)<- c("Estimate","95% LB","95% UB")

  Sim_sigma2.1 <- as.matrix(freqFit,pars="sigma")  #sigma(pupil): sqrt of residual
  Sim_sigma2.2 <- as.matrix(freqFit0,pars="sigma")
  sigma2.1 <- mean(Sim_sigma2.1^2)
  sigma2.2 <- mean(Sim_sigma2.2^2)
  sim_ES1 <- Sim_Beta_t/Sim_sigma2.1
  sim_ES2 <- Sim_Beta_t/Sim_sigma2.2


  ModelSmry <- list(Beta=Beta, Sim_Beta_t=Sim_Beta_t, sigma2.1=sigma2.1, Sim_sigma2.1=Sim_sigma2.1, sigma2.2=sigma2.2, Sim_sigma2.2=Sim_sigma2.2, sim_ES1=sim_ES1,sim_ES2=sim_ES2)
  return(ModelSmry)
}


