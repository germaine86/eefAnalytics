
###################################################################################################
############# Bayesian Multilevel Analysis of Cluster Randomised Education Trials
##################################################################################################

#' Bayesian analysis of cluster randomised education trials using Vague Priors.
#'
#' \code{crtBayes} performs analysis of cluster randomised education trials using a multilevel model under a Bayesian setting,
#' assuming vague priors.
#'
#' @export
#' @param formula the model to be analysed is of the form y ~ x1+x2+.... Where y is the outcome variable and Xs are the independent variables.
#' @param random a string variable specifying the "clustering variable" as contained in the data. See example below.
#' @param intervention a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.
#' @param baseln A string variable allowing the user to specify the reference category for intervention variable. When not specified, the first level will be used as a reference.
#' @param adaptD As this function uses rstanarm, this term provides the target average proposal acceptance probability during Stan’s adaptation period. Default is NULL.
#' @param nsim number of MCMC iterations per chain. Default is 2000.
#' @param condopt additional arguments of \code{\link[rstanarm]{stan_glm}} to be passed only to the conditional model specification (for example, defining priors only for the conditional model, etc.).
#' @param uncopt additional arguments of \code{\link[rstanarm]{stan_glm}} to be passed only to the unconditional model specification (for example, defining priors only for the unconditional model, etc.).
#' @param threshold a scalar or vector of pre-specified threshold(s) for estimating Bayesian posterior probability such that the observed effect size is greater than or equal to the threshold(s).
#' @param ... additional arguments of \code{\link[rstanarm:stan_glmer]{stan_lmer}} to be passed both to the conditional and unconditional model specifications.
#' @param data data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and credible intervals for variables specified in the model. Use \code{summary.eefAnalytics} to get Rhat and effective sample size for each estimate.
#' \item \code{ES}: Conditional Hedges' g effect size and its 95% credible intervals.
#' \item \code{covParm}: A vector of variance decomposition into between cluster variance (Schools) and within cluster variance (Pupils). It also contains intra-cluster correlation (ICC).
#' \item \code{SchEffects}: A vector of the estimated deviation of each school from the intercept.
#' \item \code{ProbES}: A matrix of Bayesian Posterior Probabilities such that the observed effect size is greater than or equal to a pre-specified threshold(s).
#' \item \code{Model}: A stan_glm object used in ES computation, this object can be used for convergence diagnostic.
#' \item \code{Unconditional}: A list of unconditional effect sizes, covParm and ProbES obtained based on between and within cluster variances from the unconditional model (model with only the intercept as a fixed effect).
#' }
#' @example inst/examples/crtBExample.R
crtBayes <- function(formula,random,intervention,baseln,adaptD,nsim=2000,condopt,uncopt,data,threshold=1:10/10,...)UseMethod("crtBayes")

#' @export
crtBayes.default <- function(formula,random,intervention,baseln,adaptD,nsim=2000,condopt,uncopt,data,threshold=1:10/10,...){stop("No correct formula input given.")}

#' @export
crtBayes.formula <- function(formula,random,intervention,baseln,adaptD,nsim=2000,condopt,uncopt,data,threshold=1:10/10,...){

  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),random, intervention))])
  data <- data[order(data[,which(colnames(data)==random)]),]

  intervention <- intervention
  trt <- data[,which(colnames(data)==intervention)]
  tmp2 <- which(colnames(data)==random)
  cluster2 <-  data[,tmp2]

  chk <- sum(rowSums(table(cluster2,trt)!=0)>1)
  if(chk >0){stop("This is not a CRT design")}
  stp <- as.character(row.names(table(cluster2,trt)))
  stp2 <- (apply(table(cluster2,trt),1,function(x)colnames(table(cluster2,trt))[x!=0]))


  tmp3 <- which(colnames(data)==intervention)
  #data[,tmp3] <- as.factor(data[,tmp3])
  if(missing(baseln)){data[,tmp3] <- as.factor(data[,tmp3])}
  if(!missing(baseln)){data[,tmp3] <- relevel(as.factor(data[,tmp3]),baseln)}
  if(missing(adaptD)){adaptD=NULL}
  mf <- model.frame(formula=formula, data=data)
  mf <- mf[order(cluster2),]
  cluster <- cluster2[order(cluster2)]
  trt <- trt[order(cluster2)]
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  posttest <- model.response(mf)
  new.data <- data.frame(post=posttest, fixedDesignMatrix[,-1],cluster=cluster)
  btp  <- which(tmp %in% paste0(intervention, unique(data[, intervention]) ))
  if(missing(condopt)){condopt=NULL}
  if(missing(uncopt)){uncopt=NULL}

  LHS.formula <- "post"
  RHS.formula <-paste(tmp[-1], collapse = "+")
  RHS.cluster <-"+(1|cluster)"
  new.formula <-  as.formula(paste0(LHS.formula,"~",RHS.formula,RHS.cluster))
  new.formula0 <- as.formula(paste0(LHS.formula,"~",RHS.cluster))


  nsim=nsim
  if(length(tmp2)!= 1){stop("Cluster variable misspecified")}
  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}
  #if(nsim < 2000){stop("nsim >= 10000 is recommended")}

  output <-  erantBAYES(formula=new.formula, unc.formula=new.formula0,intervention=intervention,data1=new.data,cluster=random,threshold=threshold,btp=btp,adaptD=adaptD,nsim=nsim,condopt=condopt,uncopt=uncopt,...)
  output$Method <- "MLM"
  output$Function <- "crtBayes"
  class(output) <- "eefAnalytics"
  invisible(output)
}

#############################################################################################################################
################################################################################################################################






## - internal

erantBAYES<- function(formula,unc.formula,intervention,data1,cluster,btp,threshold,adaptD=adaptD,nsim=nsim,condopt,uncopt,...){

  freqFit <- do.call(stan_lmer, c(list(formula = formula, seed=1234,adapt_delta=adaptD,iter=nsim, data = data1, ...), condopt))
  freqFit0 <-do.call(stan_lmer, c(list(formula = unc.formula, seed=1234,adapt_delta=adaptD,iter=nsim, data = data1, ...), uncopt))

  Sim_Beta <- data.frame(as.matrix(freqFit))
  # btp <- which(substring(names(Sim_Beta),1,nchar(intervention))==intervention &
  #                nchar(names(Sim_Beta))==(nchar(intervention)+1))
  #btp <- which(names(Sim_Beta) %in% Int.lev)

  output1 <- list()
  output2 <- list()
  ProbES <- list()
  for( i in 1:length(btp )){

    Output0=errantSummary(freqFit=freqFit,freqFit0=freqFit0,Sim_Beta= Sim_Beta, btp= btp[i])
    Beta <-  Output0$Beta
    ICC1 <- Output0$ICC1
    Sim_Beta_t <-  Output0$Sim_Beta_t
    var.B1 <- Output0$var.B1
    var.W1 <- Output0$var.W1
    var.tt1 <- Output0$var.tt1

    sim_ES.tt1 <- Sim_Beta_t/sqrt(var.tt1)
    sim_ES.W1 <- Sim_Beta_t/sqrt(var.W1)

    ICC2 <- Output0$ICC2
    var.W2 <- Output0$var.W2
    var.tt2 <- Output0$var.tt2
    sim_ES.tt2 <- Sim_Beta_t/sqrt(var.tt2)
    sim_ES.W2 <- Sim_Beta_t/sqrt(var.W2)


    ProbES0 <- data.frame(cbind(Total1=sapply(threshold, function(j) mean(sim_ES.tt1 > j)),
                                Within1 =sapply(threshold, function(j) mean(sim_ES.W1 > j)),
                                Total2=sapply(threshold, function(j) mean(sim_ES.tt2 > j)),
                                Within2 =sapply(threshold, function(j) mean(sim_ES.W2 > j))))

    row.names(ProbES0) <- paste0("P(ES>",threshold,")")
    ProbES[[i]] <- round(ProbES0,3)


    gtmp.W1 <- round(c("ES"=mean(sim_ES.W1), quantile(sim_ES.W1,probs=c(0.025,0.975))),2)
    gtmp.W2 <- round(c("ES"=mean(sim_ES.W2), quantile(sim_ES.W2,probs=c(0.025,0.975))),2)
    gtmp.tt1 <- round(c("ES"=mean(sim_ES.tt1), quantile(sim_ES.tt1,probs=c(0.025,0.975))),2)
    gtmp.tt2 <- round(c("ES"=mean(sim_ES.tt2), quantile(sim_ES.tt2,probs=c(0.025,0.975))),2)
    output.1 <- rbind(Within=gtmp.W1,Total=gtmp.tt1)
    output.2 <- rbind(Within=gtmp.W2,Total=gtmp.tt2)
    colnames(output.1) <- c("Estimate","95% LB","95% UB")
    colnames(output.2) <- c("Estimate","95% LB","95% UB")
    output1[[i]] <- round(output.1,2)
    output2[[i]] <- round(output.2,2)
  }

  names(output1) <- names(Sim_Beta)[btp]
  names(output2) <- names(Sim_Beta)[btp]
  Output <- list(Beta=round(Beta,2), ES=output1, covParm=Output0$sigmaBE1, SchEffects=Output0$SchEffects, ProbES=lapply(ProbES, function(x) x[, c("Within1","Total1")]),Model=freqFit,
                 Unconditional= list(ES=output2, covParm=Output0$sigmaBE2, ProbES=lapply(ProbES, function(x) x[, c("Within2","Total2")])))
  return(Output)
}



## - internal

errantSummary <- function(freqFit,freqFit0,Sim_Beta, btp){

  Sim_Beta_t <- as.matrix(freqFit,pars=names(Sim_Beta[btp])) #treatment effect
  Beta <- data.frame( summary(freqFit,pars=c("alpha","beta"), probs = c(0.025, 0.975), digits=2 ))[c("mean","X2.5.","X97.5.")]
  names(Beta)<- c("Estimate","95% LB","95% UB")

  #conditional ouptut
  Sim_var.W1 <- as.matrix(freqFit,pars="sigma")^2  #sigma(pupil)
  Sim_var.B1 <- as.matrix(freqFit,pars="Sigma[cluster:(Intercept),(Intercept)]")
  Sim_var.tt1 <- Sim_var.B1+Sim_var.W1
  Sim_ICC1 <- Sim_var.B1/Sim_var.tt1
  sigmaBE1 <- c(mean(Sim_var.B1),mean(Sim_var.W1),mean(Sim_var.tt1),mean(Sim_ICC1))
  names(sigmaBE1)<- c("Schools","Pupils","Total","ICC")

  #conditional ouptut
  Sim_var.W2 <- as.matrix(freqFit0,pars="sigma")^2  #sigma(pupil)
  Sim_var.B2 <- as.matrix(freqFit0,pars="Sigma[cluster:(Intercept),(Intercept)]")
  Sim_var.tt2 <- Sim_var.B2+Sim_var.W2
  Sim_ICC2 <- Sim_var.B2/Sim_var.tt2
  sigmaBE2 <- c(mean(Sim_var.B2),mean(Sim_var.W2),mean(Sim_var.tt2),mean(Sim_ICC2))
  names(sigmaBE2)<- c("Schools","Pupils","Total","ICC")

  #random effects
  SchEffects0 <- summary(freqFit, regex_pars = "b\\[\\(Intercept\\) cluster\\:", probs = c(0.025, 0.975),digits = 2)
  SchEffects1 <- data.frame(SchEffects0)[,"mean"]
  schl_ID <-  as.numeric(as.character(substr(rownames(SchEffects0), nchar("b[(Intercept) cluster\\:"), nchar(rownames(SchEffects0))-1)))
  SchEffects <- round(data.frame(cbind(Schools=schl_ID,Intercept=SchEffects1)),2)



  #all output
  ModelSmry <- list(Beta=Beta, Sim_Beta_t=Sim_Beta_t,
                    var.B1=Sim_var.B1, var.W1=Sim_var.W1, var.tt1=Sim_var.tt1, ICC1=Sim_ICC1, sigmaBE1=round(sigmaBE1,2),
                    var.B2=Sim_var.B2, var.W2=Sim_var.W2, var.tt2=Sim_var.tt2, ICC2=Sim_ICC2, sigmaBE2=round(sigmaBE2,2),
                    SchEffects =SchEffects)
  return(ModelSmry)
}






