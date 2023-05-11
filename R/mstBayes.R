
###################################################################################################
############# Bayesian Multilevel Analysis of Multisite Randomised Education Trials
##################################################################################################

#' Bayesian analysis of Multisite Randomised Education Trials using Vague Priors.
#'
#' \code{mstBayes} performs analysis of multisite randomised education trials using a multilevel model under a Bayesian setting
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
#' @param threshold a scalar or vector of pre-specified threshold(s) for estimating Bayesian posterior probability that the observed effect size is greater than or equal to the threshold(s).
#' @param ... additional arguments of \code{\link[rstanarm:stan_glmer]{stan_lmer}} to be passed both to the conditional and unconditional model specifications.
#' @param data data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and credible intervals for variables specified in the model. Use \code{summary.eefAnalytics} to get Rhat and effective sample size for each estimate.
#' \item \code{ES}: Conditional Hedges' g effect size and its 95% credible intervals.
#' \item \code{covParm}: A list of variance decomposition into between cluster variance-covariance matrix (schools and school by intervention) and within cluster variance (Pupils). It also contains intra-cluster correlation (ICC).
#' \item \code{SchEffects}: A vector of the estimated deviation of each school from the intercept and intervention slope.
#' \item \code{ProbES}: A matrix of Bayesian posterior probabilities such that the observed effect size is greater than or equal to a pre-specified threshold(s).
#' \item \code{Model}: A stan_glm object used in ES computation, this object can be used for convergence diagnostic.
#' \item \code{Unconditional}: A list of unconditional effect sizes, covParm and ProbES obtained based on between and within cluster variances from the unconditional model (model with only the intercept as a fixed effect).
#' }
#' @example inst/examples/mstBExample.R
mstBayes <- function(formula,random,intervention,baseln,adaptD,nsim=2000,condopt,uncopt,data,threshold=1:10/10,...)UseMethod("mstBayes")

#' @export
mstBayes.default <- function(formula,random,intervention,baseln,adaptD,nsim=2000,condopt,uncopt,data,threshold=1:10/10,...){stop("No correct formula input given.")}

#' @export
mstBayes.formula <- function(formula,random,intervention,baseln,adaptD,nsim=2000,condopt,uncopt,data,threshold=1:10/10,...){

  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),random, intervention))])
  data <- data[order(data[,which(colnames(data)==random)]),]

  intervention <- intervention
  trt <- data[,which(colnames(data)==intervention)]
  tmp2 <- which(colnames(data)==random)
  cluster2 <-  data[,tmp2]

  chk <- sum(rowSums(table(cluster2,trt)!=0)>1)
  if(chk ==0){stop("This is not a MST design")}
  stp <- as.character(row.names(table(cluster2,trt)))
  #stp2 <- as.numeric(apply(table(cluster2,trt),1,function(x)colnames(table(cluster2,trt))[x!=0]))


  tmp3 <- which(colnames(data)==intervention)
  #data[,tmp3] <- as.factor(data[,tmp3])
  if( missing(baseln)){data[,tmp3] <- as.factor(data[,tmp3])}
  if(!missing(baseln)){data[,tmp3] <- relevel(as.factor(data[,tmp3]),baseln)}
  mf <- model.frame(formula=formula, data=data)
  mf <- mf[order(cluster2),]
  cluster <- cluster2[order(cluster2)]
  trt <- trt[order(cluster2)]
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  #btp  <- which(substring(tmp,1,nchar(intervention))==intervention & nchar(tmp)==(nchar(intervention)+1))
  if(missing(condopt)){condopt=NULL}
  if(missing(uncopt)){uncopt=NULL}
  btp  <- which(tmp %in% paste0(intervention, unique(data[, intervention]) ))
  btp1 <- tmp[btp]
  btp2 <- gsub(intervention,"trt", btp1)
  btp3 <- data.frame(fixedDesignMatrix[,btp1])
  colnames(btp3) <- btp2
  if(missing(adaptD)){adaptD=NULL}
  posttest <- model.response(mf)
  new.data <- data.frame(post=posttest, fixedDesignMatrix[,-1],btp3,cluster=cluster)

  LHS.formula <- "post"
  RHS.formula <-paste(tmp[-1], collapse = "+")
  RHS.cluster0 <-"+(1|cluster)"
  RHS.cluster <-paste0("+(1+", paste0(btp2, collapse = "+"),"|cluster)")
  new.formula <-  as.formula(paste0(LHS.formula,"~",RHS.formula,RHS.cluster))
  new.formula0 <- as.formula(paste0(LHS.formula,"~",RHS.cluster0))


  nsim=nsim
  if(length(tmp2)!= 1){stop("Cluster variable misspecified")}
  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}
  #if(nsim < 2000){stop("nsim >= 10000 is recommended")}

  output <-  erantMstBAYES(formula=new.formula, unc.formula=new.formula0,intervention=intervention,data1=new.data,cluster=random,threshold=threshold,btp=btp,btp2=btp2,adaptD=adaptD,nsim=nsim,condopt=condopt,uncopt=uncopt,...)
  output$Method <- "MLM"
  output$Function <- "mstBayes"
  class(output) <- "eefAnalytics"
  return(output)
}


################################################################################################################################
################################################################################################################################







## - internal

erantMstBAYES<- function(formula,unc.formula,intervention,data1,cluster,threshold,btp,btp2,adaptD=adaptD,nsim=nsim,condopt,uncopt,...){

  freqFit <- do.call(stan_lmer, c(list(formula = formula, seed=1234,adapt_delta=adaptD,iter=nsim, data = data1, ...), condopt))
  freqFit0 <-do.call(stan_lmer, c(list(formula = unc.formula, seed=1234,adapt_delta=adaptD,iter=nsim, data = data1, ...), uncopt))

  Sim_Beta <- data.frame(as.matrix(freqFit))

  output1 <- list()
  output2 <- list()
  ProbES <- list()
  for( i in 1:length(btp )){

    Output0=errantMStSummary(freqFits=freqFit,Cond=F,Sim_Beta= Sim_Beta, btp= btp[i],btp1=btp,btp2=btp2 )
    Beta <-  Output0$Beta
    Sim_Beta_t <-  Output0$Sim_Beta_t
    var.B1 <- Output0$var.B1
    var.W1 <- Output0$var.W1
    var.tt1 <- Output0$var.tt1
    sim_ES.tt1 <- Sim_Beta_t/sqrt(var.tt1)
    sim_ES.W1 <- Sim_Beta_t/sqrt(var.W1)

    Output01=errantMStSummary(freqFits=freqFit0,Cond=T,Sim_Beta= Sim_Beta, btp= btp[i],btp1=btp,btp2=btp2 )

    var.W2 <- Output01$var.W1
    var.tt2 <- Output01$var.tt1
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

  names(ProbES) <- names(Sim_Beta)[btp]
  names(output1) <- names(Sim_Beta)[btp]
  names(output2) <- names(Sim_Beta)[btp]
  Output <- list(Beta=round(Beta,2), ES=output1, covParm=Output0$sigmaBE1, SchEffects=Output0$SchEffects, ProbES=lapply(ProbES, function(x) x[, c("Within1","Total1")]),Model=freqFit,
                 Unconditional= list(ES=output2, covParm=Output01$sigmaBE1, ProbES=lapply(ProbES, function(x) x[, c("Within2","Total2")])))
  return(Output)
}



## - internal


errantMStSummary <- function(freqFits, Cond=F,Sim_Beta, btp, btp1,btp2){

  data1 <- freqFits$data
  N <- length(freqFits$residuals)
  N.t <- sapply(btp1, function(i) sum(data1[,i])) # need to check how to add this

  if(Cond==T){
    # Unconditional results
    Beta <- NULL
    Sim_Beta_t <- NULL
    SchEffects <- NULL
    Sim_var.B1 <- as.matrix(freqFits,pars="Sigma[cluster:(Intercept),(Intercept)]")
    Sim_var.W1 <- as.matrix(freqFits,pars="sigma")^2  #sigma(pupil)
    Sim_var.tt1 <- Sim_var.B1+Sim_var.W1
    Sim_ICC1 <- Sim_var.B1/Sim_var.tt1
    sigmaBE1 <- round(c(mean(Sim_var.B1),mean(Sim_var.W1),mean(Sim_var.tt1),mean(Sim_ICC1)),2)
    names(sigmaBE1)<- c("Schools","Pupils","Total","ICC")
  }




  if(Cond==F){
    #conditional ouptut
    Sim_Beta_t <- as.matrix(freqFits,pars=names(Sim_Beta[btp])) #treatment effect
    Beta <- data.frame( summary(freqFits,pars=c("alpha","beta"), probs = c(0.025, 0.975), digits=2 ))[c("mean","X2.5.","X97.5.")]
    names(Beta)<- c("Estimate","95% LB","95% UB")

    Sim_var.B   <- as.matrix(freqFits, regex_pars ="Sigma\\[\\cluster\\:trt")

    name.varW <- "sigma"
    name.Vsch <- "Sigma[cluster:(Intercept),(Intercept)]"
    name.diag <- grep(paste(paste0(btp2,",",btp2),collapse = "|"), colnames(Sim_var.B), value = T)
    name.diag0<- grep(paste(paste0(btp2,",",btp2),collapse = "|"), colnames(Sim_var.B), value = T, invert = T)
    name.cov  <- grep(paste(paste0(btp2,",","\\(Intercept\\)"),collapse = "|"), colnames(Sim_var.B), value = T)


    Sim_var.W1 <- as.matrix(freqFits, pars=name.varW)^2  #sigma(pupil)
    Sim_var.sch1 <- as.matrix(freqFits, pars=name.Vsch)
    Sim_var.B1  <- Sim_var.B[, name.diag]
    Sim_cov.B1  <- Sim_var.B[, name.cov]
    Sim_Term1.1 <-  Sim_var.W1+Sim_var.sch1
    Sim_Term2.1 <- 2*Sim_cov.B1+Sim_var.B1
    sim_Nt_B1  <- rowSums(matrix(Sim_Term2.1, ncol = length(btp1))%*%matrix(N.t))
    Sim_var.tt1 <- (N*Sim_Term1.1 + sim_Nt_B1)/N
    Sim_ICC1 <- rowSums(matrix(Sim_var.B1, ncol = length(btp1)))/Sim_var.tt1


    ## variance covariance matrix
    varAll <- data.frame(summary(freqFits, regex_pars = "Sigma\\[cluster\\:", probs = 0.025,digits = 2))
    diag <- varAll[c(name.Vsch,name.diag),"mean"]
    offdiag <- varAll[name.diag0,"mean"]
    Vcov1 <- matrix(NA, ncol = length(diag), nrow = length(diag))
    rownames(Vcov1) <- c("Intercept",btp2)
    colnames(Vcov1) <- c("Intercept",btp2)
    Vcov1[lower.tri(Vcov1)] <- offdiag
    Vcov1[upper.tri(Vcov1)] <- t(Vcov1)[upper.tri(t(Vcov1))]
    diag(Vcov1) <- diag
    sigmaBE1 <- list(round(Vcov1,2),round(mean(Sim_var.W1),2),round(mean(Sim_var.tt1),2),round(mean(Sim_ICC1),2))
    names(sigmaBE1)<- c("School:trt","Pupils","Total","ICC")

    #random effects
    SchEffects00 <- data.frame(summary(freqFits, regex_pars = "b\\[\\(Intercept\\) cluster\\:", probs = c(0.025, 0.975),digits = 2))
    SchEffects0 <- SchEffects00[,"mean"]
    SchEffects1 <- sapply(btp2,function(x) data.frame(summary(freqFits, regex_pars = paste0("b\\[",x, " cluster\\:"),digits = 2))[,"mean"])

    schl_ID <-  as.numeric(as.character(substr(rownames(SchEffects00), nchar("b[(Intercept) cluster\\:"), nchar(rownames(SchEffects00))-1)))
    SchEffects <- round(data.frame(cbind(Schools=schl_ID,Intercept=SchEffects0,SchEffects1)),2)
  }


  #all output
  ModelSmry <- list(Beta=Beta, Sim_Beta_t=Sim_Beta_t,
                    var.B1=Sim_var.B1, var.W1=Sim_var.W1, var.tt1=Sim_var.tt1, ICC1=Sim_ICC1, sigmaBE1=sigmaBE1,
                    SchEffects =SchEffects)
  return(ModelSmry)
}
