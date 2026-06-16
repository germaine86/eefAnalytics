###################################################################################################
############# Bayesian Multilevel Analysis of Simple Randomised Education Trials (SRT) ###############
##################################################################################################

#' Bayesian Analysis of Simple Randomised Education Trials (SRT) using Bayesian Linear Regression Model with Vague Priors.
#'
#' \code{srtBayes} performs Bayesian multilevel analysis of Simple Randomised Education Trials (SRT), utilising vague priors
#' and JAGS language to fit the model.
#' This can also be used with schools as fixed effects.
#'
#' @export
#' @param formula The model to be analysed is of the form y~x1+x2+.... Where y is the outcome variable and Xs are the independent variables.
#' @param intervention A string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.
#' @param baseln A string variable allowing the user to specify the reference category for intervention variable. When not specified, the first level will be used as a reference.
#' @param nsim number of MCMC iterations per chain. Default is 2000.
#' @param threshold a scalar or vector of pre-specified threshold(s) for estimating Bayesian posterior probability such that the observed effect size is greater than or equal to the threshold(s).
#' @param data Data frame containing the data to be analysed.
#' @param condopt additional arguments of \code{\link[R2jags]{jags}} to be passed only to the conditional model specification (for example, defining priors only for the conditional model, etc.).
#' @param uncopt additional arguments of \code{\link[R2jags]{jags}} to be passed only to the unconditional model specification (for example, defining priors only for the unconditional model, etc.).
#' @param alpha significant level, default alpha = 0.05.
#' @param digits number of decimal places, by default digits=3
#' @param ... Common additional arguments of \code{\link[R2jags]{jags}} to be passed to both the conditional and unconditional model specifications
#'
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and credible intervals for the variables specified in the model. Use \code{summary.eefAnalytics} to get Rhat and effective sample size for each estimate.
#' \item \code{ES}: Conditional Hedges' g effect size and its 95% credible intervals.
#' \item \code{sigma2}: Residual variance.
#' \item \code{ProbES}: A matrix of Bayesian posterior probabilities such that the observed effect size is greater than or equal to a pre-specified threshold(s).
#' \item \code{Model}: A model object from \code{\link[R2jags]{jags}} and an \code{\link[MCMCvis]{MCMCsummary}} object containing only the mean and credible intervals (CIs) as columns.
#' \item \code{Unconditional}: A list of unconditional effect sizes, sigma2 and ProbES obtained based on residual variance from the unconditional model (model with only the intercept as a fixed effect).
#' }
#'
#' @example inst/examples/srtBExample.R
srtBayes  <- function(formula,intervention,baseln,nsim=2000,data,alpha=0.05,digits=3,threshold=1:10/10,condopt,uncopt,...){UseMethod("srtBayes")}

#' @export
srtBayes.default <- function(formula,intervention,baseln,nsim=2000,data,alpha=0.05,digits=3,threshold=1:10/10,condopt,uncopt,...){stop("No correct formula input given.")}

#' @export
srtBayes.formula <- function(formula,intervention,baseln,nsim=2000,data,alpha=0.05,digits=3,threshold=1:10/10,condopt,uncopt,...){

  #data and #stop and warning rules
  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),intervention))])
  tmp3 <- which(colnames(data)==intervention)
  if(!missing(baseln)){data[,tmp3] <- relevel(as.factor(data[,tmp3]),as.character(baseln))}
  if(missing(baseln)){data[,tmp3] <- as.factor(data[,tmp3])}

  if(nsim < 10000){warning("nsim >= 10000 is recommended")}
  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}


  if(missing(condopt)){condopt=NULL}
  if(missing(uncopt)){uncopt=NULL}
  output <- OLS.function(data=data,
                         formula=formula,
                         intervention=intervention,
                         nsim=nsim,
                         alpha=alpha,
                         digits=digits,
                         threshold=threshold,
                         condopt=condopt,
                         uncopt=uncopt)#BayesOutput


  output$Method <- "LM"
  output$Function <- "srtBayes"
  class(output) <- "eefAnalytics"
  return(output)
}



################################################################################################################################
################################################################################################################################
OLS.function <- function(data, formula, intervention, nsim,alpha, digits,threshold, condopt, uncopt,...){

  #### load required packages ###
  requireNamespace("R2jags", quietly = TRUE) || stop("Please install the 'R2jags' package.")
  #require(R2jags)
  requireNamespace("MCMCvis", quietly = TRUE) || stop("Please install the 'MCMCvis' package.")
  #require(MCMCvis)
  requireNamespace("coda", quietly = TRUE) || stop("Please install the 'coda' package.")
  #require(coda)
  #function(data, formula, intervention,nsim){

  #### data preparation
  Pdata <- na.omit(data[,all.vars(formula)])
  outcome <- all.vars(formula)[1] ## Y: Posttest
  dummies<- data.frame(model.matrix(formula, data=Pdata)) ## X0-Intercept, X1-Prettest, X2-Intervention
  Pdata0 <- data.frame(na.omit(dummies[,!(names(dummies) %in% "X.Intercept.")]))
  names(Pdata0) <- setdiff(names(dummies), "X.Intercept.")
  Pdata1 <- na.omit(cbind(post=Pdata[,outcome],Pdata0)) #all covariates and school dummies


  # Jags data
  var<-names(Pdata1)
  N <- nrow(Pdata1)


  #Initial values for chain 1
  lm.m <- lm(formula,Pdata)
  betaB0=as.vector(lm.m$coefficients)
  betaB1=as.vector(betaB0)
  tau1 <-1/summary(lm.m)$sigma

  #initial for other chains
  betaB2=betaB1-2; tau2 <-tau1+20
  betaB3=betaB1+5;tau3 <-tau1+3


  #list of jags initial values
  UNCjags.inits <- function(){
    list("beta0"=betaB1[1], "tau"=tau1)#Initial values for chain1
    list("beta0"=betaB2[1],"tau"=tau2)#Initial values for chain2
    list("beta0"=betaB3[1],"tau"=tau3)#Initial values for chain3
  }


  ## 1. Uncondtional model
  #========================
  UNCdata <- list(N=N, post=Pdata1$post)
  jags.UNCparams <- c("sigma")

  # model
  filenames_OLS_UNC <-  system.file("jags", "OLS_UNC.txt", package = "eefAnalytics")#file.path("inst/jags/OLS_UNC.txt")
  cat(paste("
  model {
    # Likelihood
    for (i in 1:N) {
      post[i] ~ dnorm(mu[i], tau)
      mu[i] <- beta0
    }

    # Priors
    beta0 ~ dnorm(0, 0.0001) # Overall intercept
    tau ~ dgamma(0.001, 0.001) # Precision for within-cluster variation (residual variance)
    sigma <- 1 / sqrt(tau) # Within-cluster standard deviation (Residual SD)

  }
"), file = filenames_OLS_UNC)


  # Summarise UNCONDITIONAL JAGS output
  UNC.ols<-do.call(jags, c(list(model.file=filenames_OLS_UNC,
                                data = UNCdata,
                                n.iter= nsim,
                                n.burnin = nsim/2,
                                inits=UNCjags.inits,
                                parameters.to.save=jags.UNCparams,
                                n.thin = 10,  ...), uncopt))

  UNC.ols.upd <- autojags(UNC.ols)
  UNC.sigma <-MCMCsummary(UNC.ols.upd,round = 2)["sigma",c("mean")]




  ## 2. Condtional model
  #=======================
  for(i in 1:length(var)){assign(var[i],Pdata1[,i])}
  p <- length(var) #number of covariate + intercept
  threshold1 <- stats::setNames(1:length(threshold), threshold)
  Post.T0 <- intersect(c(intervention,paste0(intervention, unique(Pdata[, intervention]))),var )
  Post.T <- stats::setNames(which(var %in% Post.T0),Post.T0)#position of intervention



  data <- as.list(Pdata1)
  data[c("N", "p", "UNC.sigma")]<- c(N,p,UNC.sigma) #jags data
  data[["Post.T"]] <- Post.T
  data[["threshold1"]] <- threshold1
  data[["threshold"]] <- as.numeric(names(threshold1))




  # COND Jags parameters to monitor
  jags.params <- c("beta","COND.sigma","UNC.sigma","COND.ES","UNC.ES","COND.ProbES","UNC.ProbES")


  # 2. COND Jags model -----------
  filenames_OLS_CON <- system.file("jags", "SRT.txt", package = "eefAnalytics")#file.path("inst/jags/SRT.txt")
  cat(paste("

              model{
              for(i in 1:N){
              post[i] ~ dnorm(mu[i],tau)
              mu[i] <-", paste0("beta[1]+",paste0("beta","[",1:length(var[-1])+1,"]","*",var[-1],"[i]",collapse ="+")),"
              }

              for(k in 1:p){beta[k]~dnorm(0.0,1.0E-06)}

              #PRIORS
              tau~dgamma(0.001,0.0001)
              sigma<-1/tau #CONVERT PRECISION TO VARIANCES
              COND.sigma <- sigma


              # EFFECT SIZE
                for(pp  in round(Post.T)){
                COND.d[pp] <- beta[pp]/sqrt(COND.sigma)#Cohen'd
                UNC.d[pp]  <- beta[pp]/sqrt(UNC.sigma)#Cohen'd
                COND.ES[pp]<- COND.d[pp]* (1- (3/(4*(N-2)-1)))#hedges'g
                UNC.ES[pp] <- UNC.d[pp] * (1- (3/(4*(N-2)-1)))#hedges' g

                #Posterior probabilities
                for(thd  in round(threshold1)){
                COND.ProbES[pp, thd]<- step(COND.ES[pp] - threshold[thd] )#conditional
                UNC.ProbES[pp, thd] <- step(UNC.ES[pp] - threshold[thd] )#unconditional
                }
                }



              }
              ")
      ,file=filenames_OLS_CON
  )

  #list of jags initial values
  jags.inits <- function(){
    list("beta"=betaB1, "tau"=tau1) #Initial values for chain1
    list("beta"=betaB2,"tau"=tau2)#Initial values for chain2
    list("beta"=betaB3,"tau"=tau3)#Initial values for chain3
  }


  #Model and Summarise CONDITIONAL JAGS output
  jag.ols <- do.call(jags, c(list(model.file=filenames_OLS_CON,
                                  data = data,
                                  n.iter= nsim,
                                  n.burnin = nsim/2,
                                  inits=jags.inits,
                                  parameters.to.save=jags.params,
                                  n.thin = 10, ...), condopt))

  jag.ols.upd <- autojags(jag.ols)
  jag.ols.sum0 <- MCMCsummary(jag.ols.upd,round = digits, probs= c(alpha/2, 1-alpha/2))
  jag.ols.sum  <- jag.ols.sum0[,c("mean", paste0(alpha/2*100,"%" ), paste0((1-alpha/2)*100,"%" ))]
  row.names(jag.ols.sum)[row.names(jag.ols.sum)%in% paste0("beta[",1:p,"]")] <- c("Intercept",var[-1])


  # Posterior probabilities
  #-------------------------
  ProbES <-  lapply(Post.T, function(i) {
    ProbESP <-  data.frame(cbind(Cond=jag.ols.sum[grep(paste0("COND.ProbES\\[",i),row.names(jag.ols.sum)),"mean"],
                                 Uncond =jag.ols.sum[grep(paste0("UNC.ProbES\\[",i),row.names(jag.ols.sum)),"mean"]))
    row.names(ProbESP) <- paste0("P(ES>",threshold,")")
    round(ProbESP,digits)
  })


  # Effect sizes
  #-----------------
  UNCOND_ols_fn <-  function(trt, COND="COND", jag.ols.sum){
    ES <-  jag.ols.sum[grep(paste0(COND,".ES\\[",trt),row.names(jag.ols.sum)),]
    colnames(ES) <- c("Estimate",paste0((1-alpha)*100, c("% LB","% UB")))
    ES
  }




  #Unconditional output
  #---------------------
  colnames(jag.ols.sum) <- c("Estimate",paste0((1-alpha)*100, c("% LB","% UB")))
  Unconditional= list(ES=t(sapply(Post.T, function(i) UNCOND_ols_fn(trt=i,COND="UNC", jag.ols.sum=jag.ols.sum))),
                      sigma2=jag.ols.sum["UNC.sigma", "Estimate"],
                      ProbES=lapply(ProbES, function(x) x[, "Uncond",drop=FALSE]))


  #Final all output
  #-----------------
  ObjectOUT <- list(
    Beta=jag.ols.sum[c("Intercept",var[-1]),],
    ES= t(sapply(Post.T, function(i) ES1=UNCOND_ols_fn(trt=i,COND="COND", jag.ols.sum=jag.ols.sum))),
    sigma2=jag.ols.sum["COND.sigma", "Estimate"],
    ProbES=lapply(ProbES, function(x) x[, "Cond",drop=FALSE]),
    Model=list(jag.ols=jag.ols,jag.ols.upd=jag.ols.upd, jag.ols.sum=jag.ols.sum),
    Unconditional=Unconditional)


  return(ObjectOUT)

}

############################################# END #######################################################


