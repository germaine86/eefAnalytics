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
#' @param nsim number of MCMC iterations per chain. Default is 2000.
#' @param data Data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and credible intervals for the variables specified in the model. Use \code{summary.eefAnalytics} to get Rhat and effective sample size for each estimate.
#' \item \code{ES}: Conditional Hedges' g effect size and its 95% credible intervals.
#' \item \code{sigma}: Residual variance.
#' \item \code{ProbES}: A matrix of Bayesian posterior probabilities such that the observed effect size is greater than or equal to a pre-specified threshold(s).
#' \item \code{Unconditional}: A list of unconditional effect sizes, sigma2 and ProbES obtained based on residual variance from the unconditional model (model with only the intercept as a fixed effect).
#'  }
#' @example inst/examples/srtBExample.R
srtBayes <- function(formula,intervention,nsim=10000,data)UseMethod("srtBayes")

#' @export
srtBayes.default <- function(formula,intervention,nsim=10000,data){stop("No correct formula input given.")}

#' @export
srtBayes.formula <- function(formula,intervention,nsim=10000,data){

  #if(nsim<2000){stop("nsim must be at least 2000")}

  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),intervention))])
  tmp3 <- which(colnames(data)==intervention)
  mf <- model.frame(formula=formula, data=data)
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  posttest <- model.response(mf)
  intervention <- intervention
  trt <- data[,which(colnames(data)==intervention)]
  btp  <- which(tmp %in% paste0(intervention, unique(data[, intervention]) ))
  LHS.formula <- "post"
  RHS.formula <-paste(tmp[-1], collapse = "+")
  new.formula0 <- as.formula(post~1)
  new.formula  <- as.formula(paste0(LHS.formula,"~",RHS.formula ))
  new.data <- data.frame(post=posttest, fixedDesignMatrix[,-1])

  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}


  BayesOutput <- OLS.function(data=data, formula=formula,intervention=intervention, nsim=nsim)
  output  <- errantSummary(bayesObject=BayesOutput)

  output$Method <- "LM"
  output$Function <- "srtBayes"
  class(output) <- "eefAnalytics"
  return(output)
}



################################################################################################################################
################################################################################################################################
OLS.function <- function(data, formula, intervention,nsim){

  #### data preparation
  Pdata <- na.omit(data[,all.vars(formula)])
  outcome <- all.vars(formula)[1] ## Y: Posttest
  dummies<- data.frame(model.matrix(formula, data=Pdata)) ## X0-Intercept, X1-Prettest, X2-Intervention
  Pdata0 <- na.omit(dummies[,!(names(dummies) %in% "X.Intercept.")])
  Pdata1<-na.omit(cbind(post=Pdata[,outcome],Pdata0)) #all covariates and school dummies

  # Jags data
  var<-names(Pdata1)
  N <- nrow(Pdata1)

  #Initial values preparation --------------------------
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
  #************************
  UNCdata <- list(N=N, post=Pdata1$post)
  jags.UNCparams <- c("sigma")
  # 1. UNC Jags model -----------
  filenames_OLS_UNC <- file.path("inst/jags/OLS_UNC.txt")
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

  # 1. UNC Jags model -----------
  #### Summarise UNCONDITIONAL JAGS output ####
  UNC.ols<-jags(model.file=filenames_OLS_UNC, data = UNCdata, n.iter= nsim, n.burnin = nsim/2, inits=UNCjags.inits,n.thin = 10, parameters.to.save=jags.UNCparams)
  UNC.ols.upd <- autojags(UNC.ols)
  UNC.sigma <-MCMCsummary(UNC.ols.upd,round = 2)["sigma",c("mean")]


  ## 2. Condtional model
  #************************
  p <- length(var)# number of Xs+1
  for(i in 1:length(var)){assign(var[i],Pdata1[,i])}
  data <- as.list(Pdata1);data[c("N", "p", "UNC.sigma")]<- c(N,p,UNC.sigma) #jags data
  Post.T <- match(intervention, var) #position of intervention

  # COND Jags parameters to monitor
  jags.params <- c("beta","COND.sigma","UNC.sigma","COND.ES","COND.g","UNC.ES","UNC.g")

  # 2. COND Jags model -----------
  filenames_OLS_CON <- file.path("inst/jags/SRT.txt")
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
              COND.ES <- ",paste0("beta[",Post.T,"]"),"/sqrt(sigma)
              UNC.ES <- ",paste0("beta[",Post.T,"]"),"/sqrt(UNC.sigma)
              COND.g <- COND.ES* (1-  (3/(4*(N-2)-1)))
              UNC.g <- UNC.ES* (1-  (3/(4*(N-2)-1)))

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

  #### Summarise CONDITIONAL JAGS output ####
  jag.ols<-jags(model.file=filenames_OLS_CON, data = data, n.iter= nsim, n.burnin = nsim/2, inits=jags.inits,n.thin = 10, parameters.to.save=jags.params)
  jag.ols.upd <- autojags(jag.ols)
  jag.ols.sum <-MCMCsummary(jag.ols.upd,round = 2)[,c("mean","2.5%" , "97.5%")]
  row.names(jag.ols.sum)[row.names(jag.ols.sum)%in% paste0("beta[",1:p,"]")] <- c("Intercept",var[-1])

  # 3. Jags Output -----------
  #### 3.1 Final Summarised output ####
  Unconditional= list(ES=data.frame(jag.ols.sum["UNC.ES",]), sigma=jag.ols.sum["COND.sigma",1])

  ObjectOUT <- list(jags.model=list(jag.ols=jag.ols,jag.ols.upd=jag.ols.upd, jag.ols.sum=jag.ols.sum),
                    beta=jag.ols.sum[c("Intercept",var[-1]),],
                    ES=data.frame(jag.ols.sum["COND.ES",]),
                    sigma=jag.ols.sum["COND.sigma",1],
                    unconditional=Unconditional)
  #return(ObjectOUT)
  #### 3.2 mcmc output ####
  mcmc_sample <- as.mcmc(jag.ols.upd)

  # Convert mcmc_sample to a matrix for each chain
  # If mcmc_sample is of class mcmc.list (multiple chains), you can use as.matrix
  mcmc_matrix <- as.matrix(mcmc_sample)
  # Extract the number of chains
  num_chains <- nchain(mcmc_sample)
  # Create a data frame with the three chains
  # Add a column indicating which chain the row belongs to
  mcmc_df <- data.frame(mcmc_matrix)
  mcmc_df$Chain <- rep(1:num_chains, each = nrow(mcmc_matrix) / num_chains)

  ###### 3.2.1 Conditional mcmc output ####
  # Beta estimates ---- from mcmc samples
  mcmc_Beta <- mcmc_df[, grep("beta", colnames(mcmc_df))]
  colnames(mcmc_Beta)<- colnames(dummies)

  # Effect Size estimates ---- from mcmc samples
  mcmc_ES <- mcmc_df[, grep("COND.ES", colnames(mcmc_df))]

  # sigma estimates ---- from mcmc samples
  mcmc_sigma <- mcmc_df[, grep("COND.sigma", colnames(mcmc_df))]


  ###### 3.2.2 UNConditional mcmc output ####
  # Effect Size estimates ---- from mcmc samples
  mcmc_UNC_ES <- mcmc_df[, grep("UNC.ES", colnames(mcmc_df))]

  # sigma estimates ---- from mcmc samples
  mcmc_UNC_sigma <- mcmc_df[, grep("UNC.sigma", colnames(mcmc_df))]


  unconditional= list(ES=mcmc_UNC_ES,
                      sigma=mcmc_UNC_sigma)

  ###### 3.2.3 Summarised mcmc output ####
  Output <- list(sigma=mcmc_sigma,
                 Beta=mcmc_Beta,
                 ES=mcmc_ES,
                 Unconditional = unconditional)
  return(Output)
}

#### I. summarise covariance parameters - internal #####
sigmaSummary <- function(bayesObject){
  covParm <- bayesObject$sigma ### MCMC
  covParm2 <- mean(covParm)
  covParm3 <- quantile(covParm,prob=c(0.025,0.975))
  covParm4 <- data.frame("Variance" = covParm2,
                         "95% LB" = covParm3[1],
                         "95% UB" = covParm3[2])
  colnames(covParm4) <- c("Variance","95% LB","95% UB")
  rownames(covParm4) <- c("estimate")
  return(covParm4)
}



#### II. summarise beta parameters - internal #####
betaSummary <- function(bayesObject){

  betaParm <- bayesObject$Beta
  betaParm2 <- colMeans(betaParm)
  betaParm3 <- t(apply(betaParm,2,function(x)quantile(x,prob=c(0.025,0.975))))
  betaParm4 <- data.frame(cbind(betaParm2,betaParm3))
  colnames(betaParm4) <- c("Estimate","95% LB","95% UB")
  rownames(betaParm4) <- colnames(betaParm)
  return(betaParm4)
}


#### III. summarise Effect Sizes - internal #####
esSummary <- function(bayesObject){

  esParm <- bayesObject$ES
  esParm2 <- mean(esParm)
  esParm3 <- quantile(esParm,prob=c(0.025,0.975))
  esParm4 <- data.frame("Estimate" = esParm2,
                        "95% LB" = esParm3[1],
                        "95% UB" = esParm3[2])
  colnames(esParm4) <- c("Estimate","95% LB","95% UB")
  rownames(esParm4) <- "ES"
  return(esParm4)
}

#### IV. summarise all Bayesian parameters - internal #####
errantSummary <- function(bayesObject){
  sigmaValues <- sigmaSummary(bayesObject=bayesObject)
  betaValues <- betaSummary(bayesObject=bayesObject)
  esValues <- esSummary(bayesObject)

  unconditional= list(ES=round(esSummary(bayesObject$Unconditional),2),
                      sigma=round(sigmaSummary(bayesObject$Unconditional),2))

  output <- list(Beta=round(betaValues,2),
                 sigma= round(sigmaValues,2),
                 ES=round(esValues,2),
                 Unconditional = unconditional)
}

############################################# END #######################################################


