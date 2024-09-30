

###################################################################################################
############# Bayesian Multilevel Analysis of Cluster Randomised Education Trials ###############
##################################################################################################

#' Bayesian analysis of Cluster Randomised Education Trials (CRT) using Vague Priors.
#'
#' \code{crtBayes} performs Bayesian multilevel analysis of cluster randomised education trials, utilising vague priors
#' and JAGS language to fit the model. It assumes hierarchical clustering, such as students within schools, and estimates
#' treatment effects while accounting for this structure.
#'
#' The function provides posterior estimates for fixed effects (predictors) and random effects (clustering) under a Bayesian framework.
#' Effect sizes are computed using Hedges' g, and variance components are decomposed into between-cluster and within-cluster variances.
#'
#' @export
#' @param formula The model to be analysed. It should be of the form y ~ x1 + x2 + ..., where y is the outcome variable and Xs are the predictors.
#' @param random A string specifying the "clustering variable" (e.g., schools or sites) as found in the dataset.
#' @param intervention A string specifying the "intervention variable" as it appears in the formula.
#' @param nsim Number of MCMC iterations to be performed. A minimum of 10,000 is recommended to ensure convergence.
#' @param data A data frame containing the variables referenced in the formula, including predictors, the clustering variable, and the intervention.
#' @return S3 object; a list consisting of:
#' \itemize{
#'   \item \code{Beta}: Estimates and credible intervals for the predictors specified in the model (posterior distributions).
#'   \item \code{ES}: Hedges' g effect size for the intervention(s). If bootstrapping is not used, 95% credible intervals are computed based on MCMC sampling.
#'   \item \code{covParm}: Variance components broken down into between-cluster variance (e.g., between schools), within-cluster variance (e.g., within pupils), and intra-cluster correlation (ICC)..
#'   \item \code{ProbES}: A matrix showing the probability of observing an effect size larger than various thresholds (0, 0.05, 0.10, ...).
#'   \item \code{Unconditional}: A list containing the unconditional effect size and variance decomposition.
#' }
#'
#' @example inst/examples/crtBExample.R
crtBayes <- function(formula,random,intervention,nsim=10000,data)UseMethod("crtBayes")

#' @export
crtBayes.default <- function(formula,random,intervention,nsim=10000,data){stop("No correct formula input given.")}

#' @export
crtBayes.formula <- function(formula,random,intervention,nsim=10000,data){

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


  LHS.formula <- "post"
  RHS.formula <-paste(tmp[-1], collapse = "+")
  RHS.cluster <-"+(1|cluster)"
  new.formula <-  as.formula(paste0(LHS.formula,"~",RHS.formula,RHS.cluster))
  new.formula0 <- as.formula(paste0(LHS.formula,"~",RHS.cluster))


  nsim=nsim
  if(length(tmp2)!= 1){stop("Cluster variable misspecified")}
  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}
  if(nsim < 2000){stop("nsim >= 10000 is recommended")}

  BayesOutput <- CRT.function(data=data, formula=formula,random=random, intervention=intervention, nsim=nsim)
  output  <- errantSummary(bayesObject=BayesOutput,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention)

  output$Method <- "MLM"
  output$Function <- "crtBayes"
  class(output) <- "eefAnalytics"
  invisible(output)
}

#############################################################################################################################
################################################################################################################################
#### I. perform Bayesian multilevel linear modeling for Multi-Stage Trials (MST) using JAGS language - internal #####
CRT.function <- function(data, formula,random, intervention, nsim){

  #### load required packages ###
  requireNamespace("R2jags", quietly = TRUE) || stop("Please install the 'R2jags' package.")
  require(R2jags)
  requireNamespace("lme4", quietly = TRUE) || stop("Please install the 'lme4' package.")
  require(lme4)
  requireNamespace("MCMCvis", quietly = TRUE) || stop("Please install the 'MCMCvis' package.")
  require(MCMCvis)
  requireNamespace("coda", quietly = TRUE) || stop("Please install the 'coda' package.")
  require(coda)

  #### check if it is a CRT design
  Pdata <- na.omit(data[,c(all.vars(formula),random)])
  chk <- sum(rowSums(table(Pdata[,random],Pdata[,intervention])!=0)>1)
  if(chk !=0){stop("This is not a CRT design, try 'MST.function' instead")}

  #### data preparation
  outcome <- all.vars(formula)[1] ## Y: Posttest
  dummies<- data.frame(model.matrix(formula, data=Pdata)) ## X0-Intercept, X1-Prettest, X2-Intervention
  Pdata0 <- na.omit(dummies[,!(names(dummies) %in% "X.Intercept.")])
  Pdata1<-na.omit(cbind(post=Pdata[,outcome],Pdata0)) #all covariates and school dummies
  Pdata1[,random] <- as.numeric(as.factor(Pdata[,random]))

  tmp2 <- which(colnames(Pdata)==random)
  cluster2 = Pdata[,tmp2]
  mf <- model.frame(formula=formula, data=Pdata)
  mf <- mf[order(cluster2),]
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=Pdata)))
  tmp <- colnames(fixedDesignMatrix)
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp

  # Jags data
  var<-names(Pdata1)
  N <- nrow(Pdata1)
  M <- length(unique(Pdata1[,random]))


  #Initial values preparation --------------------------
  #Initial values for chain 1
  formula <- update(formula,paste0("~ .+ (1|",random,")"))
  lm.m1 <- lmer(formula ,Pdata)
  betaB=lm.m1@beta
  vsigma<-as.data.frame(VarCorr(lm.m1));tau <-1/(vsigma[1,4]+0.1);tau.s <- 1/(vsigma[1,4]+0.1) #one was added to school and school:t to avoid zero in denominator
  #initial for other chains
  betaB2=betaB+15;tau2 <-tau+2;tau.s2 <-tau.s+4
  betaB3=betaB-2;tau3 <-tau+20;tau.s3 <-tau.s+10
  #list of jags initial values
  UNCjags.inits <- function(){
    list("beta0"=betaB[1], "tau"=tau, "tau_u"=tau.s, "u"=rnorm(0,1))#Initial values for chain1
    list("beta0"=betaB2[1],"tau"=tau2,"tau_u"=tau.s2,"u"=rnorm(0,1))#Initial values for chain2
    list("beta0"=betaB3[1],"tau"=tau3,"tau_u"=tau.s3, "u"=rnorm(0,1))#Initial values for chain3
  }


  ## 1. Uncondtional model
  #************************
  UNCdata <- list(N=N,M=M, school=Pdata1[,random], post=Pdata1$post)
  jags.UNCparams <- c("sigma","sigma.tt","icc")
  filenames_MLM_UNC <- file.path("/Users/qingzhang/Desktop/Durham/WP8/eefAnalytics/inst/jags/MLM_UNC.txt")
  # 1. UNC Jags model -----------
  #### Summarise UNCONDITIONAL JAGS output ####
  UNC.ols<-jags(model.file=filenames_MLM_UNC, data = UNCdata, n.iter= nsim, n.burnin = nsim/2, inits=UNCjags.inits,n.thin = 10, parameters.to.save=jags.UNCparams)
  UNC.ols.upd <- autojags(UNC.ols)
  UNC.sigma <-MCMCsummary(UNC.ols.upd,round = 2)["sigma",c("mean")]
  UNC.sigma.tt <-MCMCsummary(UNC.ols.upd,round = 2)["sigma.tt",c("mean")]
  UNC.icc<-MCMCsummary(UNC.ols.upd,round = 2)["icc",c("mean")]


  ## 2. Conditional model
  #************************
  for(i in 1:length(var)){assign(var[i],Pdata1[,i])}
  var1 <- var[!var%in% c("post"  ,random)] #list of covariates
  p <- length(var1)+1 #number of covariate + intercept
  data <- as.list(Pdata1);data[c("N","M", "p","UNC.sigma","UNC.sigma.tt","UNC.icc")]<- c(N,M,p,UNC.sigma,UNC.sigma.tt,UNC.icc) #jags data
  Post.T <- match(intervention, var) #position of intervention

  # COND Jags parameters to monitor
  jags.params <- c("COND.ES.Within","COND.ES.Total","COND.g.with","COND.g.Total","COND.sigma.Within","COND.sigma.Total","COND.icc",
                   "UNC.ES.Within","UNC.ES.Total","UNC.g.with","UNC.g.Total","UNC.sigma.Within","UNC.sigma.Total","beta","UNC.ICC")

  # 2. COND Jags model -----------
  filenames_CRT <- file.path("/Users/qingzhang/Desktop/Durham/WP8/eefAnalytics/inst/jags/CRT.txt")
  cat(paste("
                model{
                for(i in 1:N){
                post[i] ~ dnorm(mu[i],tau)
                mu[i] <-", paste0("beta[1]+",paste0("beta","[",1:length(var1)+1,"]","*",var1,"[i]",collapse ="+"),"+b1[",random,"[i]]"),"
                }

                for(j in 1:M){
                b1[j]~dnorm(0.0,tau.b1)
                }

                tau.b1~dgamma(0.001,0.0001)
                sigma.b1<-1/tau.b1
                tau~dgamma(0.001,0.0001)
                sigma<-1/tau

                for(k in 1:p){beta[k]~dnorm(0.0,1.0E-06)}

                # ICC and TOTAL VARIANCE
                sigma.Total <-sigma + sigma.b1
                COND.icc <- sigma.b1 * pow(sigma.Total ,-1)
                UNC.ICC <- UNC.icc

                #sigmas
                COND.sigma.Total <- sigma.Total
                COND.sigma.Within <- sigma.b1
                UNC.sigma.Total <- UNC.sigma.tt
                UNC.sigma.Within <- UNC.sigma

                # EFFECT SIZE
                COND.ES.Within<- ",paste0("beta[",Post.T,"]"),"/sqrt(COND.sigma.Within)#conditional
                COND.ES.Total <- ",paste0("beta[",Post.T,"]"),"/sqrt(COND.sigma.Total)#conditional
                UNC.ES.Within<- ",paste0("beta[",Post.T,"]"),"/sqrt(UNC.sigma.Within)#unconditional
                UNC.ES.Total <- ",paste0("beta[",Post.T,"]"),"/sqrt(UNC.sigma.Total)#unconditional

                COND.g.with<- COND.ES.Within* (1- (3/(4*(N-2)-1)))#conditional
                COND.g.Total <-  COND.ES.Total* (1- (3/(4*(N-2)-1)))#conditional
                UNC.g.with<- UNC.ES.Within* (1- (3/(4*(N-2)-1)))#unconditional
                UNC.g.Total <-  UNC.ES.Total* (1- (3/(4*(N-2)-1)))#unconditional
                }
                ")
      ,file=filenames_CRT
  )

  #list of jags initial values
  jags.inits <- function(){
    list("beta"=betaB, "tau"=tau, "tau.b1"=tau.s, "b1"=rnorm(0,1))#Initial values for chain1
    list("beta"=betaB2,"tau"=tau2,"tau.b1"=tau.s2,"b1"=rnorm(0,1))#Initial values for chain2
    list("beta"=betaB3,"tau"=tau3,"tau.b1"=tau.s3, "b1"=rnorm(0,1))#Initial values for chain3
  }

  #### Summarise CONDITIONAL JAGS output ####
  jag.model<-jags(model.file=filenames_CRT, data = data, n.iter= nsim, n.burnin = nsim/2, inits=jags.inits,n.thin = 10, parameters.to.save=jags.params)
  jag.model.upd <- autojags(jag.model)
  jag.model.sum <-MCMCsummary(jag.model.upd,round = 2)[,c("mean","2.5%" , "97.5%")]
  row.names(jag.model.sum)[row.names(jag.model.sum)%in% paste0("beta[",1:p,"]")] <- c("Intercept",var1)

  # 3. Jags Output -----------
  #### 3.1 Final Summarised output ####
  Unconditional= list(ES=data.frame(rbind(Within=jag.model.sum["UNC.ES.Within",],Total=jag.model.sum["UNC.ES.Total",])),
                      covParm=c(Within=jag.model.sum["UNC.sigma.Within",1],Total=jag.model.sum["UNC.sigma.Total",1], ICC= jag.model.sum["UNC.ICC",1]))

  ObjectOUT <- list(jags.model=list(jag.model=jag.model,jag.model.upd=jag.model.upd, jag.model.sum=jag.model.sum),
                    beta=jag.model.sum[c("Intercept",var1),],
                    ES=data.frame(rbind(Within=jag.model.sum["COND.ES.Within",],Total=jag.model.sum["COND.ES.Total",])),
                    covParm=c(Within=jag.model.sum["COND.sigma.Within",1],Total=jag.model.sum["COND.sigma.Total",1],ICC= jag.model.sum["COND.icc",1]),
                    unconditional=Unconditional)
  #return(ObjectOUT)


  #### 3.2 mcmc output ####
  mcmc_sample <- as.mcmc(jag.model.upd)

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
  # Find column indices for columns containing 'COND.sigma' and 'COND.icc'
  cond_sigma_cols <- grep("COND.sigma", colnames(mcmc_df))
  cond_icc_col <- grep("COND.icc", colnames(mcmc_df))
  # Covariance estimates ---- from mcmc samples
  mcmc_covParm <- mcmc_df[, c(cond_sigma_cols, cond_icc_col)]
  mcmc_covParm <- mcmc_covParm[, c("COND.sigma.Within", "COND.sigma.Total", "COND.icc")]
  colnames(mcmc_covParm) <- c("Within", "Total","ICC")

  # Beta estimates ---- from mcmc samples
  mcmc_Beta <- mcmc_df[, grep("beta", colnames(mcmc_df))]
  colnames(mcmc_Beta)<- colnames(fixedDesignMatrix)

  # Effect Size estimates ---- from mcmc samples
  mcmc_ES <- mcmc_df[, grep("COND.ES", colnames(mcmc_df))]
  mcmc_ES <- mcmc_ES[, c("COND.ES.Within", "COND.ES.Total")]
  colnames(mcmc_ES) <- c("Within","Total")

  ###### 3.2.2 UNConditional mcmc output ####
  # Find column indices for columns containing 'UNC.sigma' and 'UNC.ICC'
  UNC_sigma_cols <- grep("UNC.sigma", colnames(mcmc_df))
  UNC_icc_col <- grep("UNC.ICC", colnames(mcmc_df))
  # Covariance estimates ---- from mcmc samples
  mcmc_UNC_covParm <- mcmc_df[, c(UNC_sigma_cols, UNC_icc_col)]
  mcmc_UNC_covParm <- mcmc_UNC_covParm[, c("UNC.sigma.Within", "UNC.sigma.Total", "UNC.ICC")]
  colnames(mcmc_UNC_covParm) <- c("Within","Total","ICC")

  # Effect Size estimates ---- from mcmc samples
  mcmc_UNC_ES <- mcmc_df[, grep("UNC.ES", colnames(mcmc_df))]
  mcmc_UNC_ES <- mcmc_UNC_ES[, c("UNC.ES.Within", "UNC.ES.Total")]
  colnames(mcmc_UNC_ES) <- c("Within","Total")

  unconditional= list(ES=mcmc_UNC_ES,
                      covParm=mcmc_UNC_covParm)

  ###### 3.2.3 Summarised mcmc output ####
  Output <- list(covParm=mcmc_covParm,
                 Beta=mcmc_Beta,
                 ES=mcmc_ES,
                 Unconditional = unconditional)
  return(Output)
}





#### II. summarise covariance parameters - internal #####
covSummary <- function(bayesObject){

  covParm <- bayesObject$covParm ### MCMC
  covParm2 <- colMeans(covParm) # remove all roundings that appear in the middle of functions
  covParm3 <- t(apply(covParm,2,function(x)quantile(x,prob=c(0.025,0.975))))
  covParm4 <- data.frame(cbind(covParm2,covParm3))
  covParm5 <- sqrt(covParm4)
  colnames(covParm4) <- c("Variance","95% LB","95% UB")
  rownames(covParm4) <- c("Pupils","Total","ICC")
  covParm5 <- covParm4[,1]
  return(covParm5)
}

UNCcovSummary <- function(bayesObject){

  covParm <- bayesObject$Unconditional$covParm ### MCMC
  covParm2 <- colMeans(covParm) # remove all roundings that appear in the middle of functions
  covParm3 <- t(apply(covParm,2,function(x)quantile(x,prob=c(0.025,0.975))))
  covParm4 <- data.frame(cbind(covParm2,covParm3))
  covParm5 <- sqrt(covParm4)
  colnames(covParm4) <- c("Variance","95% LB","95% UB")
  rownames(covParm4) <- c("Pupils","Total","ICC")
  covParm5 <- covParm4[,1]
  return(covParm5)
}


#### III. summarise beta parameters - internal #####
betaSummary <- function(bayesObject){

  betaParm <- bayesObject$Beta
  betaParm2 <- colMeans(betaParm)
  betaParm3 <- t(apply(betaParm,2,function(x)quantile(x,prob=c(0.025,0.975))))
  betaParm4 <- data.frame(cbind(betaParm2,betaParm3))
  colnames(betaParm4) <- c("Estimate","95% LB","95% UB")
  rownames(betaParm4) <- colnames(betaParm)
  return(betaParm4)
}


#### IV. summarise Effect Sizes - internal #####
esSummary <- function(bayesObject,fixedDesignMatrix,intervention){

  esParm <- bayesObject$ES
  esParm2 <- colMeans(esParm)
  esParm3 <- t(apply(esParm,2,function(x)quantile(x,prob=c(0.025,0.975))))
  esParm4 <- data.frame(cbind(esParm2,esParm3))
  colnames(esParm4) <- c("Estimate","95% LB","95% UB")
  btp <- nrow(esParm4)
  if(btp <=3){return(round(esParm4,2))}

  if(btp >3){
    btp2 <- seq(3,btp,3)
    btp3 <- seq(1,btp,3)
    ouptut <- list()
    for(i in 1:length(btp2)){
      ouptut[[i]] <- round(esParm4[btp3[i]:btp2[i],],2)
    }
    rname <- colnames(fixedDesignMatrix)[substring(colnames(fixedDesignMatrix),1,nchar(intervention))==intervention]
    names(ouptut)<- rname
    return(ouptut)
  }

}


#### V. summarise minimum expected effect size - internal #####
esProb <- function(bayesObject,esOutput){
  es <- c(0,0.05,seq(0.1,1,0.1))
  esParm <- bayesObject$ES
  esParm2 <- sapply(es,function(x)colMeans(esParm>=x))
  esParm4 <- data.frame(cbind(ES=es,t(esParm2)))
  rownames(esParm4) <- NULL
  return(esParm4)

}


#### VI. summarise all Bayesian parameters - internal #####
errantSummary <- function(bayesObject,fixedDesignMatrix,intervention){
  covValues <- covSummary(bayesObject=bayesObject)
  covValues <- data.frame(covValues)
  row.names(covValues) <- c("Pupils","Total","ICC")
  covValues <- t(covValues)
  row.names(covValues ) <- NULL
  betaValues <- betaSummary(bayesObject=bayesObject)
  esValues <- esSummary(bayesObject,fixedDesignMatrix,intervention)
  es.prob <- esProb(bayesObject=bayesObject,esOutput=esValues)

  UNCcovValues <- UNCcovSummary(bayesObject=bayesObject)
  UNCcovValues <- data.frame(UNCcovValues)
  row.names(UNCcovValues) <- c("Pupils","Total","ICC")
  UNCcovValues <- t(UNCcovValues)
  row.names(UNCcovValues) <- NULL
  unconditional= list(ES=round(esSummary(bayesObject$Unconditional,fixedDesignMatrix,intervention),2),
                      covParm=round(UNCcovValues,2))

  output <- list(Beta=round(betaValues,2),
                 covParm= round(covValues,2),
                 ES=esValues,
                 ProbES=round(data.frame(es.prob),2),
                 Unconditional = unconditional)
}

############################################# END #######################################################

