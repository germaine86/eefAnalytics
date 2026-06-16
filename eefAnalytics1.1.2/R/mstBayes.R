
###################################################################################################
############# Bayesian Multilevel Analysis of Multisite Randomised Education Trials ###############
##################################################################################################

#' Bayesian analysis of Multisite Randomised Education Trials (MST) using Vague Priors.
#'
#' \code{mstBayes} performs Bayesian multilevel analysis of multisite randomised education trials, utilising vague priors
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
#' @param nSim Number of MCMC iterations to be performed. A minimum of 10,000 is recommended to ensure convergence.
#' @param data A data frame containing the variables referenced in the formula, including predictors, the clustering variable, and the intervention.
#' @return S3 object; a list consisting of:
#' \itemize{
#'   \item \code{Beta}: Estimates and credible intervals for the predictors specified in the model (posterior distributions).
#'   \item \code{ES}: Hedges' g effect size for the intervention(s). If bootstrapping is not used, 95% credible intervals are computed based on MCMC sampling.
#'   \item \code{covParm}: Variance components broken down into between-cluster variance (e.g., between schools), within-cluster variance (e.g., within pupils), and intra-cluster correlation (ICC)..
#'   \item \code{randomEffects}: Posterior estimates of random intercepts for each cluster (e.g., schools).
#'   \item \code{ProbES}: A matrix showing the probability of observing an effect size larger than various thresholds (0, 0.05, 0.10, ...).
#'   \item \code{Unconditional}: A list containing the unconditional effect size and variance decomposition.
#' }
#'
#' @example inst/examples/mstBExample.R
#'
mstBayes <- function(formula,random,intervention,nSim,data)UseMethod("mstBayes")

#' @export
mstBayes.default <- function(formula,random,intervention,nSim=nSim,data){stop("No correct formula input given.")}

#' @export
mstBayes.formula <- function(formula,random,intervention,nSim=nSim,data){
  data1 <- data[order(data[,which(colnames(data)==random)],data[,which(colnames(data)==intervention)]),]
  tmp3 <- which(colnames(data1)==intervention)
  data1[,tmp3] <- as.factor(data1[,tmp3])
  tmp2 <- which(colnames(data1)==random)
  cluster2 = data1[,tmp2]

  mf <- model.frame(formula=formula, data=data1)
  mf <- mf[order(cluster2),]
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data1)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  intervention <- intervention
  trt <- data1[,which(colnames(data1)==intervention)]
  tmp2 <- which(colnames(data1)==random)
  nsim=nSim

  if(length(tmp2)!= 1){stop("Cluster variable misspecified")}
  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}

  if(nsim < 10000){stop("nsim >= 10000 is recommended")}


  BayesOutput <- MST.function(data=data, formula=formula,random=random, intervention=intervention, nsim=nSim)
  output  <- MSTerrantSummary(bayesObject=BayesOutput,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention)


  output$Method <- "MLM"
  output$Design <- "MST"
  output$Approach <- "Bayesian"
  output$Function <- "mstBayes"
  class(output) <- "eefAnalytics"
  return(output)
}


#### I. perform Bayesian multilevel linear modeling for Multi-Stage Trials (MST) using JAGS language - internal #####
MST.function <- function(data, formula,random, intervention, nsim){

  #### load required packages ###
  requireNamespace("R2jags", quietly = TRUE) || stop("Please install the 'R2jags' package.")
  #require(R2jags)
  requireNamespace("lme4", quietly = TRUE) || stop("Please install the 'lme4' package.")
  #require(lme4)
  requireNamespace("MCMCvis", quietly = TRUE) || stop("Please install the 'MCMCvis' package.")
  #require(MCMCvis)
  requireNamespace("coda", quietly = TRUE) || stop("Please install the 'coda' package.")
  #require(coda)

  #### check if it is a MST design
  Pdata <- na.omit(data[,c(all.vars(formula),random)])
  chk <- sum(rowSums(table(Pdata[,c(random, intervention)])!=0)>1)
  if(chk ==0){stop("This is not a MST design, try 'CRT.function' instead")}

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
  tmp <- colnames(fixedDesignMatrix )
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
  # jags.UNCparams <- c("sigma", "sigma.tt", "icc", "UNC.ES.Within", "UNC.ES.Total", "UNC.g.with", "UNC.g.Total")
  # 1. UNC Jags model -----------
  filenames_MLM_UNC <- file.path("inst/jags/MLM_UNC.txt")
  cat(paste("
                model {
                # Likelihood
                for (i in 1:N) {
                post[i] ~ dnorm(mu[i], tau)
                mu[i] <- beta0 + u[school[i]]
                }

                # Random intercepts for each cluster
                for (j in 1:M) {
                u[j] ~ dnorm(0, tau_u) # Random effect for each cluster
                }

                # Priors
                beta0 ~ dnorm(0, 0.0001) # Overall intercept
                tau ~ dgamma(0.001, 0.001) # Precision for within-cluster variation (residual variance)
                sigma <- 1 / sqrt(tau) # Within-cluster standard deviation (Residual SD)

                tau_u ~ dgamma(0.001, 0.001) # Precision for between-cluster variation
                sigma_u <- 1 / sqrt(tau_u) # Between-cluster standard deviation

                # Total variance
                sigma.tt <- sigma_u^2 + sigma^2 # Total variance (between + within variance)
                # ICC calculation
                icc <- sigma_u^2 / sigma.tt # Intraclass correlation coefficient

                # ICC and TOTAL VARIANCE
                UNC.icc <- icc

                #sigmas
                UNC.sigma.Total <- sigma.tt
                UNC.sigma.Within <- sigma

                # Effect size calculation
                # UNC.ES.Within: Effect size based on within-cluster variance
                #UNC.ES.Within <- beta0 / sigma # beta0 divided by within-cluster standard deviation

                # UNC.ES.Total: Effect size based on total variance
                #UNC.ES.Total <- beta0 / sqrt(sigma.tt) # beta0 divided by total standard deviation

                # Hedges' g correction for effect size (within-cluster)
                #UNC.g.with <- UNC.ES.Within * (1 - (3 / (4 * (N - 2) - 1)))
                # Hedges' g correction for effect size (total)
                #UNC.g.Total <- UNC.ES.Total * (1 - (3 / (4 * (N - 2) - 1)))
                }
                ")
      ,file=filenames_MLM_UNC
  )

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
  CONdata <- as.list(Pdata1);CONdata[c("N","M", "p","UNC.sigma","UNC.sigma.tt","UNC.icc")]<- c(N,M,p,UNC.sigma,UNC.sigma.tt,UNC.icc) #jags data
  CONdata[["tt"]] <- Pdata[,intervention]+1 #for random part,intervention=1:2 instead of 0:1
  Post.T <- match(intervention, var) # position of intervention

  # COND Jags parameters to monitor
  jags.params <- c("COND.ES.Within","COND.ES.Total","COND.g.with","COND.g.Total","COND.sigma.Within","COND.sigma.Total","COND.sigma.between","COND.Trt.schl","COND.icc",
                   "UNC.ES.Within","UNC.ES.Total","UNC.g.with","UNC.g.Total","UNC.sigma.Within","UNC.sigma.Total","UNC.ICC","beta", "b2diff")

  # 2. COND Jags model -----------
  filenames_MST <- file.path("inst/jags/MST.txt")
  cat(paste("
                model{
                for(i in 1:N){
                post[i] ~ dnorm(mu[i],tau)
                mu[i] <-", paste0("beta[1]+",paste0("beta","[",1:length(var1)+1,"]","*",var1,"[i]",collapse ="+"),"+b1[",random,"[i]]+b2[",random,"[i],tt[i]]"),"
                }

                for(j in 1:M){
                b1[j]~dnorm(0.0,tau.b1)
                #b2[j]~dnorm(0.0,tau.b2)
                for(k in 1:2){b2[j,k]~dnorm(0.0,tau.b2)}
                b2diff[j] <- b2[j,2] - b2[j,1]
                }

                tau~dgamma(0.001,0.0001)
                tau.b1~dgamma(0.001,0.0001)
                tau.b2~dgamma(0.001,0.0001)
                sigma<-1/tau
                sigma.b1<-1/tau.b1
                sigma.b2<-1/tau.b2

                for(k in 1:p){beta[k]~dnorm(0.0,1.0E-06)}

                # ICC and TOTAL VARIANCE
                sigma.Total <-sigma + sigma.b1 + sigma.b2
                COND.icc <- (sigma.b1+sigma.b2) * pow(sigma.Total ,-1)
                UNC.ICC <- UNC.icc

                #sigmas
                COND.sigma.Total <- sigma.Total
                COND.sigma.Within <- sigma.b1
                COND.sigma.between <- sigma.b2
                COND.Trt.schl <-sigma.b2
                UNC.sigma.Total <- UNC.sigma.tt
                UNC.sigma.Within <- UNC.sigma

                # Effect size calculation
                COND.ES.Within<- ",paste0("beta[",Post.T,"]"),"/sqrt(COND.sigma.Within)#conditional
                COND.ES.Total <- ",paste0("beta[",Post.T,"]"),"/sqrt(COND.sigma.Total)#conditional
                UNC.ES.Within<- ",paste0("beta[",Post.T,"]"),"/sqrt(UNC.sigma.Within)#unconditional
                UNC.ES.Total <- ",paste0("beta[",Post.T,"]"),"/sqrt(UNC.sigma.Total)#unconditional

                # Hedges' g correction for effect size
                COND.g.with<- COND.ES.Within* (1- (3/(4*(N-2)-1)))#conditional
                COND.g.Total <-  COND.ES.Total* (1- (3/(4*(N-2)-1)))#conditional
                UNC.g.with<- UNC.ES.Within* (1- (3/(4*(N-2)-1)))#unconditional
                UNC.g.Total <-  UNC.ES.Total* (1- (3/(4*(N-2)-1)))#unconditional
                }
                ")
      ,file=filenames_MST
  )

  #list of jags initial values
  jags.inits <- function(){
    list("beta"=betaB, "tau"=tau, "tau.b1"=tau.s, "b1"=rnorm(0,1),"b2"=rnorm(0,1))#Initial values for chain1
    list("beta"=betaB2,"tau"=tau2,"tau.b1"=tau.s2,"b1"=rnorm(0,1), "b2"=rnorm(0,1))#Initial values for chain2
    list("beta"=betaB3,"tau"=tau3,"tau.b1"=tau.s3, "b1"=rnorm(0,1),"b2"=rnorm(0,1))#Initial values for chain3
  }

  #### Summarise CONDITIONAL JAGS output ####
  jag.model<-jags(model.file=filenames_MST,
                  data = CONdata,
                  n.iter= nsim,
                  n.burnin = nsim/2,
                  inits=jags.inits,
                  n.thin = 10,
                  parameters.to.save=jags.params)
  jag.model.upd <- autojags(jag.model)
  jag.model.sum <- MCMCsummary(jag.model.upd,round = 2)[,c("mean","2.5%" , "97.5%")]
  row.names(jag.model.sum)[row.names(jag.model.sum)%in% paste0("beta[",1:p,"]")] <- c("Intercept",var1)

  #### b2 should be diff between group 1 and group 2
  # 3. Jags Output -----------
  #### 3.1 Final Summarised output ####
  b2diff <- grep("b2diff", rownames(jag.model.sum))
  Unconditional= list(ES=data.frame(rbind(Within=jag.model.sum["UNC.ES.Within",],Total=jag.model.sum["UNC.ES.Total",])),
                      covParm=c(Within=jag.model.sum["UNC.sigma.Within",1],Total=jag.model.sum["UNC.sigma.Total",1], ICC= jag.model.sum["UNC.ICC",1]))

  ObjectOUT <- list(jags.model=list(jag.model=jag.model,jag.model.upd=jag.model.upd, jag.model.sum=jag.model.sum),
                    beta=jag.model.sum[c("Intercept",var1),],
                    ES=data.frame(rbind(Within=jag.model.sum["COND.ES.Within",],Total=jag.model.sum["COND.ES.Total",])),
                    covParm=c(Within=jag.model.sum["COND.sigma.Within",1], Between=jag.model.sum["COND.sigma.between",1], Trt.schl=jag.model.sum["COND.Trt.schl",1],
                              Total=jag.model.sum["COND.sigma.Total",1],ICC= jag.model.sum["COND.icc",1]),
                    randomEffects=jag.model.sum[c(b2diff),1],
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
  # Random Effects estimates ---- from mcmc samples
  mcmc_b2 <- mcmc_df[, grep("b2diff", colnames(mcmc_df))]
  colnames(mcmc_b2) <- unique(Pdata[,random])

  # Find column indices for columns containing 'COND.sigma' and 'COND.icc'
  cond_sigma_cols <- grep("COND.sigma", colnames(mcmc_df))
  cond_icc_col <- grep("COND.icc", colnames(mcmc_df))
  # Covariance estimates ---- from mcmc samples
  mcmc_covParm <- mcmc_df[, c(cond_sigma_cols, cond_icc_col)]
  mcmc_covParm <- mcmc_covParm[, c("COND.sigma.Within", "COND.sigma.between", "COND.sigma.Total", "COND.icc")]
  colnames(mcmc_covParm) <- c("Within","Between","Total","ICC")

  # Beta estimates ---- from mcmc samples
  mcmc_Beta <- mcmc_df[, grep("beta", colnames(mcmc_df))]
  colnames(mcmc_Beta)<- colnames(fixedDesignMatrix)

  # Effect Size estimates ---- from mcmc samples
  mcmc_ES <- mcmc_df[, grep("COND.ES", colnames(mcmc_df))]
  mcmc_ES <- mcmc_ES[, c("COND.ES.Within", "COND.ES.Total")]
  colnames(mcmc_ES) <- c("Within","Total")

  ###### 3.2.2 UNConditional mcmc output ####
  # Find column indices for columns containing 'COND.sigma' and 'COND.icc'
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
  Output <- list(randomEffects=mcmc_b2,
                 covParm=mcmc_covParm,
                 Beta=mcmc_Beta,
                 ES=mcmc_ES,
                 Unconditional = unconditional)
  return(Output)
}





#### II. summarise covariance parameters - internal #####
MSTcovSummary <- function(bayesObject){

  covParm <- bayesObject$covParm ### MCMC
  covParm2 <- colMeans(covParm) # remove all roundings that appear in the middle of functions
  covParm3 <- t(apply(covParm,2,function(x)quantile(x,prob=c(0.025,0.975))))
  covParm4 <- data.frame(cbind(covParm2,covParm3))
  covParm5 <- sqrt(covParm4)
  colnames(covParm4) <- c("Variance","95% LB","95% UB")
  rownames(covParm4) <- c("Pupils","Schools","Total","ICC")
  covParm5 <- covParm4[,1]
  return(covParm5)
}

MSTUNCcovSummary <- function(bayesObject){

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
MSTbetaSummary <- function(bayesObject){

  betaParm <- bayesObject$Beta
  betaParm2 <- colMeans(betaParm)
  betaParm3 <- t(apply(betaParm,2,function(x)quantile(x,prob=c(0.025,0.975))))
  betaParm4 <- data.frame(cbind(betaParm2,betaParm3))
  colnames(betaParm4) <- c("Estimate","95% LB","95% UB")
  rownames(betaParm4) <- colnames(betaParm)
  return(betaParm4)
}


#### IV. summarise Effect Sizes - internal #####
MSTesSummary <- function(bayesObject,fixedDesignMatrix,intervention){

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

#### V. summarise random effects - internal #####
MSTschSummary <- function(bayesObject){

  schParm <- bayesObject$randomEffects
  schParm2 <- colMeans(schParm)
  schParm3 <- t(apply(schParm,2,function(x)quantile(x,prob=c(0.025,0.975))))
  schParm4 <- data.frame(School = 1:ncol(bayesObject$randomEffects), cbind(schParm2,schParm3))
  colnames(schParm4) <- c("School", "Estimate","95% LB","95% UB")
  schParm5<- schParm4
  return(schParm5)

}


#### VI. summarise minimum expected effect size - internal #####
MSTesProb <- function(bayesObject,esOutput){
  es <- c(0,0.05,seq(0.1,1,0.1))
  esParm <- bayesObject$ES
  esParm2 <- sapply(es,function(x)colMeans(esParm>=x))
  esParm4 <- data.frame(t(esParm2))
  rownames(esParm4) <- paste0("P(ES>", sprintf("%.2f", es), ")")
  return(esParm4)

}


#### VII. summarise all Bayesian parameters - internal #####
MSTerrantSummary <- function(bayesObject,fixedDesignMatrix,intervention){
  covValues <- MSTcovSummary(bayesObject=bayesObject)
  covValues <- data.frame(covValues[c(2,1,3,4)])
  row.names(covValues) <- c("Schools","Pupils","Total","ICC")
  covValues <- t(covValues)
  row.names(covValues ) <- NULL
  betaValues <- MSTbetaSummary(bayesObject=bayesObject)
  esValues <- MSTesSummary(bayesObject,fixedDesignMatrix,intervention)
  schValues <-MSTschSummary(bayesObject=bayesObject)
  es.prob <- MSTesProb(bayesObject=bayesObject,esOutput=esValues)

  UNCcovValues <- MSTUNCcovSummary(bayesObject=bayesObject)
  UNCcovValues <- data.frame(UNCcovValues)
  row.names(UNCcovValues) <- c("Pupils","Total","ICC")
  UNCcovValues <- t(UNCcovValues)
  row.names(UNCcovValues) <- NULL
  unconditional= list(ES=round(MSTesSummary(bayesObject$Unconditional,fixedDesignMatrix,intervention),2),
                      covParm=round(UNCcovValues,2))

  output <- list(Beta=round(betaValues,2),
                 covParm= round(covValues,2),
                 ES=esValues,
                 ProbES=round(data.frame(es.prob),2),
                 SchEffects=round(schValues,2),
                 Unconditional = unconditional)
}

############################################# END #######################################################

