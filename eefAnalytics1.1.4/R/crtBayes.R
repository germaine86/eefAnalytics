
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
#' @param baseln A string specifying the reference category for the intervention variable. If not provided, the first level will be used as the reference (e.g., baseln = "0" for an intervention with levels 0 and 1).
#' @param nsim Number of MCMC iterations to be performed. A minimum of 10,000 is recommended to ensure convergence.
#' @param threshold a scalar or vector of pre-specified threshold(s) for estimating Bayesian posterior probability such that the observed effect size is greater than or equal to the threshold(s).
#' @param condopt additional arguments of \code{\link[R2jags]{jags}} to be passed exclusively to the conditional model (e.g., defining n.chains only for the conditional model, etc.).
#' @param uncopt additional arguments of \code{\link[R2jags]{jags}} to be passed exclusively to the unconditional model (e.g., defining n.chains only for the unconditional model, etc.).
#' @param data A data frame containing the variables referenced in the formula, including predictors, the clustering variable, and the intervention.
#' @param alpha significant level, default alpha = 0.05.
#' @param digits number of decimal places, by default digits=3
#' @param ... Common additional arguments of \code{\link[R2jags]{jags}} to be passed to both the conditional and unconditional model specifications
#'
#' @return S3 object; a list consisting of:
#' \itemize{
#'   \item \code{Beta}: Estimates and credible intervals for the predictors specified in the model (posterior distributions).
#'   \item \code{ES}: Hedges' g effect size for the intervention(s). If bootstrapping is not used, 95% credible intervals are computed based on MCMC sampling.
#'   \item \code{covParm}: Variance components broken down into between-cluster variance (e.g., between schools), within-cluster variance (e.g., within pupils), and intra-cluster correlation (ICC)..
#'   \item \code{SchEffects}: A vector of the estimated deviation of each school from the intercept.
#'   \item \code{ProbES}: A matrix showing the probability of observing an effect size larger than various thresholds (0, 0.05, 0.10, ...).
#'   \item \code{Model}: A model object from \code{\link[R2jags]{jags}} and an \code{\link[MCMCvis]{MCMCsummary}} object containing only the mean and credible intervals (CIs) as columns.
#'   \item \code{Unconditional}: A list containing the unconditional effect size and variance decomposition.
#' }
#'
#' @example inst/examples/crtBExample.R
crtBayes <- function(formula,random,intervention,baseln,nsim=10000,data,alpha=0.05, digits=3, threshold=c(0,0.05,seq(0.1,1,0.1)),condopt,uncopt,...){UseMethod("crtBayes")}

#' @export
crtBayes.default <- function(formula,random,intervention,baseln,nsim=10000,data,alpha=0.05, digits=3, threshold=c(0,0.05,seq(0.1,1,0.1)),condopt,uncopt,...){stop("No correct formula input given.")}

#' @export
crtBayes.formula <- function(formula,random,intervention,baseln,nsim=10000,data,alpha=0.05, digits=3, threshold=c(0,0.05,seq(0.1,1,0.1)),condopt,uncopt,...){

  #Data preparation
  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),random, intervention))])
  data <- data[order(data[,which(colnames(data)==random)]),]
  tmp2 <- which(colnames(data)==random)
  tmp3 <- which(colnames(data)==intervention)
  chk <- sum(rowSums(table(data[,tmp2], data[,tmp3])!=0)>1)#check crt or mst data: correct for two arms TBD

  #warning and stop rules
  if(nsim < 10000){warning("nsim >= 10000 is recommended")}
  if(length(tmp2)!= 1){stop("Cluster variable misspecified")}
  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}
  if(chk > 0){stop("This is not a CRT design, try 'MST.function' instead")}


  #intervention and set reference level
  if(!missing(baseln)){data[,tmp3] <- relevel(as.factor(data[,tmp3]),as.character(baseln))}
  if(missing(baseln)){data[,tmp3] <- as.factor(data[,tmp3])}


#Design matrix
  fixedDesignMatrix <- data.frame(model.matrix(object=formula, data=data))#as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  fixedDesignMatrix <- fixedDesignMatrix[order(data[,tmp2]),]
  colnames(fixedDesignMatrix)[1]  <- "Intercept"

  if(missing(condopt)){condopt=NULL}
  if(missing(uncopt)){uncopt=NULL}

  #Use if else depending on whether alpha=0.05
  BayesOutput <- CRT.function(data=data,
                              formula=formula,
                              random=random,
                              intervention=intervention,
                              nsim=nsim,
                              alpha=alpha,
                              digits=digits,
                              threshold=threshold,
                              condopt=condopt,
                              uncopt=uncopt)
  #output  <- CRTerrantSummary(BayesOutput,fixedDesignMatrix,intervention,threshold=threshold)
  #output  <- CRTerrantSummary(BayesOutput,fixedDesignMatrix,intervention, threshold, alpha)# GU added aplha
  output <- BayesOutput
  output$Method <- "MLM"
  output$Design <- "CRT"
  output$Approach <- "Bayesian"
  output$Function <- "crtBayes"
  class(output) <- "eefAnalytics"


  return(output)
}

#############################################################################################################################
################################################################################################################################
#### I. perform Bayesian multilevel linear modeling for Multi-Stage Trials (MST) using JAGS language - internal #####
CRT.function <- function(data, formula,random, intervention, nsim,alpha, digits,threshold, condopt,uncopt,...){

  #### load required packages ###
  requireNamespace("R2jags", quietly = TRUE) || stop("Please install the 'R2jags' package.")
  #require(R2jags)
  requireNamespace("lme4", quietly = TRUE) || stop("Please install the 'lme4' package.")
  #require(lme4)
  requireNamespace("MCMCvis", quietly = TRUE) || stop("Please install the 'MCMCvis' package.")
  #require(MCMCvis)
  requireNamespace("coda", quietly = TRUE) || stop("Please install the 'coda' package.")
  #require(coda)

  #### check if it is a CRT design
  Pdata <- na.omit(data)[,c(all.vars(formula),random)]

  #### data preparation
  outcome <- all.vars(formula)[1] ## Y: Posttest
  dummies<- data.frame(model.matrix(formula, data=Pdata)) ## X0-Intercept, X1-Prettest, X2-Intervention
  Pdata0 <- data.frame(na.omit(dummies[,!(names(dummies) %in% "X.Intercept.")]))
  names(Pdata0) <- setdiff(names(dummies), "X.Intercept.")

  Pdata1<-na.omit(cbind(post=Pdata[,outcome],Pdata0)) #all covariates and school dummies
  Pdata1[,random] <- as.numeric(as.factor(Pdata[,random]))


  # Jags data
  var<-names(Pdata1)
  N <- nrow(Pdata1)
  M <- length(unique(Pdata1[,random]))


  #Initial values preparation --------------------------
  #Initial values for chain 1


  formula <- update(formula,paste0("~ .+ (1|",random,")"))
  lm.m1 <- lmer(formula ,Pdata)
  betaB=lm.m1@beta
  vsigma<-as.data.frame(VarCorr(lm.m1));
  tau <-1/(vsigma[1,4]+0.001); #precision
  tau.s <- 1/(vsigma[1,4]+0.001) #precision:one was added to school and school:t to avoid zero in denominator

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
  #=======================
  UNCdata <- list(N=N,M=M, school=Pdata1[,random], post=Pdata1$post)
  jags.UNCparams <- c("sigma","sigma.tt","icc")
  #filenames_MLM_UNC <- file.path("inst/jags/MLM_UNC.txt")#system.file("jags", "MLM_UNC.txt", package = "eefAnalytics")

  # model
  filenames_MLM_UNC <- system.file("jags", "MLM_UNC.txt", package = "eefAnalytics")#file.path("inst/jags/MLM_UNC.txt")#file.path("inst/jags/MLM_UNC.txt")#
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



  # 1. UNC Jags model -----------
  #### Summarise UNCONDITIONAL JAGS output ####

  UNC.ols <-do.call(jags, c(list(model.file=filenames_MLM_UNC, data = UNCdata, n.iter= nsim, n.burnin = nsim/2, inits=UNCjags.inits,n.thin = 10, parameters.to.save=jags.UNCparams, ...), uncopt))
  UNC.ols.upd <- autojags(UNC.ols)
  UNC.sigma <-MCMCsummary(UNC.ols.upd,round = 2)["sigma",c("mean")]
  UNC.sigma.tt <-MCMCsummary(UNC.ols.upd,round = 2)["sigma.tt",c("mean")]
  UNC.icc<-MCMCsummary(UNC.ols.upd,round = 2)["icc",c("mean")]


  ## 2. Conditional model
  #======================
  for(i in 1:length(var)){assign(var[i],Pdata1[,i])}
  threshold1 <- stats::setNames(1:length(threshold), threshold)
  var1 <- var[!var%in% c("post"  ,random)] #list of covariates
  p <- length(var1)+1 #number of covariate + intercept
  Post.T0 <- intersect(c(intervention,paste0(intervention, unique(Pdata[, intervention]))),var )
  Post.T <- stats::setNames(which(var %in% Post.T0),Post.T0)#position of intervention
  data <- as.list(Pdata1);  data[c("N","M", "p","UNC.sigma","UNC.sigma.tt","UNC.icc")]<- c(N,M, p,UNC.sigma,UNC.sigma.tt,UNC.icc)
  data[["Post.T"]] <- Post.T
  data[["threshold1"]] <- threshold1
  data[["threshold"]] <- as.numeric(names(threshold1))
  #jags data

  ######## this gives NA as we have intervention nae has changed into dymmies######
  ### Post.T <- match(intervention, var) #position of intervention
  #################################################################################


  # COND Jags parameters to monitor
  jags.params <- c("COND.ES.Within","COND.ES.Total","COND.ProbES.Within","COND.ProbES.Total","COND.sigma.Within", "COND.sigma.Between" ,"COND.sigma.Total","COND.icc",
                   "UNC.ES.Within","UNC.ES.Total","UNC.ProbES.Within","UNC.Prob.Total","UNC.sigma.Within","UNC.sigma.Total","beta","UNC.ICC", "b1")

  # 2. COND Jags model -----------
  filenames_CRT <-system.file("jags", "CRT.txt", package = "eefAnalytics")# file.path("inst/jags/CRT.txt")#file.path("inst/jags/CRT.txt")
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
                COND.sigma.Within <- sigma
                COND.sigma.Between <- sigma.b1
                UNC.sigma.Total <- UNC.sigma.tt
                UNC.sigma.Within <- UNC.sigma

                # EFFECT SIZE
                for(pp  in round(Post.T)){
                COND.d.Within[pp] <- beta[pp]/sqrt(COND.sigma.Within)#conditional
                COND.d.Total[pp]  <- beta[pp]/sqrt(COND.sigma.Total)#conditional
                UNC.d.Within[pp]  <- beta[pp]/sqrt(UNC.sigma.Within)#unconditional
                UNC.d.Total[pp]  <- beta[pp]/sqrt(UNC.sigma.Total)#unconditional

                #hedges g ES
                COND.ES.Within[pp]<- COND.d.Within[pp]* (1- (3/(4*(N-2)-1)))#conditional
                COND.ES.Total[pp] <-  COND.d.Total[pp]* (1- (3/(4*(N-2)-1)))#conditional
                UNC.ES.Within[pp]<- UNC.d.Within[pp]* (1- (3/(4*(N-2)-1)))#unconditional
                UNC.ES.Total[pp] <-  UNC.d.Total[pp]* (1- (3/(4*(N-2)-1)))#unconditional

                #Posterior probabilities
                for(thd  in round(threshold1)){
                COND.ProbES.Within[pp, thd]<- step(COND.ES.Within[pp] - threshold[thd] )#conditional
                COND.ProbES.Total[pp, thd] <-  step(COND.ES.Total[pp] - threshold[thd] )#conditional
                UNC.ProbES.Within[pp, thd]<- step(UNC.ES.Within[pp] - threshold[thd] )#unconditional
                UNC.Prob.Total[pp, thd] <-  step(UNC.ES.Total[pp] - threshold[thd] )#unconditional
                }
                }
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
  #--------------------------------------------
  #jag.model <- jags(model.file=filenames_CRT, data = data, n.iter= nsim, n.burnin = nsim/2, inits=jags.inits,n.thin = 10, parameters.to.save=jags.params)
  jag.model <- do.call(jags, c(list(model.file=filenames_CRT, data = data, n.iter= nsim, n.burnin = nsim/2, inits=jags.inits,n.thin = 10, parameters.to.save=jags.params, ...), condopt))
  jag.model.upd <- autojags(jag.model)
  jag.model.sum0 <-MCMCsummary(jag.model.upd,round = digits, probs= c(alpha/2, 1-alpha/2))
  jag.model.sum <-jag.model.sum0[,c("mean", paste0(alpha/2*100,"%" ), paste0((1-alpha/2)*100,"%" ))]
  row.names(jag.model.sum)[row.names(jag.model.sum)%in% paste0("beta[",1:p,"]")] <- c("Intercept",var1)



  # Posterior probabilities
  #-------------------------
  ProbES <-  lapply(Post.T, function(i) {
    ProbESP <-  data.frame(cbind(Total1=jag.model.sum[grep(paste0("COND.ProbES.Total\\[",i),row.names(jag.model.sum)),"mean"],
                     Within1 =jag.model.sum[grep(paste0("COND.ProbES.Within\\[",i),row.names(jag.model.sum)),"mean"],
                     Total2=jag.model.sum[grep(paste0("UNC.Prob.Total\\[",i),row.names(jag.model.sum)),"mean"],
                     Within2 =jag.model.sum[grep(paste0("UNC.ProbES.Within\\[",i),row.names(jag.model.sum)),"mean"]))
    row.names(ProbESP) <- paste0("P(ES>",threshold,")")
    round(ProbESP,digits)
    })


  # Random effects
  #-----------------
  SchEffects <- round(data.frame(cbind(
    Schools=unique(Pdata1[,random]),
    Intercept=jag.model.sum[grep("b1",row.names(jag.model.sum)),"mean"])),digits)




  # Effect sizes
  #-----------------
  UNCOND_fn <-  function(trt, COND="COND",  Post.Ti=Post.T, jag.model.sum){
    Within <- if(length(Post.Ti)==1){"Within"}else {paste0("Within[",trt,"]")}
    Total <- if(length(Post.Ti)==1){"Total"}else {paste0("Total[",trt,"]")}
    ES=data.frame(rbind(Within=jag.model.sum[paste0(COND,".ES.",Within),],Total=jag.model.sum[paste0(COND,".ES.", Total),]))
    colnames(ES) <- c("Estimate",paste0((1-alpha)*100, c("% LB","% UB")))
    ES
  }


  # Unconditional output
  #---------------------
  colnames(jag.model.sum) <- c("Estimate",paste0((1-alpha)*100, c("% LB","% UB")))
  Unconditional= list(ES=lapply(Post.T, function(i) UNCOND_fn(trt=i,COND="UNC", jag.model.sum=jag.model.sum)),
                      covParm=c(Within=jag.model.sum["UNC.sigma.Within",1],Total=jag.model.sum["UNC.sigma.Total",1], ICC= jag.model.sum["UNC.ICC",1]),
                      ProbES=lapply(ProbES, function(x) x[, c("Within2","Total2")]))


  #Final all output
  #-----------------
  ObjectOUT <- list(
                    Beta=jag.model.sum[c("Intercept",var1),],
                    ES=lapply(Post.T, function(i) UNCOND_fn(trt=i, jag.model.sum=jag.model.sum)),
                    covParm=c(Schools=jag.model.sum["COND.sigma.Between",1], Pupils=jag.model.sum["COND.sigma.Within",1],Total=jag.model.sum["COND.sigma.Total",1],ICC= jag.model.sum["COND.icc",1]),
                    SchEffects =SchEffects,
                    ProbES=lapply(ProbES, function(x) x[, c("Within1","Total1")]),
                    Model=list(jag.model=jag.model,jag.model.upd=jag.model.upd, jag.model.sum=jag.model.sum),
                    Unconditional=Unconditional)
  return(ObjectOUT)
}

############################################# END #######################################################

