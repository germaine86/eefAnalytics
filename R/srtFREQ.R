
#' Multisite Trial Data.
#'
#' A multisite trial dataset containing 54 schools. This data contains a random sample of test data of pupils and not actual trial data.
#'
#' \itemize{
#'   \item Posttest: posttest scores
#'   \item Prettest: prettest scores
#'   \item Intervention: the indicator for the intervention groups in a two arm trial,
#' coded as 1 for intervention group and 0 for control group.
#'   \item Intervention2: a simulated indicator for intervention groups in a three arm trial.
#'   \item School: numeric school identifier
#' }
#'
#' @format A data frame with 210 rows and 5 variables
#' @name mstData
NULL

###iwq

#' Cluster Randomised Trial Data.
#'
#' A cluster randomised trial dataset containing 22 schools. The data contains a random sample of test data of pupils and not actual trial data.
#'
#' \itemize{
#'   \item Posttest: posttest scores
#'   \item Prettest: prettest scores
#'   \item Intervention: the indicator for intervention groups in a two arm trial,
#'coded as 1 for intervention group and 0 for control group.
#'   \item Intervention2: a simulated indicator for intervention groups in a three arm trial.

#'   \item School: numeric school identifier
#' }
#'
#' @format A data frame with 265 rows and 5 variables
#' @name crtData
NULL


## IMPORTS ##
#' @importFrom lme4 lmer ranef VarCorr
#' @importFrom methods is
#' @importFrom mvtnorm rmvnorm
#' @importFrom graphics abline barplot hist legend lines par plot points text mtext title
#' @importFrom stats residuals confint lm model.frame model.matrix model.response quantile rnorm na.omit update na.omit as.formula relevel resid predict
#' @importFrom ggplot2 ggplot aes element_text facet_grid geom_errorbarh geom_point geom_text geom_vline scale_x_continuous scale_y_continuous theme theme_bw unit ylab coord_cartesian
#' @importFrom R2jags jags autojags
#' @importFrom MCMCvis MCMCsummary
#' @importFrom coda nchain as.mcmc
NULL

#############################################################################
############# SRT main functions ################################################


#' Analysis of Simple Randomised Education Trial using Linear Regression Model.
#'
#' \code{srtFREQ} performs analysis of educational trials under the assumption of independent errors among pupils.
#' This can also be used with schools as fixed effects.
#'
#' @export
#' @param formula the model to be analysed is of the form y~x1+x2+.... Where y is the outcome variable and Xs are the independent variables.
#' @param intervention a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.
#' @param baseln A string variable allowing the user to specify the reference category for intervention variable. When not specified, the first level will be used as a reference.
#' @param nBoot number of bootstraps required to generate bootstrap confidence intervals.
#' @param nPerm number of permutations required to generate permutated p-value.
#' @param ci method for bootstrap confidence interval calculations; options are the Basic (Hall's) confidence interval "basic" or the simple percentile confidence interval "percentile". If not provided default will be percentile.
#' @param seed seed required for bootstrapping and permutation procedure, if not provided default seed will be used.
#' @param data data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and confidence intervals for the variables specified in the model.
#' \item \code{ES}: Conditional Hedges'g effect size and its 95% confidence intervals. If nBoot is not specified, 95% confidence intervals are based on standard errors. If nBoot is specified, they are non-parametric bootstrapped confidence intervals.
#' \item \code{sigma2}: Residual variance.
#' \item \code{Perm}: A "nPerm x w" matrix containing permutated effect sizes using residual variance. "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is produced only if \code{nPerm} is specified.
#' \item \code{Bootstrap}: A "nBoot x w" matrix containing the bootstrapped effect sizes using residual variance. "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is produced only if \code{nBoot} is specified.
#' \item \code{Unconditional}: A list of unconditional effect size, sigma2, Perm and Bootstrap obtained based on variances from the unconditional model (model with only intercept as fixed effect).
#' }
#' @example inst/examples/srtExample.R
srtFREQ <- function(formula,intervention,baseln,nBoot,nPerm,ci,seed,data)UseMethod("srtFREQ")

#' @export
srtFREQ.default <- function(formula,intervention,baseln,nBoot,nPerm,ci,seed,data){stop("No correct formula input given.")}

#' @export
srtFREQ.formula <- function(formula,intervention,baseln,nBoot,nPerm,ci,seed,data=data){

  if(!missing(nPerm) & !missing(nBoot)){stop("Either nPerm or nBoot must be specified")}
  if(missing(nPerm)){nPerm <-0}
  if(missing(nBoot)){nBoot <-0}
  data <- na.omit(data.frame(data)[ ,unique(c(all.vars(formula),intervention))])
  tmp3 <- which(colnames(data)==intervention)
  if(missing(baseln)){data[,tmp3] <- as.factor(data[,tmp3])}
  if(!missing(baseln)){data[,tmp3] <- relevel(as.factor(data[,tmp3]),baseln)}
  #data[,tmp3] <- as.factor(data[,tmp3])
  mf <- model.frame(formula=formula, data=data)
  fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
  tmp <- colnames(fixedDesignMatrix )
  tmp[1]  <- "Intercept"
  colnames(fixedDesignMatrix)<- tmp
  posttest <- model.response(mf)
  intervention <- intervention
  trt <- data[,which(colnames(data)==intervention)]
  #if( is.null(baseln)){trt <- as.factor(trt)}
  #if(!is.null(baseln)){trt <- relevel(as.factor(trt),baseln)}
  btp  <- which(tmp %in% paste0(intervention, unique(data[, intervention])))

  if(!missing(ci) & nBoot==0){stop("Please specify number of bootstraps")}
  if(missing(ci)){ci<-"percentile"}
  optc <- c("basic","percentile")
  if(!ci %in% optc ) {stop("Please specify an allowed bootstrap option")}

  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}
  output <- srt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt, btp=btp)
  output$ES <- sapply(output$ES, function(x) round(x,2), simplify = F)

  if(nPerm>0){

    if(nPerm<1000){warning("Users should specify a higher number of iterations for valid results (nPerm>=1000)")}

    permES1<- matrix(NA,nPerm,(length(unique(trt))-1))
    permES2<- matrix(NA,nPerm,(length(unique(trt))-1))
    #if(missing(seed)){set.seed(1020252)}
    if(!missing(seed)){set.seed(seed)}
    for (i in 1:nPerm){

      data[,which(colnames(data)==intervention)]<- sample(trt)
      fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))

      p2CRTFREQ <-srt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt,btp=btp)

      permES1[i,]  <-  p2CRTFREQ$ES$Conditional[,1]
      permES2[i,]  <-  p2CRTFREQ$ES$Unconditional[,1]
    }
    permES1=data.frame(permES1)
    permES2=data.frame(permES2)
    names(permES1) <- row.names(output$ES$Conditional)
    names(permES2) <- row.names(output$ES$Conditional)

    Perm<- data.frame(cbind(permES_cond=permES1, permES_uncond=permES2))
    perm.names <- paste0( rep(c("cond", "uncond"),each=dim(permES1)[2]),"_" ,rep(row.names(output$ES$Conditional),dim(permES1)[2]))
    names(Perm) <- perm.names
    output$Perm<-Perm
  }

  if(nBoot > 0){

    if(nBoot<1000){warning("Users should specify a higher number of iterations for valid results (nBoot>=1000)")}


    tid <- c(1:nrow(fixedDesignMatrix))
    #if(missing(seed)){set.seed(1020252)}
    if(!missing(seed)){set.seed(seed)}
    bootSamples <- sapply(c(1:nBoot),function(x)sample(tid,replace=TRUE))

    bootResults <- sapply(1:dim(bootSamples)[2],function(bt)srt.srt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,btp=btp,bt=bootSamples[,bt]), simplify = F)
    bootResults2 <- data.frame(do.call(rbind, bootResults))

    bootES <- apply(bootResults2,2,function(x)quantile(x,prob=c(0.025,0.975),na.rm=TRUE))
    all.ES <- as.data.frame(do.call(rbind, output$ES))[,1]

    if(ci=="basic") {
        tempb<-list()
        g=1
        for(i in 1:(length(btp))) {
          tempb[[1]]<-2*as.numeric(all.ES[g])-as.numeric(bootES["2.5%",g])
          tempb[[2]]<-2*as.numeric(all.ES[(g+1)])-as.numeric(bootES["97.5%",(g+1)])
          bootES["97.5%",g]<-tempb[[1]]
          bootES["2.5%",(g+1)]<-tempb[[2]]
          g=g+2
        }
      }

    bootES2.1 <- round(data.frame(t(rbind(all.ES,bootES))),2)
    names(bootES2.1)=colnames(output$ES$Conditional)

    condName <- grep("uncond", names(bootResults2), invert = T)
    UncondName <- grep("uncond", names(bootResults2), invert = F)
    bootES2  <- list(Conditional=bootES2.1[condName,], Unconditional=bootES2.1[UncondName,])
    bootES2 <- lapply(bootES2, function(x){ row.names(x)<-row.names(output$ES$Conditional); x})

    bootR1 <- data.frame(bootResults2[,condName])
    bootR2  <- data.frame(bootResults2[,UncondName])
    names(bootR1)<- row.names(output$ES$Conditional)
    names(bootR2)<- row.names(output$ES$Conditional)

    output$ES <-bootES2
    output$Bootstrap <- bootResults2
  }


  output1 <- list()
  output1$Beta   <- output$Beta
  output1$ES     <- output$ES$Conditional
  output1$sigma2 <- output$sigma2["Conditional"]
  if(nPerm > 0){output1$permES    <- round(permES1,2)}
  if(nBoot > 0){output1$Bootstrap <- round(bootR1,2)}

  output1$Unconditional$ES     <- output$ES$Unconditional
  output1$Unconditional$sigma2 <- output$sigma2["Unconditional"]
  if(nPerm > 0){output1$Unconditional$permES <- round(permES2,2)}
  if(nBoot > 0){output1$Unconditional$Bootstrap <- round(bootR2,2)}


  output1$Method <- "LM"
  if(nBoot > 0){output1$CI <- ci}
  output1$Function <- "srtFREQ"
  class(output1) <- "eefAnalytics"
  return(output1)
}

########################################################################################

## - internal SRT functions
srt <- function(posttest,fixedDesignMatrix,intervention,trt,btp){

  freqFit <- lm(posttest~ fixedDesignMatrix-1)
  cit <- confint(freqFit)
  citt <- rowSums(is.na(cit))
  betaB <- data.frame(cbind(summary(freqFit)$coefficients[which(citt==0),1],cit[which(citt==0),]))
  row.names(betaB)<- colnames(fixedDesignMatrix)[which(citt==0)]
  colnames(betaB) <- c("Estimate","95% LB ","95% UB")
  betaB <- betaB

  btp3 <- substring(colnames(fixedDesignMatrix)[btp],nchar(intervention)+1,nchar(colnames(fixedDesignMatrix)[btp]))
  tmpTRT <- table(trt)
  baseln2 <- names(tmpTRT)[!names(tmpTRT)%in%btp3 ]
  sigma1 <- c(Conditional=summary(freqFit)$sigma,
              Unconditional=summary(lm(posttest~1))$sigma)

  output2.0 <- matrix(NA,length(btp),3 )
  colnames(output2.0)<- c("Estimate","95% LB","95% UB")
  row.names(output2.0) <- row.names(betaB)[btp]
  # changed n.c n.t and baseln2
  output2 <- list()
  for( j in names(sigma1)){
    sd.pool <- sigma1[j]
    for( i in 1:length(btp)){
      beta <- betaB[btp[i],1]
      cd <- (beta/sd.pool)
      trt2 <- c(baseln2,btp3[i])
      n.c <- tmpTRT[names(tmpTRT)==trt2[1]]
      n.t <- tmpTRT[names(tmpTRT)==trt2[2]]
      var.cd <- ((n.t+n.c)/(n.t*n.c)+cd^2/(2*(n.t+n.c)))
      se.cd <- sqrt(var.cd)
      cd.lb <- (cd - 1.96*se.cd)
      cd.ub <- (cd + 1.96*se.cd)
      j.df <- (1 - (3/(4*(n.t+n.c-2)-1)))
      g <- (j.df*cd)
      var.g <- (j.df^2 * var.cd)
      se.g <- sqrt(var.g)
      g.lb <- (g - 1.96*se.g)
      g.ub <- (g + 1.96*se.g)

      output2.0[i,] <- c(g, g.lb, g.ub)
    }
    output2[[j]] <- output2.0
  }

  sigma2=round(sigma1^2,2)
  output <- list(Beta=round(betaB,2),ES=output2,sigma2=sigma2)
  return(output)

}




## - internal
srt.srt<- function(posttest,fixedDesignMatrix,intervention,bt,btp,variances=variances){

  posttest2 <- posttest[bt]
  fixedDesignMatrix2 <- fixedDesignMatrix[bt,]

  freqFit <- try(lm(posttest2~ fixedDesignMatrix2-1),silent=TRUE)
  if(attr(freqFit,"class")!="try-error"){

    betaB <- data.frame(summary(freqFit)$coefficients[,1])
    ntpp <- as.character(sapply(as.character(rownames(betaB)),function(x)substring(x,(nchar("fixedDesignMatrix2")+1),nchar(x))))
    row.names(betaB)<- ntpp
    betaB <- betaB

    sigma1 <- c(Conditional=summary(freqFit)$sigma,
                Unconditional=summary(lm(posttest~1))$sigma)

    output1<- sapply(names(sigma1), function(x) sapply(1:length(btp ), function (i) betaB[btp[i],1]/sigma1[x] ))
  }

  boot.names <- paste0( rep(c("cond", "uncond"),each=length(btp)),"_" ,rep(row.names(betaB)[btp],length(btp)))
  output1<- matrix(output1,1,(length(btp)*2))
  colnames(output1)<- boot.names
  return(output1)
}




