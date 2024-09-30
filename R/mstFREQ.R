
#############################################################################
############# MST main functions ################################################

#############################################################################
############# MST main functions ################################################
#' Analysis of Multisite Randomised Education Trials using Multilevel Model under a Frequentist Setting.
#'
#' \code{mstFREQ} performs analysis of multisite randomised education trials using a multilevel model under a frequentist setting.
#'
#' @export
#' @param formula the model to be analysed is of the form y ~ x1+x2+.... Where y is the outcome variable and Xs are the independent variables.
#' @param random a string variable specifying the "clustering variable" as contained in the data. See example below.
#' @param intervention a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.
#' @param baseln A string variable allowing the user to specify the reference category for intervention variable. When not specified, the first level will be used as a reference.
#' @param nBoot number of bootstraps required to generate bootstrap confidence intervals.
#' @param nPerm number of permutations required to generate permutated p-value.
#' @param type method of bootstrapping including case re-sampling at student level "case(1)", case re-sampling at school level "case(2)", case re-sampling at both levels "case(1,2)" and residual bootstrapping using "residual". If not provided, default will be case re-sampling at student level.
#' @param ci method for bootstrap confidence interval calculations; options are the Basic (Hall's) confidence interval "basic" or the simple percentile confidence interval "percentile". If not provided default will be percentile.
#' @param seed seed required for bootstrapping and permutation procedure, if not provided default seed will be used.
#' @param data data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}: Estimates and confidence intervals for variables specified in the model.
#' \item \code{ES}: Conditional Hedge's g effect size (ES) and its 95% confidence intervals. If nBoot is not specified, 95% confidence intervals are based on standard errors. If nBoot is specified, they are non-parametric bootstrapped confidence intervals.
#' \item \code{covParm}: A list of variance decomposition into between cluster variance-covariance matrix (schools and school by intervention) and within cluster variance (Pupils). It also contains intra-cluster correlation (ICC).
#' \item \code{SchEffects}: A vector of the estimated deviation of each school from the intercept and intervention slope.
#' \item \code{Perm}: A "nPerm x 2w" matrix containing permutated effect sizes using residual variance and total variance. "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is produced only when \code{nPerm} is specified.
#' \item \code{Bootstrap}: A "nBoot x 2w" matrix containing the bootstrapped effect sizes using residual variance (Within) and total variance (Total). "w" denotes number of intervention. "w=1" for two arm trial and "w=2" for three arm trial excluding the control group. It is only produced when \code{nBoot} is specified.
#' \item \code{Unconditional}: A list of unconditional effect sizes, covParm, Perm and Bootstrap obtained based on variances from the unconditional model (model with only the intercept as a fixed effect).
#' }
#' @example inst/examples/mstExample.R
mstFREQ<- function(formula,random,intervention,baseln,nPerm,data,type,ci,seed,nBoot)UseMethod("mstFREQ")

#' @export
mstFREQ.default <- function(formula,random,intervention,baseln,nPerm,data,type,ci,seed,nBoot){stop("No correct formula input given.")}

#' @export
mstFREQ.formula <- function(formula,random,intervention,baseln,nPerm,data,type,ci,seed,nBoot){

  data <- na.omit(data[ ,unique(c(all.vars(formula),random, intervention))])
  data <- data[order(data.frame(data)[,which(colnames(data)==random)],data[,which(colnames(data)==intervention)]),]
  trt <- data[,which(colnames(data)==intervention)]
  #trt <- as.factor(trt)
  if(missing(baseln)){trt <- as.factor(trt)}
  if(!missing(baseln)){trt <- relevel(as.factor(trt),baseln)}

  tmp2 <- which(colnames(data)==random)
  cluster2 <- data[,tmp2]

  chk <- sum(rowSums(table(cluster2,trt)!=0)>1)
  if(chk ==0){stop("This is not a MST design")}

  if(!missing(nPerm) & !missing(nBoot)){stop("Either nPerm or nBoot must be specified")}
  if(missing(nPerm)){nPerm <-0}
  if(missing(nBoot)){nBoot <-0}
  if(!missing(type) & nBoot==0 | !missing(ci) & nBoot==0){stop("Please specify number of bootstraps")}
  if(missing(type)){type<-"case(1)"}
  if(missing(ci)){ci<-"percentile"}
  optc <- c("basic","percentile")
  optt <- c("case(1)","residual","case(2)","case(1,2)")
  if(!type %in% optt | !ci %in% optc ) {stop("Please specify an allowed bootstrap option")}
  tmp3 <- which(colnames(data)==intervention)
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
  posttest <- model.response(mf)
  intervention <- intervention
  btp  <- which(tmp %in% paste0(intervention, unique(data[, intervention]) ))

  if(length(tmp2)!= 1){stop("Cluster variable misspecified")}
  if(length(tmp3)!= 1){stop("Intervention variable misspecified")}

  output <- rbd(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt,cluster=cluster,btp=btp)

  if(nPerm > 0){
    if(nPerm<999){warning("Users should specify a higher number of iterations for valid results (nPerm>=1000)")}
    output$Perm <- mst.perm(formula,data,trt,intervention,nPerm,random,cluster,btp,seed)
    output$Conditional$Perm <- round(data.frame(output$Perm$Conditional),2)
    output$Unconditional$Perm <- round(data.frame(output$Perm$Unconditional),2)
  }

  if(nBoot>0){
    if(nBoot<999){warning("Users should specify a higher number of iterations for valid results (nBoot>=1000)")} #moved up from line 99
    if(type=="residual"){
      cov<-list()
      ran<-list()
      res<-list()
      pred<-list()
      Fit<-list()
      sigma2=c(Cond=NA, Uncond=NA)#matrix( nrow = 1, ncol = 2)

      # Fit[[2]] <- lmer(posttest~ 1+(1|cluster))
      # cov[[2]]<- as.matrix(VarCorr(Fit[[2]])$cluster)[,1]
      # sigma2[,2]<- as.numeric(summary(Fit[[2]])$sigma^2)
      #
      # Fit[[1]] <- lmer(posttest~ fixedDesignMatrix-1+(1+trt|cluster))
      # cov[[1]]<- as.matrix(VarCorr(Fit[[1]])$cluster)[,1:(length(btp)+1)]
      # sigma2[,1]<- as.numeric(summary(Fit[[1]])$sigma^2)

      for(j in names(sigma2) ){#1:2
        #(un)conditional models
        if(j=="Uncond"){Fit[[j]] <- lmer(posttest~ 1+(1|cluster))}#added
        if(j=="Cond"){Fit[[j]] <- lmer(posttest~ fixedDesignMatrix-1+(1+trt|cluster))} #added
        btp1 <- ifelse(j=="Cond",(length(btp)+1),1) #added
        cov[[j]] <- as.matrix(VarCorr(Fit[[j]])$cluster)[,1:btp1]
        sigma2[j] <- as.numeric(summary(Fit[[j]])$sigma^2)
        ran[[j]] <- ranef(Fit[[j]])$cluster
        #res[[j]]<-resid(Fit[[j]])
        pred[[j]] <- predict(Fit[[j]], re.form=NA)
        res[[j]] <- resid(Fit[[j]]) - mean(resid(Fit[[j]]))#res[[j]]-mean(res[[j]])
        for(i in 1:btp1){
          ran[[j]][,i] <- ran[[j]][,i]-mean(ran[[j]][,i]) #ran[[1]][,i]<- ran[[1]][,i]-mean(ran[[1]][,i])
        }
      }
      #ran[[2]][,1]<-ran[[2]][,1]-mean(ran[[2]][,1])

      tryCatch(
        for(i in names(sigma2) ){#1:2
          res[[i]]<-reinfl(varcor=sigma2[i],J=length(res[[i]]),res=as.matrix(res[[i]]))
          ran[[i]]<-reinfl(varcor=cov[[i]],J=length(ran[[i]][,1]),res=as.matrix(ran[[i]]))#dose not run
        },
        error = function(e) {message(e,"Reflating of residuals has failed; variance estimates of residuals may be shrunk towards zero.")})

      clusters <- list(as.numeric(rownames(ran[[1]])), as.numeric(rownames(ran[[2]])))
      ran<-Map(cbind, ran, cluster = clusters)
      fixedDesignMatrix<-cbind(fixedDesignMatrix,cluster)

      merged<-sapply(names(sigma2), function(x) {tryCatch(merge(cbind(cluster,pred[[x]],res[[x]]),ran[[x]],all.x=TRUE,by="cluster"))})

      colnames(merged[["Cond"]])[2:3]<-c("feprd","e.rsd")
      colnames(merged[["Uncond"]])[2:4]<-c("Ufeprd","Ue.rsd","Urand.ef")

      fixedDesignMatrix<-cbind(fixedDesignMatrix,cbind(merged[[1]][,-which(names(merged[[1]]) %in% "cluster")],merged[[2]][,-which(names(merged[[2]]) %in% "cluster")]))
    }

    tid <- c(1:nrow(fixedDesignMatrix))

    #set.seed(1020252)
    if(!missing(seed)){set.seed(seed)}
    bootSamples <- NULL

    for(ii in 1:length(unique(cluster))){
      selID <- tid[cluster==unique(cluster)[ii]]
      if(length(selID)>0){
        selID2<- sapply(c(1:nBoot),function(x)selID [sample(1:length(selID), length(selID),replace=TRUE)])
        bootSamples <- rbind(bootSamples ,selID2)
      }
    }

    bootResults <- apply(bootSamples ,2,function(bt)rbd.rbd(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt,cluster=cluster,type=type,bt=bt, btp=btp))
    #bootSamples<-as.matrix(c(1:nBoot))
    #bootResults <- apply(bootSamples ,1,function(bt)rbd.rbd(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt,cluster=cluster,type=type,bt=bt, btp=btp))

    bootES <- bootCompile(output=output,trt=trt,bootResults=bootResults,intervention=intervention,ci=ci)
    output$ES <- bootES
    output$Bootstrap <- bootResults
    output$Conditional$Bootstrap <- data.frame(t(sapply(1:nBoot, function(i)unlist(bootResults[[i]]$Conditional))))
    output$Unconditional$Bootstrap <- data.frame(t(sapply(1:nBoot, function(i)unlist(bootResults[[i]]$Unconditional))))
    Bnames <- gsub("Estimate1","Within", gsub("Estimate2", "Total",names(output$Conditional$Bootstrap)))
    names(output$Conditional$Bootstrap) <- Bnames
    names(output$Unconditional$Bootstrap) <- Bnames
  }
  output1 <- list()
  output1$Beta   <- output$Beta
  output1$covParm <- output$covParm$Conditional
  output1$ES     <- output$ES$Conditional
  output1$SchEffects     <- output$SchEffects
  if(nPerm > 0){output1$permES    <- output$Conditional$Perm}
  if(nBoot > 0){output1$Bootstrap <- output$Conditional$Bootstrap}

  output1$Unconditional$ES     <- output$ES$Unconditional
  output1$Unconditional$covParm <- output$covParm$Unconditional
  if(nPerm > 0){output1$Unconditional$permES <- output$Unconditional$Perm}
  if(nBoot > 0){output1$Unconditional$Bootstrap <-  output$Unconditional$Bootstrap}


  output1$Method <- "MLM"
  if(nBoot > 0){output1$Type <- type}
  if(nBoot > 0){output1$CI <- ci}
  output1$Function <- "mstFREQ"
  class(output1) <- "eefAnalytics"
  return(output1)
}



###########################################################################################

## - internal

rbd <- function(posttest,fixedDesignMatrix,intervention,trt,cluster,btp){

  freqFit <- lmer(posttest~ fixedDesignMatrix-1+(1+trt|cluster))
  np<- row.names(summary(freqFit)$coef)
  cit <- confint(freqFit,np)
  betaB <- data.frame(cbind(summary(freqFit)$coefficients[,1],cit))
  row.names(betaB)<- colnames(fixedDesignMatrix)
  colnames(betaB) <- c("Estimate","95% LB ","95% UB")
  betaB <- betaB
  var.B2<- as.matrix(VarCorr(freqFit)$cluster)
  var.B2_1<-as.data.frame(VarCorr(freqFit)$cluster)
  var.B2_2<-as.data.frame(VarCorr(freqFit))
  var.B3 <- diag(var.B2)[-1]
  var.sch <- var.B2[1,1]
  vcov.schTrt <- t(as.matrix(na.omit(var.B2_2$vcov[var.B2_2$var1=="(Intercept)"][-1])))
  var.E <- var.B3
  var.W<- summary(freqFit)$sigma^2
  N <- as.numeric(length(summary(freqFit)$res))
  N.t <- as.matrix(summary(trt))[-1];
  var.tt <- var.sch+var.W+sum(N.t/N*(var.B3+2*vcov.schTrt))
  ICC1 <- sum(var.B2)/var.tt
  sigmaBE1 <- list(round(var.B2_1,2),round(var.W,2),round(var.tt,2),round(ICC1,2))
  names(sigmaBE1)<-c('School:trt', 'Pupils', 'Total', 'ICC')
  freqFit1<- lmer(posttest~1+(1|cluster))
  var.B21<- as.matrix(VarCorr(freqFit1)$cluster)
  var.B21_1<-as.data.frame(VarCorr(freqFit1)$cluster)
  var.sch1 <- var.B21[1,1]
  var.W1<- summary(freqFit1)$sigma^2
  var.tt1 <- var.sch1+var.W1
  ICC2 <- sum(var.B21)/var.tt1
  sigmaBE2 <- list(round(var.B21_1,2),round(var.W1,2),round(var.tt1,2),round(ICC2,2))
  names(sigmaBE2)<-c('School', 'Pupils', 'Total', 'ICC')
  ICC <-(c(conditional=ICC1,unconditional=ICC2))
  sigmaBE <- list(sigmaBE1,sigmaBE2)
  names(sigmaBE)<-c('Conditional','Unconditional')
  sigmaBE <- sigmaBE
  schRand <- data.frame(unique(cluster),ranef(freqFit)$cluster)
  names(schRand)[1]<- "Schools"
  names(schRand)[2]<- "Intercept"
  #names(schRand)<- c("Schools","Intercept", paste0('trt', 1:(ncol(ranef(freqFit)$cluster)-1)))
  #btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)

  output2 <- list()
  for( i in 1:length(btp)){

    beta <- betaB[btp[i],1]
    group <- fixedDesignMatrix[,btp[i]]
    var.e<- var.E[i]
    esWithin <- g.within.mst(var.w=var.W, var.e=var.e, beta=beta, group=group, schoolID=cluster)
    esTotal <- g.total.mst(var.w=var.W, var.e=var.e, var.tt=var.tt, beta=beta, group=group, schoolID=cluster)
    outputc <- data.frame(rbind(esWithin,esTotal))
    colnames(outputc) <- c("Estimate","95% LB","95% UB")
    rownames(outputc) <- c("Within","Total")
    output2[[i]] <- round(outputc,2)
  }
  names(output2) <- row.names(betaB)[btp]
  output3<-list()
  for( i in 1:length(btp)){

    beta <- betaB[btp[i],1]
    group <- fixedDesignMatrix[,btp[i]]
    esWithin1 <- g.within(var.w=var.W1, beta=beta, icc=ICC2, group=group, schoolID=cluster)
    esTotal1 <- g.total(var.tt=var.tt1, beta=beta, icc=ICC2, group=group, schoolID=cluster)
    outputu <- data.frame(rbind(esWithin1,esTotal1))
    colnames(outputu) <- c("Estimate","95% LB","95% UB")
    rownames(outputu) <- c("Within","Total")
    output3[[i]] <- round(outputu,2)
  }
  names(output3) <- row.names(betaB)[btp]
  output2.3<-list(Conditional=output2,Unconditional=output3)
  output <- list(Beta=round(betaB,2),covParm=sigmaBE,ES=output2.3,SchEffects=round(schRand,2))

  return(output)
}


## - internal

rbdP <- function(posttest,fixedDesignMatrix,intervention,trt,cluster,btp){

  freqFit <- lmer(posttest~ fixedDesignMatrix-1+(1+trt|cluster))
  np<- row.names(summary(freqFit)$coef)
  cit <- confint(freqFit,np)
  betaB <- data.frame(cbind(summary(freqFit)$coefficients[,1],cit))
  row.names(betaB)<- colnames(fixedDesignMatrix)
  #colnames(betaB) <- c("Estimate","95% LB ","95% UB")
  betaB <- betaB
  var.B2<- as.matrix(VarCorr(freqFit)$cluster)
  var.B2_1<-as.data.frame(VarCorr(freqFit)$cluster)
  var.B2_2<-as.data.frame(VarCorr(freqFit))
  var.B3 <- diag(var.B2)[-1]
  var.sch <- var.B2[1,1]
  vcov.schTrt <- t(as.matrix(na.omit(var.B2_2$vcov[var.B2_2$var1=="(Intercept)"][-1])))
  var.E <- var.B3
  var.W<- summary(freqFit)$sigma^2
  N <- as.numeric(length(summary(freqFit)$res))
  N.t <- as.matrix(summary(trt))[-1];
  var.tt <- var.sch+var.W+sum(N.t/N*(var.B3+2*vcov.schTrt))
  ICC1 <- sum(var.B2)/var.tt
  sigmaBE1 <- list(round(var.B2_1,2),round(var.W,2),round(var.tt,2),round(ICC1,2))
  names(sigmaBE1)<-c('School:trt', 'Pupils', 'Total', 'ICC')
  freqFit1<- lmer(posttest~1+(1|cluster))
  var.B21<- as.matrix(VarCorr(freqFit1)$cluster)
  var.B21_1<-as.data.frame(VarCorr(freqFit1)$cluster)
  var.sch1 <- var.B21[1,1]
  var.W1<- summary(freqFit1)$sigma^2
  var.tt1 <- var.sch1+var.W1
  ICC2 <- sum(var.B21)/var.tt1
  sigmaBE2 <- list(round(var.B21_1,2),round(var.W1,2),round(var.tt1,2),round(ICC2,2))
  names(sigmaBE2)<-c('School', 'Pupils', 'Total', 'ICC')
  ICC <-(c(conditional=ICC1,unconditional=ICC2))
  sigmaBE <- list(sigmaBE1,sigmaBE2)
  names(sigmaBE)<-c('Conditional','Unconditional')
  sigmaBE <- sigmaBE
  schRand <- data.frame(unique(cluster),ranef(freqFit)$cluster)
  names(schRand)<- c("Schools","Intercept", paste0('trt', 1:(ncol(ranef(freqFit)$cluster)-1)))
  #btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)

  output2 <- list()
  for( i in 1:length(btp)){

    beta <- betaB[btp[i],1]
    group <- fixedDesignMatrix[,btp[i]]
    var.e<- var.E[i]
    esWithin <- g.within.mst(var.w=var.W, var.e=var.e, beta=beta, group=group, schoolID=cluster)
    esTotal <- g.total.mst(var.w=var.W, var.e=var.e, var.tt=var.tt, beta=beta, group=group, schoolID=cluster)
    outputc <- round(data.frame(rbind(esWithin,esTotal)),2)
    colnames(outputc) <- c("Estimate","95% LB","95% UB")
    rownames(outputc) <- c("Within","Total")
    output2[[i]] <- round(outputc,2)
  }
  names(output2) <- row.names(betaB)[btp]
  output3<-list()
  for( i in 1:length(btp)){

    beta <- betaB[btp[i],1]
    group <- fixedDesignMatrix[,btp[i]]
    esWithin1 <- g.within(var.w=var.W1, beta=beta, icc=ICC2, group=group, schoolID=cluster)
    esTotal1 <- g.total(var.tt=var.tt1, beta=beta, icc=ICC2, group=group, schoolID=cluster)
    outputu <- data.frame(rbind(esWithin1,esTotal1))
    colnames(outputu) <- c("Estimate","95% LB","95% UB")
    rownames(outputu) <- c("Within","Total")
    output3[[i]] <- round(outputu,2)
  }
  names(output3) <- row.names(betaB)[btp]
  output2.3<-list(Conditional=output2,Unconditional=output3)
  output <- list(ES=output2.3)

  return(output)
}

## - internal

rbd.rbd <- function(posttest,fixedDesignMatrix,intervention,trt,cluster,type,bt,btp){

  if(type=="case(1)") {
    posttest2 <- posttest[bt]
    fixedDesignMatrix2 <- fixedDesignMatrix[bt,]
    trt2 <- trt[bt]
    cluster2 <- cluster[bt]

  }
  if(type=="case(2)") {
    chk<-0
    q<-0
    trt2<-as.data.frame(trt)
    colnames(trt2)[1]<-"trt2"
    binded<-as.data.frame(cbind(posttest,fixedDesignMatrix,cluster,trt2))
    while(chk %in% c(0,1)){
      q<-q+1
      if(q>100) {stop("Re-sampling of cases has failed.")}
      sdata <- split(binded, binded[,"cluster"])
      schsamp <- sdata[sample(length(sdata), replace = TRUE)]
      rbinded <- do.call(rbind, schsamp)

      fixedDesignMatrix2<- as.matrix(rbinded[,which(names(rbinded) %in% colnames(fixedDesignMatrix))])
      rownames(fixedDesignMatrix2) <- 1:nrow(fixedDesignMatrix2)
      posttest2<-as.matrix(rbinded[,1])
      posttest<-posttest2
      cluster2<-as.matrix(rbinded[,"cluster"]) #maybe change the way the cluster is identified?
      cluster <- cluster2
      trt2<-as.matrix(rbinded[,"trt2"])
      chk<-colMeans(as.matrix(fixedDesignMatrix2[,btp]))
    }
  }
  if(type=="case(1,2)") {
    chk<-0
    q<-0
    trt2<-as.data.frame(trt)
    colnames(trt2)[1]<-"trt2"
    binded<-as.data.frame(cbind(posttest,fixedDesignMatrix,cluster,trt2))
    while(chk %in% c(0,1)){
      q<-q+1
      if(q>100) {stop("Re-sampling of cases has failed.")}
      sdata <- split(binded, binded[,"cluster"])
      schsamp <- sdata[sample(length(sdata), replace = TRUE)]
      samp <- lapply(schsamp, function(x) x[sample(nrow(x), replace = TRUE), ])
      rbinded <- do.call(rbind, samp)

      posttest2<-as.matrix(rbinded[,1])
      posttest<-posttest2
      fixedDesignMatrix2<-as.matrix(rbinded[,which(names(rbinded) %in% colnames(fixedDesignMatrix))])
      rownames(fixedDesignMatrix2) <- 1:nrow(fixedDesignMatrix2)
      cluster2<-as.matrix(rbinded[,"cluster"])
      cluster<-cluster2
      trt2<-as.matrix(rbinded[,"trt2"])
      chk<-colMeans(as.matrix(fixedDesignMatrix2[,btp]))
    }
  }

  if(type=="residual") {
    k<-ncol(fixedDesignMatrix)-3
    raneffcols<-colnames(fixedDesignMatrix[,c((k-length(btp)):k,(k+3))])
    condrancols<-colnames(fixedDesignMatrix[,c((k-length(btp)):k)])

    fixedDesignMatrix<-as.data.frame(fixedDesignMatrix)
    res<-as.data.frame(fixedDesignMatrix[ ,c("e.rsd","Ue.rsd")])
    ran<-as.data.frame(fixedDesignMatrix[ ,c("cluster",raneffcols)])
    pred<-as.data.frame(fixedDesignMatrix[ ,c("feprd","Ufeprd")])

    ressamp <- as.data.frame(res[sample(length(res[,1]), replace = TRUE),]) #or res[bt]
    sdata <- split(ran,ran[,"cluster"])
    sdata<-lapply(sdata, function(x) x[1,]) #collapse random effects by group

    ransamp<-cbind(do.call(rbind, sdata))
    rownames(ransamp) <- 1:nrow(ransamp)

    ransamp[,raneffcols]<-ransamp[sample(length(ransamp[,1]), replace = TRUE),raneffcols]
    binded<-cbind(fixedDesignMatrix[,"cluster"],pred,ressamp)
    colnames(binded)[1]<-"cluster"
    binded<-merge(binded,ransamp,all.x=TRUE,by="cluster")
    posttest2<-as.matrix(rowSums(binded[,c("feprd","e.rsd",condrancols)]))
    Uposttest2<-as.matrix(rowSums(binded[,c("Ufeprd","Ue.rsd","Urand.ef")]))
    fixedDesignMatrix<- as.matrix(fixedDesignMatrix[,-which(names(fixedDesignMatrix) %in% c("cluster","feprd","e.rsd","Ufeprd","Ue.rsd",raneffcols))])
    rownames(fixedDesignMatrix) <- 1:nrow(fixedDesignMatrix)
    fixedDesignMatrix2<-fixedDesignMatrix
    cluster2 <- as.matrix(cluster)
    trt2<-as.matrix(as.numeric((trt)))
  }

  freqFit <- try(lmer(posttest2~ fixedDesignMatrix2-1+(1+trt2|cluster2)),silent=TRUE)
  output2 <- NULL

  if(!is(freqFit, "try-error")){
    betaB <- data.frame(summary(freqFit)$coefficients[,1])
    row.names(betaB)<- colnames(fixedDesignMatrix)
    betaB <- betaB
    var.B2<- as.matrix(VarCorr(freqFit)$cluster)
    var.B2_1<-as.data.frame(VarCorr(freqFit)$cluster)
    var.B2_2<-as.data.frame(VarCorr(freqFit))
    var.B3 <- diag(var.B2)[-1]
    var.sch <- var.B2[1,1]
    vcov.schTrt <- t(as.matrix(na.omit(var.B2_2$vcov[var.B2_2$var1=="(Intercept)"][-1])))
    var.E <- var.B3
    var.W<- summary(freqFit)$sigma^2
    N <- as.numeric(length(summary(freqFit)$res))
    N.t <- as.matrix(summary(trt))[-1];
    var.tt <- var.sch+var.W+sum(N.t/N*(var.B3+2*vcov.schTrt))
    ICC1 <- sum(var.B2)/var.tt
    sigmaBE1 <- list(round(var.B2_1,2),round(var.W,2),round(var.tt,2),round(ICC1,2))
    names(sigmaBE1)<-c('School:trt', 'Pupils', 'Total', 'ICC')
    if(type=="residual") {posttest2<-Uposttest2}
    freqFit1<- lmer(posttest2~1+(1|cluster2))
    var.B21<- as.matrix(VarCorr(freqFit1)$cluster)
    var.B21_1<-as.data.frame(VarCorr(freqFit1)$cluster)
    var.sch1 <- var.B21[1,1]
    var.W1<- summary(freqFit1)$sigma^2
    var.tt1 <- var.sch1+var.W1
    ICC2 <- sum(var.B21)/var.tt1
    sigmaBE2 <- list(round(var.B21_1,2),round(var.W1,2),round(var.tt1,2),round(ICC2,2))
    names(sigmaBE2)<-c('School', 'Pupils', 'Total', 'ICC')
    ICC <-(c(conditional=ICC1,unconditional=ICC2))
    sigmaBE <- list(sigmaBE1,sigmaBE2)
    names(sigmaBE)<-c('Conditional','Unconditional')
    sigmaBE <- sigmaBE
    schRand <- data.frame(unique(cluster),ranef(freqFit)$cluster)
    names(schRand)<- c("Schools","Intercept", paste0('trt', 1:(ncol(ranef(freqFit)$cluster)-1)))
    #btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)


    output3 <- list()
    for( i in 1:length(btp)){

      beta <- betaB[btp[i],1]
      group <- fixedDesignMatrix[,btp[i]]
      var.e<- var.E[i]
      esWithin <- beta/sqrt(var.W)
      esTotal <- beta/sqrt(var.tt)

      outputc <- data.frame(rbind(esWithin,esTotal))
      names(outputc) <- c("Estimate")
      rownames(outputc) <- c("Within","Total")
      output3[[i]] <- round(outputc,2)
    }
    names(output3) <- row.names(betaB)[btp]
    output4<-list()
    for( i in 1:length(btp)){

      beta <- betaB[btp[i],1]
      group <- fixedDesignMatrix[,btp[i]]
      esWithin1 <- beta/sqrt(var.W1)
      esTotal1 <- beta/sqrt(var.tt1)

      outputu <- data.frame(rbind(esWithin1,esTotal1))
      names(outputu) <- c("Estimate")
      rownames(outputu) <- c("Within","Total")
      output4[[i]] <- round(outputu,2)
    }
    names(output4) <- row.names(betaB)[btp]
    output2<-list(Conditional=output3,Unconditional=output4)

  }
  return(output2)
}



## - internal

g.within.mst <- function(var.w, var.e, beta, group, schoolID){
  t <- group; id <- schoolID
  d.w <- (beta/sqrt(var.w))
  n.itc <- as.data.frame.matrix(table(id,t))
  n.it <- n.itc[n.itc$`0`!=0 & n.itc$`1`!=0, "1"]
  n.ic <- n.itc[n.itc$`0`!=0 & n.itc$`1`!=0, "0"]
  M <- length(unique(id))
  N.t <- sum(table(id[t==1])); N.c <- sum(table(id[t==0]))
  N <- (N.t + N.c)
  nin<-(n.it*n.ic)/(n.it+n.ic)
  ni <- n.it+n.ic
  vterm1 <- 1/(sum(var.w/(var.e+var.w/nin)))
  vterm2 <- ((d.w^2)/((2*N-4*M)))
  se <- sqrt(vterm1+vterm2)
  LB <- (d.w-1.96*se); UB <- (d.w+1.96*se)
  output <- data.frame(d.w, LB, UB)
  names(output) <- c("g", "LB", "UB")
  return(output)
}

## - internal

g.total.mst <- function(var.w, var.e, var.tt, beta, group, schoolID){
  t <- group; id <- schoolID
  n.itc <- as.data.frame.matrix(table(id,t))
  n.it <- n.itc[n.itc$`0`!=0 & n.itc$`1`!=0, "1"]
  n.ic <- n.itc[n.itc$`0`!=0 & n.itc$`1`!=0, "0"]
  M <- length(unique(id))
  N.t <- sum(table(id[t==1])); N.c <- sum(table(id[t==0]))
  N <- (N.t + N.c)
  nin<-(n.it*n.ic)/(n.it+n.ic)
  ni <- n.it+n.ic
  d.t <- (beta/sqrt(var.tt))
  vterm1 <- 1/(sum(var.tt/(var.e+var.w/nin)))
  vterm2 <- ((d.t^2)/((2*N-4*M)))
  se <- sqrt(vterm1+vterm2)
  LB <- (d.t-1.96*se); UB <- (d.t+1.96*se)
  output <- data.frame(d.t, LB, UB)
  names(output)<- c("g", "LB", "UB")
  return(output)
}


## - internal

mst.perm <- function(formula,data,trt,intervention,nPerm,random,cluster,btp,seed){

  data2 <- data
  g <- matrix(NA,nPerm,2*(length(unique(trt))-1))
  g.unc <- matrix(NA,nPerm,2*(length(unique(trt))-1))

  for(i in 1:nPerm){
    #set.seed(12890*i+1)
    if(!missing(seed)){set.seed(seed*i+1)}
    tryCatch({data2[,which(colnames(data)==intervention)]<-unlist(tapply(trt,cluster,function(x)sample(x)))
    data3 <- data2[order(data2[,which(colnames(data2)==random)],data2[,which(colnames(data2)==intervention)]),]
    cluster = data3[,which(colnames(data3)==random)]
    mf <- model.frame(formula=formula, data=data3)
    fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data3)))
    tmp <- colnames(fixedDesignMatrix )
    tmp[1]  <- "Intercept"
    colnames(fixedDesignMatrix)<- tmp
    posttest <- model.response(mf)
    intervention <- intervention
    trt2 <- data3[,which(colnames(data3)==intervention)]
    p2CRTFREQ <-rbdP(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt2,cluster=cluster, btp=btp)

    chkppp <- data.frame(cond=unlist(p2CRTFREQ$ES$Conditional),uncond=unlist(p2CRTFREQ$ES$Unconditional))
    chkppp2 <- c(seq(1,6*(length(unique(trt))-1),6),seq(2,6*(length(unique(trt))-1),6))
    chkppp3 <- chkppp2[order(chkppp2)]
    g[i,]  <-  chkppp[chkppp3,"cond"]
    g.unc[i,]  <-  chkppp[chkppp3,"uncond"]}, error=function(e){})
  }
  ntpp <- rep(names(p2CRTFREQ$ES$Conditional),2)
  ntpp <- ntpp[order(ntpp)]
  wt <- rep(c("Within","Total"),length(names(p2CRTFREQ$ES$Conditional)))
  colnames(g) <- paste(ntpp ,wt,sep="")
  colnames(g.unc) <- paste(ntpp ,wt,sep="")
  g1 <- list(Conditional=g,Unconditional=g.unc)
  return(g1)
}


#function that centers and reflates residuals
reinfl<-function(varcor,J,res){
  Lr<-t(chol(varcor))
  S<-crossprod(res)/J#crossprod(res)= t(res)%*%res
  Ls<-t(chol(S))# if you want to consider the "non-negative definite matrix" as well
  # you could cosider the use of pivote: t(chol(S, pivot = TRUE))
  A<-t(Lr%*%solve(Ls))
  out<-res%*%A
  return(out)
}
