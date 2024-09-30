#' A plot method for an eefAnalytics S3 object obtained from the eefAnalytics package.
#'
#' Plots different figures based on output from eefAnalytics package.
#'
#' @export
#' @param x an output object from the eefAnalytics package.
#' @param group a string/scalar value indicating which intervention to plot.
#' This must be one of the values of intervention variable excluding the control group.
#' For a two arm trial, the maximum number of values to consider is 1 and 2 for three arm trial.
#' @param Conditional a logical value to indicate whether to plot the conditional effect size.
#' The default is Conditional=TRUE, otherwise Conditional=FALSE should be specified for plot based on the unconditional effect size.
#'  Conditional variance is total or residual variance from a multilevel model with fixed effects, whilst unconditional variance is total variance or residual variance from a multilevel model with only intercept as fixed effect.
#' @param ES_Total A logical value indicating whether to plot the effect size based on total variance or within school variance.
#' The default is ES_Total=TRUE, to plot the effect size using total variance.
#' ES_Total=FALSE should be specified for the effect size based on within school or residuals variance.
#' @param slope A logical value indicating whether to return the plot of random intercept (default is slope=FALSE).
#' return other school-by-intervention interaction random slope (s) is slope=TRUE.
#' This argument is suitable only for mstBayes and mstFREQ functions.
#' @param ... arguments passed to \code{\link[graphics]{plot.default}}
#' @details Plot produces a graphical visualisation depending on which model is fitted:
#' \itemize{
#' \item For \code{srtFREQ()}, plot can only be used when \code{nBoot} or \code{nPerm} is specified to visualise the distribution of bootstrapped or permutated values.
#' \item For \code{crtFREQ()} or \code{mstFREQ()}, plot shows the distribution of random intercepts when \code{group=NULL}.
#' It produces histogram of permutated or bootstrapped values when \code{group} is specified and either \code{nBoot} or \code{nPerm} is also specified.
#' }
#' @return Returns relevant plots for each model.
#' @example inst/examples/plotExample.R
plot.eefAnalytics <- function(x,group, Conditional=TRUE,ES_Total=TRUE,slope=FALSE,...){
  if(missing(group)){group=NULL}
    plotObject(analyticObject=x,group=group, Conditional=Conditional,ES_Total=ES_Total,slope=slope, compare=FALSE,modelNames=FALSE,...)

}



#' A plot function to compare different eefAnalytics S3 objects from the eefAnalytics package.
#'
#' @description It generates bar plot that compares the effect size from eefAnalytics' methods.
#'
#' @export
#' @param eefAnalyticsList A list of eefAnalytics S3 objects from eefAnalytics package.
#' @param group a string/scalar value indicating which intervention to plot.
#' This must be one of the values of intervention variable excluding the control group.
#' For a two arm trial, the maximum number of values to consider is 1 and 2 for three arm trial.
#' @param Conditional  a logical value to indicate whether to plot conditional effect size.
#' The default is Conditional=TRUE, otherwise Conditional=FALSE should be specified for plot based on unconditional effect size.
#'  Conditional variance is total or residual variance a multilevel model with fixed effects, whilst unconditional variance is total variance or residual variance from a multilevel model with only intercept as fixed effect.
#' @param ES_Total A logical value indicating whether to plot the effect size based on total variance or within school variance.
#' The default is ES_Total=TRUE, to plot effect size using total variance.
#' ES_Total=FALSE should be specified for effect size based on within school or residuals variance.
#' @param modelNames a string factor containing the names of model to compare. See examples below.
#' @details \code{ComparePlot} produces a bar plot which compares the effect sizes and the associated confidence intervals from the different models.
#' For a multilevel model, it shows the effect size based on residual variance and total variance.
#'
#' @return Returns a bar plot to compare the different methods. The returned figure can be further modified as any \code{\link[ggplot2]{ggplot}}
#' @example inst/examples/compareExample.R
ComparePlot <- function(eefAnalyticsList,group, Conditional=TRUE,ES_Total=TRUE,modelNames){
  if(!is(eefAnalyticsList,"list")){stop("eefAnalyticsList is not a list.")}

  if(!all(unlist(lapply(eefAnalyticsList,function(x) is(x,"eefAnalytics"))))){stop("Not all list objects are a eefAnalytics class object.")}

  if(missing(modelNames)){stop("modelNames must be specified.")}
  if(missing(group)){stop("group must be specified.")}

  plotObject(analyticObject=eefAnalyticsList,group=group, Conditional=Conditional,ES_Total=ES_Total, compare=TRUE,modelNames=modelNames)

}




##################################
#      Internal plot function    #
##################################


plotObject <- function(analyticObject,group, Conditional,ES_Total,slope, compare,modelNames,...){

  if(Conditional ==TRUE){analyticObject2=analyticObject; Condname="Conditional"}
  if(Conditional ==FALSE){analyticObject2=analyticObject$Unconditional; Condname="Unconditional"}
  if(ES_Total ==TRUE){ES_TW<-"Total"}
  if(ES_Total ==FALSE){ES_TW<-"Within"}

  if(compare==TRUE & !is.null(names(analyticObject))){stop("Specify the list of objects to compare")}
  if(!is.null(group)){
  trtname <-rownames(analyticObject2$ES)
  if(is.null(trtname)){trtname <-names(analyticObject2$ES)}
  trtpos <-substr(trtname, (nchar(trtname)-nchar(group)+1),  nchar(trtname))==group
  trt <- trtname[trtpos]
  }
  #bootstrap, permutation and Pprobability plot for SRT model
  #---------------------------------------------------------
  if(sum(analyticObject$Method=="LM") ==1) {
    if(is.null(group)){stop("Group must be specified.")}
    if(!is.null(group)& sum(names(analyticObject)=="ProbES")== 0){
      if(sum(names(analyticObject)=="Bootstrap"|
             names(analyticObject)=="permES")==0){stop("Only relevant for bootstrapped or permutated values")}
      if(sum(names(analyticObject)=="Bootstrap")==1){
        ntp <- nrow(as.matrix(analyticObject$ES))
        if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}

        obs.est <- analyticObject2$ES[trt,1]
        tmp2 <- as.numeric(analyticObject2$Bootstrap[,trt])
        xlabs=paste0(Condname," Bootstrap estimates")
        hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main="")
        abline(v=obs.est,col="red",lwd=2,lty=1)
        abline(v=0,col="grey48",lwd=2,lty=1)
        legend("topright",c("Observed Estimate","Zero-Line"),col=c("red","grey48"),bty="n",lty=1,lwd=2)
      }

      if(sum(names(analyticObject2)=="permES")==1){
        ntp <- nrow(as.matrix(analyticObject2$ES))
        if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}


        Perm.names <-names(analyticObject2$permES)

        obs.est <- analyticObject2$ES[trt,1]
        tmp2 <- as.numeric(analyticObject2$permES[,trt])
        pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
        xlabs=paste0("Permutation values (PermES) based on ",Condname, " ES")

        hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main=paste("P(|PermES| > |ES|)=",pvalue,sep=""))
        abline(v=obs.est,col="red",lwd=2,lty=2)
        abline(v=-obs.est,col="red",lwd=2,lty=2)
        legend("topright",c("(-) Observed Estimate"),col=c("red"),bty="n",lty=2,lwd=2)
      }
    }
    if(sum(names(analyticObject)=="ProbES")> 0 ){
      if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
      tmp2 <- analyticObject2$ProbES[[which(trtpos==TRUE)]]

      thd0<- regmatches(rownames(tmp2),  gregexpr("[[:digit:]]+\\.*[[:digit:]]*",rownames(tmp2)))
      thd <- as.numeric(unlist(thd0))

      par_original <- par()[c("mar","xpd")]
      par_original0<- par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(thd, tmp2[,1], col="blue", type="b",ylim=c(0,max(tmp2)), xlab="Threshold", ylab="Posterior probability")
      on.exit(par(par_original0))
      on.exit(par(par_original))
    }

  }


  #bootstrap, permutation, and Pprobability plot for CRT and MST model
  #--------------------------------------------------------------------
  if(sum(analyticObject$Method=="MLM")==1){

    if(is.null(group)){

      tmp000  <- data.frame(analyticObject$SchEffects)#use analyticObject since both (un)condition has the same SchEffects object.
      if(slope==FALSE){
        tmp00 <- tmp000[,grep("Schools|Intercept|Estimate",names(tmp000))]
        mar11 <-  c(5, 4, 4, 2) + 0.1
      }
      if(slope==TRUE & dim(tmp000)[2] ==2){stop("x must be mstFREQ or mstBAyes object")}
      if(slope==TRUE & dim(tmp000)[2] >2){
        tmp00 <- tmp000[,!(names(tmp000) %in% "Intercept")]
        if(dim(tmp00)[2]==2){mar11 <- c(5, 4, 4, 2) + 0.1}
        if(dim(tmp00)[2]==3){mar11 <- c(5, 2, 4, 0) + 1.0}
        if(dim(tmp00)[2] >3){mar11 <- c(3, 2, 0, 0) + 1.0}}

      op <- par(mfrow = c(floor(dim(tmp00)[2]/2),round(dim(tmp00)[2]/2)),
                mar = mar11)


      for(i in 2:dim(tmp00)[2]){
        tmp <- data.frame(y=tmp00[,i],x=c(1:length(tmp00[,i])))
        tmp2 <- tmp[order(tmp$y),]
        ylabs=gsub("trt", "Intervention ",gsub("Estimate","Intercept", names(tmp00)[i]))
        barplot(tmp2$y,names.arg=tmp2$x,las=2,col="cornflowerblue",border="cornflowerblue")
        if(dim(tmp00)[2]<=2){mtext(ylabs, side = 2.5, line = 2)}
        if(dim(tmp00)[2] >2){mtext(ylabs, side = 2, line = 1.7, cex = 0.8)}
      }
      lines1=-2.5
      if(dim(tmp00)[2] >3){lines1=-1}
      title(xlab="School labels", outer = TRUE, line = lines1,cex.lab = 1.2)
      on.exit(par(op))
    }





    if( !is.null(group) & sum(names(analyticObject)=="Bootstrap")>0){
      ntp <- length(analyticObject2$ES)
      if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
      Boot.names <-names(analyticObject2$Bootstrap)
      obs.est <- analyticObject2$ES[[trt]][ES_TW,1]
      tmp2 <- as.numeric(analyticObject2$Bootstrap[,grep(ES_TW,grep(trt,Boot.names, ignore.case = T, value = T))])
      pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
      xlabs=paste0("Bootstrap estimates for ",Condname, " ES_",ES_TW)
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main="")
      abline(v=obs.est,col="red",lwd=2,lty=1)
      abline(v=0,col="grey48",lwd=2,lty=1)
      legend("topright",c("Observed Estimate","Zero-Line"),col=c("red","grey48"),bty="n",lty=1,lwd=2)
    }


    if( !is.null(group) & sum(names(analyticObject)=="permES")>0){
      ntp <- ifelse(is.list(analyticObject$ES),length(analyticObject$ES),1)
      if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
      Perm.names <-names(analyticObject2$permES)
      obs.est <- analyticObject2$ES[[trt]][ES_TW,1]
      tmp2 <- as.numeric(analyticObject2$permES[,grep(ES_TW,grep(trt,Perm.names, ignore.case = T, value = T))])
      pvalue <- round(mean(abs(tmp2)> abs(obs.est)),3)
      xlabs=paste0("Permutation values(PermES) based on ",Condname, " ES_",ES_TW)
      hist(tmp2,breaks=30,col="white",border="cornflowerblue",xlab=xlabs,main=paste("P(|PermES| > |ES|)=",pvalue,sep=""))
      abline(v=obs.est,col="red",lwd=2,lty=2)
      abline(v=-obs.est,col="red",lwd=2,lty=2)
      legend("topright",c("Observed Estimate"),col=c("red"),bty="n",lty=2,lwd=2)
    }


    if( !is.null(group) &sum(names(analyticObject)=="ProbES")> 0 ){

      if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
      tmp2 <- analyticObject2$ProbES[[which(trtpos==TRUE)]]
      tmp2.within<- tmp2[, grep("with",names(tmp2), ignore.case = TRUE)]
      tmp2.total <- tmp2[, grep("total",names(tmp2), ignore.case = TRUE)]
      thd0<- regmatches(rownames(tmp2),  gregexpr("[[:digit:]]+\\.*[[:digit:]]*",rownames(tmp2)))
      thd <- as.numeric(unlist(thd0))

      par_original <- par()[c("mar","xpd")]
      op<- par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(thd, tmp2.within, col="blue", type="b",ylim=c(0,max(tmp2)), xlab="Threshold", ylab="Posterior probability",...)
      lines(thd, tmp2.total, col="red", type="b", lty=2)
      legend("topright", legend=c("within", "total"), col=c("blue", "red"), lty=1:2, cex=0.8)
      on.exit(par(op))
      on.exit(par(par_original))

    }

  }

  # error bar for model comparing models
  #-------------------------------------
  if(is.null(names(analyticObject))){

    ltp <- names(analyticObject)
    if(!is.null(ltp)){stop("Specify list of eefAnalytics objects for comparison")}
    if(is.null(group)){stop("Group number must be defined")}
    ntp <- length(analyticObject)
    if(length(modelNames)!= ntp){stop("Names must be equal to the number of eefAnalytics objects")}

    es.mean <- es.lower <- es.upper <- p.name <- var.name <- NULL
    for(k in 1:ntp){
      tmp <- analyticObject[[k]]

      if(tmp$Method=="LM"){
        trtname <-rownames(tmp$ES)
        trtpos <-substr(trtname, (nchar(trtname)-nchar(group)+1),  nchar(trtname))==group
        trt <- trtname[trtpos]

        if(Conditional==TRUE){tmp2 <- as.matrix(tmp$ES)}
        if(Conditional==FALSE){tmp2 <- as.matrix(tmp$Unconditional$ES)}
        trtname <-rownames(tmp$ES)
        if(is.null(trtname)){trtname <-names(tmp$ES)}
        trtpos <-substr(trtname, (nchar(trtname)-nchar(group)+1),  nchar(trtname))==group
        trt <- trtname[trtpos]
        if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
        es.mean1 <-tmp2[trt,1]
        es.lower1 <-tmp2[trt,2]
        es.upper1 <-tmp2[trt,3]
        p.name1 <- rep(modelNames[k],length(es.mean1))
        var.name1 <- rep("Within",length(es.mean1))

      }

      if(tmp$Method=="MLM"){
        trtname <-names(tmp$ES)
        if(is.null(trtname)){trtname <-names(tmp$ES)}
        trtpos <-substr(trtname, (nchar(trtname)-nchar(group)+1),  nchar(trtname))==group
        trt <- trtname[trtpos]
        if( sum((trtname %in%trt))==0){stop("Group must be one of the intervention values")}
        if(Conditional==TRUE) {tmp2 <- tmp$ES[[trt]]}
        if(Conditional==FALSE){tmp2 <- tmp$Unconditional$ES[[trt]]}
        es.mean1 <- tmp2[,1]
        es.lower1 <-tmp2[,2]
        es.upper1 <-tmp2[,3]
        p.name1 <- rep(modelNames[k],length(es.mean1))
        var.name1 <- rownames(tmp2)


      }


      es.mean <- c(es.mean,es.mean1)
      es.lower <- c(es.lower,es.lower1)
      es.upper <- c(es.upper,es.upper1)
      p.name <- c(p.name,p.name1)
      var.name <- c(var.name,var.name1)

    }


    MyData1 <- data.frame(ES=es.mean,LB.95=es.lower,UB.95=es.upper,Variance=var.name,Name=p.name)
    MyData1<- MyData1[order(MyData1$Variance,decreasing = T),]
    MyData1$Anot <- paste0(MyData1$ES, " [", MyData1$LB.95,", ", MyData1$UB.95,"]")
    MyData1$Index <- 1:dim(MyData1)[1]
    MyData1$Xaxis <- (max(MyData1$UB.95))+0.05
    Mybreaks <-  round(c(min(MyData1$LB.95),(min(MyData1$LB.95)+(max(MyData1$UB.95)))/2,max(MyData1$UB.95)),2)
    xlimits <- c(min(min(MyData1$LB.95),0),(max(MyData1$UB.95))+0.4)

    Ann_text <- data.frame(Index = length(MyData1$Variance[MyData1$Variance=="Within"])+0.5,
                           ES = MyData1$Xaxis[1],LB.95=0, UB.95=0,lab = "Text",
                           Variance = factor("Total",levels = c("Within", "Total")))

    #ggplot
    p <- ggplot(data=MyData1, aes(x=ES, y=Name, xmin=LB.95, xmax=UB.95))
    p <- p + geom_point()
    p <- p + geom_errorbarh(height=.1)
        p <- p + scale_x_continuous(limits=xlimits ,breaks=Mybreaks, name=expression(paste("Hedge's ", italic("g"))))
    p <- p + geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)
    if(sum(unique(MyData1$Variance)%in% "Total")>0){p <- p + facet_grid(Variance~., scales= "free", space="free")}
    p <- p + ylab("Models")
    p <- p + theme_bw()
    p <- p + theme(axis.text.y =element_text(color="black"))
    p <- p + theme(text=element_text(size=16, color="black"))
    p <- p + theme(panel.spacing = unit(1, "lines"))
    p <- p + geom_text(aes(x = Xaxis,y = Name, label = Anot),hjust = 0)
    p <- p + geom_text(data = Ann_text, y=Inf,label = "95% CI",hjust = -1.1,vjust = -0.5,size=4,fontface = "bold")
    p <- p + coord_cartesian(clip = "off")
    p <- p + theme(plot.margin = unit(c(30,5,5,5), "point"))
    p

  }

}

