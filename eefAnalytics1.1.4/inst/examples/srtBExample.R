if(interactive()){

  data(mstData)

  ########################################################
  ## Bayesian analysis of simple randomised trials      ##
  ########################################################

  output <- srtBayes(Posttest~ Intervention+Prettest,
                     alpha = 0.2,
                     digits=4,
                     intervention="Intervention",
                     nsim=10000,
                     data=mstData)

  ### Fixed effects
  beta <- output$Beta
  beta

  ### Effect size
  ES1 <- output$ES
  ES1

  ### Effect size
  ES2 <- output$Unconditional$ES
  ES2

  ## Covariance matrix
  covParm1 <- output$sigma2
  covParm1


  ## Unconditional Covariance matrix
  covParm2 <- output$Unconditional$sigma2
  covParm2


  ## Prob ES
  ProbES1 <- output$ProbES
  ProbES1


  ## Prob  based on Unconditional ES
  ProbES2 <- output$Unconditional$ProbES
  ProbES2


  ### plot posterior probability of an effect size to be bigger than a pre-specified threshold

  plot(output,group=1)


}
