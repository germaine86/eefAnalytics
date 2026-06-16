if(interactive()){

  data(mstData)

  ########################################################
  ## Bayesian analysis of multisite randomised trials   ##
  ########################################################

  output <- mstBayes(formula = Posttest ~ Prettest + Intervention,
                     random = "School",
                     intervention = "Intervention",
                     nSim = 10000,
                     data = mstData)
  output

  ### Fixed effects
  beta <- output$Beta
  beta

  ### Effect size
  ES1 <- output$ES
  ES1

  ## Covariance matrix
  covParm <- output$covParm
  covParm

  ## Prob ES
  ProbES <- output$ProbES
  ProbES

  ## Unconditional
  Unconditional <- output$Unconditional
  Unconditional

  ## Random Effect
  randomEffects <- output$SchEffects
  randomEffects


  ### plot random effects for schools

  plot(output)

  ### plot posterior probability of an effect size to be bigger than a pre-specified threshold

  plot(output,group=1)


  #############################################################################################
  ## Bayesian analysis of multisite randomised trials using informative priors for treatment ##
  #############################################################################################

  ### define priors for explanatory variables

  my_prior <- normal(location = c(0,6), scale = c(10,1))

  ### specify the priors for the conditional model only

  output2 <- mstBayes(Posttest~ Prettest+Intervention,random="School",
                      intervention="Intervention",nsim=2000,data=mstData,
                      condopt=list(prior=my_prior))

  ### Fixed effects
  beta2 <- output2$Beta
  beta2

  ### Effect size
  ES2 <- output2$ES
  ES2
}
