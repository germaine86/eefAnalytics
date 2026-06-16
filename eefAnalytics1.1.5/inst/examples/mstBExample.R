if(interactive()){

  data(mstData)

  ########################################################
  ## Bayesian analysis of multisite randomised trials   ##
  ########################################################

  output <- mstBayes(formula = Posttest ~ Prettest + Intervention,
                     random = "School",
                     intervention = "Intervention",
                     alpha = 0.05,
                     digits = 3,
                     nsim = 10000,
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

}
