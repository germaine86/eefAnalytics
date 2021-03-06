if(interactive()){

  data(crtData)

  ########################################################
  ## Bayesian analysis of cluster randomised trials     ##
  ########################################################

  output <- crtBayes(Posttest~ Intervention+Prettest,random="School",
                     intervention="Intervention",nsim=2000,data=crtData)

  ### Fixed effects
  beta <- output$Beta
  beta

  ### Effect size
  ES1 <- output$ES
  ES1

  ## Covariance matrix
  covParm <- output$covParm
  covParm

  ### plot random effects for schools

  plot(output)

  ### plot posterior probability of an effect size to be bigger than a pre-specified threshold

  plot(output,group=1)
}
