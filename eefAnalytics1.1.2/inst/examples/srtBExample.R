if(interactive()){

data(mstData)

########################################################
## Bayesian analysis of simple randomised trials      ##
########################################################

output <- srtBayes(Posttest~ Intervention+Prettest,
		intervention="Intervention",nsim=10000,data=mstData)

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

###########################################################################################
## Bayesian analysis of simple randomised trials using informative priors for treatment  ##
###########################################################################################

### define priors for explanatory variables

my_prior <- normal(location = c(0,6), scale = c(10,1))

### specify the priors for the conditional model only

output2 <- srtBayes(Posttest~ Prettest+Intervention,
                    intervention="Intervention",
                    nsim=2000,data=mstData,
                    condopt=list(prior=my_prior))

### Fixed effects
beta2 <- output2$Beta
beta2

### Effect size
ES2 <- output2$ES
ES2
}
