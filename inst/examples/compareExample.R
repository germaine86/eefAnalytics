if(interactive()){

data(mstData)
###############
##### SRT #####
###############

outputSRT <- srtFREQ(Posttest~ Intervention + Prettest,
                     intervention = "Intervention", data = mstData)

outputSRTBoot <- srtFREQ(Posttest~ Intervention + Prettest,
                         intervention = "Intervention",nBoot=1000, data = mstData)

###############
##### MST #####
###############

outputMST <- mstFREQ(Posttest~ Intervention + Prettest,
                     random = "School", intervention = "Intervention", data = mstData)

outputMSTBoot <- mstFREQ(Posttest~ Intervention + Prettest,
                         random = "School", intervention = "Intervention",
                         nBoot = 1000, data = mstData)

##################
##### Bayesian #####
##################

outputSRTbayes <- srtBayes(Posttest~ Intervention + Prettest,
                           intervention = "Intervention",
                           nsim = 2000, data = mstData)

## comparing different results

ComparePlot(list(outputSRT,outputSRTBoot,outputMST,outputMSTBoot,outputSRTbayes),
            modelNames =c("ols", "olsBoot","MLM","MLMBoot","OLSBayes"),group=1)


}
