if(interactive()){

#### read data
data(mstData)
data(crtData)


###############
##### SRT #####
###############

##### Bootstrapped

outputSRTBoot <- srtFREQ(Posttest~ Intervention + Prettest,
                         intervention = "Intervention",nBoot=1000, data = mstData)
plot(outputSRTBoot,group=1)

##### Permutation
outputSRTPerm <- srtFREQ(Posttest~ Intervention + Prettest,
                         intervention = "Intervention",nPerm=1000, data = mstData)

plot(outputSRTPerm,group=1)


###############
##### MST #####
###############


#### Random intercepts
outputMST <- mstFREQ(Posttest~ Intervention + Prettest,
                     random = "School", intervention = "Intervention", data = mstData)
plot(outputMST)


#### Bootstrapped
outputMSTBoot <- mstFREQ(Posttest~ Intervention + Prettest,
                         random = "School", intervention = "Intervention",
                         nBoot = 1000, data = mstData)

plot(outputMSTBoot)
plot(outputMSTBoot,group=1)

#### Permutation
outputMSTPerm <- mstFREQ(Posttest~ Intervention + Prettest,
                         random = "School", intervention = "Intervention",
                         nPerm = 1000, data = mstData)
plot(outputMSTPerm)
plot(outputMSTPerm,group=1)



###############
##### CRT #####
###############

#### Random intercepts
outputCRT <- crtFREQ(Posttest~ Intervention + Prettest, random = "School",
                     intervention = "Intervention", data = crtData)
plot(outputCRT)


## Bootstrapped
outputCRTBoot <- crtFREQ(Posttest~ Intervention + Prettest, random = "School",
                         intervention = "Intervention", nBoot = 1000, data = crtData)

plot(outputCRTBoot,group=1)


##Permutation
outputCRTPerm <- crtFREQ(Posttest~ Intervention + Prettest, random = "School",
                         intervention = "Intervention", nPerm = 1000, data = crtData)

plot(outputCRTPerm,group=1)
}
