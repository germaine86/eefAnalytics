if(interactive()){

  #### read data
  data(mstData)
  data(crtData)


  ###############
  ##### SRT #####
  ###############

  outputSRT <- srtFREQ(Posttest~ Intervention + Prettest,
                       intervention = "Intervention",data = mstData)
  summary(outputSRT)


  outputSRTB <- srtBayes(Posttest~ Intervention + Prettest,
                         intervention = "Intervention", data = mstData)
  summary(outputSRTB)

  ###############
  ##### MST #####
  ###############


  #### Random intercepts
  outputMST <- mstFREQ(Posttest~ Intervention + Prettest,
                       random = "School", intervention = "Intervention", data = mstData)
  summary(outputMST)



  ###############
  ##### CRT #####
  ###############

  #### Random intercepts
  outputCRT <- crtFREQ(Posttest~ Intervention + Prettest,
                       random = "School", intervention = "Intervention", data = crtData)
  summary(outputCRT)


}

