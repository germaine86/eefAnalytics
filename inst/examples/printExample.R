if(interactive()){

  #### read data
  data(mstData)
  data(crtData)


  ###############
  ##### SRT #####
  ###############

  outputSRT <- srtFREQ(Posttest~ Intervention + Prettest,
                       intervention = "Intervention",data = mstData)
  print(outputSRT)


  outputSRTB <- srtBayes(Posttest~ Intervention + Prettest,
                         intervention = "Intervention", data = mstData)
  print(outputSRTB)

  ###############
  ##### MST #####
  ###############


  #### Random intercepts
  outputMST <- mstFREQ(Posttest~ Intervention + Prettest,
                       random = "School", intervention = "Intervention", data = mstData)
  print(outputMST)



  ###############
  ##### CRT #####
  ###############

  #### Random intercepts
  outputCRT <- crtFREQ(Posttest~ Intervention + Prettest,
                       random = "School", intervention = "Intervention", data = crtData)
  print(outputCRT)

}


