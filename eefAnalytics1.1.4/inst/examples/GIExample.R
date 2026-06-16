######### EXAMPLE ONE: crtData #########
\dontrun{
  data(crtData)
  output1 <- GainIndex(data = crtData, formula = Posttest~Prettest, random = "School",n.iter = 200,
                       intervention = "Intervention", NA.omit = T, alpha = 0.05)
  output1


  ########## EXAMPLE TWO: mstData ######
  data(mstData)
  output1 <- GainIndex(data = mstData, formula = Posttest~Prettest, random = "School",n.iter = 200,
                       intervention = "Intervention", NA.omit = T, alpha = 0.05)
  output1
}
