######### EXAMPLE 1: crtData #########
load("data/crtData.rda")
output1 <- GainIndex(data = crtData, formula = Posttest~Prettest, random = "School",
                     intervention = "Intervention", NA.omit = T, alpha = 0.05)
output1


########## EXAMPLE 2: mstData ######
load("data/mstData.rda")
output1 <- GainIndex(data = mstData, formula = Posttest~Prettest, random = "School",
                     intervention = "Intervention", NA.omit = T, alpha = 0.05)
output1

