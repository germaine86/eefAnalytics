# ============================================================
# R Script: All code examples from the eefAnalytics paper
# "The eefAnalytics R package for the Analysis of
#  Randomised Controlled Trials in Education"
# Uwimpuhwe, Singh, Zhang, Einbeck et al.
# ============================================================

# ---- Setup --------------------------------------------------
# Install eefAnalytics from CRAN if not already installed
# install.packages("eefAnalytics")

sink("output_log.txt")        # start capturing
library(eefAnalytics)


# ============================================================
# SECTION: Frequentist Functions
# ============================================================

# ---- srtFREQ ------------------------------------------------
# Simple Randomised Trial (SRT) - analysing mstData as if SRT
output1 <- srtFREQ(
  Posttest ~ Intervention + Prettest,
  intervention = "Intervention",
  data = mstData
)

# Print model summary (conditional and unconditional effect sizes)
output1

# Access specific components
output1$Beta        # beta coefficients
output1$ES          # conditional effect size
output1$Unconditional$ES  # unconditional effect size


# ---- crtFREQ ------------------------------------------------
# Cluster Randomised Trial (CRT)
output2 <- crtFREQ(
  Posttest ~ Intervention + Prettest,
  random = "School",
  intervention = "Intervention",
  data = crtData
)

# Print model summary
output2


# Extract beta coefficients and variance components
output2[c("Beta", "covParm")]

# Extract predicted random intercepts for schools (first 6)
head(output2$SchEffects)


# ---- mstFREQ ------------------------------------------------
# Multisite Trial (MST)
output3 <- mstFREQ(
  Posttest ~ Intervention + Prettest,
  random = "School",
  intervention = "Intervention",
  data = mstData
)


# Extract unconditional effect sizes
output3$Unconditional$ES


# ============================================================
# SECTION: Permutation and Bootstrap Options
# ============================================================

# ---- Bootstrap (srtFREQ with nBoot) -------------------------
outputb <- srtFREQ(
  Posttest ~ Intervention + Prettest,
  intervention = "Intervention",
  nBoot = 1000,
  data = mstData
)

# ---- Permutation (crtFREQ with nPerm) -----------------------
outputp <- crtFREQ(
  Posttest ~ Intervention + Prettest,
  random = "School",
  intervention = "Intervention",
  nPerm = 1000,
  data = crtData
)


# ============================================================
# SECTION: Bayesian Functions
# ============================================================

nsim <- 10000
thd  <- c(0.1, 0.2)

# ---- srtBayes -----------------------------------------------
# SRT Bayesian: model without school variable
output4 <- srtBayes(
  Posttest ~ Intervention + Prettest,
  intervention = "Intervention",
  nsim = nsim,
  data = mstData,
  threshold = thd,
  digits = 2
)

output4$ES                 # conditional ES
output4$Unconditional$ES   # unconditional ES

# SRT Bayesian: model with school as fixed factor
output5 <- srtBayes(
  Posttest ~ Intervention + Prettest + factor(School),
  intervention = "Intervention",
  nsim = nsim,
  data = mstData,
  threshold = thd
)

output5$ES
output5$Unconditional$ES


# ---- crtBayes -----------------------------------------------
# CRT Bayesian: two-arm intervention
thd  <- c(0.1, 0.2)
nsim <- 2000

output6 <- crtBayes(
  Posttest ~ Prettest + Intervention,
  intervention = "Intervention",
  random = "School",
  nsim = nsim,
  data = crtData,
  threshold = thd,
  digits = 2
)

output6$ES   # conditional ES (within and total)

# CRT Bayesian: three-arm intervention
# Check the distribution of the three-arm intervention variable
table(crtData$Intervention2)

thd  <- 0.1
nsim <- 1000

output7 <- crtBayes(
  Posttest ~ Prettest + Intervention2,
  intervention = "Intervention2",
  random = "School",
  nsim = nsim,
  data = crtData,
  threshold = thd,
  digits = 2
)

output7$ES       # conditional ES for each intervention arm

# Posterior probabilities of exceeding threshold
output7$ProbES


# ---- mstBayes -----------------------------------------------
# MST Bayesian: two-arm intervention
thd  <- c(0.1, 0.2)
nsim <- 10000

output8 <- mstBayes(
  Posttest ~ Prettest + Intervention,
  intervention = "Intervention",
  random = "School",
  nsim = nsim,
  data = mstData,
  threshold = thd
)

output8$ES                    # conditional ES
output8$Unconditional         # unconditional ES, variance, posterior probs

# MST Bayesian: three-arm intervention
output9 <- mstBayes(
  Posttest ~ Prettest + Intervention2,
  intervention = "Intervention2",
  random = "School",
  nsim = nsim,
  data = mstData,
  threshold = thd
)

output9$ES                    # conditional ES (both arms)
output9$Unconditional$ES      # unconditional ES (both arms)


# ============================================================
# SECTION: Gain Index
# ============================================================

output10 <- GainIndex(
  data = mstData,
  formula = Posttest ~ Prettest,
  random = "School",
  intervention = "Intervention",
  NA.omit = TRUE,
  alpha = 0.05
)

output10          # GI, proportions table, timing


# ============================================================
# SECTION: Plots
# ============================================================
# Set working directory to folder where figures will be saved
# Replace the path below with your actual folder path.
Working_folder <-  "WHERE_TO_SAVE_FIGURES"
setwd(Working_folder)


# ---- Random intercepts plot ---------------------------------
# Requires a CRT or MST model object
png("Picture 1.png", width = 12, height = 9, units = "cm", res = 300)
plot(output2, group = NULL)
dev.off()

# ---- Bootstrap plot -----------------------------------------
png("Picture 2.png", width = 12, height = 9, units = "cm", res = 300)
plot(outputb, group = 1)
dev.off()

# ---- Permutation plot (p-value readable from plot) ----------
png("Picture 3.png", width = 12, height = 9, units = "cm", res = 300)
plot(outputp, group = 1)
dev.off()

# ---- Posterior probabilities plot ---------------------------
# Re-fit output9 with a finer threshold grid for the plot
thd  <- c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
nsim <- 2000

output9 <- mstBayes(
  Posttest ~ Prettest + Intervention,
  intervention = "Intervention",
  random = "School",
  nsim = nsim,
  data = mstData,
  threshold = thd
)


png("Picture 4.png", width = 12, height = 9, units = "cm", res = 300)
plot(output9, group = 1, Conditional = TRUE)
dev.off()

# ---- Forest plot / sensitivity analysis ---------------------


#original figure 
#-----------------------
Picture5 <- ComparePlot(
  list(output1, output3, output8),
  modelNames = c("ols", "MST", "MSTBayes"),
  group = 1
)
#save Picture 5
ggplot2::ggsave("Picture 5.png", 
                plot = Picture5,   
                width = 10, 
                height = 6, 
                dpi = 300)



# Customised forest plot
#-----------------------
library(ggplot2)
Picture6 <- ComparePlot(
  list(output1, output3, output8),
  modelNames = c("ols", "MST", "MSTBayes"),
  group = 1
) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(name = "Hedges' g") +
  ylab("Models")

#save Picture 6
ggplot2::ggsave("Picture 6.png",
                plot = Picture6, 
                width = 10, 
                height = 6, 
                dpi = 300)
