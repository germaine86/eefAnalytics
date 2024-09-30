#' Calculate the Gain Index (GI) using JAGS
#'
#' This function computes the Gain Index and other related statistics for educational trials.
#' Gain index provides a proportion of pupils who would not have make good progress without intervention.
#' This function supports flexible configurations for JAGS modeling.
#'
#' @param data A list containing the data for the JAGS model which must include
#'        columns: School, Posttest, Pretest, Intervention. Data should not have
#'        any missing values in these columns.
#' @param formula the model to be analysed is of the form y ~ x1+x2+.... Where y is the outcome variable and Xs are the independent variables. Formula does not need to include `Intervention` variable.
#' @param random a string variable specifying the "clustering variable" as contained in the data. See example below.
#' @param intervention a string variable specifying the "intervention variable" as appearing in the formula and the data. See example below.
#' @param NA.omit Optional; a logic to check if omitting missing value.
#'        If NA.omit = TRUE, results will output the percentage of missing value in the four required columns and then JAGS results.
#'        If NA.omit = FALSE, will give a warning "Please handle missing values before using GainIndex()."
#'        If not provided, the function uses default TRUE.
#' @param inits Optional; a list of initial values for the JAGS model. If NULL,
#'        the function generates default initial values.
#' @param n.iter Total number of iterations for the MCMC simulation.
#' @param n.burnin Number of burn-in iterations to be discarded before analysis.
#' @param n.chains Number of chains to run in the MCMC simulation.
#' @param model.file Optional; a custom path to the JAGS model file.
#'        If not provided, the function uses default path.
#' @param alpha significant level, default alpha = 0.05.
#'
#' @return An S3 object containing the following components:
#' \describe{
#'   \item{GI}{A data frame containing the Gain Index and its 95\% confidence intervals,
#'             as well as the Progress Index and its 95\% confidence intervals.}
#'   \item{Proportions}{A data frame showing the proportion of participants achieving
#'                      each level of gain (low and high) for both control and
#'                      intervention groups.}
#'   \item{Timing}{A vector with execution time details, including user and elapsed
#'                 time in seconds.}
#' }
#'
#' @example inst/examples/GIExample.R
#'
#' @export
#'
# GainIndex() function ----
# Include validation and preparation steps at the beginning of the GainIndex() function
GainIndex <- function(data, formula, random, intervention, NA.omit=TRUE, n.iter = 20000, n.chains = 3, inits = NULL, model.file = NULL, alpha = 0.05) {
  requireNamespace("R2jags", quietly = TRUE) || stop("Please install the 'R2jags' package.")
  require(R2jags)

  group = 2


  # Prepare data
  if (group == 2) {
    prepared_data <- prepare_data1(data, formula, random, intervention, NA.omit=NA.omit)
    model.file <- "inst/jags/modelJags_2groups.txt"

  } else if (group == 3) {
    prepared_data <- prepare_data2(data, formula, random, intervention, NA.omit=NA.omit)

    model.file <- "inst/jags/modelJags_3groups.txt"
  } else {
    stop("Invalid group number. Please select either 2 or 3.")
  }

  # Default initial values if none provided
  if (is.null(inits)) {
    # Adapt initial values based on Group
    if (group == 2) {
      jags.inits <- list(
        list(beta=c(-0.1,0.1), tau=1, tau.b=1, pi=c(0.8,0.2), b=rnorm(prepared_data$ns, 0, 1), .RNG.seed=1),
        list(beta=c(-0.5,0.70), tau=0.2, tau.b=0.4, pi=c(0.7,0.3), b=rnorm(prepared_data$ns, -0.5, 1.4), .RNG.seed=2),
        list(beta=c(-0.2,0.35), tau=0.8, tau.b=0.3, pi=c(0.4,0.5), b=rnorm(prepared_data$ns, 0.2, 5), .RNG.seed=3)
      )
    } else {
      jags.inits <- list(
        list(beta=c(-0.4,0,0.4), tau=1, tau.b=1, pi=c(0.6,0.2,0.2), b=rnorm(prepared_data$ns, 0, 1), .RNG.seed=1),
        list(beta=c(-0.5,0,0.70), tau=0.2, tau.b=0.4, pi=c(0.5,0.1,0.4), b=rnorm(prepared_data$ns, -0.5, 1.4), .RNG.seed=2),
        list(beta=c(-0.2,0,0.35), tau=0.8, tau.b=0.3, pi=c(0.4,0.1,0.4), b=rnorm(prepared_data$ns, 0.2, 5), .RNG.seed=3)
      )
    }
  }

  # Monitor the variables of interest
  vars.to.monitor <- if (group == 2) {
    c("beta", "sigma", "sigma.b", "pi", "b", "Ti", "T2", "PGI", "GI")
  } else {
    c("beta", "sigma", "sigma.b", "pi", "b", "Ti", "T3", "PGI", "GI")
  }

  timing <- system.time({
    jags.samples <- jags(data=prepared_data,
                         parameters.to.save=vars.to.monitor,
                         n.burnin=n.iter/2,
                         inits=jags.inits,
                         n.chains=n.chains,
                         n.iter=n.iter,
                         model.file=model.file)
  })

  results <- jags.samples$BUGSoutput$summary
  Result <- data.frame(results)


  # Calculate Gain Index & Progress Index

  GI <- Result[grep("GI", rownames(Result)), c("mean","X2.5.","X97.5.")][1,]
  percent <- paste0((1 - alpha) * 100, "%")
  colnames(GI) <- c("Estimate", paste0(percent, " LB"), paste0(percent, " UB"))
  #colnames(GI) <- c("Estimate","(1-alpha)% LB","(1-alpha)% UB")

  mcmc <- as.mcmc(jags.samples)
  GI_iteration <- rowMeans(sapply(mcmc, function(x) data.frame(x)$`PGI.2.` - data.frame(x)$`PGI.1.`)) # Means of 3 chains in each iteration (in total of 1,000)
  GI[1,2] <- quantile(GI_iteration, probs = alpha/2) ## 2.5%
  GI[1,3] <- quantile(GI_iteration, probs = 1 - alpha/2) ## 97.5%
  ## 95% CI
  ### give options to calculate CI ---- option: significant level alpha
  ## probs = alpha/2 OR 1 - alpha/2

  #GI[2,1] <- (Result[grep("PGI", rownames(Result)), c("mean")][2] / Result[grep("PGI", rownames(Result)), c("mean")][1] - 1)*100

  PI_iteration <- rowMeans(sapply(mcmc, function(x) data.frame(x)$`PGI.2.`/data.frame(x)$`PGI.1.` - 1)) # Means of 3 chains in each iteration (in total of 1,000)
  # paired PI
  #
  GI[2,1] <- mean(PI_iteration)*100
  GI[2,2] <- quantile(PI_iteration, probs = alpha/2)*100 ## 2.5%
  GI[2,3] <- quantile(PI_iteration, probs = 1 - alpha/2)*100 ## 97.5%
  rownames(GI) <- c("Gain Index","Progress Index")
  GI <- round(GI,2)
  GI[2,] <- paste(GI[2,],"%",sep = "")

  # Create Table_mcmc based on Group
  if (group == 2) {
    # To calculate and format Table_mcmc for 2 groups
    mcmc <- as.mcmc(jags.samples)[[1]]
    mcmc0 <- mcmc[,grep("Ti",rownames(Result),value = T)]
    Table_mcmc0 <- apply(mcmc0, 1, function(x) prop.table(table(x, prepared_data$t), 2))
    Prob_mcmc <- rowMeans(Table_mcmc0)
    names(Prob_mcmc) <- c("01", "02", "11", "12")
    Table_mcmc_1 <- data.frame(round(rbind(Prob_mcmc[1:2] * prepared_data$nc, Prob_mcmc[3:4] * prepared_data$nt)))
    rownames(Table_mcmc_1) <- c("Control", "Intervention")
    names(Table_mcmc_1) <- c("-ve gain", "+ve gain")
    Table_mcmc_2a <- Table_mcmc_1["Intervention", "+ve gain"] / sum(Table_mcmc_1["Intervention", ])
    Table_mcmc_2b <- Table_mcmc_1["Control", "+ve gain"] / sum(Table_mcmc_1["Control", ])
    Table_mcmc <- cbind(Table_mcmc_1, round(rbind(Table_mcmc_2b, Table_mcmc_2a), 2))
    colnames(Table_mcmc) <- c("-ve gain", "+ve gain", " proportion of +ve gain")

  } else if (group == 3) {
    # To calculate and format Table_mcmc for 3 groups
    mcmc0 <- mcmc[,grep("Ti",rownames(Result),value = T)]
    Table_mcmc0 <- apply(mcmc0[[1]], 1, function(x) prop.table(table(x, prepared_data$t), 2))
    # ROW 1-3 : Control for each group (T1-T3)
    # ROW 4-6: Intervention for each group (T1-T3)
    Prob_mcmc <- rowMeans(Table_mcmc0)
    names(Prob_mcmc) <- c("01", "02", "03", "11", "12","13")
    Table_mcmc_1 <- data.frame(round(rbind(Prob_mcmc[1:3] * prepared_data$nc, Prob_mcmc[4:6] * prepared_data$nt)))
    rownames(Table_mcmc_1) <- c("Control", "Intervention")
    names(Table_mcmc_1) <- c("low gain","middle gain", "good gain")
    Table_mcmc_2c <- Table_mcmc_1["Intervention", "middle gain"] / sum(Table_mcmc_1["Intervention", ])
    Table_mcmc_2d <- Table_mcmc_1["Control", "middle gain"] / sum(Table_mcmc_1["Control", ])
    Table_mcmc <- cbind(Table_mcmc_1, rbind(Table_mcmc_2d, Table_mcmc_2c))
    Table_mcmc_2a <- Table_mcmc_1["Intervention", "good gain"] / sum(Table_mcmc_1["Intervention", ])
    Table_mcmc_2b <- Table_mcmc_1["Control", "good gain"] / sum(Table_mcmc_1["Control", ])
    Table_mcmc <- cbind(Table_mcmc, rbind(Table_mcmc_2b, Table_mcmc_2a))
    colnames(Table_mcmc) <- c("low gain","middle gain", "good gain", " proportion of middle gain", " proportion of good gain")
    Table_mcmc$` proportion of middle and good gain` <- rowSums(Table_mcmc_1[, c("middle gain","good gain")]) / rowSums(Table_mcmc_1)
  }


  # Prepare output
  output <- list(GI = GI, Proportions = round(Table_mcmc, 2), Timing = timing)
  return(output)
}


# Validate input data ----
validate_data <- function(data, formula, random, intervention, NA.omit=TRUE) {
  required_columns <- c("School", all.vars(formula)[1], all.vars(formula)[-1], "Intervention")
  data <- data.frame(data)
  data$School <- data[,random]
  data$Intervention <- data[,intervention]

  missing_cols <- setdiff(required_columns, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse=", "))
  }
  if (!all(sapply(data[,all.vars(formula)], is.numeric))) {
    stop("Columns 'Posttest' and 'Pretest' must be numeric")
  }
  data1 <- data[,required_columns]
  if (anyNA(data1)) {
    print(paste0("Data contains ", round((nrow(data1)-nrow(na.omit(data1)))/nrow(data1)*100,2),  "% of rows with NA values in required columns."))

    if (NA.omit==TRUE) { data1 <- na.omit(data1)}
    else {stop("Please handle missing values before using GainIndex().")}
  }
  return(data1)
}

# Data preparation function ----
# Convert prepared data into a format suitable for JAGS
# Prepared for 2 groups setting ----
prepare_data1 <- function(data, formula, random, intervention, NA.omit = T, ...) {
  data=validate_data(data, formula, random, intervention, NA.omit)
  data$School <- as.numeric(factor(data$School, levels = unique(data$School)))
  data <- data[order(data$School), ]
  data$res <- residuals(lm(formula, data = data))

  # Preparing Ti vector
  # initialise the first and last positions specially
  Ti <- rep(NA, length(data$res))
  Ti[1] <- 1
  Ti[length(Ti)] <- 2

  # Return the prepared data as a list
  return(list(
    School = data$School,
    t = data$Intervention,
    y = data$res,  # Use data$res as y for consistency
    N = length(data$School),
    ns = length(unique(data$School)),
    nt = sum(data$Intervention == 1),
    nc = sum(data$Intervention == 0),
    LB = min(data$res),
    UB = max(data$res),
    Ti = Ti,
    eta = rep(1, 2)  # Assuming eta needs to be a vector of length 2
  ))
}

# Prepared for 3 groups setting ----
prepare_data2 <- function(data, formula,NA.omit=T,...) {
  data=validate_data(data, formula, random, intervention, NA.omit)
  data$School <- as.numeric(factor(data$School, levels = unique(data$School)))
  data <- data[order(data$School), ]
  data$res <- residuals(lm(formula, data = data))

  # Preparing Ti vector
  # initialize the first, second and last positions specially (3 groups)
  Ti <- rep(NA, length(data$res))
  Ti[1] <- 1
  Ti[2] <- 2
  Ti[length(Ti)] <- 3


  # Return the prepared data as a list
  return(list(
    School = data$School,
    t = data$Intervention,
    y = data$res,  # Use data$res as y for consistency
    N = length(data$School),
    ns = length(unique(data$School)),
    nt = sum(data$Intervention == 1),
    nc = sum(data$Intervention == 0),
    LB = min(data$res),
    UB = max(data$res),
    Ti = Ti,
    eta = rep(1, 3)  # Assuming eta needs to be a vector of length 3
  ))
}
