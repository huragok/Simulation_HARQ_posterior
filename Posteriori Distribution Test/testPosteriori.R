library(MVN)
library(XML)

# The filenames of the datasets
filenames <- c("samples_1_400.xml", "samples_2_400.xml", "samples_3_400.xml", "samples_4_400.xml",
               "samples_1_100.xml", "samples_2_100.xml", "samples_3_100.xml", "samples_4_100.xml",
               "samples_1_10.xml", "samples_2_10.xml", "samples_3_10.xml", "samples_4_10.xml",
               "samples_1_1.xml", "samples_2_1.xml", "samples_3_1.xml", "samples_4_1.xml")

filename.output.settings <- "settings.csv"
filename.output.results.failure <- "results_failure.csv"
filename.output.results.success <- "results_success.csv"

n.files <- length(filenames) # Number of datasets to test

# Initialization
index <- seq(n.files)
n.independent.channels <- rep(0, n.files)
m <- rep(0, n.files)
sigma2.actual <- rep(0, n.files)

n.failure <- rep(0, n.files)
p.value.mardia.skew.failure <- rep(0, n.files)
p.value.mardia.kurt.failure <- rep(0, n.files)
p.value.hz.failure <- rep(0, n.files)
p.value.royston.failure <- rep(0, n.files)
p.value.param.failure <- rep(0, n.files)
p.value.param.relaxed.failure <- rep(0, n.files)
beta.ML.failure <- rep(0, n.files)
K.ML.failure <- rep(0, n.files)
theta.ML.failure <- rep(0, n.files)
sigma2.ML.failure <- rep(0, n.files)

n.success <- rep(0, n.files)
p.value.mardia.skew.success <- rep(0, n.files)
p.value.mardia.kurt.success <- rep(0, n.files)
p.value.hz.success <- rep(0, n.files)
p.value.royston.success <- rep(0, n.files)
p.value.param.success <- rep(0, n.files)
p.value.param.relaxed.success <- rep(0, n.files)
beta.ML.success <- rep(0, n.files)
K.ML.success <- rep(0, n.files)
theta.ML.success <- rep(0, n.files)
sigma2.ML.success <- rep(0, n.files)

# Log-likelihood test on the parameters based on Wilks's theorem
paramTest <- function(data, meanTarget, covTarget, sl = 0.01) {
    n <- nrow(data)
    p <- ncol(data)
    
    mu <- colMeans(data) # ML estimation of the mean
    Sigma <- cov(data) * (n - 1) / n # ML estimation of the covariance
    
    data.minusMeanTarget <- data - matrix(rep(1, n) %x% meanTarget, ncol = p, byrow = TRUE) # Subtract meanTarget from all samples
    W <- -n * (log(det(Sigma)) - log(det(covTarget))) - n * p + sum(diag(solve(covTarget, t(data.minusMeanTarget) %*% data.minusMeanTarget)))# The test statistic -2 * log(Lambda) 
    
    result <- list();
    result$name <- "Log-likelihood mean and covariance test"
    result$p.value <- 1 - pchisq(W, p * (p + 1) / 2 + p)
    if (result$p.value < sl) result$Result <- "H0 rejected: mean and covariance of the samples do not match the target." # Reject
    else result$Result <- "H0 accepted: mean and covariance of the samples match the target." # Accept
    #result$mean.ML <- mu
    #result$cov.ML <- Sigma
    
    return (result)
}

# Test whether the posteriori channel is of the same type as the prior distribution but with different parameters (beta, K, theta, sigma2) using Wilks's theorem
paramRelaxedTest <- function(data, sl = 0.01) {
    n <- nrow(data)
    p <- ncol(data)
    
    Sigma <- cov(data) * (n - 1) / n # ML estimation of the covariance
    
    mean.real.h <- 4 * sum(data[,1:(p/4)]) / n / p
    mean.imag.h <- 4 * sum(data[,(p/4+1):(p/2)]) / n / p
    real.h.minusMean <-data[,1:(p/4)] - matrix(rep(1, n * p / 4) %x% mean.real.h, ncol = p / 4, byrow = TRUE)   
    imag.h.minusMean <-data[,(p/4+1):(p/2)] - matrix(rep(1, n * p / 4) %x% mean.imag.h, ncol = p / 4, byrow = TRUE)
    var.half.h <- 2 * (sum(real.h.minusMean * real.h.minusMean) + sum(imag.h.minusMean * imag.h.minusMean)) / n / p;
    var.noise <- 4 * sum(data[,(p/2+1):p] * data[,(p/2+1):p]) / n / p 

    W <- -n * log(det(Sigma)) + n * p / 2 * (log(var.half.h) + log(var.noise / 2)) # The test statistic -2 * log(Lambda) 
    
    result <- list();
    result$name <- "Log-likelihood relaxed Gaussian test"
    result$p.value <- 1 - pchisq(W, p * (p + 1) / 2 + p - 4)
    if (result$p.value > sl) result$Result <- "H0 accepted: channel and noise follow Rician and CSCG distribution respectively." # Accept
    else result$Result <- "H0 rejected: channel and noise do not follows Rician and CSCG distribution respectively." # Reject
    
    result$K.ML <- (mean.real.h ^ 2 + mean.imag.h ^ 2) / (2 * var.half.h)
    result$beta.ML <- mean.real.h ^ 2 + mean.imag.h ^ 2 + 2 * var.half.h
    result$theta.ML <- atan(mean.imag.h / mean.real.h)
    result$sigma2.ML <- var.noise
    
    return (result)
}

for (i in seq(n.files)) {
    # Import the data
    data <- xmlParse(filenames[i])
    xml.data <- xmlToList(data)

    N <- as.numeric(xml.data$.attrs[["N"]])
    beta <- as.numeric(xml.data$.attrs[["beta"]])
    K <- as.numeric(xml.data$.attrs[["K"]])
    theta <- as.numeric(xml.data$.attrs[["theta"]])
    n.s <- as.numeric(xml.data$.attrs[["nSuccess"]])
    n.f <- as.numeric(xml.data$.attrs[["nFailure"]])
    n.ic <- as.numeric(xml.data$.attrs[["Nprb"]])
    sigma2 <- as.numeric(xml.data$.attrs[["sigma2"]])

    samples.success <- matrix(unlist(lapply(xml.data$Success, function(entry) as.numeric(unlist(strsplit(entry, "\\s+"))))), ncol = 4 * N, byrow = TRUE)
    samples.failure <- matrix(unlist(lapply(xml.data$Failure, function(entry) as.numeric(unlist(strsplit(entry, "\\s+"))))), ncol = 4 * N, byrow = TRUE)

    # MVN Tests (Korkmaz, Selcuk, Dincer Goksuluk, and Gokmen Zararsiz. "MVN: An R Package for Assessing Multivariate Normality." A peer-reviewed, open-access publication of the R Foundation for Statistical Computing (2014): 151.)
    ## Mardia's test
    result.failure.mardia <- mardiaTest(samples.failure)
    result.success.mardia <- mardiaTest(samples.success)

    ## Henze-Zirkler's MVN test
    result.failure.hz <- hzTest(samples.failure)
    result.success.hz <- hzTest(samples.success)
    
    ## Royston's MVN test
    if (n.s >= 2000 && n.f >= 2000) {
        result.failure.royston <- roystonTest(samples.failure[1:2000,])
        result.success.royston <- roystonTest(samples.success[1:2000,])
    }

    ## Parameter match test
    meanTarget <- c(rep(sqrt(K * beta / (K + 1)) * cos(theta), N), rep(sqrt(K * beta / (K + 1)) * sin(theta), N), rep(0, 2 * N))
    covTarget <-  diag(c(rep(beta / 2 / (K + 1), times = 2 * N), rep(sigma2 / 2, times = 2 * N)))
    result.failure.param <- paramTest(samples.failure, meanTarget, covTarget)
    result.success.param <- paramTest(samples.success, meanTarget, covTarget)

    ## Relaxed parameter match test
    result.failure.param.relaxed <- paramRelaxedTest(samples.failure)
    result.success.param.relaxed <- paramRelaxedTest(samples.success)

    # Save the test results
    n.independent.channels[i] <- n.ic
    m[i] <- N - 1
    sigma2.actual[i] <- sigma2
    
    n.failure[i] <- n.f
    p.value.mardia.skew.failure[i] <- result.failure.mardia@p.value.skew
    p.value.mardia.kurt.failure[i] <- result.failure.mardia@p.value.kurt
    p.value.hz.failure[i] <-  result.failure.hz@p.value
    p.value.royston.failure[i] <- result.failure.royston@p.value
    p.value.param.failure[i] <- result.failure.param$p.value
    p.value.param.relaxed.failure[i] <- result.failure.param.relaxed$p.value
    beta.ML.failure[i] <- result.failure.param.relaxed$beta.ML
    K.ML.failure[i] <- result.failure.param.relaxed$K.ML
    theta.ML.failure[i] <- result.failure.param.relaxed$theta.ML
    sigma2.ML.failure[i] <- result.failure.param.relaxed$sigma2.ML
    
    n.success[i] <- n.s
    p.value.mardia.skew.success[i] <- result.success.mardia@p.value.skew
    p.value.mardia.kurt.success[i] <- result.success.mardia@p.value.kurt
    p.value.hz.success[i] <-  result.success.hz@p.value
    p.value.royston.success[i] <- result.success.royston@p.value
    p.value.param.success[i] <- result.success.param$p.value
    p.value.param.relaxed.success[i] <- result.success.param.relaxed$p.value
    beta.ML.success[i] <- result.success.param.relaxed$beta.ML
    K.ML.success[i] <- result.success.param.relaxed$K.ML
    theta.ML.success[i] <- result.success.param.relaxed$theta.ML
    sigma2.ML.success[i] <- result.success.param.relaxed$sigma2.ML
       
    print(paste(filenames[i], "tested..."))
    flush.console()
}

# Write the simulation results to csv files
## Simulation settings
print("Writing results to .csv files...")
flush.console()
settings <- data.frame(index, n.independent.channels, m, n.failure, n.success, sigma2.actual)
write.table(settings, sep=",", row.names = F, file = filename.output.settings)
print(paste("Simulation settings written to", filename.output.settings))
flush.console()

## Hypothesis tests on the failure samples
results.failure <- data.frame(index, p.value.mardia.skew.failure, p.value.mardia.kurt.failure, p.value.hz.failure, p.value.royston.failure, p.value.param.failure, p.value.param.relaxed.failure, beta.ML.failure, K.ML.failure, theta.ML.failure, sigma2.ML.failure)
write.table(results.failure, sep=",", row.names = F, file = filename.output.results.failure)
print(paste("Results on failure samples written to", filename.output.results.failure))
flush.console()

## Hypothesis tests on the success samples
results.success <- data.frame(index, p.value.mardia.skew.success, p.value.mardia.kurt.success, p.value.hz.success, p.value.royston.success, p.value.param.success, p.value.param.relaxed.success, beta.ML.success, K.ML.success, theta.ML.success, sigma2.ML.success)
write.table(results.success, sep=",", row.names = F, file = filename.output.results.success)
print(paste("Results on success samples written to", filename.output.results.success))

print("Done.")
