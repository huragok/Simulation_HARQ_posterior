library(MVN)
library(XML)

# Import the data
data <- xmlParse("samples_4_1.xml")
xml.data <- xmlToList(data)

N <- as.numeric(xml.data$.attrs[["N"]])
beta <- as.numeric(xml.data$.attrs[["beta"]])
n.success <- as.numeric(xml.data$.attrs[["nSuccess"]])
n.failure <- as.numeric(xml.data$.attrs[["nFailure"]])
n.independent.channel <- as.numeric(xml.data$.attrs[["Nprb"]])
sigma2 <- as.numeric(xml.data$.attrs[["sigma2"]])

samples.success <- matrix(unlist(lapply(xml.data$Success, function(entry) as.numeric(unlist(strsplit(entry, "\\s+"))))), ncol = 8 * N, byrow = TRUE)
samples.failure <- matrix(unlist(lapply(xml.data$Failure, function(entry) as.numeric(unlist(strsplit(entry, "\\s+"))))), ncol = 8 * N, byrow = TRUE)
#samples.complete <- rbind(samples.success, samples.failure)

# MVN Tests (Korkmaz, Selcuk, Dincer Goksuluk, and Gokmen Zararsiz. "MVN: An R Package for Assessing Multivariate Normality." A peer-reviewed, open-access publication of the R Foundation for Statistical Computing (2014): 151.)
## Mardia's test
pdf(paste("qqplot/qqplotMardia_", N, "_", n.independent.channel, ".pdf", sep = ""))
result.failure.mardia <- mardiaTest(samples.failure, qqplot = TRUE)
dev.off()
#result.complete.mardia <- mardiaTest(samples.complete, qqplot = TRUE) # For comparison

## Henze-Zirkler's MVN test
pdf(paste("qqplot/qqplotHz_", N, "_", n.independent.channel, ".pdf", sep = ""))
result.failure.hz <- hzTest(samples.failure, qqplot = TRUE)
dev.off()
#result.complete.hz <- hzTest(samples.complete, qqplot = TRUE)  # For comparison

## Royston's MVN test
if (n.success >= 2000 && n.failure >= 2000) {
    pdf(paste("qqplot/qqplotRoyston_", N, "_", n.independent.channel,".pdf", sep = ""))
    result.failure.royston <- roystonTest(samples.failure[1:2000,], qqplot = TRUE)
    dev.off()
    #result.complete.royston <- roystonTest(samples.complete[1:2000, ], qqplot = TRUE)  # For comparison
}

# Log-likelihood test on the parameters based on Wilks's theorem
paramTest <- function(data, meanTarget, covTarget, sl = 0.05) {
    n <- nrow(data)
    p <- ncol(data)
    
    mu <- colMeans(data) # ML estimation of the mean
    Sigma <- cov(data) * (n - 1) / n # ML estimation of the covariance
    
    data.minusMeanTarget <- data - matrix(rep(1, n) %x% meanTarget, ncol = p, byrow = TRUE) # Subtract meanTarget from all samples
    W <- -n * (log(det(Sigma)) - log(det(covTarget))) - n * p + sum(diag(solve(covTarget, t(data.minusMeanTarget) %*% data.minusMeanTarget)))# The test statistic -2 * log(Lambda) 
    
    result <- list();
    result$name <- "Log-likelihood mean and covariance test"
    result$p.value <- 1 - pchisq(W, p * (p + 1) / 2 + p)
    if (result$p.value < sl) result$Result <- "Mean and covariance of the samples do not match the target." # Reject
    else result$Result <- "Mean and covariance of the samples match the target." # Accept
    result$mean.ML <- mu
    result$cov.ML <- Sigma
    
    return (result)
}

meanTarget <- rep(0, 8 * N)
covTarget <-  diag(c(rep(beta / 2, times = 4 * N), rep(sigma2 / 2, times = 4 * N)))
result.failure.param <- paramTest(samples.failure, meanTarget, covTarget)
#result.complete.param <- paramTest(samples.complete, meanTarget, covTarget)

# Test whether the posteriori channel is a scaled version of the prior distribution using Wilks's theorem
paramScaledTest <- function(data, sl = 0.05) {
    n <- nrow(data)
    p <- ncol(data)
    
    Sigma <- cov(data) * (n - 1) / n # ML estimation of the covariance
    a.h <- 4 * sum(data[,1:(p/4)] * data[,1:(p/4)]) / n / p
    a.g <- 4 * sum(data[,(p/4+1):(p/2)] * data[,(p/4+1):(p/2)]) / n / p
    a.R <- 4 * sum(data[,(p/2+1):(3*p/4)] * data[,(p/2+1):(3*p/4)]) / n / p
    a.D <- 4 * sum(data[,(3*p/4+1):p] * data[,(3*p/4+1):p]) / n / p
    
    W <- -n * log(det(Sigma)) + n * p / 4 * (log(a.h) + log(a.g) + log(a.R) + log(a.D)) # The test statistic -2 * log(Lambda) 
    
    result <- list();
    result$name <- "Log-likelihood scaled Gaussian test"
    result$p.value <- 1 - pchisq(W, p * (p + 1) / 2 + p - 4)
    if (result$p.value > sl) result$Result <- "Means are 0 and covariance is a scaled identity matrix." # Reject
    else result$Result <- "Means are not 0 or covariance is not a scaled identity matrix." # Accept
    result$beta.h.ML <- 2 * a.h
    result$beta.g.ML <- 2 * a.g
    result$sigma2.n.R.ML <- 2 * a.R
    result$sigma2.n.D.ML <- 2 * a.D
    
    return (result)
}
result.failure.paramScaled <- paramScaledTest(samples.failure)
#result.complete.paramScaled <- paramScaledTest(samples.complete, meanTarget, covTarget)

# Output
cat("N =", N, "n.independent.channel =", n.independent.channel, ", n.failure =", n.failure, ", n.success =", n.success, ", beta =", beta, ", sigma2 =", sigma2, "\n")
print("MVN distribution tests ********************************************")
print(result.failure.mardia)
print(result.failure.hz)
print(result.failure.royston)
print("Parameter tests ***************************************************")
print(result.failure.param)
print("Parameter scaled tests ********************************************")
print(result.failure.paramScaled)
