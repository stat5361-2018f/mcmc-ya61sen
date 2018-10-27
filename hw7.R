## generate a random sample from a mixture normal with true parameter values
delta <- 0.7 # true value to be estimated based on the data
n <- 100
set.seed(123)
u <- rbinom(n, prob = delta, size = 1)
x <- rnorm(n, ifelse(u == 1, 7, 10) , 0.5)
mydata <- data.frame(x)

## posterior density
logpost <- function(theta,data) {
  mu_1 <- theta[1]; mu_2 <- theta[2]; 
  sigma2_1 <- theta[3]; sigma2_2 <- theta[4]; delta <- theta[5]
  x <- data$x
  library("invgamma")
  return(sum(log(delta*dnorm(x,mu_1,sigma2_1^0.5)+(1-delta)*dnorm(x,mu_2,sigma2_2^0.5))) +
           dnorm(mu_1,0,10,log = T) + dnorm(mu_2,0,10,log = T) +
           dinvgamma(sigma2_1,shape=0.5,scale=10,log = T) +
           dinvgamma(sigma2_2,shape=0.5,scale=10,log = T))
}

## MCMC
mymcmc <- function(niter, thetaInit, data, nburn) {
  p <- length(thetaInit)
  thetaCurrent <- thetaInit
  ## define a function for full conditional sampling  
  logFC <- function(th, idx) {
    theta <- thetaCurrent
    theta[idx] <- th
    logpost(theta, data)
  }
  out <- matrix(thetaInit, niter, p, byrow = TRUE)
  ## Gibbs sampling
  for (i in 2:niter) {
    for (j in 1:p) {
      ## general-purpose arms algorithm
      # Indicator function
      indF <- function(x, idx) {
        if (idx==1 | idx==2) {(x<30)*(x>-30)}
        else if (idx==3 | idx==4) {(x<10)*(x>0)}
        else {(x<1)*(x>=0)}
      }
      out[i, j] <- thetaCurrent[j] <-
        HI::arms(thetaCurrent[j], logFC, indF, 1, idx = j)
    }
  }
  out[-(1:nburn), ]
}

## simulation
niter <- 5000; nburn <- 1000
thetaInit <- c(0,0,0.2,0.2,0.5)
sim <- mymcmc(niter, thetaInit, mydata, nburn)
name_para <- c("mu_1","mu_2","sigma2_1","sigma2_2","delta")
for (i in 1:5) {
  plot(ts(sim[,i]),main=(paste(name_para[i])))
  hist(sim[,i],main=(paste(name_para[i])))
}

