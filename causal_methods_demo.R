### CBPS uses exact condition (in contrast to the previous one using overidentified condition)

rm(list=ls())

library(nnet)
library(CBPS)
library(multilevelMatching)
library(Matching)
source("~/function_weighting.R")

n <- 1500
ntrt <- 3 # number of treatment groups
npair <- choose(ntrt, 2)

### coefficients for propensity score model
beta1 <- matrix(rep(0,4))
beta2 <- matrix(c(log(0.6), 0.1, 0.1, 0.2)) # intercept determines treatment prevalence
beta3 <- matrix(c(log(0.8), 0.1, 0.2, 0.1))

# coefficients for outcome model
alpha1 <- matrix(c(log(0.05), 0.5, 0.2, -0.3)) # change the intercept will result in the change of the prevalence of event
alpha2 <- matrix(c(log(0.2), 0.5, 0.3, -0.3))
alpha3 <- matrix(c(log(0.3), 0.4, 0.2, 0.2))


set.seed(1)

X1 <- rnorm(n,0,1)
X2 <- rnorm(n,0,1)
X3 <- rbinom(n,1,0.5)
X <- cbind(1,X1,X2,X3)
exb1 <- exp(X%*%beta1)
exb2 <- exp(X%*%beta2)
exb3 <- exp(X%*%beta3)
pi1 <- as.numeric(exb1/(exb1+exb2+exb3)) # true prob of receiving treatment 1
pi2 <- as.numeric(exb2/(exb1+exb2+exb3))
pi3 <- as.numeric(exb3/(exb1+exb2+exb3))
pi <- cbind(pi1,pi2,pi3)

# multinomial draws for treatment assignment
Z.mat <- t(apply(pi, 1, rmultinom, n = 1, size = 1)) # treatment
Z <- apply(Z.mat, 1, function(x) which(x==1))
Z <- as.factor(Z)

# potential outcomes
pY1 <- 1/(1+exp(-X %*% alpha1))
pY2 <- 1/(1+exp(-X %*% alpha2))
pY3 <- 1/(1+exp(-X %*% alpha3))
Y1 <- rbinom(n, 1, pY1)
Y2 <- rbinom(n, 1, pY2)
Y3 <- rbinom(n, 1, pY3)
Y <- as.numeric(Y1*(Z==1) + Y2*(Z==2) + Y3*(Z==3)) # observed outcome

  
sim.data <- data.frame(Y=Y, Z=Z, X)
sim.data$V1 <- NULL

allcombn <- combn(ntrt,2) # all pairwise combinations

### Estimation of Propensity Score (2 Methods) ###
  
## Multinomial logistic regression
multinom.fit <- nnet::multinom(Z~X1+X2+X3, data=sim.data, trace=F)
gps.glm <- fitted(multinom.fit)
wts.glm <- 1/gps.glm

## Covariate balancing propensity score
# CBPS function from the package CBPS
CBPS.fit <- CBPS::CBPS(Z~X1+X2+X3, ATT=0, method="exact")
gps.CBPS <- CBPS.fit$fitted.values
wts.CBPS <- 1/gps.CBPS
  
### 1.Outcome Regression estimator (G-formula)
# Maybe we can put them into a function such that the input is Y, X, and Z (or a data set that contains there variables)
# and the output is the average treatment effect tau.
y1.fit <- glm(Y~X1+X2+X3, data=sim.data, subset=(Z==1), family="binomial")
y2.fit <- glm(Y~X1+X2+X3, data=sim.data, subset=(Z==2), family="binomial")
y3.fit <- glm(Y~X1+X2+X3, data=sim.data, subset=(Z==3), family="binomial")
  
h1 <- predict(y1.fit, newdata=sim.data, type="response") # estimated average response
h2 <- predict(y2.fit, newdata=sim.data, type="response")
h3 <- predict(y3.fit, newdata=sim.data, type="response")
h <- cbind(h1, h2, h3)

tau.OR <- c(mean(h2-h1), mean(h3-h1), mean(h3-h2))
  
### 2. Horvitz Thompson estimator
tau.HT.glm <- apply(allcombn, 2, function(x) mean(Y*(Z==x[2])*wts.glm[,x[2]]) - mean(Y*(Z==x[1])*wts.glm[,x[1]]))
tau.HT.CBPS <- apply(allcombn, 2, function(x) mean(Y*(Z==x[2])*wts.CBPS[,x[2]]) - mean(Y*(Z==x[1])*wts.CBPS[,x[1]]))

### 3. Inverse Probability Weighting estimator
tau.IPW.glm <- estTau.weighting(Y=Y, Z=Z, gps=gps.glm, hx=matrix(0, nrow=n, ncol=ntrt), target="all")
tau.IPW.CBPS <- estTau.weighting(Y=Y, Z=Z, gps=gps.CBPS, hx=matrix(0, nrow=n, ncol=ntrt), target="all")

### 4. Matching Weights estimator (Yoshida et al. 2017)
tau.MW.glm <- estTau.weighting(Y=Y, Z=Z, gps=gps.glm, hx=matrix(0, nrow=n, ncol=ntrt), target="matching")
tau.MW.CBPS <- estTau.weighting(Y=Y, Z=Z, gps=gps.CBPS, hx=matrix(0, nrow=n, ncol=ntrt), target="matching")

### 5. Overlap Weights estimator (Li and Li, 2019)
tau.OW.glm <- estTau.weighting(Y=Y, Z=Z, gps=gps.glm, hx=matrix(0, nrow=n, ncol=ntrt), target="overlap")
tau.OW.CBPS <- estTau.weighting(Y=Y, Z=Z, gps=gps.CBPS, hx=matrix(0, nrow=n, ncol=ntrt), target="overlap")
  
### 6. Augmented IPW (AIPW)
tau.DR.glm <- tau.DR.CBPS <- NULL
for(k in 1:ncol(allcombn)){
  t1 <- allcombn[1,k]; t2 <- allcombn[2,k]
  tau.DR.glm[k] <- tau.HT.glm[k] - mean(((Z==t2)-gps.glm[,t2])*h[,t2]/gps.glm[,t2]) + mean(((Z==t1)-gps.glm[,t1])*h[,t1]/gps.glm[,t1])
  tau.DR.CBPS[k] <- tau.HT.CBPS[k] - mean(((Z==t2)-gps.CBPS[,t2])*h[,t2]/gps.CBPS[,t2]) + mean(((Z==t1)-gps.CBPS[,t1])*h[,t1]/gps.CBPS[,t1])
}
 
### 7. Augmented MW 
tau.AMW.glm <- estTau.weighting(Y=Y, Z=Z, gps=gps.glm, hx=h, target="matching")
tau.AMW.CBPS <- estTau.weighting(Y=Y, Z=Z, gps=gps.CBPS, hx=h, target="matching")

### 8. Augmented OW
tau.AOW.glm <- estTau.weighting(Y=Y, Z=Z, gps=gps.glm, hx=h, target="overlap")
tau.AOW.CBPS <- estTau.weighting(Y=Y, Z=Z, gps=gps.CBPS, hx=h, target="overlap")

### 9. Matching on covariates
# Yang et al. 2016 implemented this method in their package multilevelMatching
res.MCOV <- suppressMessages(multilevelMatching::multiMatch(Y=Y, W=as.numeric(Z), X=as.matrix(cbind(X1,X2,X3)), match_on="covariates"))
tau.MCOV <- res.MCOV$results$Estimate

### 10. Matching on a scalar function of GPS
# From the package multilevelMatching
res.MGPSS.glm <- suppressMessages(multilevelMatching::multiMatch(Y=Y, W=as.numeric(Z), X=as.matrix(cbind(X1,X2,X3)), match_on="multinom", trimming=0))
tau.MGPSS.glm <- res.MGPSS.glm$results$Estimate

# Yang et al. did not implement the matching procedure for CBPS, so I modified their codes that can accomodate matching on CBPS.
# Note: the function GPSMatch can be applied to MGPSS.glm and MGPSV as well, but it only gives point estimates,
# since there is no variance calculation formula available for MGPSS.CBPS.

res.MGPSS.CBPS <- GPSMatch(Y=Y, W=as.numeric(Z), X=gps.CBPS.c, method="scalar")
tau.MGPSS.CBPS <- res.MGPSS.CBPS$tauestimate

### 11. Matching on the vector of GPS
res.MGPSV.glm <- suppressMessages(multilevelMatching::multiMatch(Y=Y, W=as.numeric(Z), X=gps.glm[,1:(ncol(gps.glm)-1)], match_on="covariates"))
tau.MGPSV.glm <- res.MGPSV.glm$results$Estimate

res.MGPSV.CBPS <- suppressMessages(multilevelMatching::multiMatch(Y=Y, W=as.numeric(Z), X=gps.CBPS.c[,1:(ncol(gps.CBPS)-1)], match_on="covariates"))
tau.MGPSV.CBPS <- res.MGPSV.CBPS$results$Estimate

### 12. PENCOMP (extension of Zhou et al. 2019 JASA)

pVarList <- c("X1","X2","X3") # variables that are associated with treatment
oVarList <- c("X1","X2","X3") # variables that are associated with outcome

nboot <- 100 # PENCOMP can be very slow due to the fitting of spline models.
tau.GAM.boot <- var.GAM.boot <- matrix(NA, nrow=nboot, ncol=npair)

for(k in 1:nboot){ # create nboot imputed data sets
  res.gamcomp <- gamcomp(dataOR=sim.data, propenVarList=pVarList, outcomeVarList=oVarList,
                         treat.varname="Z", outcome.varname="Y", original=0, numKnot=10)

  if(is.list(res.gamcomp)){ # the gamcomp function may return NA due to GAM fitting issues
    tau.GAM.boot[k,] <- res.gamcomp$ATE
    var.GAM.boot[k,] <- res.gamcomp$var
  }
}

# pool the imputed results 
# probably put them into a function for which the input is the bootstrap results tau.GAM.boot and var.GAM.boot
# and the output is the pooled ATE (tau.GAM.glm), pooled variance of ATE, and upper and lower bound of 95% CI.
tau.GAM.glm <- apply(tau.GAM.boot, 2, function(x) mean(x, na.rm=T)) # ATE

Wd <- apply(var.GAM.boot,2,function(x) mean(x, na.rm=T))  # within imputation variance
m <- apply(tau.GAM.boot,2,function(x) sum(complete.cases(x))) # number of complete cases for each pairwise comparison
Bd <- sapply(1:ncol(tau.GAM.boot), function(a) (1/(m[a]-1))*(sum((tau.GAM.boot[,a]-mean(tau.GAM.boot[,a], na.rm=T))^2, na.rm=T))) # between imputation variance
Td <- Wd + (1+1/m)*Bd # total variability associated with mean
df <- (m-1)*(1+(1/(1/m+1))*(Wd/Bd))^2 # degree of freedom

se.GAM.glm <- sqrt(Td) # standard error of ATE

tau.GAM.upper <- tau.GAM.glm + qt(0.975,df=df)*se.GAM.glm # upper bound of 95% CI
tau.GAM.lower <- tau.GAM.glm - qt(0.975,df=df)*se.GAM.glm # lower bound of 95% CI



