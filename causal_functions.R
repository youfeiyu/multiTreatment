

estTau.weighting <- function(Y, Z, gps, hx, target){
  # estimate pairwise average treatment effects for all groups
  # Y: outcome variable
  # Z: treatment variable
  # gps: matrix of generalized propensity scores
  # hx: matrix of function of X in the augmented term
  # target: target population, can be 'all', 'matching', 'overlap'
  
  if(target=="all"){
    fw <- rep(1, length(Z))
  } else if(target=="matching"){
    fw <- apply(gps, 1, min)
  } else if(target=="overlap"){
    fw <- 1/apply(1/gps, 1, sum)
  }
  
  ntrt <- length(table(Z))
  allcombn <- combn(ntrt, 2)
  
  Z.mat <- matrix(rep(Z,ntrt), ncol=ntrt)
  Z.dummy <- t(apply(Z.mat, 1, function(x) x==1:ntrt))
  
  bw <- fw/rowSums(Z.dummy*gps) # balancing weights
  
  mu.j <- NULL
  for(j in 1:length(table(Z))){
    #mu.j[j] <- mean(bw*(Z==j)*(Y-hx[,j]))/sum(fw) + weighted.mean(hx[,j], fw)
    mu.j[j] <- weighted.mean(Y-hx[,j], bw*(Z==j)) + weighted.mean(hx[,j], fw)
  }
  
  apply(allcombn, 2, function(a) mu.j[a[2]] - mu.j[a[1]])
  
}

GPSMatch <- function(Y, W, X, method){
  
  # Y: vector of outcome
  # W: vector of treatment
  # X: propensity scores
  
  ## order the treatment increasingly
  if(1-is.unsorted(W)){
    sorted <- sort(W,index.return=TRUE)
    sorted <- list(x=sorted)
    sorted$ix <- 1:length(W)
  }
  
  if(is.unsorted(W)){
    sorted <- sort(W, index.return=TRUE)
  }
  
  W <- W[sorted$ix]
  N <- length(Y) # number of observations
  X <- as.matrix(X)
  X <- X[sorted$ix,]
  Y <- Y[sorted$ix]
  
  trtnumber <- length(unique(W)) # number of treatment levels
  trtlevels <- unique(W) # all treatment levels
  pertrtlevelnumber <- table(W) # number of observations by treatment level
  taunumber <- choose(trtnumber,2)  # number of pairwise treatment effects
  
  tauestimate <- rep(NA,taunumber)
  
  Yiw <- matrix(NA,N,trtnumber) # Yiw is the full imputed data set
  
  for(kk in 1:trtnumber){
    thistrt <- trtlevels[kk]
    if(kk==1){fromto <- 1:pertrtlevelnumber[1]}
    if(kk>1){fromto <- (1:pertrtlevelnumber[kk])+sum(pertrtlevelnumber[1:(kk-1)])}
    W1 <- (W!=thistrt)
    if(method=="scalar"){
      out <- Matching::Match(Y=Y, Tr=W1, X=X[,kk], distance.tolerance=0, replace=T, ties=F, Weight=2)
    } else if(method=="vector"){
      Weight.matrix <- diag(ncol(X)-1)
      out <- Matching::Match(Y=Y, Tr=W1, X=X[,1:(ncol(X)-1)], distance.tolerance=0, replace=T, ties=F, Weight=3, Weight.matrix=Weight.matrix)
    } 
    
    mdata = out$mdata
    
    Yiw[which(W==thistrt), kk] <- Y[which(W==thistrt)] # observed potential outcome Y(thistrt)
    Yiw[out$index.treated, kk] <- mdata$Y[which(mdata$Tr==0)] # imputed potential outcome Y(thistrt)
  }
  
  allcombn = combn(1:trtnumber,2)
  
  cname1<-NULL
  for(i in 1:ncol(allcombn)){
    thistrt = trtlevels[allcombn[1,i]]
    thattrt = trtlevels[allcombn[2,i]]
    cname1[i] = paste("EY(", thattrt, ")-EY(", thistrt, ")", sep="")
    tauestimate[i] = mean(Yiw[,thattrt]- Yiw[,thistrt], na.rm=T)
  }
  
  names(tauestimate) <- cname1
  
  # Sort the indices back to the original order
  imputed <- Yiw
  imputed[sorted$ix,] <- Yiw
  
  return(list(tauestimate = tauestimate,
              impute_mat = imputed))
}

formulaF <- function(varList, y.name){
  return ( as.formula(paste(y.name, "~ ", paste(c(varList), collapse = "+"))) )
}


gamcomp <- function(dataOR, propenVarList, outcomeVarList, treat.varname, outcome.varname, original=0, numKnot) {
  
  tryCatch ( 
    {
      ntrt <- length(table(dataOR[, treat.varname])) # number of treatment groups
      npair <- choose(ntrt, 2)
      data <- NULL
      if (original==1){
        data <- dataOR
      } else {
        
        # Create bootstrap samples stratified on treatment
        sampledID <- NULL
        for(z in 1:ntrt){
          trtID <- which(dataOR[, treat.varname]==z) # subjects in treatment group z
          sampledID <- c(sampledID, sample(trtID, replace=T))
        }
        data <- dataOR[sampledID, ]  ##in the simulation studies, we did stratified bootstraps, results were similar with random bootstraps 
      }
      
      treatInd <- dataOR[, treat.varname]  ###indices of treatment in the original dataset
      treatBoot <- data[, treat.varname]  ###indices of treatment in the bootstrap sample
      Yobs <- dataOR[, outcome.varname]
      
      ######## propensity score model ########################################################### 
      propen.formula <- formulaF(varList=propenVarList, y.name=treat.varname)
      
      model2a <- nnet::multinom(propen.formula, data=data, trace=F)
      prob.boot <- model2a$fitted.values
      
      data.pslogit <- log(prob.boot/(1-prob.boot))   ####spline on the logit of propensity
      data.pslogit <- data.pslogit[,1:(ncol(prob.boot)-1)]
      
      
      ####fit the prediction model for all the models
      ###for the outcome model
      ##############################################################
      ##############################################################
      ###use equally spaced fixed knots assuming K knots
      
      imputedY <- matrix(NA, nrow(data), ntrt)
      
      for(j in 1:ntrt){
        # Fit prediction model for Y
        data.j <- data[data$Z==j,]
        propen.score <- data.pslogit[data$Z==j,]
        space <- (apply(propen.score,2,max)-apply(propen.score,2,min))/(numKnot+1)
        knots <- sweep(sapply(space, function(x) x*(1:numKnot)), 2, apply(propen.score,2,min), "+")
        
        linearB <- NULL
        for(l in 1:(ntrt-1)){
          linearB.l <- outer(propen.score[,l], knots[,l], "-")
          linearB.l <- linearB.l * (linearB.l > 0)
          linearB <- cbind(linearB, linearB.l)
        }
        linearB <- as.matrix(linearB)
        
        response <- data.j[, outcome.varname]
        covariateX <- NULL
        for(kk in 1:length(outcomeVarList)){
          covariateX <- cbind(covariateX, data.j[, outcomeVarList[kk]])
        }
        covariateX <- cbind(rep(1,nrow(data.j)), covariateX, propen.score)
        
        for(jj in 1:ncol(linearB)){
          eval(parse(text=(paste0("linearB",jj,"=linearB[,",jj,"]"))))
        }
        
        pspp.j <- mgcv::gam(as.formula(paste("response~covariateX-1+", paste("s(linearB",1:ncol(linearB),",bs='re')", collapse="+", sep=""))), family=binomial)
        
        # Impute the missing potential outcomes
        newData <- dataOR
        newData[, treat.varname] <- j
        pred2a <- predict(model2a, newdata=newData, type="probs")
        
        newData.pslogit <- log(pred2a/(1-pred2a))
        newData.pslogit <- newData.pslogit[,1:(ncol(pred2a)-1)]
        propen.score.new <- newData.pslogit
        linearB <- NULL
        for(l in 1:(ntrt-1)){
          linearB.l <- outer(propen.score.new[,l], knots[,l], "-")
          linearB.l <- linearB.l * (linearB.l > 0)
          linearB <- cbind(linearB, linearB.l)
        }
        linearB <- as.matrix(linearB)
        
        covariateX <- NULL
        for(kk in 1:length(outcomeVarList)){
          covariateX <- cbind(covariateX, newData[, outcomeVarList[kk]])
        }
        covariateX <- cbind(rep(1,nrow(newData)), covariateX, propen.score.new)
        
        pYlogit <- as.numeric(cbind(covariateX, linearB) %*% as.matrix(coef(pspp.j)))
        imputedY[,j] <- rbinom(nrow(dataOR), 1, exp(pYlogit)/(1+exp(pYlogit)))
        # keep the observed outcome the same
        imputedY[treatInd==j,j] = Yobs[treatInd==j]
      }
      
      #######include everyone
      allcombn <- combn(ntrt, 2)
      yhat <- apply(imputedY,2,mean)
      ATE <- NULL; var.ATE <- NULL
      for(ii in 1:npair){
        p1 <- yhat[allcombn[1,ii]]; p2 <- yhat[allcombn[2,ii]]
        
        ATE[ii] <- p2 - p1
        var.ATE[ii] <- p2*(1-p2)/nrow(dataOR) + p1*(1-p1)/nrow(dataOR)
        
        #var.ATE[k] <- var(imputedY[,allcombn[2,k]]-imputedY[,allcombn[1,k]])/nrow(dataOR)
      }
      
      return(list(ATE=ATE, var=var.ATE)) ###estimate and variance of ATE
    }, error=function(e) return(NA) )
}



