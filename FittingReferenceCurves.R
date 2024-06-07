########################################
## Fitting Reference Viability Curves ##
########################################

library(drc)
library(msm)


##Brain-Cousens function
my.bc <- function(x, b, c, d, e, f){
  summand <- exp(b*(log(x)-log(e)))
  return(c+((d-c+f*x)/(1+summand)))
}

##AIC for a constant function
aic.const.mod <- function(resp){
  nn <- length(resp) - sum(is.na(resp))
  mu <- mean(resp, na.rm=TRUE)
  sig2 <- (nn-1)/(nn)*var(resp, na.rm=TRUE)
  
  loglik <- -(nn/2)*log(2*pi*sig2)-(1/(2*sig2))*sum((resp-mu)^2, na.rm=TRUE)
  
  ## Analogously to the calculation of the AIC for the 4pLL and the BC model, in addition the the model parameter (in the constant
  ##  model given by the overall mean), the estimated variance is one parameter in the sense of the penalty term in the AIC
  aic <- -2*loglik+2*2
  return(aic)
}

# Curve fitting function with model selection step
curvefit.singlecomp <- function(daten, alpha, main.title, don.title, lambda=80){
  aic.4pll <- tryCatch({
    obj.4pll <- drm(daten[,"resp"] ~ daten[,"conc"], fct=LL2.4())
    mselect(obj.4pll)[2]
  },
  error=function(e){
    return(Inf)
  })
  
  aic.bc <- tryCatch({
    obj.bc <- drm(daten[,"resp"] ~ daten[,"conc"], fct=BC.5())
    aic <- mselect(obj.bc)[2]
    if(coef(obj.bc)[1] < 1 | coef(obj.bc)[5] < 0){
      aic <- Inf
    }
    aic
  },
  error=function(e){
    return(Inf)
  })
  
  aic.const <- aic.const.mod(daten[,"resp"])
  
  ##model selection step
  if(aic.const < min(aic.4pll, aic.bc)) opt.modell <- "const"
  if(aic.4pll < min(aic.const, aic.bc)) opt.modell <- "LL2.4"
  if(aic.bc < min(aic.const, aic.4pll)) opt.modell <- "BC.5"
  
  ##Constant model
  if(opt.modell == "const"){
    fit <- mean(daten[,"resp"])
    resp.norm <- (daten[,"resp"]/fit)*100
    daten <- cbind(daten, resp.norm)
    
    gof <- NA
    est <- paste0(">", max(daten[,"conc"]))
    
    est.return <- c("Const", est, NA, NA, NA, NA)
    names(est.return) <- c("Model", "Estimate", "StandardError", "CI.Lower", "CI.Upper", "GoF")
    return(list(Est = est.return,
                Model = list(Type = "Const",
                             GoF = NA,
                             est = NA,
                             Pars = mean(resp.norm),
                             data = daten)))
  }


  ##4pLL model
  if(opt.modell == "LL2.4"){
    # first fit
    obj <- drm(daten[,"resp"] ~ daten[,"conc"], fct=LL2.4())
    bb <- obj$coefficients["b:(Intercept)"]
    cc <- obj$coefficients["c:(Intercept)"]
    dd <- obj$coefficients["d:(Intercept)"]
    ee <- obj$coefficients["e:(Intercept)"]
    
    # normalized fit
    norm <- ifelse(bb > 0, dd, cc)
    resp.norm <- (daten[,"resp"]/norm)*100
    daten <- cbind(daten, resp.norm)
    obj.norm <- drm(resp.norm ~ daten[,"conc"], fct=LL2.4())
    
    gof <- 1-obj.norm$fit$ovalue/((length(resp.norm)-1)*var(resp.norm, na.rm=TRUE))
    
    ## The EC value is not calculated if the mean of the highest concentration is larger than 90%
    if(mean(resp.norm[which(daten[,"conc"] == max(daten[,"conc"]))]) > 90){
      est <- paste0(">", max(daten[,"conc"]))
      se <- NA
      est.lo <- NA
      est.up <- NA
    }
    if(mean(resp.norm[which(daten[,"conc"] == max(daten[,"conc"]))]) <= 90){
      ## inverse
      est <- exp(coef(obj.norm)["e:(Intercept)"])*((coef(obj.norm)["d:(Intercept)"]-lambda)/(lambda-coef(obj.norm)["c:(Intercept)"]))^(1/coef(obj.norm)["b:(Intercept)"])
      
      nu <- df.residual(obj.norm)
      seformula <- sprintf("~ log(exp(x4)*((x3-%f)/(%f-x2))^(1/x1))",lambda,lambda)
      se <- deltamethod(as.formula(seformula),coef(obj.norm),vcov(obj.norm))
      est.lo <- exp(log(est)-qt(1-alpha/2,df.residual(obj.norm))*se)
      est.up <- exp(log(est)+qt(1-alpha/2,df.residual(obj.norm))*se)
      
      if(est > max(daten[,"conc"]) | is.na(est)){
        est <- paste0(">", max(daten[,"conc"]))
        se <- NA
        est.lo <- NA
        est.up <- NA
      }
    }
    
    est.return <- c("4pLL", ifelse(is.numeric(est), round(est, 6), est), round(se, 6), round(est.lo, 6), round(est.up, 6), round(gof, 3))
    names(est.return) <- c("Model", "Estimate", "StandardError", "CI.Lower", "CI.Upper", "GoF")
    return(list(Est = est.return,
                Model = list(Type = "4pLL",
                             GoF = gof,
                             est = est,
                             Pars = coef(obj.norm),
                             data = daten)))
  }
  
  ##BC model
  if(opt.modell == "BC.5"){
    
    # first fit
    obj <- drm(daten[,"resp"] ~ daten[,"conc"], fct=BC.5())
    
    norm <- coef(obj)[3]
    resp.norm <- (daten[,"resp"]/norm)*100
    daten <- cbind(daten, resp.norm)
    
    # second/normalized fit numerically not possible, so parameters are hard-coded
    pars <- coef(obj)
    pars[c(2,3,5)] <- (pars[c(2,3,5)]/norm)*100
    
    ovalue <- sum((my.bc(daten[,"conc"], pars[1], pars[2], pars[3], pars[4], pars[5])-resp.norm)^2, na.rm=TRUE)
    gof <- 1-ovalue/((length(resp.norm)-1)*var(resp.norm, na.rm=TRUE))
    
    if(mean(resp.norm[which(daten[,"conc"] == max(daten[,"conc"]))]) > 90 | is.na(mean(resp.norm[which(daten[,"conc"] == max(daten[,"conc"]))]))){
      est <- paste0(">", max(daten[,"conc"]))
    }
    if(mean(resp.norm[which(daten[,"conc"] == max(daten[,"conc"]))], na.rm=TRUE) <= 90 & !is.na(mean(resp.norm[which(daten[,"conc"] == max(daten[,"conc"]))]))){
      
      # calculate EC value via grid search
      grid <- exp(seq(from=log(sort(unique(daten[,"conc"]))[2]), to= log(max(daten[,"conc"])), length.out=100000))
      # take care of wrong results due to the non-monotonous nature of the curve - if several intersections are found, use the last
      ind.iq <- which(abs(my.bc(grid, pars[1], pars[2], pars[3], pars[4], pars[5]) - lambda) < 0.05)
      if(length(ind.iq) == 0) est <- NA
      if(length(ind.iq) > 0) est <- grid[ind.iq[length(ind.iq)]]

      est.lo <- NA
      est.up <- NA
      se <- NA
      
      if(est > max(daten[,"conc"]) | is.na(est)){
        est <- paste0(">", max(daten[,"conc"]))
      }
    }
    
    est.return <- c("BC", ifelse(is.numeric(est), round(est, 6), est), NA, NA, NA, round(gof, 3))
    names(est.return) <- c("Model", "Estimate", "StandardError", "CI.Lower", "CI.Upper", "GoF")
    return(list(Est = est.return,
                Model = list(Type ="BC",
                             GoF = gof,
                             est = est,
                             Pars = pars,
                             data = daten)))
  }
}



### Fitting the curves
load("Cytotox-HepG2.RData")


## Extract the EC values only
lambda <- 80
ECvalues <- lapply(names(list.concresponse), function(comp){
  print(comp)
  data <- list.concresponse[[comp]]
  zwerg <- lapply(sort(names(data)), function(don){
    daten <- data[[don]]
    curvefit.singlecomp(daten, alpha=0.05, lambda=lambda)
  })
  out.mat <- t(sapply(zwerg, function(set){
    set$Est
  }))

  colnames(out.mat) <- c("Model", "Estimate", "StandardError", "CI.Lower", "CI.Upper", "GoF")
  rownames(out.mat) <- sort(names(data))
  out.mat
})
names(ECvalues) <- names(list.concresponse)

# print the ECvalues
ECvalues


## Save the entire fitted models
lambda <- 80
models <- lapply(names(list.concresponse), function(comp){
  print(comp)
  data <- list.concresponse[[comp]]
  zwerg <- lapply(sort(names(data)), function(don){
    daten <- data[[don]]
    curvefit.singlecomp(daten, alpha=0.05, lambda=lambda)
  })
  out <- lapply(zwerg, function(set){
    set$Model
  })
  out
})
names(models) <- names(list.concresponse)

# save the models
save(models, file = "Models-CTB-HepG2.RData")





