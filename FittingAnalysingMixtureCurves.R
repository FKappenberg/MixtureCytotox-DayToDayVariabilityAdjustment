####################################
## Fitting the mixture curves     ##
##  and adjusting for day-to-day  ##
##  variability                   ##
####################################

library(drc)

##Brain-Cousens function
my.bc <- function(x, b, c, d, e, f){
  summand <- exp(b*(log(x)-log(e)))
  return(c+((d-c+f*x)/(1+summand)))
}

##4pLL function 
my.4pll <- function(x, b, c, d, etilde){
  summand <- exp(b*(log(x)-etilde))
  return(c+(d-c)/(1+summand))
}



# modeldata: list of appropriate numbers of models (depending on approach / in the responsible approach,
# on how many models the median is based)
# rtilde: viability in the current measurement
# approach: either "responsible" (only the responsible curve(s) is/are considered) or "median"
ctilde_function <- function(modeldata, rtilde, approach){

  if(approach == "responsible"){
    # differentiate between only one responsible or two responsible curves
    if(length(modeldata) == 1){
      # BC Curve
      if(modeldata[[1]]$Type == "BC"){
        # gridsearch to determine the concentration where the curve attains the value rtilde
        grid <- exp(seq(from=log(sort(unique(modeldata[[1]]$data[,2]))[2]/10),
                        to = log(max(modeldata[[1]]$data[,2])),
                        length.out=100000))
        ctilde <- grid[min(which(my.bc(grid, 
                                       modeldata[[1]]$Pars[1],
                                       modeldata[[1]]$Pars[2],
                                       modeldata[[1]]$Pars[3],
                                       modeldata[[1]]$Pars[4],
                                       modeldata[[1]]$Pars[5])<=rtilde))]
        if(is.na(ctilde)) ctilde <- max(modeldata[[1]]$data[,2])
        return(ctilde)    
      }  
      
      # 4pLL Curve
      if(modeldata[[1]]$Type == "4pLL"){
        # gridsearch to determine the concentration where the curve attains the value rtilde
        grid <- exp(seq(from=log(sort(unique(modeldata[[1]]$data[,2]))[2]/10),
                        to = log(max(modeldata[[1]]$data[,2])),
                        length.out=100000))
        ctilde <- grid[min(which(my.4pll(grid, 
                                         modeldata[[1]]$Pars[1],
                                         modeldata[[1]]$Pars[2],
                                         modeldata[[1]]$Pars[3],
                                         modeldata[[1]]$Pars[4])<=rtilde))]
        if(is.na(ctilde)) ctilde <- max(modeldata[[1]]$data[,2])
        return(ctilde)
      }  
    }
    if(length(modeldata) == 2){
      # deal with the first curve
      if(modeldata[[1]]$Type == "BC"){
        # gridsearch to determine the concentration where the curve attains the value rtilde
        grid <- exp(seq(from=log(sort(unique(modeldata[[1]]$data[,2]))[2]/10),
                        to = log(max(modeldata[[1]]$data[,2])),
                        length.out=100000))
        ctilde_1 <- grid[min(which(my.bc(grid, 
                                         modeldata[[1]]$Pars[1],
                                         modeldata[[1]]$Pars[2],
                                         modeldata[[1]]$Pars[3],
                                         modeldata[[1]]$Pars[4],
                                         modeldata[[1]]$Pars[5])<=rtilde))]
        if(is.na(ctilde_1)) ctilde_1 <- max(modeldata[[1]]$data[,2])
      }  
      
      # 4pLL Curve
      if(modeldata[[1]]$Type == "4pLL"){
        # gridsearch to determine the concentration where the curve attains the value rtilde
        grid <- exp(seq(from=log(sort(unique(modeldata[[1]]$data[,2]))[2]/10),
                        to = log(max(modeldata[[1]]$data[,2])),
                        length.out=100000))
        ctilde_1 <- grid[min(which(my.4pll(grid, 
                                           modeldata[[1]]$Pars[1],
                                           modeldata[[1]]$Pars[2],
                                           modeldata[[1]]$Pars[3],
                                           modeldata[[1]]$Pars[4])<=rtilde))]
        if(is.na(ctilde_1)) ctilde_1 <- max(modeldata[[1]]$data[,2])
      }  
      # deal with the second curve
      if(modeldata[[2]]$Type == "BC"){
        # gridsearch to determine the concentration where the curve attains the value rtilde
        grid <- exp(seq(from=log(sort(unique(modeldata[[2]]$data[,2]))[2]/10),
                        to = log(max(modeldata[[2]]$data[,2])),
                        length.out=100000))
        ctilde_2 <- grid[min(which(my.bc(grid, 
                                         modeldata[[2]]$Pars[1],
                                         modeldata[[2]]$Pars[2],
                                         modeldata[[2]]$Pars[3],
                                         modeldata[[2]]$Pars[4],
                                         modeldata[[2]]$Pars[5])<=rtilde))]
        if(is.na(ctilde_2)) ctilde_2 <- max(modeldata[[2]]$data[,2])
      }  
      
      # 4pLL Curve
      if(modeldata[[2]]$Type == "4pLL"){
        # gridsearch to determine the concentration where the curve attains the value rtilde
        grid <- exp(seq(from=log(sort(unique(modeldata[[1]]$data[,2]))[2]),
                        to = log(max(modeldata[[1]]$data[,2])),
                        length.out=100000))
        ctilde_2 <- grid[min(which(my.4pll(grid, 
                                           modeldata[[2]]$Pars[1],
                                           modeldata[[2]]$Pars[2],
                                           modeldata[[2]]$Pars[3],
                                           modeldata[[2]]$Pars[4])<=rtilde))]
        if(is.na(ctilde_2)) ctilde_2 <- max(modeldata[[2]]$data[,2])
      }  
      return(mean(c(ctilde_1, ctilde_2)))
    }
  }

  if(approach == "median"){
    ctildes <- sapply(1:length(modeldata), function(i){
      if(modeldata[[i]]$Type == "BC"){
        # gridsearch to determine the concentration where the curve attains the value rtilde
        grid <- exp(seq(from=log(sort(unique(modeldata[[i]]$data[,2]))[2]/10),
                        to = log(max(modeldata[[i]]$data[,2])),
                        length.out=100000))
        ctilde <- grid[min(which(my.bc(grid, 
                                       modeldata[[i]]$Pars[1],
                                       modeldata[[i]]$Pars[2],
                                       modeldata[[i]]$Pars[3],
                                       modeldata[[i]]$Pars[4],
                                       modeldata[[i]]$Pars[5])<=rtilde))]
        return(ctilde)    
      }  
      
      # 4pLL Curve
      if(modeldata[[i]]$Type == "4pLL"){
        # gridsearch to determine the concentration where the curve attains the value rtilde
        grid <- exp(seq(from=log(sort(unique(modeldata[[i]]$data[,2]))[2]/10),
                        to = log(max(modeldata[[i]]$data[,2])),
                        length.out=100000))
        ctilde <- grid[min(which(my.4pll(grid, 
                                         modeldata[[i]]$Pars[1],
                                         modeldata[[i]]$Pars[2],
                                         modeldata[[i]]$Pars[3],
                                         modeldata[[i]]$Pars[4])<=rtilde))]
        return(ctilde)  
      }
      # Flat curve
      if(modeldata[[i]]$Type == "Const"){
        return(NA)
      }
    })
    if(any(is.na(ctildes))) ctildes[which(is.na(ctildes))] <- Inf
    out <- median(ctildes)
    if(out == Inf) out <- max(sapply(1:length(modeldata), function(i) max(modeldata[[i]]$data[,2])))
    return(out)
  }
}





load("MixtureData.RData")
load("Models-CTB-HepG2.RData")

set3_ctrlmean <- tapply(set3_drcdata$Ctrl, set3_drcdata$biological.replicate, mean)

set3_drcdata_forfitting <- data.frame(
  Prop = rep(c(0, 1/128, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1), 4),
  Resp1 = (unlist(c(set3_drcdata[1, 2:10], set3_drcdata[2, 2:10], set3_drcdata[3, 2:10], set3_drcdata[4, 2:10]))/set3_ctrlmean[1])*100,
  Resp2 = (unlist(c(set3_drcdata[5, 2:10], set3_drcdata[6, 2:10], set3_drcdata[7, 2:10], set3_drcdata[8, 2:10]))/set3_ctrlmean[2])*100,
  Resp3 = (unlist(c(set3_drcdata[9, 2:10], set3_drcdata[10, 2:10], set3_drcdata[11, 2:10], set3_drcdata[12, 2:10]))/set3_ctrlmean[3])*100 
)


## calculate the replicate-wise means of the negative controls
ctrl_set3_1 <- mean(c(set3_einzel$technical.replicate.1[1],
                      set3_einzel$technical.replicate.2[1],
                      set3_einzel$technical.replicate.3[1],
                      set3_einzel$technical.replicate.4[1]))

ctrl_set3_2 <- mean(c(set3_einzel$technical.replicate.1[12],
                      set3_einzel$technical.replicate.2[12],
                      set3_einzel$technical.replicate.3[12],
                      set3_einzel$technical.replicate.4[12]))

ctrl_set3_3 <- mean(c(set3_einzel$technical.replicate.1[23],
                      set3_einzel$technical.replicate.2[23],
                      set3_einzel$technical.replicate.3[23],
                      set3_einzel$technical.replicate.4[23]))


## Biol_Rep 1
apap_set3_1 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "APAP" & 
                                                    set3_einzel$biological.replicate == 1),2:3]))/ctrl_set3_1 *100
clon_set3_1 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "CLON" & 
                                                    set3_einzel$biological.replicate == 1),2:3]))/ctrl_set3_1 *100
cyp_set3_1 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "CYP" & 
                                                   set3_einzel$biological.replicate == 1),2:3]))/ctrl_set3_1 *100
inah_set3_1 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "INAH" & 
                                                    set3_einzel$biological.replicate == 1),2:3]))/ctrl_set3_1 *100
lab_set3_1 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "LAB" & 
                                                   set3_einzel$biological.replicate == 1),2:3]))/ctrl_set3_1 *100
lev_set3_1 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "LEV" & 
                                                   set3_einzel$biological.replicate == 1),2:3]))/ctrl_set3_1 *100
oxc_set3_1 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "OXC" & 
                                                   set3_einzel$biological.replicate == 1),2:3]))/ctrl_set3_1 *100
spb_set3_1 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "SPB" & 
                                                   set3_einzel$biological.replicate == 1),2:3]))/ctrl_set3_1 *100
vitc_set3_1 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "VITC" & 
                                                    set3_einzel$biological.replicate == 1),2:3]))/ctrl_set3_1 *100
vpa_set3_1 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "VPA" & 
                                                   set3_einzel$biological.replicate == 1),2:3]))/ctrl_set3_1 *100


## Biol_Rep 2
apap_set3_2 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "APAP" & 
                                                    set3_einzel$biological.replicate == 2),2:3]))/ctrl_set3_2 *100
clon_set3_2 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "CLON" & 
                                                    set3_einzel$biological.replicate == 2),2:3]))/ctrl_set3_2 *100
cyp_set3_2 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "CYP" & 
                                                   set3_einzel$biological.replicate == 2),2:3]))/ctrl_set3_2 *100
inah_set3_2 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "INAH" & 
                                                    set3_einzel$biological.replicate == 2),2:3]))/ctrl_set3_2 *100
lab_set3_2 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "LAB" & 
                                                   set3_einzel$biological.replicate == 2),2:3]))/ctrl_set3_2 *100
lev_set3_2 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "LEV" & 
                                                   set3_einzel$biological.replicate == 2),2:3]))/ctrl_set3_2 *100
oxc_set3_2 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "OXC" & 
                                                   set3_einzel$biological.replicate == 2),2:3]))/ctrl_set3_2 *100
spb_set3_2 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "SPB" & 
                                                   set3_einzel$biological.replicate == 2),2:3]))/ctrl_set3_2 *100
vitc_set3_2 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "VITC" & 
                                                    set3_einzel$biological.replicate == 2),2:3]))/ctrl_set3_2 *100
vpa_set3_2 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "VPA" & 
                                                   set3_einzel$biological.replicate == 2),2:3]))/ctrl_set3_2 *100



## Biol_Rep 2
apap_set3_3 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "APAP" & 
                                                    set3_einzel$biological.replicate == 3),2:3]))/ctrl_set3_3 *100
clon_set3_3 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "CLON" & 
                                                    set3_einzel$biological.replicate == 3),2:3]))/ctrl_set3_3 *100
cyp_set3_3 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "CYP" & 
                                                   set3_einzel$biological.replicate == 3),2:3]))/ctrl_set3_3 *100
inah_set3_3 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "INAH" & 
                                                    set3_einzel$biological.replicate == 3),2:3]))/ctrl_set3_3 *100
lab_set3_3 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "LAB" & 
                                                   set3_einzel$biological.replicate == 3),2:3]))/ctrl_set3_3 *100
lev_set3_3 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "LEV" & 
                                                   set3_einzel$biological.replicate == 3),2:3]))/ctrl_set3_3 *100
oxc_set3_3 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "OXC" & 
                                                   set3_einzel$biological.replicate == 3),2:3]))/ctrl_set3_3 *100
spb_set3_3 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "SPB" & 
                                                   set3_einzel$biological.replicate == 3),2:3]))/ctrl_set3_3 *100
vitc_set3_3 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "VITC" & 
                                                    set3_einzel$biological.replicate == 3),2:3]))/ctrl_set3_3 *100
vpa_set3_3 <- unlist(as.vector(set3_einzel[which(set3_einzel$Abbreviation == "VPA" & 
                                                   set3_einzel$biological.replicate == 3),2:3]))/ctrl_set3_3 *100



# this contains the "rtilde" values, i.e. the viabilities at the reference concentration c_i
set3_tagesaktuell_df <- data.frame(
  Compound = factor(rep(rep(c("APAP", "CLON", "CYP", "INAH", "LAB", "LEV", "OXC", "SPB", "VITC", "VPA"), each=2), 3)),
  TechRep = rep(c(1,2), 30),
  BiolRep = factor(rep(c(1, 2, 3), each=20)),
  Value = c(apap_set3_1, clon_set3_1, cyp_set3_1, inah_set3_1, lab_set3_1, lev_set3_1, oxc_set3_1, spb_set3_1, vitc_set3_1, vpa_set3_1,
            apap_set3_2, clon_set3_2, cyp_set3_2, inah_set3_2, lab_set3_2, lev_set3_2, oxc_set3_2, spb_set3_2, vitc_set3_2, vpa_set3_2,
            apap_set3_3, clon_set3_3, cyp_set3_3, inah_set3_3, lab_set3_3, lev_set3_3, oxc_set3_3, spb_set3_3, vitc_set3_3, vpa_set3_3)
)



comp_set3 <- c("APAP", "CLON", "CYP", "INAH", "LAB", "LEV", "OXC", "SPB", "VITC", "VPA")

set3_ctilde <- lapply(comp_set3, function(comp){
  print(comp)
  
  ## account for the situation that some of the compounds have several corresponding measurements 
  # (e.g. LAB was measured as LAB, LAB2 in the original dataset)
  compindata <- grep(comp, names(models))
  
  if(length(compindata) == 1){
    modelsdata <- models[[compindata]]
  }
  if(length(compindata) == 2){
    modelsdata <- c(models[[compindata[1]]], models[[compindata[2]]])
  }
  if(length(compindata) == 3){
    modelsdata <- c(models[[compindata[1]]], models[[compindata[2]]], models[[compindata[3]]])
  }
  
  orec20 <- originalec20_set3[comp]
  
  #find out which of the replicates for the compounds in the original data
  # was responsible for the EC20, median
  ec20allmodels <- sapply(1:length(modelsdata), function(i) modelsdata[[i]]$est)
  orec20exact <- median(ec20allmodels)
  print(round(orec20exact, 2))
  if(length(ec20allmodels) %% 2 == 1){
    print(any(ec20allmodels == orec20exact))
    whichrep <- which(ec20allmodels == orec20exact)
  }
  if(length(ec20allmodels) %% 2 == 0){
    print(mean(sort(ec20allmodels)[c(length(ec20allmodels)/2, length(ec20allmodels)/2+1)]) == orec20exact)
    
    whichrep <- which(ec20allmodels %in% sort(ec20allmodels)[c(length(ec20allmodels)/2, length(ec20allmodels)/2+1)])
  }
  
  
  rtilde_rep1 <- tapply(set3_tagesaktuell_df$Value[which(set3_tagesaktuell_df$BiolRep == 1)],
                        set3_tagesaktuell_df$Compound[which(set3_tagesaktuell_df$BiolRep == 1)],
                        mean)[comp]
  rtilde_rep2 <- tapply(set3_tagesaktuell_df$Value[which(set3_tagesaktuell_df$BiolRep == 2)],
                        set3_tagesaktuell_df$Compound[which(set3_tagesaktuell_df$BiolRep == 2)],
                        mean)[comp]
  rtilde_rep3 <- tapply(set3_tagesaktuell_df$Value[which(set3_tagesaktuell_df$BiolRep == 3)],
                        set3_tagesaktuell_df$Compound[which(set3_tagesaktuell_df$BiolRep == 3)],
                        mean)[comp]
  
  ## approach a) determine the tilde(c_i) based on whatever curve was responsible for c as the median
  if(length(whichrep)==1){
    respcurve <- list(modelsdata[[whichrep]])
  }
  if(length(whichrep)==2){
    respcurve <- list(modelsdata[[whichrep[1]]], modelsdata[[whichrep[2]]])
  }
  
  # determine tilde(c) as the concentration where the curve attains the value tilde(r)
  ctilde_resp_rep1 <- ctilde_function(respcurve, rtilde_rep1, approach = "responsible")
  ctilde_resp_rep2 <- ctilde_function(respcurve, rtilde_rep2, approach = "responsible")
  ctilde_resp_rep3 <- ctilde_function(respcurve, rtilde_rep3, approach = "responsible")
  
  ## approach b) determine the tilde(c_i) based on all curves and take the median afterwards
  ctilde_median_rep1 <- ctilde_function(modelsdata, rtilde_rep1, approach = "median")
  ctilde_median_rep2 <- ctilde_function(modelsdata, rtilde_rep2, approach = "median")
  ctilde_median_rep3 <- ctilde_function(modelsdata, rtilde_rep3, approach = "median")
  
  
  out <- c(ctilde_resp_rep1, ctilde_resp_rep2, ctilde_resp_rep3,
           ctilde_median_rep1, ctilde_median_rep2, ctilde_median_rep3)
  names(out) <- c("Resp-Rep1", "Resp-Rep2", "Resp-Rep3",
                  "Med-Rep1", "Med-Rep2", "Med-Rep3")
  return(out)
})
names(set3_ctilde) <- comp_set3



##Adjusted Mean Budget

# Consider responsible curve only:
adjusted_budget_set3_resp_rep1 <- mean(sapply(names(set3_ctilde), function(comp){
  set3_ctilde[[comp]]["Resp-Rep1"]/originalec20_set3[comp]*100
}))
print(adjusted_budget_set3_resp_rep1)
adjusted_budget_set3_resp_rep2 <- mean(sapply(names(set3_ctilde), function(comp){
  set3_ctilde[[comp]]["Resp-Rep2"]/originalec20_set3[comp]*100
}))
print(adjusted_budget_set3_resp_rep2)
adjusted_budget_set3_resp_rep3 <- mean(sapply(names(set3_ctilde), function(comp){
  set3_ctilde[[comp]]["Resp-Rep3"]/originalec20_set3[comp]*100
}))
print(adjusted_budget_set3_resp_rep3)

# Consider the median:
adjusted_budget_set3_med_rep1 <- mean(sapply(names(set3_ctilde), function(comp){
  set3_ctilde[[comp]]["Med-Rep1"]/originalec20_set3[comp]*100
}))
print(adjusted_budget_set3_med_rep1)
adjusted_budget_set3_med_rep2 <- mean(sapply(names(set3_ctilde), function(comp){
  set3_ctilde[[comp]]["Med-Rep2"]/originalec20_set3[comp]*100
}))
print(adjusted_budget_set3_med_rep2)
adjusted_budget_set3_med_rep3 <- mean(sapply(names(set3_ctilde), function(comp){
  set3_ctilde[[comp]]["Med-Rep3"]/originalec20_set3[comp]*100
}))
print(adjusted_budget_set3_med_rep3)



## Fit the mixture curves and normalize them
fit_set3_rep1 <- drm(Resp1~Prop, data=set3_drcdata_forfitting, fct=LL2.4())
fit_set3_rep2 <- drm(Resp2~Prop, data=set3_drcdata_forfitting, fct=LL2.4())
fit_set3_rep3 <- drm(Resp3~Prop, data=set3_drcdata_forfitting, fct=LL2.4())

set3_drcdata_forfitting$Resp1Ren <- (set3_drcdata_forfitting$Resp1/coef(fit_set3_rep1)[3])*100
set3_drcdata_forfitting$Resp2Ren <- (set3_drcdata_forfitting$Resp2/coef(fit_set3_rep2)[3])*100
set3_drcdata_forfitting$Resp3Ren <- (set3_drcdata_forfitting$Resp3/coef(fit_set3_rep3)[3])*100

fit_set3_rep1_ren <- drm(Resp1Ren ~ Prop, data=set3_drcdata_forfitting, fct=LL2.4())
fit_set3_rep2_ren <- drm(Resp2Ren ~ Prop, data=set3_drcdata_forfitting, fct=LL2.4())
fit_set3_rep3_ren <- drm(Resp3Ren ~ Prop, data=set3_drcdata_forfitting, fct=LL2.4())




# Matrix with the adjusted Budget and the corresponding adjusted proportion
# (corresponds to the first two columns in Table 3 in the paper)
adjustedeverythingmatrix <- matrix(ncol = 2, nrow=6)
adjustedeverythingmatrix[,1] <- c(
  adjusted_budget_set3_resp_rep1, adjusted_budget_set3_med_rep1,
  adjusted_budget_set3_resp_rep2, adjusted_budget_set3_med_rep2,
  adjusted_budget_set3_resp_rep3, adjusted_budget_set3_med_rep3
)
adjustedeverythingmatrix[,2] <- 100/c(
  adjusted_budget_set3_resp_rep1, adjusted_budget_set3_med_rep1,
  adjusted_budget_set3_resp_rep2, adjusted_budget_set3_med_rep2,
  adjusted_budget_set3_resp_rep3, adjusted_budget_set3_med_rep3
)*0.1

colnames(adjustedeverythingmatrix) <- c("Adjusted Budget", "Adjusted Proportion")
rownames(adjustedeverythingmatrix) <- c("Rep1-Responsible", "Rep1-Median",
                                        "Rep2-Responsible", "Rep2-Median",
                                        "Rep3-Responsible", "Rep3-Median")
print(round(adjustedeverythingmatrix, 3))



# Matrix with the final predicted viability value 
# (corresponds to the last column in Table 3 in the paper)

pred_viabilitymatrix <- matrix(nrow = 2, ncol = 3)
pred_viabilitymatrix[,1] <- predict(fit_set3_rep1_ren, 
                                    newdata=data.frame(Prop = 100/c(adjusted_budget_set3_resp_rep1, 
                                                                    adjusted_budget_set3_med_rep1) *0.1))
pred_viabilitymatrix[,2] <- predict(fit_set3_rep2_ren, 
                                    newdata=data.frame(Prop = 100/c(adjusted_budget_set3_resp_rep2, 
                                                                    adjusted_budget_set3_med_rep2) *0.1))
pred_viabilitymatrix[,3] <- predict(fit_set3_rep3_ren, 
                                    newdata=data.frame(Prop = 100/c(adjusted_budget_set3_resp_rep3, 
                                                                    adjusted_budget_set3_med_rep3) *0.1))
rownames(pred_viabilitymatrix) <- c("Responsible", "Median")
colnames(pred_viabilitymatrix) <- c("Rep1", "Rep2", "Rep3")
print(round(pred_viabilitymatrix, 1))








