#Functions for calculating lake metabolism (MLE)
#Input: DateTime_UTC, doobs, dosat, kgas, zmix, lux, wtr, dummy

#NLL function with linear GPP-light and respiration-water temperature relationships
nllfn <- function (pars, datain) {
  
  Pmax <- exp(pars[1])
  Rmax <- exp(pars[2])
  init_conc <- exp(pars[3])
  
  #Define variables from datain
  nobs <- dim(datain)[1]
  irr <- datain$lux
  doobs <- datain$doobs
  k.gas <- datain$kgas
  Rwtr <- datain$wtr
  zmix <- datain$zmix
  dummy <- datain$dummy
  dosat <- datain$dosat
  
  #Set up output
  dohat <- rep(NA,nobs)
  atmflux <- rep(NA,nobs)
  dohat[1] <- init_conc
  
  #Metabolism model
  for (i in 1:(nobs-1)) {
    atmflux[i] <- dummy[i] * -k.gas[i] * (dohat[i] - dosat[i]) / zmix[i]
    dohat[i+1] <- dohat[i] + (Pmax*irr[i]) - (Rmax*1.073^(Rwtr[i]-20)) + atmflux[i]
  }
  
  #Calculation of residuals, sigma2 and negative log likelihood
  res <- doobs - dohat
  nres <- length(res)
  SSE <- sum(res^2)
  sigma2 <- SSE/nres
  NLL <- -sum(dnorm(doobs, dohat, sd=sqrt(sigma2), log=TRUE)) 
  return(NLL)
}

#Function for obtaining predicted oxygen
domodel <- function(pars, datain) {
  
  nobs <- dim(datain)[1]
  irr <- datain$lux
  doobs <- datain$doobs
  k.gas <- datain$kgas
  Rwtr <- datain$wtr
  zmix <- datain$zmix
  dummy <- datain$dummy
  dosat <- datain$dosat
  
  dopred <- rep(NA,nobs)
  atmflux <- rep(NA,nobs)
  dopred[1] <- pars$init_conc
  
  for (i in 1:(nobs-1)) {
    atmflux[i] <- dummy[i] * -k.gas[i] * (dopred[i] - dosat[i]) / zmix[i]  
    dopred[i+1] <- dopred[i] + (pars$gppcoef*irr[i]) - (pars$rcoef*1.073^(Rwtr[i]-20)) + atmflux[i]
  }
  
  return(dopred)
}

#Function collecting metabolism tools
metab_calc <- function(df){
  datain <- df
  
  parguess <- log(c(5E-5, 2E-3, datain$doobs[1])) #1.3E-4 #7E-6
  
  fit <- tryCatch(optim(parguess, nllfn, datain = datain, method = "Nelder-Mead"), error = function(err){NULL}) #BFGS
  if(is.null(fit)){return(NA)}
  
  gppcoef <- exp(fit$par[1])
  rcoef <- exp(fit$par[2])
  init_conc <- exp(fit$par[3])
  convergence <- fit$convergence
  
  GPP <- mean(gppcoef*datain$lux)*144
  R <- mean(rcoef*1.073^(datain$wtr-20))*144
  NEP <- GPP - R
  
  pars <- list(gppcoef = gppcoef, rcoef = rcoef, init_conc = init_conc) 
  dopred <- domodel(pars = pars, datain = datain)
  r_spear <- cor(dopred, datain$doobs, method = "spearman")
  rmse <- sqrt(mean((datain$doobs-dopred)^2))
  
  daily <- data.frame(DateTime_UTC_min = min(datain$DateTime_UTC), DateTime_UTC_max = max(datain$DateTime_UTC),
                      gppcoef = gppcoef, rcoef = rcoef, init_conc = init_conc, convergence = convergence,
                      GPP = GPP, R = R, NEP = NEP, r_spear = r_spear, rmse = rmse, 
                      wtr_mean = mean(datain$wtr), lux_mean = mean(datain$lux))
  
  obs_pred <- data.frame(DateTime_UTC = datain$DateTime_UTC, dopred = dopred, doobs = datain$doobs)
  
  return(list("daily" = daily, "obs_pred" = obs_pred))
}

#Function for calculating gas exchange velocity as the mean of three empirical models
k_gas_ensemble <- function(wnd, wtr, area){
  k600_cole <- k.cole.base(wnd)
  k600_crucius <- k.crusius.base(wnd)
  k600_vachon <- k.vachon.base(wnd, area)
  k_cole <- k600.2.kGAS.base(k600_cole, wtr, gas="O2")
  k_crucius <- k600.2.kGAS.base(k600_crucius, wtr, gas="O2")
  k_vachon <- k600.2.kGAS.base(k600_vachon, wtr, gas="O2")
  k_mean <- mean(k_cole, k_crucius, k_vachon)
  return(k_mean)
}
k_gas_ensemble_vec <- Vectorize(k_gas_ensemble)

#Function for calculating chemical enhancement factor
k_enhance_factor #<- 