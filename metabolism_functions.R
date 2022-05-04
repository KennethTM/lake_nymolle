#Functions for calculating lake metabolism (MLE)
#Input for oxygen: datetime, doobs, dosat, kgas, zmix, lux, wtr, dummy
#Input for oxygen: datetime, kgas, zmix, lux, wtr_dic, dummy, dic, pH

#Oxygen metabolism model
#NLL function with linear GPP-light and respiration-water temperature relationships
oxygen_nll <- function (pars, datain) {
  
  Pmax <- exp(pars[1])
  Rmax <- exp(pars[2])
  doinit <- exp(pars[3])
  
  #Define variables from datain
  nobs <- nrow(datain)
  irr <- datain$lux
  doobs <- datain$doobs
  k.gas <- datain$kgas_o2
  Rwtr <- datain$wtr
  zmix <- datain$zmix
  dummy <- datain$dummy
  dosat <- datain$dosat
  
  #Set up output
  dohat <- rep(NA,nobs)
  atmflux <- rep(NA,nobs)
  dohat[1] <- doinit
  
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
oxygen_predict <- function(pars, datain) {
  
  nobs <- nrow(datain)
  irr <- datain$lux
  doobs <- datain$doobs
  k.gas <- datain$kgas_o2
  Rwtr <- datain$wtr
  zmix <- datain$zmix
  dummy <- datain$dummy
  dosat <- datain$dosat
  
  dopred <- rep(NA,nobs)
  atmflux <- rep(NA,nobs)
  dopred[1] <- pars$doinit
  
  for (i in 1:(nobs-1)) {
    atmflux[i] <- dummy[i] * -k.gas[i] * (dopred[i] - dosat[i]) / zmix[i]  
    dopred[i+1] <- dopred[i] + (pars$gppcoef*irr[i]) - (pars$rcoef*1.073^(Rwtr[i]-20)) + atmflux[i]
  }
  
  return(dopred)
}

#Function calculating oxygen metabolism
oxygen_metab <- function(df){
  datain <- df
  
  parguess <- log(c(5E-5, 2E-3, datain$doobs[1])) #1.3E-4 #7E-6
  
  fit <- tryCatch(optim(parguess, oxygen_nll, datain = datain, method = "Nelder-Mead"), error = function(err){NULL}) #BFGS
  if(is.null(fit)){return(NA)}
  
  gppcoef <- exp(fit$par[1])
  rcoef <- exp(fit$par[2])
  doinit <- exp(fit$par[3])
  convergence <- fit$convergence
  
  GPP <- mean(gppcoef*datain$lux)*144
  R <- mean(rcoef*1.073^(datain$wtr-20))*144
  NEP <- GPP - R
  
  pars <- list(gppcoef = gppcoef, rcoef = rcoef, doinit = doinit) 
  dopred <- oxygen_predict(pars = pars, datain = datain)
  r_spear <- cor(dopred, datain$doobs, method = "spearman")
  rmse <- sqrt(mean((datain$doobs-dopred)^2))
  
  daily <- data.frame(datetime_min = min(datain$datetime), datetime_max = max(datain$datetime),
                      gppcoef = gppcoef, rcoef = rcoef, doinit = doinit, convergence = convergence,
                      GPP = GPP, R = R, NEP = NEP, r_spear = r_spear, rmse = rmse, 
                      wtr_mean = mean(datain$wtr), lux_mean = mean(datain$lux))
  
  obs_pred <- data.frame(datetime = datain$datetime, dopred = dopred, doobs = datain$doobs)
  
  return(list("oxygen_daily" = daily, "oxygen_predict" = obs_pred))
}

#DIC metabolism model
dic_nll <- function (pars, datain) {
  
  Pmax <- exp(pars[1])
  Rmax <- exp(pars[2])
  dic_init <- exp(pars[3])
  
  nobs <- nrow(datain)
  irr <- datain$lux
  k.gas <- datain$kgas_co2
  Rwtr <- datain$wtr_dic
  zmix <- datain$zmix
  dummy <- datain$dummy
  dic <- datain$dic
  ph <- datain$ph
  calc <- datain$calcification

  #Set up output
  dic_hat <- rep(NA,nobs)
  atmflux <- rep(NA,nobs)
  dic_hat[1] <- dic_init

  #Metabolism model
  for (i in 1:(nobs-1)){
    carb_sys <- aquaenv(S=0, t=Rwtr[i], SumCO2 = dic_hat[i], pH = ph[i])
    atmflux[i] <- dummy[i] * -k.gas[i] * (carb_sys$CO2 - carb_sys$CO2_sat) / zmix[i]
    dic_hat[i+1] <- dic_hat[i] - (Pmax*irr[i]) + (Rmax*1.073^(Rwtr[i]-20)) + atmflux[i] + calc[i]
  }
  
  #Calculation of residuals, sigma2 and negative log likelihood
  res <- dic - dic_hat
  nres <- length(res)
  SSE <- sum(res^2)
  sigma2 <- SSE/nres
  NLL <- -sum(dnorm(dic, dic_hat, sd=sqrt(sigma2), log=TRUE)) 
  return(NLL)
}

#Function for obtaining predicted oxygen
dic_predict <- function(pars, datain) {
  
  nobs <- nrow(datain)
  irr <- datain$lux
  k.gas <- datain$kgas_co2
  Rwtr <- datain$wtr_dic
  zmix <- datain$zmix
  dummy <- datain$dummy
  dic <- datain$dic
  ph <- datain$ph
  calc <- datain$calcification
  
  dic_pred <- rep(NA,nobs)
  atmflux <- rep(NA,nobs)
  dic_pred[1] <- pars$dic_init
  
  for (i in 1:(nobs-1)){
    carb_sys <- aquaenv(S=0, t=Rwtr[i], SumCO2 = dic_pred[i], pH = ph[i])
    atmflux[i] <- dummy[i] * -k.gas[i] * (carb_sys$CO2 - carb_sys$CO2_sat) / zmix[i]
    dic_pred[i+1] <- dic_pred[i] - (pars$gppcoef*irr[i]) + (pars$rcoef*1.073^(Rwtr[i]-20)) + atmflux[i] + calc[i]
  }
  
  return(dic_pred)
}

#Function collecting DIC metabolism functions
dic_metab <- function(df){
  datain <- df
  
  parguess <- log(c(5E-6, 3E-6, datain$dic[1]))
  
  fit <- tryCatch(optim(parguess, dic_nll, datain = datain, method = "Nelder-Mead"), error = function(err){NULL}) #BFGS
  if(is.null(fit)){return(NA)}
  
  gppcoef <- exp(fit$par[1])
  rcoef <- exp(fit$par[2])
  dic_init <- exp(fit$par[3])
  convergence <- fit$convergence
  
  GPP <- mean(gppcoef*datain$lux)*144
  R <- mean(rcoef*1.073^(datain$wtr_dic-20))*144
  NEP <- GPP - R
  
  pars <- list(gppcoef = gppcoef, rcoef = rcoef, dic_init = dic_init) 
  dic_pred <- dic_predict(pars = pars, datain = datain)
  r_spear <- cor(dic_pred, datain$dic, method = "spearman")
  rmse <- sqrt(mean((datain$dic-dic_pred)^2))
  
  daily <- data.frame(datetime_min = min(datain$datetime), datetime_max = max(datain$datetime),
                      gppcoef = gppcoef, rcoef = rcoef, dic_init = dic_init, convergence = convergence,
                      GPP = GPP, R = R, NEP = NEP, r_spear = r_spear, rmse = rmse, 
                      wtr_mean = mean(datain$wtr_dic), lux_mean = mean(datain$lux))
  
  obs_pred <- data.frame(datetime = datain$datetime, dic_pred = dic_pred, dic = datain$dic)
  
  return(list("dic_daily" = daily, "dic_predict" = obs_pred))
}

#Function for calculating gas exchange velocity as the mean of three empirical models
k_gas_ensemble <- function(wnd, wtr, area, gas){
  k600_cole <- k.cole.base(wnd)
  k600_crucius <- k.crusius.base(wnd)
  k600_vachon <- k.vachon.base(wnd, area)
  k_cole <- k600.2.kGAS.base(k600_cole, wtr, gas=gas)
  k_crucius <- k600.2.kGAS.base(k600_crucius, wtr, gas=gas)
  k_vachon <- k600.2.kGAS.base(k600_vachon, wtr, gas=gas)
  k_mean <- mean(k_cole, k_crucius, k_vachon)
  return(k_mean)
}
k_gas_ensemble_vec <- Vectorize(k_gas_ensemble)

#Calculation of chemical enhancement (Hoover and Berkshire mode, Wanninkhof & Knox 1996)
k_gas_enchance <- function(kco2, wtr, ph, S=0){
  
  wtr_k <- wtr + 273.15
  
  #Constants kco2 (r1, s-1) and kd (r2, l/(mol*s)) Johnson L&O 1982
  r1 <- exp(1246.98 + 0*S^0.5 + (-6.19*10^4)/wtr_k + -183.0*log(wtr_k))
  r2 <- exp(1346.24 + -0.126*S^0.5 + (-6.44*10^4)/wtr_k + -196.4*log(wtr_k))
  
  ah <- 10^(-ph)
  
  #Dickson & millero 1987, S between 0 and 40
  pk1 <- (-840.39/wtr_k + 19.894 - 3.0189*log(wtr_k))*S^0.5 + 0.00668*S + 6320.81/wtr_k - 126.3405 + 19.568*log(wtr_k)
  k1 <- 10^(-pk1)

  pk2 <- (-690.59/wtr_k + 17.176 - 2.6719*log(wtr_k))*S^0.5 + 0.0217*S + 5143.69/wtr_k - 90.1833 + 14.613*log(wtr_k)
  k2 <- 10^(-pk2)
  
  #Kw from seacarb r package, S between 0 and 45 and temperature between 0 and 45
  lnKw_1 <- -13847.26/wtr_k + 148.9802 - 23.6521 * log(wtr_k)
  lnKw_2 <- (118.67/wtr_k - 5.977 + 1.0495 * log(wtr_k)) * sqrt(S) - 0.01615 * S
  lnKw <- lnKw_1 + lnKw_2
  kw <- exp(lnKw)

  oh <- kw/ah
  
  r <- r1 + r2*oh  
  
  t <- 1 + ah^2/(k1*k2 + k1*ah)
  
  #Jahne 1987, co2 diffusion coeffcient
  d <- 5019 * exp(-19.51/(0.00831451*wtr_k)) * 10^-5 #unit cm2/s
  
  kco2_cm_s <- kco2*100/(24*60*60) #convert unit m/d to cm/s
  
  z <- d/kco2_cm_s
  
  q <- (r*t/d)^0.5
  
  alpha <- t/((t-1) + (tanh(q*z)/(q*z)))
  
  return(alpha)
  
}


#Make units of daily metabolism estimates as mmol/m3/d (*depth to get m2) 
#(dic is in mol/l/d or *10^6 mmol/m3/d - convert before metabolism calc, o2 is g/m3/d convert by /32*1000 )
