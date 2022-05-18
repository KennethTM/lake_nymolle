#Functions for calculating lake metabolism (MLE)

#Jassby-Platt GPP-light relationship
GPP <- function(Pmax, Palpha, light){
  return(Pmax*tanh(Palpha*light/Pmax))
}

#Arrhenius temperature dependant respiration with Q10 ~ 2
R <- function(Rmax, wtr, temp_const = 1.073){
  return(Rmax*temp_const^(wtr-20))
}

#NLL function
NLL <- function(obs, hat){
  res <- obs - hat
  nres <- length(res)
  SSE <- sum(res^2)
  sigma2 <- SSE/nres
  NLL <- -sum(dnorm(obs, hat, sd=sqrt(sigma2), log=TRUE))
  return(NLL)
}

#Oxygen metabolism model
#Input for oxygen: datetime, doobs, dosat, kgas, zmix, lux, wtr, dummy
oxygen_nll <- function (pars, datain) {
  
  Pmax <- exp(pars[1])
  Palpha <- exp(pars[2])
  Rmax <- exp(pars[3])

  #Define variables from datain
  nobs <- nrow(datain)
  irr <- datain$lux
  doobs <- datain$doobs
  k.gas <- datain$kgas_o2
  wtr <- datain$wtr
  zmix <- datain$zmix
  dummy <- datain$oxygen_dummy
  dosat <- datain$dosat
  
  #Set up output
  dohat <- rep(NA,nobs)
  atmflux <- rep(NA,nobs)
  dohat[1] <- doobs[1]
  
  #Metabolism model
  for (i in 1:(nobs-1)) {
    atmflux[i] <- dummy[i] * -k.gas[i] * (dohat[i] - dosat[i]) / zmix[i]
    dohat[i+1] <- dohat[i] + GPP(Pmax, Palpha, irr[i]) - R(Rmax, wtr[i]) + atmflux[i]
  }
  
  #Calculation of residuals, sigma2 and negative log likelihood
  nll <- NLL(doobs, dohat)
  return(nll)
}

#Function for obtaining predicted oxygen
oxygen_predict <- function(pars, datain) {
  
  nobs <- nrow(datain)
  irr <- datain$lux
  doobs <- datain$doobs
  k.gas <- datain$kgas_o2
  wtr <- datain$wtr
  zmix <- datain$zmix
  dummy <- datain$oxygen_dummy
  dosat <- datain$dosat
  
  dopred <- rep(NA,nobs)
  atmflux <- rep(NA,nobs)
  dopred[1] <- doobs[1]
  
  for (i in 1:(nobs-1)) {
    atmflux[i] <- dummy[i] * -k.gas[i] * (dopred[i] - dosat[i]) / zmix[i]  
    dopred[i+1] <- dopred[i] + GPP(pars$Pmax, pars$Palpha, irr[i]) - R(pars$Rmax, wtr[i]) + atmflux[i]
  }
  
  return(dopred)
}

#Function calculating oxygen metabolism
oxygen_metab <- function(datain){

  if(year(datain$date[1]) == 2019){
    parguess <- log(c(0.25, 5e-6, 0.1))
  }else{
    parguess <- log(c(0.02, 2e-6, 0.01))
  }
  
  fit <- tryCatch(optim(parguess, oxygen_nll, datain = datain, method="Nelder-Mead", control = list(maxit=2500)), error = function(err){NULL})
  if(is.null(fit)){return(NA)}
  
  Pmax <- exp(fit$par[1])
  Palpha <- exp(fit$par[2])
  Rmax <- exp(fit$par[3])
  convergence <- fit$convergence
  
  GPP <- mean(GPP(Pmax, Palpha, datain$lux))*144
  R <- mean(R(Rmax, datain$wtr))*144
  NEP <- GPP - R
  
  pars <- list(Pmax = Pmax, Palpha = Palpha, Rmax = Rmax)
  dopred <- oxygen_predict(pars = pars, datain = datain)
  r_spear <- cor(dopred, datain$doobs, method = "spearman")
  rmse <- sqrt(mean((datain$doobs-dopred)^2))
  
  daily <- data.frame(datetime_min = min(datain$datetime), datetime_max = max(datain$datetime),
                      Pmax = Pmax, Palpha = Palpha, Rmax = Rmax, convergence = convergence, 
                      GPP = GPP/32*1000, R = R/32*1000, NEP = NEP/32*1000, #convert g/m3/d to mmol/m3/day
                      r_spear = r_spear, rmse = rmse, 
                      wtr_mean = mean(datain$wtr), lux_mean = mean(datain$lux))
  
  predictions <- data.frame(datetime = datain$datetime, dopred = dopred, doobs = datain$doobs)
  
  return(list("oxygen_daily" = daily, "oxygen_predict" = predictions))
}

#DIC metabolism model
#Input for dic: datetime, kgas, zmix, lux, wtr_dic, dummy, dic, pH
dic_nll <- function(pars, datain){
  
  Pmax <- exp(pars[1])
  Palpha <- exp(pars[2])
  Rmax <- exp(pars[3])

  nobs <- nrow(datain)
  irr <- datain$lux
  k.gas <- datain$kgas_co2
  wtr <- datain$wtr_dic
  zmix <- datain$zmix
  dummy <- datain$dic_dummy
  dic <- datain$dic
  ph <- datain$ph
  calc <- datain$calcification

  dic_hat <- rep(NA,nobs)
  atmflux <- rep(NA,nobs)
  dic_hat[1] <- dic[1]

  for (i in 1:(nobs-1)){
    carb_sys <- aquaenv(S=0, t=wtr[i], SumCO2 = dic_hat[i], pH = ph[i])
    atmflux[i] <- dummy[i] * -k.gas[i] * (carb_sys$CO2 - carb_sys$CO2_sat) / zmix[i]
    dic_hat[i+1] <- dic_hat[i] - GPP(Pmax, Palpha, irr[i]) + R(Rmax, wtr[i]) + atmflux[i] + calc[i]
    
  }
  
  nll <- NLL(dic, dic_hat)
  return(nll)
}

#Function for obtaining predicted dic
dic_predict <- function(pars, datain) {
  
  nobs <- nrow(datain)
  irr <- datain$lux
  k.gas <- datain$kgas_co2
  wtr <- datain$wtr_dic
  zmix <- datain$zmix
  dummy <- datain$dic_dummy
  dic <- datain$dic
  ph <- datain$ph
  calc <- datain$calcification
  
  dic_pred <- rep(NA,nobs)
  atmflux <- rep(NA,nobs)
  dic_pred[1] <- dic[1]
  
  for (i in 1:(nobs-1)){
    carb_sys <- aquaenv(S=0, t=wtr[i], SumCO2 = dic_pred[i], pH = ph[i])
    atmflux[i] <- dummy[i] * -k.gas[i] * (carb_sys$CO2 - carb_sys$CO2_sat) / zmix[i]
    dic_pred[i+1] <- dic_pred[i] - GPP(pars$Pmax, pars$Palpha, irr[i]) + R(pars$Rmax, wtr[i]) + atmflux[i] + calc[i]
  }
  
  return(dic_pred)
}

#Function collecting DIC metabolism functions
dic_metab <- function(datain){

  parguess <- log(c(7E-6, 2e-9, 2E-6))
  
  fit <- tryCatch(optim(parguess, dic_nll, datain = datain, method="Nelder-Mead", control = list(maxit=2500)), error = function(err){NULL}) #BFGS
  if(is.null(fit)){return(NA)}
  
  Pmax <- exp(fit$par[1])
  Palpha <- exp(fit$par[2])
  Rmax <- exp(fit$par[3])
  convergence <- fit$convergence
  
  GPP <- mean(GPP(Pmax, Palpha, datain$lux))*144
  R <- mean(R(Rmax, datain$wtr_dic))*144
  NEP <- GPP - R
  
  pars <- list(Pmax = Pmax, Palpha = Palpha, Rmax = Rmax) 
  dic_pred <- dic_predict(pars = pars, datain = datain)
  r_spear <- cor(dic_pred, datain$dic, method = "spearman")
  rmse <- sqrt(mean((datain$dic-dic_pred)^2))
  
  daily <- data.frame(datetime_min = min(datain$datetime), datetime_max = max(datain$datetime),
                      Pmax = Pmax, Palpha = Palpha, Rmax = Rmax, convergence = convergence,
                      GPP = GPP*10^6, R = R*10^6, NEP = NEP*10^6, r_spear = r_spear, rmse = rmse, #units to mmol/m3/day
                      wtr_mean = mean(datain$wtr_dic), lux_mean = mean(datain$lux))
  
  predictions <- data.frame(datetime = datain$datetime, dic_pred = dic_pred, dic = datain$dic)
  
  return(list("dic_daily" = daily, "dic_predict" = predictions))
}

k_schilder <- function(wnd){
  k600 <- 0.9+0.97*wnd 
  k600_m_day <- k600*24/100
  return(k600_m_day)
}

#Function for calculating gas exchange velocity (m/day) as the mean of three empirical models
k_gas_ensemble <- function(wnd, wtr, gas){
  k600_cole <- k.cole.base(wnd)
  k600_crucius <- k.crusius.base(wnd)
  k600_schilder <- k_schilder(wnd)
  k600_all <- c(k600_schilder, k600_crucius, k600_cole)
  k_all <- sapply(k600_all, function(k600){k600.2.kGAS.base(k600, wtr, gas=gas)})
  k_all_mean <- mean(k_all)
  return(k_all_mean)
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
  
  #Calculate chemical enhancement factor alpha
  alpha <- t/((t-1) + (tanh(q*z)/(q*z)))
  
  return(alpha)
  
}
