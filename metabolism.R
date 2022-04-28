source("rawdata.R")
source("metabolism_functions.R")

#Calculation lake metabolism from diel oxygen concentrations
#Oxygen metabolism
depth_2019 <- 0.66

zmix_2019 <- sensor_data_2019 |> 
  select(datetime, wtr_0.08 = Temp_8, wtr_0.21 = Temp_21, wtr_0.29 = Temp_29, wtr_0.37 = Temp_37, wtr_0.49 = Temp_49) |> 
  as.data.frame() |> 
  ts.thermo.depth() |> 
  as_tibble() |> 
  mutate(zmix = ifelse(is.na(thermo.depth), depth_2019, thermo.depth)) |> 
  select(-thermo.depth)

metab_data_2019 <- sensor_data_2019 |> 
  left_join(wnd_predict) |> 
  left_join(zmix_2019) |> 
  left_join(select(dic_all, datetime, wtr_dic, ph, dic, anc_predicted)) |> 
  mutate(dosat = o2.at.sat.base(Temp_37),
         oxygen_depth = 0.39,
         dummy = as.numeric(zmix > oxygen_depth),
         wnd_pred_10 = wind.scale.base(wnd_pred, 0.5),
         date = as_date(datetime),
         lux = Light_8/10000,
         calcification = c(0, diff(anc_predicted/1000))/2) |> 
  mutate(kgas_o2 = k_gas_ensemble_vec(wnd_pred_10, Temp_8, 71082.35, "O2")/24/6,
         kgas_co2 = k_gas_ensemble_vec(wnd_pred_10, Temp_8, 71082.35, "CO2")/24/6) |> 
  select(date, datetime, doobs = `DO_mg/L_39`, dosat, lux, 
         wtr = Temp_37, zmix, dummy, kgas_o2, kgas_co2, dic, wtr_dic, ph, calcification)

daily_metab_2019 <- metab_data_2019 |> 
  nest(data = c(datetime, doobs, dosat, lux, wtr, zmix, dummy, kgas_o2, kgas_co2, dic, wtr_dic, ph, calcification)) |> 
  mutate(n = map_int(data, nrow),
         depth = 0.66) |> 
  filter(n == 144) |> 
  slice(3) |> 
  mutate(dic_mle = map(data, ~dic_metab(.x)),
         dic_daily = map(dic_mle, ~.x$dic_daily),
         dic_predict = map(dic_mle, ~.x$dic_predict))
  # mutate(oxygen_mle = map(data, ~oxygen_metab(.x)),
  #        oxygen_daily = map(oxygen_mle, ~.x$oxygen_daily),
  #        oxygen_predict = map(oxygen_mle, ~.x$oxygen_predict))

daily_metab_2019 |> 
  unnest(dic_predict) |>  
  ggplot(aes(datetime))+
  geom_line(aes(y = dic), col = "black")+
  geom_line(aes(y = dic_pred), col = "orange")

daily_metab_2019 |> 
  unnest(oxygen_daily) |> View()
