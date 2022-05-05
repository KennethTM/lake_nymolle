source("rawdata.R")
source("metabolism_functions.R")

#Calculation lake metabolism from diel oxygen and dic concentrations
#Oxygen metabolism

metab_data_2019 <- sensor_data_2019 |> 
  left_join(wnd_predict) |> 
  left_join(zmix_2019) |> 
  left_join(select(dic_all, datetime, wtr_dic, ph, dic, anc_predicted)) |> 
  mutate(dosat = o2.at.sat.base(Temp_37),
         oxygen_depth = o2_depth_2019,
         dummy = as.numeric(zmix > oxygen_depth),
         date = as_date(datetime),
         lux = Light_8,
         calcification = (c(0, diff(anc_predicted))/2)*10^6) |> 
  mutate(kgas_o2_m_day = k_gas_ensemble_vec(wnd_pred_10, Temp_8, lake_area, "O2"),
         kgas_co2_m_day = k_gas_ensemble_vec(wnd_pred_10, Temp_8, lake_area, "CO2"),
         kgas_o2 = kgas_o2_m_day/24/6,
         kgas_co2_alpha = k_gas_enchance(kgas_co2_m_day, Temp_37, ph),
         kgas_co2 = (kgas_co2_m_day*kgas_co2_alpha)/24/6,
         dic = dic*10^6) |> 
  select(date, datetime, doobs = `DO_mg/L_39`, dosat, lux, 
         wtr = Temp_37, zmix, dummy, kgas_o2, kgas_co2, dic, wtr_dic, ph, calcification)

daily_metab_2019 <- metab_data_2019 |> 
  nest(data = c(datetime, doobs, dosat, lux, wtr, zmix, dummy, kgas_o2, kgas_co2, dic, wtr_dic, ph, calcification)) |> 
  mutate(n = map_int(data, nrow),
         depth = depth_2019) |> 
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
