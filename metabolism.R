source("rawdata.R")
source("metabolism_functions.R")

#Metabolism 2019
metab_data_2019 <- sensor_data_2019 |> 
  left_join(wnd_predict) |> 
  left_join(zmix_2019) |> 
  left_join(dic_2019[, c("datetime", "wtr_dic", "ph", "dic", "anc_predicted")]) |> 
  mutate(dosat = o2.at.sat.base(Temp_37),
         oxygen_depth = sensor_depth$oxygen$`2019`,
         oxygen_dummy = as.numeric(zmix > oxygen_depth),
         dic_depth = sensor_depth$ph$`2019`,
         dic_dummy = as.numeric(zmix > dic_depth),
         date = as_date(datetime),
         lux = Light_8,
         kgas_o2_m_day = k_gas_ensemble_vec(wnd_pred_10, Temp_8, "O2"),
         kgas_o2 = kgas_o2_m_day/24/6,
         kgas_co2_m_day = k_gas_ensemble_vec(wnd_pred_10, Temp_8, "CO2"),
         chem_enh = k_gas_enchance(kgas_co2_m_day, Temp_8, ph),
         kgas_co2 = kgas_co2_m_day*chem_enh/24/6,
         calcification = (c(0, diff(anc_predicted))/2)) |> 
  select(date, datetime, doobs = `DO_mg/L_39`, dosat, lux, wtr = Temp_37, zmix, 
         oxygen_dummy, dic_dummy, kgas_o2, kgas_co2, dic, wtr_dic, ph, calcification)

daily_metab_2019 <- metab_data_2019 |> 
  group_by(date) |> 
  nest() |> 
  mutate(n = map_int(data, nrow)) |> 
  filter(n == 144) |> 
  ungroup() |> 
  mutate(depth = sensor_depth$total$`2019`) |> 
  mutate(dic_mle = map(data, ~dic_metab(.x)),
         dic_daily = map(dic_mle, ~.x$dic_daily),
         dic_predict = map(dic_mle, ~.x$dic_predict)) |>
  mutate(oxygen_mle = map(data, ~oxygen_metab(.x)),
         oxygen_daily = map(oxygen_mle, ~.x$oxygen_daily),
         oxygen_predict = map(oxygen_mle, ~.x$oxygen_predict))

daily_metab_2019 |> 
  unnest(dic_predict) |>  
  ggplot(aes(datetime))+
  geom_line(aes(y = dic), col = "black")+
  geom_line(aes(y = dic_pred), col = "orange")

daily_metab_2019 |> 
  unnest(oxygen_predict) |>  
  ggplot(aes(datetime))+
  geom_line(aes(y = doobs), col = "black")+
  geom_line(aes(y = dopred), col = "orange")

#Metabolism 2020
metab_data_2020 <- bind_rows(sensor_data_2020_1, sensor_data_2020_2, sensor_data_2020_3) |> 
  left_join(wnd_predict) |> 
  left_join(zmix_2020) |> 
  mutate(dosat = o2.at.sat.base(Temp_10748214),
         oxygen_dummy = as.numeric(zmix > oxygen_depth),
         date = as_date(datetime),
         lux = Light_10748214,
         kgas_o2_m_day = k_gas_ensemble_vec(wnd_pred_10, Temp_10748214, "O2"),
         kgas_o2 = kgas_o2_m_day/24/6) |> 
  select(date, depth, datetime, doobs = `DO_mgL_500098`, dosat, lux, wtr = Temp_10748214, 
         zmix, oxygen_dummy, kgas_o2)

daily_metab_2020 <- metab_data_2020 |> 
  group_by(date, depth) |> 
  nest() |> 
  mutate(n = map_int(data, nrow)) |> 
  filter(n == 144) |>
  ungroup() |> 
  mutate(oxygen_mle = map(data, ~oxygen_metab(.x)),
         oxygen_daily = map(oxygen_mle, ~.x$oxygen_daily),
         oxygen_predict = map(oxygen_mle, ~.x$oxygen_predict))

daily_metab_2020 |> 
  unnest(oxygen_predict) |>  
  ggplot(aes(datetime))+
  geom_line(aes(y = doobs), col = "black")+
  geom_line(aes(y = dopred), col = "orange")

#Save daily metabolism calculations
daily_metab_all <- list("2019" = daily_metab_2019, "2020" = daily_metab_2020)

saveRDS(daily_metab_all, "data/daily_metab.rds")
