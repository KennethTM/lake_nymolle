source("rawdata.R")
source("metabolism_functions.R")

#Calculation lake metabolism from diel oxygen concentrations
#Oxygen metabolism
#Input: DateTime_UTC, doobs, dosat, kgas, zmix, lux, wtr, dummy
depth_2019 <- 0.66

zmix_2019 <- sensor_data_2019 |> 
  select(datetime = Date_time, wtr_0.08 = Temp_8, wtr_0.21 = Temp_21, wtr_0.29 = Temp_29, wtr_0.37 = Temp_37, wtr_0.49 = Temp_49) |> 
  as.data.frame() |> 
  ts.thermo.depth() |> 
  as_tibble() |> 
  mutate(zmix = ifelse(is.na(thermo.depth), depth_2019, thermo.depth)) |> 
  select(-thermo.depth)

metab_data_2019 <- sensor_data_2019 |> 
  left_join(wnd_predict, by = c("Date_time" = "datetime")) |> 
  left_join(zmix_2019, by = c("Date_time" = "datetime")) |> 
  mutate(dosat = o2.at.sat.base(Temp_37),
         oxygen_depth = 0.39,
         dummy = as.numeric(zmix > oxygen_depth),
         wnd_pred_10 = wind.scale.base(wnd_pred, 1),
         date = as_date(Date_time)) |> 
  mutate(kgas = k_gas_ensemble_vec(wnd_pred_10, Temp_8, 71082.35)/24/6) |> 
  select(date, DateTime_UTC = Date_time, doobs = `DO_mg/L_39`, dosat, lux = Light_8, 
         wtr = Temp_37, zmix, dummy, kgas)

daily_metab_2019 <- metab_data_2019 |> 
  nest(data = c(DateTime_UTC, doobs, dosat, lux, wtr, zmix, dummy, kgas)) |> 
  mutate(n = map_int(data, nrow),
         depth = 0.66) |> 
  filter(n == 144) |> 
  mutate(metab_mle = map(data, ~metab_calc(.x)),
         daily = map(metab_mle, ~.x$daily),
         obs_pred = map(metab_mle, ~.x$obs_pred))

daily_metab_2019 |> 
  unnest(obs_pred) |>  
  ggplot(aes(DateTime_UTC))+
  geom_line(aes(y = doobs), col = "black")+
  geom_line(aes(y = dopred), col = "orange")

daily_metab_2019 |> 
  unnest(daily) |> View()
