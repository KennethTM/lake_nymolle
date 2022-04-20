source("libs_and_funcs.R")

#Rawdata

#Wind data from 2020
#Sensor malfunction after few days, predict wind speed using dmi data

wnd <- read_csv("data/nymolle_wind.csv", skip=1)

wnd_clean <- wnd |> 
  select(2:3) |> 
  set_names(c("datetime_gmt2", "wnd")) |> 
  mutate(datetime = round_date(mdy_hms(datetime_gmt2), "10 mins") - 2*60*60) |> 
  select(datetime, wnd) |> 
  filter(datetime < ymd_hms("2020-05-04 00:00:00"))

dmi_wnd_1 <- read_csv("data/dmi_wind_1.csv")
dmi_wnd_2 <- read_csv("data/dmi_wind_2.csv")

dmi_wnd <- bind_rows(dmi_wnd_1, dmi_wnd_2) |> 
  rename(dmi_wnd = wnd)

#Fit model for wind speed between lake and dmi weather station
wnd_model_data <- wnd_clean |> 
  left_join(dmi_wnd)

wnd_m0 <- lm(wnd ~ 1, data = wnd_model_data)
wnd_m1 <- lm(wnd ~ dmi_wnd - 1, data = wnd_model_data)
wnd_m2 <- lm(wnd ~ dmi_wnd + I(dmi_wnd^2) - 1, data = wnd_model_data)

anova(wnd_m0, wnd_m1, wnd_m2)
summary(wnd_m2)

#Predict for all observations
wnd_predict <- dmi_wnd |> 
  mutate(wnd_pred = predict(wnd_m2, newdata = data.frame(dmi_wnd = dmi_wnd))) |> 
  select(-dmi_wnd)

#Calibration curves for alk versus specifik conductivity
calcurve_2019 <- tribble(~spec_cond, ~alk,
                         453.6,	1.527,
                         314.4,	1.036,
                         236.2,	0.779,
                         193.6,	0.6405,
                         162.3,	0.532)

calcurve_2019_model <- lm(alk ~ spec_cond, data = calcurve_2019)
summary(calcurve_2019_model)

calcurve_2020 <- tribble(~spec_cond, ~alk,
                         465.5,	1.404,
                         327.6,	0.946,
                         242.3,	0.701,
                         193.1,	0.526,
                         157.6,	0.461)
calcurve_2020_model <- lm(alk ~ spec_cond, data = calcurve_2020)
summary(calcurve_2020_model)

#Sensor data
sensor_data_2019 <- read_excel("data/sensor_data_2019.xlsx")
sensor_data_2020_1 <- read_excel("data/sensor_data_2020.xlsx", sheet = 1)
sensor_data_2020_2 <- read_excel("data/sensor_data_2020.xlsx", sheet = 2)
sensor_data_2020_3 <- read_excel("data/sensor_data_2020.xlsx", sheet = 3)

#Water temperature data
wtr_2019 <- sensor_data_2019 |> 
  select(datetime = Date_time, wtr_1 = Temp_8, wtr_2 = Temp_21, wtr_3 = Temp_29, wtr_4 = Temp_37, wtr_5 = Temp_49)

wtr_2020_1 <- sensor_data_2020_1 |> 
  select(datetime = Date_time, wtr_1 = Temp_10748214, wtr_4 = Temp_10748223, wtr_5 = Temp_10675577)

wtr_2020_2 <- sensor_data_2020_2 |> 
  select(datetime = Date_time, wtr_1 = Temp_10748214, wtr_3 = Temp_10748223, wtr_4 = Temp_10675577, wtr_5 = Temp_10748206)

wtr_2020_3 <- sensor_data_2020_3 |> 
  select(datetime = Date_time, wtr_1 = Temp_10748214, wtr_4 = Temp_10675577, wtr_5 = Temp_10748206)

wtr_all <- bind_rows(wtr_2019, wtr_2020_1, wtr_2020_2, wtr_2020_3)

#Oxygen sensor data
oxygen_2019 <- sensor_data_2019 |> 
  select(datetime = Date_time, oxygen_2 = `DO_%_39`, oxygen_3 = `DO_%_65`)

oxygen_2020_1 <- sensor_data_2020_1 |> 
  select(datetime = Date_time, oxygen_1 = DO_procent_500098, oxygen_2 = DO_procent_348923)

oxygen_2020_2 <- sensor_data_2020_2 |> 
  select(datetime = Date_time, oxygen_1 = DO_procent_500098, oxygen_2 = DO_procent_348923, oxygen_3 = DO_procent_652844)

oxygen_2020_3 <- sensor_data_2020_3 |> 
  mutate(dosat = o2.at.sat.base(Temp_10675577),
         DO_procent_652844 = DO_mgL_652844/dosat*100) |>
  select(datetime = Date_time, oxygen_1 = DO_procent_500098, oxygen_2 = DO_procent_348923, oxygen_3 = DO_procent_652844)

oxygen_all <- bind_rows(oxygen_2019, oxygen_2020_1, oxygen_2020_2, oxygen_2020_3)
