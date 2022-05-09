source("libs_and_funcs.R")

#Rawdata loading and processing

#Wind data from 2020
#Sensor malfunction after few days, predict wind speed using dmi data
wnd_raw <- read_csv("data/nymolle_wind.csv", skip=1)

wnd_clean <- wnd_raw |> 
  select(2:3) |> 
  set_names(c("datetime_gmt2", "wnd")) |> 
  mutate(datetime = round_date(mdy_hms(datetime_gmt2), "10 mins") - 2*60*60) |> 
  select(datetime, wnd) |> 
  filter(datetime < ymd_hms("2020-05-04 00:00:00"))

dmi_wnd_1 <- read_csv("data/dmi_wind_1.csv")
dmi_wnd_2 <- read_csv("data/dmi_wind_2.csv")

dmi_wnd_all <- bind_rows(dmi_wnd_1, dmi_wnd_2) |> 
  rename(wnd_dmi = wnd) |> 
  arrange(datetime)

#Global radiation data from DMI station
dmi_rad_1 <- read_csv("data/dmi_globrad_1.csv")
dmi_rad_2 <- read_csv("data/dmi_globrad_2.csv")

daily_glob_rad <- bind_rows(dmi_rad_1, dmi_rad_2) |> 
  mutate(date = as_date(datetime)) |> 
  group_by(date) |> 
  summarise(globrad = mean(globrad))

#Perform linear interpolation for small gaps
wnd_seq <- data.frame(datetime = seq(min(dmi_wnd_1$datetime), max(dmi_wnd_2$datetime), "10 min"))

dmi_wnd_interp <- left_join(wnd_seq, dmi_wnd_all) |> 
  mutate(wnd_dmi = na.approx(wnd_dmi, na.rm = FALSE, maxgap=18))

#Fit model for wind speed between lake and dmi weather station
wnd_model_data <- left_join(wnd_clean, dmi_wnd_all)

wnd_m0 <- lm(wnd ~ 1, data = wnd_model_data)
wnd_m1 <- lm(wnd ~ wnd_dmi - 1, data = wnd_model_data)
wnd_m2 <- lm(wnd ~ wnd_dmi + I(wnd_dmi^2) - 1, data = wnd_model_data)

anova(wnd_m0, wnd_m1, wnd_m2)
summary(wnd_m2)

#Predict for all observations
wnd_predict <- dmi_wnd_interp |> 
  as_tibble() |> 
  mutate(wnd_pred = predict(wnd_m2, newdata = data.frame(wnd_dmi = wnd_dmi)),
         wnd_pred_10 = wind.scale.base(wnd_pred, 0.5)) |> 
  select(-wnd_dmi, -wnd_pred)

#Calibration curves for alkalinity versus specific conductivity for both 2019 and 2020
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
#Adjust time zone and round datetime
sensor_data_2019 <- read_excel("data/sensor_data_2019.xlsx") |> 
  mutate(datetime = Date_time - 2*60*60,
         datetime = round_date(datetime, "10 mins")) |> 
  select(-Date_time)
sensor_data_2020_1 <- read_excel("data/sensor_data_2020.xlsx", sheet = 1) |> 
  mutate(datetime = Date_time - 1*60*60,
         datetime = round_date(datetime, "10 mins"))|> 
  select(-Date_time)
sensor_data_2020_2 <- read_excel("data/sensor_data_2020.xlsx", sheet = 2) |> 
  mutate(datetime = Date_time - 2*60*60,
         datetime = round_date(datetime, "10 mins"))|> 
  select(-Date_time)
sensor_data_2020_3 <- read_excel("data/sensor_data_2020.xlsx", sheet = 3) |> 
  mutate(datetime = Date_time - 2*60*60,
         datetime = round_date(datetime, "10 mins"))|> 
  select(-Date_time)

#Depth info for each measurement period
sensor_depth <- list("total" = list("2019" = 0.66, "2020_1" = 0.665, "2020_2" = 0.665, "2020_3" = 0.57),
                     "oxygen" = list("2019" = 0.39, "2020_1" = 0.115, "2020_2" = 0.06, "2020_3" = 0.085),
                     "ph" = list("2019" = 0.27, "2020_1" = 0.515, "2020_2" = 0.37, "2020_3" = 0.165))

#Water temperature data
wtr_2019 <- sensor_data_2019 |> 
  select(datetime, wtr_1 = Temp_8, wtr_2 = Temp_21, wtr_3 = Temp_29, wtr_4 = Temp_37, wtr_5 = Temp_49)

wtr_2020_1 <- sensor_data_2020_1 |> 
  select(datetime, wtr_1 = Temp_10748214, wtr_4 = Temp_10748223, wtr_5 = Temp_10675577)

wtr_2020_2 <- sensor_data_2020_2 |> 
  select(datetime, wtr_1 = Temp_10748214, wtr_3 = Temp_10748223, wtr_4 = Temp_10675577, wtr_5 = Temp_10748206)

wtr_2020_3 <- sensor_data_2020_3 |> 
  select(datetime, wtr_1 = Temp_10748214, wtr_4 = Temp_10675577, wtr_5 = Temp_10748206)

wtr_2020 <- bind_rows(wtr_2020_1, wtr_2020_2, wtr_2020_3)

wtr_all <- bind_rows(wtr_2019, wtr_2020)

#Datetime sequence used for figures
datetime_seq <- data.frame(datetime = seq(min(wtr_all$datetime), max(wtr_all$datetime), "10 mins"))

#Oxygen sensor data
oxygen_2019 <- sensor_data_2019 |> 
  select(datetime, oxygen_2 = `DO_%_39`, oxygen_3 = `DO_%_65`)

oxygen_2020_1 <- sensor_data_2020_1 |> 
  select(datetime, oxygen_1 = DO_procent_500098, oxygen_2 = DO_procent_348923)

oxygen_2020_2 <- sensor_data_2020_2 |> 
  select(datetime, oxygen_1 = DO_procent_500098, oxygen_2 = DO_procent_348923, oxygen_3 = DO_procent_652844)

oxygen_2020_3 <- sensor_data_2020_3 |> 
  mutate(dosat = o2.at.sat.base(Temp_10675577),
         DO_procent_652844 = DO_mgL_652844/dosat*100) |>
  select(datetime, oxygen_1 = DO_procent_500098, oxygen_2 = DO_procent_348923, oxygen_3 = DO_procent_652844)

oxygen_2020 <- bind_rows(oxygen_2020_1, oxygen_2020_2, oxygen_2020_3)

oxygen_all <- bind_rows(oxygen_2019, oxygen_2020)

#DIC calculations from pH and alk (predicted from spec cond)
# dic_2019 <- sensor_data_2019 |> 
#   mutate(anc_predicted = predict(calcurve_2019_model, newdata = data.frame(spec_cond=`Sp.Cond._18`)),
#          anc_predicted = anc_predicted/1000) |> #mmol/l to mol/l
#   select(datetime, wtr_dic = Temp_21, ph=pH_27, spec_cond = `Sp.Cond._18`, anc_predicted) |> 
#   mutate(aquaenv = pmap(list(wtr_dic, ph, anc_predicted), ~aquaenv(S=0, t=..1, SumCO2 = NULL, pH = ..2, TA = ..3)),
#          dic = map_dbl(aquaenv, ~.$SumCO2))

#Perform DIC calculations and cache result
# saveRDS(dic_2019, "data/dic_2019.rds")

dic_2019 <- readRDS("data/dic_2019.rds")

#zmix calculations from vertical temperature profiles
zmix_2019 <- sensor_data_2019 |> 
  select(datetime, wtr_0.08 = Temp_8, wtr_0.21 = Temp_21, wtr_0.29 = Temp_29, wtr_0.37 = Temp_37, wtr_0.49 = Temp_49) |> 
  as.data.frame() |> 
  ts.thermo.depth() |> 
  as_tibble() |> 
  mutate(zmix = ifelse(is.na(thermo.depth), sensor_depth$total$`2019`, thermo.depth)) |> 
  select(-thermo.depth)

zmix_2020_1 <- sensor_data_2020_1 |> 
  select(datetime, wtr_0.125 = Temp_10748214, wtr_0.415 = Temp_10748223, wtr_0.545 = Temp_10675577) |> 
  as.data.frame() |> 
  ts.thermo.depth() |> 
  as_tibble() |> 
  mutate(zmix = ifelse(is.na(thermo.depth), sensor_depth$total$`2020_1`, thermo.depth)) |> 
  select(-thermo.depth)

zmix_2020_2 <- sensor_data_2020_2 |> 
  select(datetime, wtr_0.005 = Temp_10748214, wtr_0.25 = Temp_10748223, wtr_0.38 = Temp_10675577, wtr_0.53 = Temp_10748206) |> 
  as.data.frame() |> 
  ts.thermo.depth() |> 
  as_tibble() |> 
  mutate(zmix = ifelse(is.na(thermo.depth), sensor_depth$total$`2020_2`, thermo.depth)) |> 
  select(-thermo.depth)

zmix_2020_3 <- sensor_data_2020_3 |> 
  select(datetime, wtr_0.02 = Temp_10748214, wtr_0.355 = Temp_10675577, wtr_0.505 = Temp_10748206) |> 
  as.data.frame() |> 
  ts.thermo.depth() |> 
  as_tibble() |> 
  mutate(zmix = ifelse(is.na(thermo.depth), sensor_depth$total$`2020_3`, thermo.depth)) |> 
  select(-thermo.depth)

zmix_2020 <- bind_rows(zmix_2020_1 |> 
                         mutate(oxygen_depth = sensor_depth$oxygen$`2020_1`,
                                dic_depth = sensor_depth$ph$`2020_1`,
                                depth = sensor_depth$total$`2020_1`), 
                       zmix_2020_2 |> 
                         mutate(oxygen_depth = sensor_depth$oxygen$`2020_2`,
                                dic_depth = sensor_depth$ph$`2020_2`,
                                depth = sensor_depth$total$`2020_2`),
                       zmix_2020_3 |> 
                         mutate(oxygen_depth = sensor_depth$oxygen$`2020_3`,
                                dic_depth = sensor_depth$ph$`2020_3`,
                                depth = sensor_depth$total$`2020_3`))

#Plant survey rawdata
lake_poly <- st_read("data/lake_nymolle.sqlite")

lake_area <- as.numeric(st_area(lake_poly))

# #Read data from plant survey
# plants <- read_excel("data/nymolle_plants_2019.xlsx") |>
#   select(pkt = Punkt, long, lat, species = `Art latin`,
#          total_cover = `Total dækningsgrad %`, depth = `Dybde [m]`,
#          species_cover = `Dækningsgrad art%`) |>
#   filter(!is.na(long)) |>
#   st_as_sf(crs = 4326, coords = c("long", "lat")) |>
#   mutate(chara_cover = ifelse(grepl("Chara*", species), species_cover, 0))
# 
# #Save to file to edit coordinates in Google earth
# st_write(plants, "data/plants.kml")
# st_write(lake_poly, "data/lake_poly.kml")

plants_edit <- st_read("data/plants_edit.kml") |> 
  select(pkt, species, total_cover, depth, species_cover, chara_cover) |> 
  st_transform(25832) |> 
  st_zm()

#Use mean x-y per pkt due to small differences cause by editting
plant_points <- bind_cols(st_coordinates(plants_edit), st_drop_geometry(plants_edit)) |> 
  group_by(pkt) |> 
  summarise(X = mean(X), Y=mean(Y),
            depth = mean(depth),
            total_cover = mean(total_cover),
            chara_cover = sum(chara_cover)) |> 
  mutate(chara_cover = ifelse(chara_cover > 100, 100, chara_cover),
         chara_pres_abs = ifelse(chara_cover > 0, "Presence", "Absence")) |> 
  st_as_sf(coords=c("X", "Y"), crs=25832)
