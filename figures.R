source("rawdata.R")

#Figures

#Figure 1 - depth and plant cover maps
raster_template <- rast(vect(lake_poly), resolution = c(1, 1))

plants_voronoi <- voronoi(vect(plant_points))

total_cover_raster <- rasterize(plants_voronoi, raster_template, field="total_cover")
total_cover_raster_mask <- mask(total_cover_raster, vect(lake_poly))

chara_cover_raster <- rasterize(plants_voronoi, raster_template, field="chara_cover")
chara_cover_raster_mask <- mask(chara_cover_raster, vect(lake_poly))

height_raster <- rasterize(plants_voronoi, raster_template, field="height")
height_raster_mask <- mask(height_raster, vect(lake_poly))

#Bathymetric map
depth_zero <- lake_poly |> 
  st_cast("MULTILINESTRING") %>% 
  st_cast("LINESTRING") %>% 
  st_line_sample(density = 0.2) %>% 
  st_cast("POINT") %>% 
  st_as_sf() %>% 
  mutate(depth = 0) %>% 
  rename(geometry = x)

depth_all <- rbind(depth_zero, plant_points[, "depth"])

tps_mod <- Tps(st_coordinates(depth_all), depth_all$depth)
depth_interp <- interpolate(raster_template, tps_mod)
depth_interp_mask <- mask(depth_interp, vect(lake_poly))
depth_interp_mask[depth_interp_mask < 0] = 0

#Lake stats, e.g. depth, biomass, and chemistry
z_mean <- mean(depth_interp_mask[], na.rm =TRUE) #0.94 m
lake_area

littoral_site_chara_height <- c(34, 35, 33, 44, 32, 35, 38) #cm
mean(littoral_site_chara_height)
sd(littoral_site_chara_height)

littoral_site_chara_biomass <- c(1039,	1287, 1367) #g dw/m2
mean(littoral_site_chara_biomass)
sd(littoral_site_chara_biomass)

summary(chemistry)

#Table S1
chemistry |> 
  select(secchi_depth, alk_mmol_l, chla_ug_l, tp_mg_p_l, tn_mg_n_l) |> 
  summarise_all(list("mean" = mean,"sd" = sd))

#Proportion of lake below 1 meter and chara cover > 75
sum((total_cover_raster_mask > 75 & depth_interp_mask < 1)[], na.rm=TRUE)/lake_area*100

#Mean depth of lake where depth >1 meter
depth_minus_veg <- depth_interp_mask - height_raster_mask
depth_minus_veg[(total_cover_raster_mask > 75 & depth_interp_mask < 1)] <- NA
mean(depth_minus_veg[], na.rm = TRUE)
sum(!is.na(depth_minus_veg[]))

#Data frames for plotting
depth_df <- as.data.frame(depth_interp_mask, xy = TRUE)
total_cover_df <- as.data.frame(total_cover_raster_mask, xy = TRUE)
chara_cover_df <- as.data.frame(chara_cover_raster_mask, xy = TRUE)

depth_map <- ggplot()+
  geom_raster(data = depth_df, aes(x=x, y=y, fill=lyr.1))+
  geom_sf(data = lake_poly, fill=NA)+
  coord_sf(datum=25832)+
  annotate("point", x = littoral_site_x, y = littoral_site_y, col="red", size = 3)+
  #annotate("point", x = open_site_x, y = open_site_y, col="red", size = 2)+
  scale_fill_viridis_c(option="mako", direction = -1, name = "Depth (m)", 
                       guide = guide_colorbar(reverse = TRUE))+
  ylab(NULL)+
  xlab(NULL)+
  theme(axis.ticks = element_blank(), axis.text = element_blank())+
  annotation_scale(location = "br")+
  annotation_north_arrow(location = "tr", height = unit(1, "cm"), width = unit(1, "cm"))

cover_map <- ggplot()+
  geom_raster(data = total_cover_df, aes(x=x, y=y, fill=total_cover))+
  geom_sf(data = lake_poly, fill=NA)+
  geom_sf(data=plant_points, aes(shape = chara_pres_abs))+
  coord_sf(datum=25832)+
  scale_fill_gradient(low = brewer.pal(5, "Greens")[1], high = brewer.pal(5, "Greens")[5], name = "Macrophyte\ncover (%)")+
  scale_shape_manual(values = c(1, 19), name = "Charophytes")+
  ylab(NULL)+
  xlab(NULL)+
  theme(axis.ticks = element_blank(), axis.text = element_blank())

figure_1 <- depth_map + cover_map + plot_annotation(tag_levels = "A")+plot_layout(ncol=1, guides = "collect")

figure_1

ggsave("figures/figure_1.png", figure_1, width = 129, height = 160, units = "mm")
ggsave("figures/figure_1.pdf", figure_1, width = 129, height = 160, units = "mm")

#Figure 2 - open water water temperature, ph and oxygen profiles

#Stats
#Gradient
profile |> 
  filter(between(month(date), 5, 9)) |> 
  group_by(date) |> 
  summarise(temp_grad = (temp[which(depth_m == 0.25)] - temp[which(depth_m == 3)])/(3-0.25)) |> 
  summary()

#Surface stats
profile |> 
  filter(between(month(date), 5, 9) & depth_m == 0.25) |> 
  summary()

#CO2 stats
profile |> 
  filter(between(month(date), 5, 9)) |> 
  group_by(date) |> 
  summarise(ph = mean(ph), temp=mean(temp)) |> 
  na.omit() |> 
  left_join(chemistry[, c("date", "alk_mmol_l")]) |> 
  mutate(aquaenv = pmap(list(temp, ph, alk_mmol_l), ~aquaenv(S=0, t=..1, SumCO2 = NULL, pH = ..2, TA = ..3/1000)),
         co2 = map_dbl(aquaenv, ~.$CO2)*10^6) |> 
  summary()

profile_long <- profile |>  
  mutate(date = dmy(date), 
         Date = strftime(date, format = "%b %Y"),
         Date = factor(Date, levels = unique(strftime(date, format = "%b %Y")))) |>  
  select(-oxygen_mg_l, -date) |> 
  gather(variable, value, -Date, -depth_m) |> 
  mutate(variable = case_when(variable == "ph" ~ "'pH'",
                              variable == "oxygen_perc" ~ "'Oxygen saturation (%)'",
                              variable == "temp" ~ "Water~temperature~'('*degree*C*')'"),
         variable = factor(variable, levels=c("Water~temperature~'('*degree*C*')'", 
                                              "'Oxygen saturation (%)'",
                                              "'pH'"))) |> 
  na.omit()

figure_2 <- profile_long |> 
  ggplot(aes(value, depth_m, col = Date))+
  geom_point()+
  geom_path()+
  scale_color_brewer(palette = "Dark2")+
  scale_y_reverse(breaks=seq(0, 3, 0.5), limits = c(3, 0))+
  facet_grid(.~variable, scales = "free", labeller = label_parsed)+
  xlab(NULL)+
  ylab("Depth (m)")+
  theme(strip.background = element_blank())

figure_2

ggsave("figures/figure_2.png", figure_2, width = 174, height = 84, units = "mm")
ggsave("figures/figure_2.pdf", figure_2, width = 174, height = 84, units = "mm")

#Figure 3 - water temperature dynamics
wtr_plot_data <- left_join(datetime_seq, wtr_all) |> 
  gather(variable, value, -datetime) |> 
  mutate(Position = factor(parse_number(variable)),
         period = year(datetime),
         datetime_hour = round_date(datetime, "hour")) |> 
  group_by(datetime_hour, Position, period) |> 
  summarise(value = mean(value)) |> 
  ungroup() |> 
  filter(between(datetime_hour, min(wtr_2019$datetime), max(wtr_2019$datetime)) | 
           between(datetime_hour, min(wtr_2020$datetime), max(wtr_2020$datetime))) 

wtr_all_plot <- wtr_plot_data |> 
  mutate(`Mean depth (cm)` = case_when(Position == "1" ~ "7",
                                    Position == "2" ~ "21",
                                    Position == "3" ~ "32",
                                    Position == "4" ~ "42",
                                    Position == "5" ~ "51"),
         `Mean depth (cm)` = factor(`Mean depth (cm)`, levels = c("7", "21", "32", "42", "51"))) |> 
  ggplot(aes(datetime_hour, value, col = `Mean depth (cm)`))+
  geom_rect(data = rect_df, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "white", linetype=2, col="black")+
  geom_text(data = rect_df, inherit.aes = FALSE, aes(x=x, y=27, label=letter))+
  geom_line()+
  facet_grid(.~period, scales="free_x", space = "free_x")+
  scale_color_viridis_d(direction=-1)+
  scale_x_datetime(date_labels = "%d %b", date_breaks = "10 days")+
  xlab("Date")+
  ylab(expression(Water~temperature~'('*degree*C*')'))+
  theme(strip.background = element_blank())

wtr_sub_2019 <- wtr_plot_data |> 
  filter(between(datetime_hour, xmin_2019, xmax_2019)) |> 
  mutate(`Depth (cm)` = case_when(Position == "1" ~ "8",
                                       Position == "2" ~ "21",
                                       Position == "3" ~ "29",
                                       Position == "4" ~ "37",
                                       Position == "5" ~ "49"),
         `Depth (cm)` = factor(`Depth (cm)`, levels = c("8", "21", "29", "37", "49"))) |> 
  ggplot(aes(datetime_hour, value, col = `Depth (cm)`))+
  geom_line()+
  scale_color_viridis_d(direction=-1)+
  scale_x_datetime(date_labels = "%d %b")+
  xlab("Date")+
  ylab(expression(Water~temperature~'('*degree*C*')'))

wtr_sub_2020 <- wtr_plot_data |> 
  filter(between(datetime_hour, xmin_2020, xmax_2020)) |>
  mutate(`Depth (cm)` = case_when(Position == "1" ~ "5",
                                  Position == "3" ~ "25",
                                  Position == "4" ~ "38",
                                  Position == "5" ~ "51"),
         `Depth (cm)` = factor(`Depth (cm)`, levels = c("5", "25", "38", "51"))) |> 
  na.omit() |> 
  ggplot(aes(datetime_hour, value, col = `Depth (cm)`))+
  geom_line()+
  scale_color_manual(values = viridis(5, direction = -1)[c(1, 3, 4, 5)])+ 
  scale_x_datetime(date_labels = "%d %b")+
  xlab("Date")+
  ylab(expression(Water~temperature~'('*degree*C*')'))

figure_3 <- wtr_all_plot/(wtr_sub_2019+wtr_sub_2020)+
  plot_annotation(tag_levels = "A")+
  plot_layout(heights = c(1,1)) &
  theme(legend.position='bottom') & 
  guides(col=guide_legend(title.position = "top", title.hjust = 0.5))

figure_3

ggsave("figures/figure_3.png", figure_3, width = 174, height = 165, units = "mm")
ggsave("figures/figure_3.pdf", figure_3, width = 174, height = 165, units = "mm")

#Figure 4 - oxygen dynamics
oxygen_pal <- brewer.pal(n = 3, name = "Dark2")[c(2, 1, 3)]

oxygen_plot_data <- left_join(datetime_seq, oxygen_all) |> 
  gather(variable, value, -datetime) |> 
  mutate(Position = factor(parse_number(variable)),
         period = year(datetime),
         datetime_hour = round_date(datetime, "hour")) |> 
  group_by(datetime_hour, Position, period) |> 
  summarise(value = mean(value)) |> 
  ungroup() |> 
  filter(between(datetime_hour, min(wtr_2019$datetime), max(wtr_2019$datetime)) | 
           between(datetime_hour, min(wtr_2020$datetime), max(wtr_2020$datetime))) 

oxygen_all_plot <- oxygen_plot_data |> 
  mutate(`Mean depth (cm)` = case_when(Position == "1" ~ "9",
                                       Position == "2" ~ "30",
                                       Position == "3" ~ "51"),
         `Mean depth (cm)` = factor(`Mean depth (cm)`, levels = c("9", "30", "51"))) |> 
  ggplot(aes(datetime_hour, value, col = `Mean depth (cm)`))+
  geom_rect(data = rect_df, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "white", linetype=2, col="black")+
  geom_text(data = rect_df, inherit.aes = FALSE, aes(x=x, y=y, label=letter))+
  geom_line()+
  facet_grid(.~period, scales="free_x", space = "free_x")+
  scale_color_manual(values=oxygen_pal)+
  scale_x_datetime(date_labels = "%d %b", date_breaks = "10 days")+
  xlab("Date")+
  ylab(expression(Oxygen~saturation~'(%)'))+
  theme(strip.background = element_blank())

oxygen_sub_2019 <- oxygen_plot_data |> 
  filter(between(datetime_hour, xmin_2019, xmax_2019)) |> 
  mutate(`Depth (cm)` = case_when(Position == "2" ~ "39",
                                  Position == "3" ~ "65"),
         `Depth (cm)` = factor(`Depth (cm)`, levels = c("39", "65"))) |> 
  na.omit() |> 
  ggplot(aes(datetime_hour, value, col = `Depth (cm)`))+
  geom_line()+
  scale_color_manual(values=oxygen_pal[c(2, 3)])+
  scale_x_datetime(date_labels = "%d %b")+
  xlab("Date")+
  ylab(expression(Oxygen~saturation~'(%)'))

oxygen_sub_2020 <- oxygen_plot_data |> 
  filter(between(datetime_hour, xmin_2020, xmax_2020)) |> 
  mutate(`Depth (cm)` = case_when(Position == "1" ~ "6",
                                  Position == "2" ~ "23",
                                  Position == "3" ~ "40"),
         `Depth (cm)` = factor(`Depth (cm)`, levels = c("6", "23", "40"))) |> 
  ggplot(aes(datetime_hour, value, col = `Depth (cm)`))+
  geom_line()+
  scale_color_manual(values=oxygen_pal)+
  scale_x_datetime(date_labels = "%d %b")+
  xlab("Date")+
  ylab(expression(Oxygen~saturation~'(%)'))

figure_4 <- oxygen_all_plot/(oxygen_sub_2019+oxygen_sub_2020)+
  plot_annotation(tag_levels = "A")+
  plot_layout(heights = c(1,1)) &
  theme(legend.position='bottom') & 
  guides(col=guide_legend(title.position = "top", title.hjust = 0.5))

figure_4

ggsave("figures/figure_4.png", figure_4, width = 174, height = 165, units = "mm")
ggsave("figures/figure_4.pdf", figure_4, width = 174, height = 165, units = "mm")

#Figure 5 - ph and dic
ph_dic_col <- c("pH" = brewer.pal(n = 4, name = "Dark2")[3],
                "DIC" = brewer.pal(n = 4, name = "Dark2")[4])

figure_5 <- dic_2019 |> 
  mutate(dic_mmol = dic*1000) |> 
  ggplot(aes(datetime))+
  geom_line(aes(y = ph, col="pH"))+
  geom_line(aes(y = dic_mmol*5, col="DIC"))+
  scale_y_continuous(name = "pH", sec.axis = sec_axis(~. / 5, name=expression("DIC (mmol L"^{-1}*")")))+
  scale_color_manual(values=ph_dic_col)+
  scale_x_datetime(date_labels = "%d %b", date_breaks = "5 days")+
  xlab("Date")+
  theme(strip.background = element_blank(), 
        legend.title = element_blank(),
        legend.position = "bottom")

figure_5

ggsave("figures/figure_5.png", figure_5, width = 129, height = 75, units = "mm")
ggsave("figures/figure_5.pdf", figure_5, width = 129, height = 75, units = "mm")

#Figure 6 - diel dic, ph and calcification dynamics
diel_data <- dic_2019 %>% 
  mutate(hco3 = map_dbl(aquaenv, ~.$HCO3),
         co3 = map_dbl(aquaenv, ~.$CO3),
         co2 = map_dbl(aquaenv, ~.$CO2),
         co2_sat = map_dbl(aquaenv, ~.$CO2_sat),
         calcification = -1*(c(0, diff(anc_predicted))/2)*6,
         date = as_date(datetime),
         hours = as.numeric(datetime - floor_date(datetime, "day"))/60/60) #calc in units mol/l/hour

mean_sun <- sensor_data_2019 %>% 
  select(datetime, Light_8) %>% 
  mutate(dark = as.numeric(Light_8 == 0), date = as_date(datetime), dark_shift = c(0, diff(dark)) ) %>% 
  group_by(date) %>% 
  summarise(rise = datetime[which.min(dark_shift)], set = datetime[which.max(dark_shift)]) %>% 
  slice(2:(n()-1)) %>% 
  summarise(rise_mean = mean(rise), set_mean = mean(set)) %>% 
  mutate(rise_mean_hour = as.numeric(rise_mean - floor_date(rise_mean, "day")),
         set_mean_hour = as.numeric(set_mean - floor_date(set_mean, "day")))

#Calc stats
diel_dic_stats <- diel_data |> 
  mutate(date = as_date(datetime)) |> 
  filter(between(hours, mean_sun$rise_mean_hour, mean_sun$set_mean_hour)) |> 
  group_by(date) |> 
  summarise(mean_calc = mean(calcification)*1000*(mean_sun$set_mean_hour-mean_sun$rise_mean_hour), #calc in mmol/l/daytime period
            drop_dic = (max(dic)-min(dic))*1000) |>  #dic drop in mmol/l/daytime period
  mutate(calc_dic_prop = mean_calc/drop_dic*100)

summary(diel_dic_stats)

ph_diel <- diel_data %>% 
  ggplot(aes(hours, y = ph))+
  annotate("rect", xmin=-Inf, xmax = mean_sun$rise_mean_hour, ymin=-Inf, ymax=Inf, fill= grey(0.8))+
  annotate("rect", xmax=Inf, xmin = mean_sun$set_mean_hour, ymin=-Inf, ymax=Inf, fill= grey(0.8))+
  geom_line(aes(group=date), alpha=0.2, col="black")+
  geom_smooth(method="gam", formula = y~s(x, bs="cc"), col=brewer.pal(n = 4, name = "Dark2")[3], fill = brewer.pal(n = 4, name = "Dark2")[3])+
  scale_x_continuous(breaks = seq(0, 24, 6))+
  ylab("pH")+
  xlab("Time of day")

dic_diel_col <- c("dic_mmol" = brewer.pal(n = 4, name = "Dark2")[4],
                  "hco3_mmol" = brewer.pal(n = 4, name = "Dark2")[1],
                  "co3_mmol" = brewer.pal(n = 4, name = "Dark2")[2])

dic_diel_labels <- c("DIC",
                     expression(HCO[3]^{"-"}),
                     expression(CO[3]^{"2-"}))

dic_diel <- diel_data %>% 
  mutate(dic_mmol = dic*1000,
         hco3_mmol = hco3*1000,
         co3_mmol = co3*1000) %>% 
  gather(variable, value, dic_mmol, hco3_mmol, co3_mmol) |> 
  ggplot(aes(hours, value, col=variable))+
  annotate("rect", xmin=-Inf, xmax = mean_sun$rise_mean_hour, ymin=-Inf, ymax=Inf, fill= grey(0.8))+
  annotate("rect", xmax=Inf, xmin = mean_sun$set_mean_hour, ymin=-Inf, ymax=Inf, fill= grey(0.8))+
  geom_line(aes(group = paste0(date, variable)), alpha=0.2)+
  geom_smooth(se=FALSE, method="gam", formula = y~s(x, bs="cc"))+
  scale_x_continuous(breaks = seq(0, 24, 6))+
  scale_color_manual(values=dic_diel_col, labels = dic_diel_labels)+
  xlab("Time of day")+
  ylab(expression("DIC, HCO"[3]^{"-"}*", and CO"[3]^{"2-"}~"(mmol L"^{-1}*")"))+
  theme(legend.title = element_blank(),
        legend.position = "bottom")

diel_co2_sat <- diel_data %>% 
  group_by(hours) %>% 
  summarise(co2_sat_mean = mean(co2_sat)*10^6)

co2_diel <- diel_data %>% 
  mutate(co2_umol = co2*10^6) %>% 
  ggplot(aes(hours, y = co2_umol))+
  annotate("rect", xmin=-Inf, xmax = mean_sun$rise_mean_hour, ymin=-Inf, ymax=Inf, fill= grey(0.8))+
  annotate("rect", xmax=Inf, xmin = mean_sun$set_mean_hour, ymin=-Inf, ymax=Inf, fill= grey(0.8))+
  geom_line(aes(group=date), alpha=0.2, col="black")+
  geom_line(data = diel_co2_sat, aes(y=co2_sat_mean), show.legend = FALSE, linetype=2)+
  geom_smooth(method="gam", formula = y~s(x, bs="cc"), col=brewer.pal(n = 4, name = "Dark2")[3], fill = brewer.pal(n = 4, name = "Dark2")[3])+
  scale_x_continuous(breaks = seq(0, 24, 6))+
  ylab(expression(CO[2]~concentration~"("*mu*mol~L^{-1}*")"))+
  xlab("Time of day")+
  coord_cartesian(ylim=c(0, 40))

calc_diel <- diel_data %>% 
  mutate(calcification_umol = calcification*10^6) |> 
  ggplot(aes(hours, y = calcification_umol))+
  annotate("rect", xmin=-Inf, xmax = mean_sun$rise_mean_hour, ymin=-Inf, ymax=Inf, fill= grey(0.8))+
  annotate("rect", xmax=Inf, xmin = mean_sun$set_mean_hour, ymin=-Inf, ymax=Inf, fill= grey(0.8))+
  geom_line(aes(group=date), alpha=0.05, col="black")+
  geom_hline(yintercept = 0, linetype = 3)+
  geom_smooth(method="gam", formula = y~s(x, bs="cc"), col=brewer.pal(n = 4, name = "Dark2")[3], fill = brewer.pal(n = 4, name = "Dark2")[3])+
  ylab(expression("Calcification ("*mu*mol~L^{-1}~h^{-1}*")"))+
  xlab("Time of day")+
  scale_x_continuous(breaks = seq(0, 24, 6))+
  coord_cartesian(ylim=c(-15, 15))

figure_6 <- ph_diel + dic_diel + co2_diel + calc_diel + plot_layout(ncol=2)+plot_annotation(tag_levels = "A")

figure_6

ggsave("figures/figure_6.png", figure_6, width = 174, height = 150, units = "mm")
ggsave("figures/figure_6.pdf", figure_6, width = 174, height = 150, units = "mm")

#Figure 7. Lake metabolism (R, GPP and NEP) based on O2 (solid line) and DIC (dashed line) with smoother to shown trend.
daily_metab <- readRDS("data/daily_metab.rds")

daily_depth <- bind_rows(daily_metab$`2019`[, c("date", "depth")], 
                         daily_metab$`2020`[, c("date", "depth")])

oxygen_daily <- bind_rows(daily_metab$`2019`$oxygen_daily,
                          daily_metab$`2020`$oxygen_daily) |> 
  filter(Rmax > 1e-3) #remove 5 days, 2 poor fits, 3 unrealistic R values

dic_daily <- bind_rows(daily_metab$`2019`$dic_daily)

all_daily <- bind_rows(add_column(oxygen_daily, method = "Oxygen"),
                       add_column(dic_daily, method = "DIC")) |> 
  mutate(period = year(datetime_min),
         date = as_date(datetime_min))

metab_colors <- c("GPP" = brewer.pal(n = 4, name = "Dark2")[1],
                  "NEP" = brewer.pal(n = 4, name = "Dark2")[3],
                  "R" = brewer.pal(n = 4, name = "Dark2")[2])

figure_7_data <- all_daily |> 
  gather(variable, value, GPP, R, NEP) |> 
  left_join(daily_depth) |> 
  mutate(value_m2 = value*depth,
         value_m2 = ifelse(variable == "R", value_m2*-1, value_m2))

figure_7 <- figure_7_data |> 
  ggplot(aes(date, value_m2, col = variable, linetype=method, shape = method))+
  geom_hline(yintercept = 0, linetype=3)+
  geom_point()+
  facet_grid(.~period, scales="free_x", space = "free_x")+
  scale_color_manual(values = metab_colors, name="Component")+
  scale_linetype_manual(values = c(2, 1))+
  scale_shape_manual(values = c(1, 16), name="Method")+
  scale_x_date(date_labels = "%d %b", date_breaks = "10 days")+
  ylab(expression("Metabolism (mmol m"^{-2}~d^{-1}*")"))+
  xlab("Date")

figure_7

ggsave("figures/figure_7.png", figure_7, width = 174, height = 75, units = "mm")
ggsave("figures/figure_7.pdf", figure_7, width = 174, height = 75, units = "mm")

#Table S4
#Metabolism summary table
figure_7_data |> 
  select(period, method, variable, value_m2) |> 
  group_by(period, method, variable) |>
  mutate(value_m2 = ifelse(variable == "R", value_m2*-1, value_m2)) |> 
  summarise(mean = mean(value_m2), sd=sd(value_m2), min= min(value_m2), max=max(value_m2), n=n()) |> 
  mutate(label = paste0(round(mean, digits = 0), " (±", round(sd, digits = 0), ", ", 
                        round(min, digits = 0), "–", round(max, digits = 0), ")")) |> 
  select(period, method, variable, label) |> 
  spread(variable, label) |> 
  write_csv("figures/table_s4.csv")

#Figure 8. 
#A) O2 vs DIC rates with 1:1 line, maybe model II regression fit. 
#Filled points O2 and empty points are DIC
#Solid lines are 02 and dashed are DIC
figure_8_a_data <- figure_7_data |> 
  select(date, variable, value_m2, method) |> 
  spread(method, value_m2) |> 
  na.omit()

lm2_o2_dic <- lmodel2(Oxygen~DIC, data = figure_8_a_data)
lm2_o2_dic

#photosynthetic quotient
pq <- figure_8_a_data[figure_8_a_data$variable == "GPP", ]$Oxygen/figure_8_a_data[figure_8_a_data$variable == "GPP", ]$DIC
pq_mean <- mean(pq)
pq_se <- sd(pq)/sqrt(length(pq))
pq_low <- pq_mean - pq_se*1.96
pq_high <- pq_mean + pq_se*1.96

#respiratory quotient
rq <- figure_8_a_data[figure_8_a_data$variable == "R", ]$DIC/figure_8_a_data[figure_8_a_data$variable == "R", ]$Oxygen
rq_mean <- mean(rq)
rq_se <- sd(rq)/sqrt(length(rq))
rq_low <- rq_mean - rq_se*1.96
rq_high <- rq_mean + rq_se*1.96

#t-tests for difference between dic and oxygen
t.test(figure_8_a_data[figure_8_a_data$variable == "GPP", ]$DIC, figure_8_a_data[figure_8_a_data$variable == "GPP", ]$Oxygen, paired = TRUE)
t.test(figure_8_a_data[figure_8_a_data$variable == "R", ]$DIC, figure_8_a_data[figure_8_a_data$variable == "R", ]$Oxygen, paired = TRUE)
t.test(figure_8_a_data[figure_8_a_data$variable == "NEP", ]$DIC, figure_8_a_data[figure_8_a_data$variable == "NEP", ]$Oxygen, paired = TRUE)

figure_8_a <- figure_8_a_data |> 
  ggplot(aes(DIC, Oxygen, col=variable))+
  geom_abline(intercept = 0, slope=1, linetype=3)+
  geom_point()+
  scale_color_manual(values = metab_colors, name="Component")+
  ylab(expression(Metabolism[DIC]~"(mmol m"^{-2}~d^{-1}*")"))+
  xlab(expression(Metabolism[Oxygen]~"(mmol m"^{-2}~d^{-1}*")"))+
  ylim(-700, 700)+
  xlim(-700, 700)

#B) R vs GPP (normalized to 20 degrees) with model II regression fit. 
figure_8_b_data <- all_daily |> 
  left_join(daily_depth) |> 
  mutate(GPP_m2 = GPP*depth,
         R_m2 = R*depth,
         GPP_m2_20 = GPP_m2/(1.073^(wtr_mean - 20)),
         R_m2_20 = R_m2/(1.073^(wtr_mean - 20)))

metab_dic_20 <- figure_8_b_data |> 
  filter(method=="DIC")
lm2_dic <- lmodel2(R_m2_20~GPP_m2_20, data = metab_dic_20)
lm2_pred_dic <- data.frame(GPP_m2_20 = seq(min(metab_dic_20$GPP_m2_20), max(metab_dic_20$GPP_m2_20)))
lm2_dic_df <- lm2_pred_dic |> 
  bind_cols(predict(lm2_dic, method="SMA", interval = "confidence", newdata=lm2_pred_dic))

metab_oxygen_20 <- figure_8_b_data |> 
  filter(method=="Oxygen")
lm2_oxygen <- lmodel2(R_m2_20~GPP_m2_20, data = metab_oxygen_20)
lm2_pred_oxygen <- data.frame(GPP_m2_20 = seq(min(metab_oxygen_20$GPP_m2_20), max(metab_oxygen_20$GPP_m2_20)))
lm2_oxygen_df <- lm2_pred_oxygen |> 
  bind_cols(predict(lm2_oxygen, method="SMA", interval = "confidence", newdata=lm2_pred_oxygen))

figure_8_b <- figure_8_b_data |> 
  ggplot()+
  geom_abline(intercept = 0, slope=1, linetype=3)+
  #geom_ribbon(data = lm2_dic_df, aes(GPP_m2_20, fit, ymin=upr, ymax=lwr), fill=grey(0.5, alpha=0.3))+
  geom_line(data = lm2_dic_df, aes(GPP_m2_20, fit), size=1, col="grey")+
  #geom_ribbon(data = lm2_oxygen_df, aes(GPP_m2_20, fit, ymin=upr, ymax=lwr), fill=grey(0.5, alpha=0.3))+
  geom_line(data = lm2_oxygen_df, aes(GPP_m2_20, fit), size=1)+
  geom_point(aes(GPP_m2_20, R_m2_20, col=method))+
  scale_color_manual(values = c("Oxygen" = "black", "DIC" = "grey"), name="Method")+
  ylab(expression(R[20]~"(mmol m"^{-2}~d^{-1}*")"))+
  xlab(expression(GPP[20]~"(mmol m"^{-2}~d^{-1}*")"))+
  ylim(0, 600)+
  xlim(0, 600)

figure_8 <- figure_8_a + figure_8_b + plot_layout(ncol=1)+plot_annotation(tag_levels = "A")

figure_8

ggsave("figures/figure_8.png", figure_8, width = 129, height = 180, units = "mm")
ggsave("figures/figure_8.pdf", figure_8, width = 129, height = 180, units = "mm")

#Table S3
#Species list
plants_edit |> 
  st_drop_geometry() |> 
  group_by(species) |> 
  summarise(count = n(),
            max_depth = max(depth),
            avg_cover = mean(species_cover)) |> 
  na.omit() |> #filter(!str_detect(species, "Chara")) |> summary()
  write_csv("figures/table_s3.csv")
