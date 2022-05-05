source("rawdata.R")

#Figures

#Figure 1 - depth and plant cover maps
raster_template <- rast(vect(lake_poly), resolution = c(1, 1))

plants_voronoi <- voronoi(vect(plant_points))

total_cover_raster <- rasterize(plants_voronoi, raster_template, field="total_cover")
total_cover_raster_mask <- mask(total_cover_raster, vect(lake_poly))

chara_cover_raster <- rasterize(plants_voronoi, raster_template, field="chara_cover")
chara_cover_raster_mask <- mask(chara_cover_raster, vect(lake_poly))

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
  scale_fill_gradient(low = brewer.pal(5, "Greens")[1], high = brewer.pal(5, "Greens")[5], name = "Cover (%)")+
  scale_shape_manual(values = c(1, 19), name = "Charophytes")+
  ylab(NULL)+
  xlab(NULL)+
  theme(axis.ticks = element_blank(), axis.text = element_blank())

#add image
chara_img <- readJPEG("data/chara_image.jpg")
chara_img_grob <- rasterGrob(chara_img)

figure_1 <- depth_map + cover_map + chara_img_grob + plot_annotation(tag_levels = "A")+plot_layout(ncol=1, heights = c(1, 1, 0.95))

ggsave("figures/figure_1.png", figure_1, width = 129, height = 234, units = "mm")

#Figure 2 - open water water temperature, ph and oxygen profiles
profile <- read.delim2("data/profile.txt")

profile_long <- profile |>  
  mutate(date = dmy(date), 
         Date = strftime(date, format = "%b. %Y"),
         Date = factor(Date, levels = unique(strftime(date, format = "%b. %Y")))) |>  
  select(-oxygen_mg_l, -date) |> 
  gather(variable, value, -Date, -depth_m) |> 
  mutate(variable = case_when(variable == "ph" ~ "'pH'",
                              variable == "oxygen_perc" ~ "'Oxygen saturation (%)'",
                              variable == "temp" ~ "Water~temperature~'('*degree*C*')'")) |> 
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
  ggplot(aes(datetime_hour, value, col = Position))+
  geom_rect(data = rect_df, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=grey(0.8))+
  geom_line()+
  facet_grid(.~period, scales="free_x", space = "free_x")+
  scale_color_viridis_d(direction=-1)+
  scale_x_datetime(date_labels = "%d %b.", date_breaks = "10 days")+
  xlab("Date")+
  ylab(expression(Water~temperature~'('*degree*C*')'))+
  theme(strip.background = element_blank())

wtr_sub_2019 <- wtr_plot_data |> 
  filter(between(datetime_hour, xmin_2019, xmax_2019)) |> 
  ggplot(aes(datetime_hour, value, col = Position))+
  geom_line(show.legend = FALSE)+
  scale_color_viridis_d(direction=-1)+
  scale_x_datetime(date_labels = "%d %b.")+
  xlab("Date")+
  ylab(expression(Temp.~'('*degree*C*')'))

wtr_sub_2020 <- wtr_plot_data |> 
  filter(between(datetime_hour, xmin_2020, xmax_2020)) |> 
  ggplot(aes(datetime_hour, value, col = Position))+
  geom_line(show.legend = FALSE)+
  scale_color_viridis_d(direction=-1)+
  scale_x_datetime(date_labels = "%d %b.")+
  xlab("Date")+
  ylab(expression(Temp.~'('*degree*C*')'))

figure_3 <- wtr_all_plot/(wtr_sub_2019+wtr_sub_2020)+
  plot_annotation(tag_levels = "A")+
  plot_layout(guides = "collect", heights = c(1,0.5)) &
  theme(legend.position='bottom')

figure_3

ggsave("figures/figure_3.png", figure_3, width = 174, height = 120, units = "mm")

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
  ggplot(aes(datetime_hour, value, col = Position))+
  geom_rect(data = rect_df, inherit.aes = FALSE, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=grey(0.8))+
  geom_line()+
  facet_grid(.~period, scales="free_x", space = "free_x")+
  scale_color_manual(values=oxygen_pal)+
  scale_x_datetime(date_labels = "%d %b.", date_breaks = "10 days")+
  xlab("Date")+
  ylab(expression(Oxygen~saturation~'(%)'))+
  theme(strip.background = element_blank())

oxygen_sub_2019 <- oxygen_plot_data |> 
  filter(between(datetime_hour, xmin_2019, xmax_2019)) |> 
  ggplot(aes(datetime_hour, value, col = Position))+
  geom_line(show.legend = FALSE)+
  scale_color_manual(values=oxygen_pal)+
  scale_x_datetime(date_labels = "%d %b.")+
  xlab("Date")+
  ylab(expression(O[2]~'(%)'))

oxygen_sub_2020 <- oxygen_plot_data |> 
  filter(between(datetime_hour, xmin_2020, xmax_2020)) |> 
  ggplot(aes(datetime_hour, value, col = Position))+
  geom_line(show.legend = FALSE)+
  scale_color_manual(values=oxygen_pal)+
  scale_x_datetime(date_labels = "%d %b.")+
  xlab("Date")+
  ylab(expression(O[2]~'(%)'))

figure_4 <- oxygen_all_plot/(oxygen_sub_2019+oxygen_sub_2020)+
  plot_annotation(tag_levels = "A")+
  plot_layout(guides = "collect", heights = c(1,0.5)) &
  theme(legend.position='bottom')

figure_4

ggsave("figures/figure_4.png", figure_4, width = 174, height = 120, units = "mm")

#Figure 5 - ph and dic
ph_dic_col <- brewer.pal(n = 4, name = "Dark2")[c(3, 4)]

ph_dic_data <- dic_all |> 
  select(datetime, ph, dic) |> 
  right_join(datetime_seq) |> 
  mutate(period = year(datetime),
         datetime_hour = round_date(datetime, "hour")) |> 
  group_by(datetime_hour, period) |> 
  summarise(ph = mean(ph), dic=mean(dic) * 1000) |> 
  ungroup() |> 
  filter(between(datetime_hour, min(wtr_2019$datetime), max(wtr_2019$datetime)) | 
           between(datetime_hour, min(wtr_2020$datetime), max(wtr_2020$datetime))) 

a.diff <- max(ph_dic_data$ph, na.rm = TRUE) - min(ph_dic_data$ph, na.rm = TRUE)
b.diff <- max(ph_dic_data$dic, na.rm = TRUE) - min(ph_dic_data$dic, na.rm = TRUE)
a.min <- min(ph_dic_data$ph, na.rm = TRUE)
b.min <- min(ph_dic_data$dic, na.rm = TRUE)

figure_5 <- ph_dic_data |> 
  ggplot(aes(datetime_hour))+
  geom_line(aes(y = ph, col="pH"))+
  geom_line(aes(y = (dic - b.min) / b.diff * a.diff + a.min, col="DIC"))+
  scale_y_continuous(name = "pH", 
                     sec.axis = sec_axis(~((. -a.min) * b.diff / a.diff) + b.min, name=expression("DIC (mmol L"^{-1}*")")))+
  scale_color_manual(values=ph_dic_col)+
  facet_grid(.~period, scales="free_x", space = "free_x")+
  scale_x_datetime(date_labels = "%d %b.", date_breaks = "10 days")+
  xlab("Date")+
  theme(strip.background = element_blank(), 
        legend.title = element_blank(),
        legend.position = "bottom")

figure_5

ggsave("figures/figure_5.png", figure_5, width = 174, height = 84, units = "mm")

#Figure 6. Lake metabolism (R, GPP and NEP) based on O2 (solid line) and DIC (dashed line). Or just points and smoothing?
metab_colors <- c("GPP" = brewer.pal(n = 4, name = "Dark2")[1],
                  "NEP" = brewer.pal(n = 4, name = "Dark2")[2],
                  "R" = brewer.pal(n = 4, name = "Dark2")[3])

#Figure 7. A) R vs GPP (normalized to 20 degrees) with model II regression fit. 
#B) GPP vs mean lux with linear or saturating curve. 
#C) O2 vs DIC rates with 1:1 line, maybe model II regression fit. 
#Filled points O2 and empty points are DIC
#Solid lines are 02 and dashed are DIC

#Table 1. Lake characteristics (area, mean depth), and chemistry (secchi depth, total P, alkalinity), site plant biomass.
