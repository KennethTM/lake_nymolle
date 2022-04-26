source("rawdata.R")

#Figures

#Figure 1 - depth and plant cover maps

#https://rspatial.org/terra/analysis/4-interpolation.html

lake_poly <- st_read("data/lake_nymolle.sqlite")

plants_raw <- read_excel("data/nymolle_plants_2019.xlsx")

plants <- plants_raw |> 
  select(pkt = Punkt, long, lat, species = `Art latin`, 
         total_cover = `Total dækningsgrad %`, depth = `Dybde [m]`,
         species_cover = `Dækningsgrad art%`) |> 
  filter(!is.na(long)) |> 
  st_as_sf(crs = 4326, coords = c("long", "lat")) |> 
  st_transform(25832) |> 
  mutate(chara_cover = ifelse(grepl("Chara*", species), species_cover, 0))

depth <- plants |> 
  select(depth) |> 
  distinct()

cover <- plants |> 
  select(total_cover) |> 
  distinct()

# cover_df <- cbind(st_coordinates(cover), cover = cover$total_cover) |> 
#   as.data.frame() |> 
#   rename(x = X, y = Y)

chara_cover <- plants |> 
  group_by(pkt) |> 
  summarise(chara_cover = max(chara_cover))

# chara_cover_df <- cbind(st_coordinates(chara_cover), pres = chara_cover$chara_cover) |> 
#   as.data.frame() |> 
#   rename(x = X, y = Y)

raster_template <- rast(vect(lake_poly), resolution = c(1, 1))

cover_voronoi <- voronoi(vect(cover))
cover_raster <- rasterize(cover_voronoi, raster_template, field="total_cover")
cover_raster_mask <- mask(cover_raster, vect(lake_poly))

chara_cover_voronoi <- voronoi(vect(chara_cover))
chara_cover_raster <- rasterize(chara_cover_voronoi, raster_template, field="chara_cover")
chara_cover_raster_mask <- mask(chara_cover_raster, vect(lake_poly))

cover_raster_mask_df <- as.data.frame(cover_raster_mask, xy = TRUE)
chara_cover_raster_mask_df <- as.data.frame(chara_cover_raster_mask > 20, xy = TRUE)

ggplot()+
  geom_raster(data = cover_raster_mask_df, aes(x=x, y=y, fill=total_cover))+
  #geom_raster(data = chara_cover_raster_mask_df, aes(x=x, y=y, alpha=chara_cover))+
  scale_alpha_manual(values = c(0, 0.5))+
  geom_sf(data = lake_poly, fill=NA)+
  geom_sf(data = depth)

#MAKE POINT SHP CHARA PRES/ABS

#Dybdekort
#Plant cover + Chara pres/abs
#Image

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

