source("rawdata.R")

#Figures

#Figure 2 - water temperature and oxygen dynamics
datetime_seq <- data.frame(datetime = seq(min(wtr_all$datetime), max(wtr_all$datetime), "10 mins"))

oxygen_pal <- brewer.pal(n = 7, name = "RdPu")[c(3, 5, 7)]
wtr_pal <- brewer.pal(n = 7, name = "YlOrBr")[3:7]

rect_df <- data.frame(xmin = c(ymd_hm("2019-08-29 00:00"), ymd_hm("2020-04-08 00:00")),
                      xmax = c(ymd_hm("2019-09-02 00:00"), ymd_hm("2020-04-12 00:00")),
                      ymin = -Inf,
                      ymax = Inf,
                      period = c(2019, 2020))

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
  scale_color_manual(values=wtr_pal)+
  scale_x_datetime(date_labels = "%d %b.", date_breaks = "10 days")+
  xlab("Date")+
  ylab(expression(Water~temperature~'('*degree*C*')'))+
  theme(strip.background = element_blank())

wtr_sub_2019 <- wtr_plot_data |> 
  filter(between(datetime_hour, ymd_hm("2019-08-29 00:00"), ymd_hm("2019-09-02 00:00"))) |> 
  ggplot(aes(datetime_hour, value, col = Position))+
  geom_line(show.legend = FALSE)+
  scale_color_manual(values=wtr_pal)+
  scale_x_datetime(date_labels = "%d %b.")+
  xlab("Date")+
  ylab(expression(Water~temp.~'('*degree*C*')'))

wtr_sub_2020 <- wtr_plot_data |> 
  filter(between(datetime_hour, ymd_hm("2020-04-08 00:00"), ymd_hm("2020-04-12 00:00"))) |> 
  ggplot(aes(datetime_hour, value, col = Position))+
  geom_line(show.legend = FALSE)+
  scale_color_manual(values=wtr_pal)+
  scale_x_datetime(date_labels = "%d %b.")+
  xlab("Date")+
  ylab(expression(Water~temp.~'('*degree*C*')'))


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
  filter(between(datetime_hour, ymd_hm("2019-08-29 00:00"), ymd_hm("2019-09-02 00:00"))) |> 
  ggplot(aes(datetime_hour, value, col = Position))+
  geom_line(show.legend = FALSE)+
  scale_color_manual(values=wtr_pal)+
  scale_x_datetime(date_labels = "%d %b.")+
  xlab("Date")+
  ylab(expression(Oxygen~sat.~'(%)'))
  
oxygen_sub_2020 <- oxygen_plot_data |> 
  filter(between(datetime_hour, ymd_hm("2020-04-08 00:00"), ymd_hm("2020-04-12 00:00"))) |> 
  ggplot(aes(datetime_hour, value, col = Position))+
  geom_line(show.legend = FALSE)+
  scale_color_manual(values=wtr_pal)+
  scale_x_datetime(date_labels = "%d %b.")+
  xlab("Date")+
  ylab(expression(Oxygen~sat.~'(%)'))


figure_2_a <- wtr_all_plot/(wtr_sub_2019+wtr_sub_2020)+oxygen_all_plot/(oxygen_sub_2019+oxygen_sub_2020)+plot_annotation(tag_levels = "A")+plot_layout(guides = "collect", heights = c(1,0.7, 1, 0.7))

figure_2

ggsave("figures/figure_2.png", figure_2, width = 174, height = 234, units = "mm")


#Figure 4 - open water water temperature, ph and oxygen profiles
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

figure_4 <- profile_long |> 
  ggplot(aes(value, depth_m, col = Date))+
  geom_point()+
  geom_path()+
  scale_color_brewer(palette = "Dark2")+
  scale_y_reverse(breaks=seq(0, 3, 0.5), limits = c(3, 0))+
  facet_grid(.~variable, scales = "free", labeller = label_parsed)+
  xlab(NULL)+
  ylab("Depth (m)")+
  theme(strip.background = element_blank())

figure_4

ggsave("figures/figure_4.png", figure_4, width = 174, height = 84, units = "mm")
