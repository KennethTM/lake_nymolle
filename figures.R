source("libs_and_funcs.R")

#Figures

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
