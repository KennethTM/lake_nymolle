#Libraries and functions

library(tidyverse);library(lubridate);library(patchwork);library(viridisLite)
library(httr);library(jsonlite);library(readxl);library(LakeMetabolizer);library(rLakeAnalyzer)
library(RColorBrewer);library(AquaEnv);library(zoo);library(sf);library(terra)
library(fields);library(jpeg);library(grid);library(ggspatial)
library(lmodel2);library(ggpmisc)

#Figure sizing. For most journals the figures should be 39 mm, 84 mm, 129 mm, or 174 mm wide and not higher than 234 mm.
#ggplot theme
theme_pub <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(colour = "black"), 
        panel.border = element_rect(fill = NA, colour = "black"),
        strip.background = element_rect(fill = "white"))
theme_set(theme_pub)

#Sys.setenv(TZ="GMT");Sys.setlocale("LC_TIME", "English")

#Open site coordinates UTM32
open_site_x <- 698630.5	
open_site_y <- 6169134

#Littoral site coordinates UTM32
littoral_site_x <- 698751.9
littoral_site_y <- 6169169

#Dates and rects used for figures
xmin_2019 <- ymd_hm("2019-08-30 00:00")
xmax_2019 <- ymd_hm("2019-09-03 00:00")
xmin_2020 <- ymd_hm("2020-04-08 00:00")
xmax_2020 <- ymd_hm("2020-04-12 00:00")

rect_df <- data.frame(xmin = c(xmin_2019, xmin_2020),
                      xmax = c(xmax_2019, xmax_2020),
                      ymin = -Inf,
                      ymax = Inf,
                      period = c(2019, 2020),
                      letter = c("B", "C"),
                      x = c(ymd_hm("2019-09-01 00:00"), ymd_hm("2020-04-10 00:00")),
                      y = 175)
