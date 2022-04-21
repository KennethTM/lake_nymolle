#Libraries and functions

library(tidyverse);library(lubridate);library(patchwork);library(viridisLite)
library(httr);library(jsonlite);library(readxl);library(LakeMetabolizer)
library(RColorBrewer)

#Figure sizing. For most journals the figures should be 39 mm, 84 mm, 129 mm, or 174 mm wide and not higher than 234 mm.
#ggplot theme
theme_pub <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(colour = "black"), 
        panel.border = element_rect(fill = NA, colour = "black"),
        strip.background = element_rect(fill = "white"))
theme_set(theme_pub)

Sys.setenv(TZ="GMT");Sys.setlocale("LC_TIME", "English")

#Open site coordinates UTM32
open_site_x <- 698630.5	
open_site_y <- 6169134

#Littoral site coordinates UTM32
littoral_site_x <- 698751.9
littoral_site_y <- 6169169
