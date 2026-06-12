#Plotting coordinates to map

#Load necessary libraries:
library(readxl)
library(readr)
library(maps)
library(ggplot2)
library(ggrepel)
library(tidyverse)

#Load acridum dataset
acridum_data <- read_excel("d01a_acridum_map_pg_metadata.xlsx")
sites1 <- as.data.frame(acridum_data)

#Start plotting
worldmap = map_data('world')
set.seed(24)

#add noise to the data to avoid overlapping points and add size variance for date_bp
sites2 <- sites1
sites2$Lat <- as.numeric(sites2$Lat)
sites2$Lat = sites2$Lat + rnorm(length(sites2$Lat),0,0.7)
sites2$Long = sites2$Long + rnorm(length(sites2$Long),0,0.7)

#ggplot map
map1 <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'lightgrey', color = 'lightgrey') + 
  coord_fixed(1.3) + 
  theme_void()+
  geom_point(data = sites2, 
             aes(x = as.numeric(Long), 
                 y = as.numeric(Lat), fill = Study), size = 5, alpha = 0.8, pch = 21, colour = "white") + 
  scale_fill_identity() +
  scale_size_area(max_size = 10) + 
  theme(legend.position = 'bottom', legend.direction = 'horizontal') + 
  theme(title = element_text(size = 12))

map1


