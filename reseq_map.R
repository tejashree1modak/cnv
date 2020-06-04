#Script to plot resequencing samples on a map

#install.packages(c("maps", "mapdata"))
library(tidyverse)
library(mapdata)
library(maps)
library(ggrepel)
states <- map_data("state")
head(states)

sampling_loc <- read.csv("/Users/tejashree/Documents/Projects/cnv/genome/Lat_Long_all_pop.csv",
                         header = TRUE, stringsAsFactors = FALSE)

sampling <- states %>%
  filter(region %in% c("maine","new hampshire", "massachusetts","connecticut", "rhode island", "vermont",
                       "new york","delaware", "virginia","maryland", "new jersey", "pennsylvania",
                       "florida","georgia","south carolina", "north carolina",
                       "texas","alabama","mississippi","louisiana"))
ggplot(data = sampling) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") + 
  geom_point(data = sampling_loc, aes(x = Long, y = Lat), color = "black", size = 1) +
  geom_text_repel(data = sampling_loc, aes(x = Long, y = Lat, label = Pop.ID.), 
                  hjust = 0, nudge_x = 1) +
  labs(x="Longitude", y="Latitude") +
  coord_quickmap() + theme_minimal()
ggsave("reseq_map.png")
