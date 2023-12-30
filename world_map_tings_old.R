# Load required libaries
library(sf); library(ggplot2); library(dplyr); library(rnaturalearth)
library(ggspatial); library(rgdal)

# Get high-quality natural earth data
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(region_un != "Antarctica") %>%
  filter(region_un != "Seven seas (open ocean)") %>%
  filter(continent != "Seven seas (open ocean)")

data <- data.frame(country = world$geounit, value = runif(n = nrow(world))) ## add in BPSV impact here
world_robinson <- st_transform(world, crs = "ESRI:54030")
merged_data <- merge(world_robinson, data, by.x = 'geounit', by.y = 'country')

# Plotting the output 
ggplot() +
  geom_sf(data = merged_data, aes(fill = continent), color = NA, border = NA) +
  coord_sf(crs = st_crs("ESRI:54030")) +
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none") +
  theme_void() +
  scale_fill_brewer(type = "qual", palette = 2, direction = -1) +
  guides(fill = "none")

