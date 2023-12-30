# Load required libaries
library(sf); library(ggplot2); library(dplyr); library(rnaturalearth)
library(ggspatial); library(rgdal)
source("functions/extract_process_generate_fits.R")
source("functions/run_sars_x.R")
source("functions/helper_functions.R")
source("main.R")

## Generating BPSV impact
iso_list <- substr(list.files(path = "outputs/Figure3_SC2_Counterfactual_Impact/raw/"), 5, 7)
deaths_df <- data.frame(iso = rep(NA_character_, length(iso_list)), 
                        empirical_deaths = rep(NA_real_, length(iso_list)), 
                        deaths_BPSV = rep(NA_real_, length(iso_list)))

for (i in 1:lenght(iso_list)) {
  temp <- readRDS(paste0("outputs/Figure3_SC2_Counterfactual_Impact/raw/raw_", iso_list[i], "_fit.rds"))
  temp_overall <- evaluate_country_impact(original_fit = temp, country_iso = iso_list[i])
  saveRDS(object = temp_overall, file = paste0("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual/BPSV_counterfactual_", iso_list[i], "_fit.rds"))
  deaths <- temp_overall %>%
    group_by(scenario) %>%
    summarise(med = sum(med)) %>%
    pivot_wider(names_from = scenario, values_from = med)
  deaths_df$iso[i] <- iso_list[i]
  deaths_df$empirical_deaths[i] <- deaths$no_bpsv_deaths
  deaths_df$deaths_BPSV[i] <- deaths$bpsv_deaths
  print(i)
}


# Get high-quality natural earth data
world <- ne_countries(scale = "medium", returnclass = "sf")
data <- data.frame(country = world$geounit, value = runif(n = nrow(world))) ## add in BPSV impact here
world_robinson <- st_transform(world, crs = "ESRI:54030")
merged_data <- merge(world_robinson, data, by.x = 'geounit', by.y = 'country')

# Getting bounding box and graticules
load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
NE_box_rob <- spTransform(NE_box, CRSobj = PROJ)
NE_graticules_rob <- spTransform(NE_graticules, CRSobj = PROJ)

# Plotting the output 
ggplot() +
  geom_sf(data = merged_data, aes(fill = value), color = "black") +
  geom_path(data=NE_graticules_rob, aes(long, lat, group=group), linetype="dashed", color="grey50", size = 0.05) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.1) + 
  coord_sf(crs = st_crs("ESRI:54030")) +
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white")) +
  theme_void() 

