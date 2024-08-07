# Load required libraries
library(tidyverse); library(rnaturalearth); library(sf); #library(rgdal); # library(squire.page); 

# Sourcing required functions
source("functions/extract_process_generate_fits.R")
source("functions/run_sars_x.R")
source("functions/helper_functions.R")

# Getting vaccination rates by country and summarising average by income strata
vacc_rates <- read.csv("data/country_vaccination_rates.csv")
vacc_rates$iso3c <- countrycode::countrycode(vacc_rates$country, "country.name", "iso3c")
wb_income_strata <- wbstats::wb_countries() %>%
  select(iso3c, iso2c, country, income_level_iso3c) %>%
  filter(!is.na(income_level_iso3c)) %>%
  left_join(vacc_rates, by = c("iso3c", "country")) %>%
  mutate(daily_percent_vaccinated = 2 * Percent.population)

vaccination_rate_summarised <- wb_income_strata %>%
  group_by(income_level_iso3c) %>%
  summarise(mean_vaccination_rate = mean(daily_percent_vaccinated, na.rm = TRUE) / 100,
            median_vaccination_rate = median(daily_percent_vaccinated, na.rm = TRUE) / 100,
            weekly_mean_vaccination_rate = 7 * mean_vaccination_rate,
            weekly_median_vaccination_rate = 7 * median_vaccination_rate) %>%
  filter(income_level_iso3c != "INX")

# Extracting unprocessed squire.page fits (required for Linux desktop to run as it struggles with download
# call in grab_fit)
fresh_run_fits <- FALSE
if (fresh_run_fits) {
  country_ISOs <- unique(squire::population$iso3c)
  for (i in 199:length(country_ISOs)) {
    # Getting squire.page country fit draws 
    excess <- grab_fit(country_ISOs[i], TRUE, TRUE)
    
    # Saving output
    saveRDS(object = excess, file = paste0("outputs/Figure3_SC2_Counterfactual_Impact/unprocessed_outputs/" , country_ISOs[i], "_fit.rds"))
    print(paste0("i = ", i, ", ISO = ", country_ISOs[i]))
  }
}

# Loading country ISOs and extracting squire.page fits
fresh_run_fits <- FALSE
if (fresh_run_fits) {
  country_ISOs <- unique(squire::population$iso3c)
  for (i in 1:length(country_ISOs)) {
    # Load squire.page again
    library(squire.page)
    
    # Selecting specific country ISO
    temp_ISO <- country_ISOs[i]
    
    # Getting squire.page country fit draws 
    temp <- get_country_draws(country_iso = temp_ISO)
    
    # Saving output
    saveRDS(object = temp, file = paste0("outputs/Figure3_SC2_Counterfactual_Impact/raw/raw_" , temp_ISO, "_fit.rds"))
    print(paste0("i = ", i, ", ISO = ", temp_ISO))
  }
}

## Generating BPSV impact estimates for SC2 whilst varying coverage, BPSV deployment date and vaccination rates
iso_list <- substr(list.files(path = "outputs/Figure3_SC2_Counterfactual_Impact/raw/"), 5, 7)

## BPSV Coverage Scenarios
low_coverage <- 0.4
mid_coverage <- 0.6
high_coverage <- 0.8
coverage_df <- data.frame(income_level_iso3c = c("LIC", "LMC", "UMC", "HIC"),
                          low_coverage = rep(low_coverage, 4),
                          mid_coverage = rep(mid_coverage, 4),
                          high_coverage = rep(high_coverage, 4),
                          variable_coverage = c(0.2, 0.4, 0.6, 0.8))

## BPSV Start Dates (pinned to when cumulative reported COVID-19 deaths reached certain threshold)
bpsv_start_date_1 <- as.Date("2020-01-09") # date of first reported deaths worldwide (https://www.bbc.co.uk/news/world-54337098)
bpsv_start_date_10 <- as.Date("2020-01-26") # date of 10 reported deaths worldwide (OWID)
bpsv_start_date_100 <- as.Date("2020-02-02") # date of 100 reported deaths worldwide (OWID)
bpsv_start_date_1000 <- as.Date("2020-02-10") # date of 1000 reported deaths worldwide (OWID)
bpsv_start_dates <- c(bpsv_start_date_1, bpsv_start_date_10, bpsv_start_date_100, bpsv_start_date_1000)
bpsv_start_dates_df <- data.frame(start_trigger = c("1Deaths", "10Deaths", "100Deaths", "1000Deaths"),
                                  bpsv_start_date = bpsv_start_dates)

## Initialising Dataframe to Store the Outputs
deaths_df <- data.frame(start_trigger = rep(NA_character_, 1),
                        coverage_scenario = rep(NA_character_, 1),
                        coverage = rep(NA_real_, 1),
                        deaths_BPSV = rep(NA_real_, 1),
                        empirical_deaths = rep(NA_real_, 1), 
                        iso = rep(NA_character_, 1))
                        
## Running All the Different Scenarios and Looping Over Country
new_run <- TRUE
for (i in 1:length(iso_list)) {
  if (new_run) {
    
    ## Manual catching of ISO3c edge cases not caught above
    if (iso_list[i] == "GUF") {
      iso_list[i] <- "GUY"
      income_strata <- wb_income_strata %>%
        filter(iso3c == iso_list[i]) 
    } else if (iso_list[i] == "TWN") {
      income_strata <- data.frame(iso3c = "TWN", income_level_iso3c = "HIC")
    } else if (iso_list[i] == "VEN") {
      income_strata <- data.frame(iso3c = "VEN", income_level_iso3c = "LIC")
    } else {
      income_strata <- wb_income_strata %>%
        filter(iso3c == iso_list[i]) 
    }
    income_group <- income_strata$income_level_iso3c
    
    ## Defining vaccination rate (based on which income group the country belongs to)
    vaccination_rate <- vaccination_rate_summarised$weekly_mean_vaccination_rate[vaccination_rate_summarised$income_level_iso3c == income_group]
    
    ## Defining coverage rate (for use in the scenario where coverage is variable by income strata)
    low_coverage <- coverage_df$low_coverage[coverage_df$income_level_iso3c == income_group]
    mid_coverage <- coverage_df$mid_coverage[coverage_df$income_level_iso3c == income_group]
    high_coverage <- coverage_df$high_coverage[coverage_df$income_level_iso3c == income_group]
    variable_coverage <- coverage_df$variable_coverage[coverage_df$income_level_iso3c == income_group]
    run_coverage_df <- data.frame(coverage = c(low_coverage, mid_coverage, high_coverage, variable_coverage),
                                  coverage_scenario = c("low", "mid", "high", "variable"))
    
    ## Loading in the raw fits (without the BPSV)
    temp <- readRDS(paste0("outputs/Figure3_SC2_Counterfactual_Impact/raw/raw_", iso_list[i], "_fit.rds"))
    
    ## Setting up parallelisation of all of the scenarios
    iso <- iso_list[i]
    num_cores <- 16
    scenarios_df <- expand.grid(bpsv_start_date = bpsv_start_dates, coverage_scenario = c("low", "mid", "high", "variable")) %>%
      left_join(bpsv_start_dates_df, "bpsv_start_date") %>%
      left_join(run_coverage_df, "coverage_scenario") 
    
    ## Running the scenarios
    x <- parallel::mclapply(1:nrow(scenarios_df), mc.cores = num_cores, function(scenario) {
      temp_evaluation <- evaluate_country_impact2(original_fit = temp, country_iso = iso, 
                                                  vaccination_rate = vaccination_rate, 
                                                  bpsv_start_date = scenarios_df$bpsv_start_date[scenario], 
                                                  coverage = scenarios_df$coverage[scenario])
      temp_evaluation$start_trigger <- scenarios_df$start_trigger[scenario]
      temp_evaluation$coverage_scenario <- scenarios_df$coverage_scenario[scenario]
      temp_evaluation$coverage <- scenarios_df$coverage[scenario]
      
      file_string <- paste0("Detect", scenarios_df$start_trigger[scenario], "_", scenarios_df$coverage_scenario[scenario], "VaccCoverage")
      saveRDS(object = temp_evaluation, 
              file = paste0("outputs/Figure3_SC2_Counterfactual_Impact/", file_string, "_", iso, "_fit.rds"))
      
      return(temp_evaluation)
    })
  } 
  all_files <- list.files(path = "outputs/Figure3_SC2_Counterfactual_Impact", full.names = TRUE)
  iso_files <- all_files[grepl(iso, all_files)]
  deaths <- iso_files %>%
    map_dfr(readRDS) %>% 
    group_by(start_trigger, coverage_scenario, coverage, scenario) %>%
    summarise(med = sum(med)) %>%
    pivot_wider(names_from = scenario, values_from = med) %>%
    rename(empirical_deaths = no_bpsv_deaths, deaths_BPSV = bpsv_deaths)
  deaths$iso <- iso_list[i]
  deaths_df <- rbind(deaths_df, deaths)
  print(i)
}

## Overall deaths averted
overall_deaths_averted <- deaths_df %>%
  group_by(start_trigger) %>%
  summarise(total_bpsv_deaths = sum(deaths_BPSV),
            total_empirical_deaths = sum(empirical_deaths))


data_list <- vector(mode = "list", length = length(iso_list))
for (i in 1:length(iso_list)) {
  bpsv <- readRDS(paste0("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual_vaccRate_firstDeathStart/BPSV_counterfactual_vaccRate_firstDeathStart_", iso_list[i], "_fit.rds")) %>%
    mutate(iso = iso_list[i])
  data_list[[i]] <- bpsv
}
impact_deaths_df <- bind_rows(data_list)

reported_deaths <- impact_deaths_df %>%
  filter(scenario == "no_bpsv_deaths") %>%
  select(iso, scenario, date, reported_deaths) %>%
  group_by(date) %>%
  summarise(total_reported_deaths = sum(reported_deaths, na.rm = TRUE)) %>%
  filter(date < as.Date("2020-11-28"))
plot(reported_deaths$date, reported_deaths$total_reported_deaths)
sum(reported_deaths$total_reported_deaths)

overall_impact_deaths <- impact_deaths_df %>%
  group_by(scenario, date) %>%
  filter(date < as.Date("2020-11-28")) %>%
  dplyr::summarise(total_deaths = sum(med, na.rm = TRUE), 
                   total_low = sum(`025`, na.rm = TRUE),
                   total_high = sum(`975`, na.rm = TRUE)) %>%
  mutate(cumulative = cumsum(total_deaths),
         cumulative_low = cumsum(total_low),
         cumulative_high = cumsum(total_high)) %>%
  pivot_wider(names_from = "scenario",
              values_from = c("total_deaths", "total_low", "total_high", "cumulative", "cumulative_low", "cumulative_high"))

a <- ggplot(overall_impact_deaths) +
  geom_line(aes(x = date, y = cumulative_no_bpsv_deaths), colour = "#748386", linewidth = 1) +
  geom_line(aes(x = date, y = cumulative_bpsv_deaths), colour = "#E9614F", linewidth = 1) +
  geom_ribbon(aes(x = date, ymin = cumulative_bpsv_deaths, ymax = cumulative_no_bpsv_deaths), 
              alpha = 0.2, fill = "#F2C7C1") +
  labs(x = "", y = "Cumulative COVID-19 Deaths") +
  theme_bw() +
  scale_y_continuous(labels = c("1M", "2M", "3M", "4M", "5M", "6M"),
                     breaks = c(1e6, 2e6, 3e6, 4e6, 5e6, 6e6))

b <- ggplot(overall_impact_deaths) +
  geom_line(aes(x = date, y = total_deaths_no_bpsv_deaths), colour = "#748386", linewidth = 1) +
  geom_line(aes(x = date, y = total_deaths_bpsv_deaths), colour = "#E9614F", linewidth = 1) +
  geom_ribbon(aes(x = date, ymin = total_deaths_bpsv_deaths, ymax = total_deaths_no_bpsv_deaths), 
              alpha = 0.2, fill = "#F2C7C1") +
  labs(x = "", y = "Daily COVID-19 Deaths") +
  theme_bw() 

c <- cowplot::plot_grid(a, b, nrow = 2)
ggsave(filename = "figures/Figure_3_BPSV_SC2_Impact/NEW_Figure3_total_plots.pdf",
       plot = c,
       width = 3.33, height = 3)

# Get high-quality natural earth data
world <- ne_countries(scale = "medium", returnclass = "sf")
data <- data.frame(country = world$geounit, value = runif(n = nrow(world))) ## add in BPSV impact here
world_robinson <- st_transform(world, crs = "ESRI:54030")
merged_data <- merge(world_robinson, deaths_df, by.x = 'iso_a3', by.y = 'iso', all.x = TRUE)
merged_data$perc_averted <- 100 * (1 - merged_data$deaths_BPSV / merged_data$empirical_deaths)

# Getting bounding box and graticules
load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
NE_box_rob <- spTransform(NE_box, CRSobj = PROJ)
NE_graticules_rob <- spTransform(NE_graticules, CRSobj = PROJ)

# Plotting the output 
world_map <- ggplot() +
  geom_sf(data = merged_data, aes(fill = perc_averted), color = "black") +
  geom_path(data=NE_graticules_rob, aes(long, lat, group=group), linetype="dashed", color="grey50", size = 0.05) +
  geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.1) + 
  coord_sf(crs = st_crs("ESRI:54030")) +
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white")) +
  theme_void() +
  theme(legend.position = "left") +
  scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(40, 100), name = "% Deaths\nAverted")

## Individual plots for Italy, Iran and India
rescaled_values <- (merged_data$perc_averted - 40) / (100 - 40)
rescaled_values <- pmin(pmax(rescaled_values, 0), 1)  # Ensure values are within [0, 1]
n_colors <- 256  # You can adjust this number
viridis_palette <- viridis::viridis(n_colors, option = "rocket", direction = -1)
country_colors <- viridis_palette[as.integer(rescaled_values * (n_colors - 1)) + 1]
color_mapping <- data.frame(Country = merged_data$iso_a3, Color = country_colors)

### Italy
age_groups <- c("65-69", "70-74", "75-79", "80+")
squire:::get_population("Italy") %>%
  filter(age_group %in% age_groups) %>%
  summarise(total = sum(n))
bpsv_ITA <- readRDS("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual_vaccRate_firstDeathStart/BPSV_counterfactual_vaccRate_firstDeathStart_ITA_fit.rds") %>%
  mutate(iso = "ITA")
ita_plot <- ggplot(data = bpsv_ITA) +
    geom_point(aes(x = date, y = reported_deaths), col = "grey") +
    geom_ribbon(aes(x = date, ymin = `025`, ymax = `975`, fill = scenario), alpha = 0.2) +
    geom_line(aes(x = date, y = med, col = scenario)) +
    theme_bw() +
    scale_colour_manual(values = c(color_mapping$Color[color_mapping$Country == "ITA"], "#001A23")) +
    scale_fill_manual(values = c(color_mapping$Color[color_mapping$Country == "ITA"], "#001A23")) +
    labs(x = "Date", y = "COVID-19 Deaths Per Day") +
    theme(legend.position = "none")

squire:::get_population("Iran") %>%
  filter(age_group %in% age_groups) %>%
  summarise(total = sum(n))
bpsv_IRN <- readRDS("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual_vaccRate_firstDeathStart/BPSV_counterfactual_vaccRate_firstDeathStart_IRN_fit.rds") %>%
  mutate(iso = "IRN")
irn_plot <- ggplot(data = bpsv_IRN) +
  geom_point(aes(x = date, y = reported_deaths), col = "grey") +
  geom_ribbon(aes(x = date, ymin = `025`, ymax = `975`, fill = scenario), alpha = 0.2) +
  geom_line(aes(x = date, y = med, col = scenario)) +
  theme_bw() +
  scale_colour_manual(values = c(color_mapping$Color[color_mapping$Country == "IRN"], "#001A23")) +
  scale_fill_manual(values = c(color_mapping$Color[color_mapping$Country == "IRN"], "#001A23")) +
  labs(x = "Date", y = "COVID-19 Deaths Per Day") +
  theme(legend.position = "none")

bpsv_COL <- readRDS("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual_vaccRate_firstDeathStart/BPSV_counterfactual_vaccRate_firstDeathStart_COL_fit.rds") %>%
  mutate(iso = "COL")
col_plot <- ggplot(data = bpsv_COL) +
  geom_point(aes(x = date, y = reported_deaths), col = "grey") +
  geom_ribbon(aes(x = date, ymin = `025`, ymax = `975`, fill = scenario), alpha = 0.2) +
  geom_line(aes(x = date, y = med, col = scenario)) +
  theme_bw() +
  scale_colour_manual(values = c(color_mapping$Color[color_mapping$Country == "COL"], "#001A23")) +
  scale_fill_manual(values = c(color_mapping$Color[color_mapping$Country == "COL"], "#001A23")) +
  labs(x = "Date", y = "COVID-19 Deaths Per Day") +
  theme(legend.position = "none")

squire:::get_population("Bangladesh") %>%
  filter(age_group %in% age_groups) %>%
  summarise(total = sum(n))
bpsv_BGD <- readRDS("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual/BPSV_counterfactual_BGD_fit.rds") %>%
  mutate(iso = "BGD")
bgd_plot <- ggplot(data = bpsv_BGD) +
  geom_point(aes(x = date, y = reported_deaths), col = "grey") +
  geom_ribbon(aes(x = date, ymin = `025`, ymax = `975`, fill = scenario), alpha = 0.2) +
  geom_line(aes(x = date, y = med, col = scenario)) +
  theme_bw() +
  scale_colour_manual(values = c(color_mapping$Color[color_mapping$Country == "BGD"], "#001A23")) +
  scale_fill_manual(values = c(color_mapping$Color[color_mapping$Country == "BGD"], "#001A23")) +
  labs(x = "Date", y = "COVID-19 Deaths Per Day") +
  theme(legend.position = "none")

row1 <- cowplot::plot_grid(world_map, ita_plot, nrow = 1, rel_widths = c(2, 1), labels = c("A", "B"))
row2 <- cowplot::plot_grid(col_plot, irn_plot, bgd_plot, ncol = 3, labels = c("C", "D", "E"))
full_plot <- cowplot::plot_grid(row1, row2, nrow = 2)
ggsave(filename = "figures/Figure_3_BPSV_SC2_Impact/NEW_Figure3_full_plots.pdf",
       plot = full_plot,
       width = 12, height = 8)

country_row1 <- cowplot::plot_grid(NULL, ita_plot, nrow = 1, rel_widths = c(2, 1), labels = c("A", "B"))
country_row2 <- cowplot::plot_grid(col_plot, irn_plot, bgd_plot, ncol = 3, labels = c("C", "D", "E"))
country_plots <- cowplot::plot_grid(country_row1, country_row2, nrow = 2)
ggsave(filename = "figures/Figure_3_BPSV_SC2_Impact/NEW_Figure3_country_plots.pdf",
       plot = country_plots,
       width = 10, height = 6)



  
# Bangladesh, Colombia, Italy and Iran

### SCRAP CODE ###
# overall_ITA <- evaluate_country_impact("ITA")
# saveRDS(object = overall_ITA, file = "ITA_BPSV_counterfactual_impact.rds")
# overall_IRN <- evaluate_country_impact("IRN")
# saveRDS(object = overall_IRN, file = "IRN_BPSV_counterfactual_impact.rds")
# 
# overall_ITA <- readRDS("ITA_BPSV_counterfactual_impact.rds")
# overall_IRN <- readRDS("IRN_BPSV_counterfactual_impact.rds")
# 
# ita_time_plot <- ggplot(data = overall_ITA) +
#   geom_point(aes(x = date, y = reported_deaths), col = "grey") +
#   geom_ribbon(aes(x = date, ymin = `025`, ymax = `975`, fill = scenario), alpha = 0.2) +
#   geom_line(aes(x = date, y = med, col = scenario)) +
#   theme_bw() +
#   scale_colour_manual(values = c("#E67552", "#001A23")) +
#   scale_fill_manual(values = c("#E67552", "#001A23")) +
#   labs(x = "Date", y = "COVID-19 Deaths Per Day") +
#   theme(legend.position = "none")
# 
# deaths_ITA <- overall_ITA %>%
#   group_by(scenario) %>%
#   summarise(med = sum(med),
#             `025` = sum(`025`),
#             `975` = sum(`975`))
# 
# ita_deaths_plot <- ggplot(deaths_ITA) +
#   geom_bar(aes(x = scenario, y = med, fill = scenario), stat = "identity") +
#   geom_errorbar(aes(x = scenario, ymin = `025`, ymax = `975`, group = scenario), width = 0.3) +
#   labs(x = "", y = "COVID-19 Deaths") +
#   scale_x_discrete(labels = c("Deaths with\nBPSV", "Actual\nDeaths")) +
#   scale_fill_manual(values = c("#E67552", "#001A23")) +
#   theme_bw() +
#   theme(legend.position = "none")
# 
# irn_time_plot <- ggplot(data = overall_IRN) +
#   geom_point(aes(x = date, y = reported_deaths), col = "grey") +
#   geom_ribbon(aes(x = date, ymin = `025`, ymax = `975`, fill = scenario), alpha = 0.2) +
#   geom_line(aes(x = date, y = med, col = scenario)) +
#   theme_bw() +
#   scale_colour_manual(values = c("#E67552", "#001A23")) +
#   scale_fill_manual(values = c("#E67552", "#001A23")) +
#   labs(x = "Date", y = "COVID-19 Deaths Per Day") +
#   theme(legend.position = "none")
# 
# deaths_IRN <- overall_IRN %>%
#   group_by(scenario) %>%
#   summarise(med = sum(med),
#             `025` = sum(`025`),
#             `975` = sum(`975`))
# 
# irn_deaths_plot <- ggplot(deaths_IRN) +
#   geom_bar(aes(x = scenario, y = med, fill = scenario), stat = "identity") +
#   geom_errorbar(aes(x = scenario, ymin = `025`, ymax = `975`, group = scenario), width = 0.3) +
#   labs(x = "", y = "COVID-19 Deaths") +
#   scale_x_discrete(labels = c("Deaths with\nBPSV", "Actual\nDeaths")) +
#   scale_fill_manual(values = c("#E67552", "#001A23")) +
#   theme_bw() +
#   theme(legend.position = "none")
# 
# italy <- cowplot::plot_grid(ita_time_plot, ita_deaths_plot, 
#                             rel_widths = c(2.5, 1), labels = c("A", "B"))
# iran <- cowplot::plot_grid(irn_time_plot, irn_deaths_plot, 
#                            rel_widths = c(2.5, 1), labels = c("C", "D"))
# 
# cowplot::plot_grid(italy, iran, nrow = 2)
# ggplot(overall_impact_deaths) +
#   geom_line(aes(x = date, y = cumulative_no_bpsv_deaths)) +
#   geom_ribbon(aes(x = date, ymin = cumulative_low_no_bpsv_deaths, ymax = cumulative_high_no_bpsv_deaths ), alpha = 0.2) +
#   geom_line(aes(x = date, y = cumulative_bpsv_deaths)) +
#   geom_ribbon(aes(x = date, ymin = cumulative_low_bpsv_deaths, ymax = cumulative_high_bpsv_deaths ), alpha = 0.2) +
#   labs(x = "", y = "COVID-19 Deaths") +
#   theme_bw() +
#   scale_y_continuous(labels = c("1M", "2M", "3M", "4M", "5M", "6M"),
#                      breaks = c(1e6, 2e6, 3e6, 4e6, 5e6, 6e6))
## 1 Death Start
# temp_overall_1 <- evaluate_country_impact2(original_fit = temp, country_iso = iso_list[i], 
#                                            vaccination_rate = vaccination_rate, 
#                                            bpsv_start_date = bpsv_start_date_1, coverage = 0.6)
# temp_overall_1$start_trigger <- "1_deaths"
# saveRDS(object = temp_overall_1, file = paste0("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual_VaccRate_1DeathsStart/BPSV_counterfactual_VaccRate_1DeathsStart_", iso_list[i], "_fit.rds"))
# 
# ## 10 Deaths Start
# temp_overall_10 <- evaluate_country_impact2(original_fit = temp, country_iso = iso_list[i], 
#                                             vaccination_rate = vaccination_rate, 
#                                             bpsv_start_date = bpsv_start_date_10, coverage = 0.6)
# temp_overall_10$start_trigger <- "10_deaths"
# saveRDS(object = temp_overall_10, file = paste0("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual_VaccRate_10DeathsStart/BPSV_counterfactual_VaccRate_10DeathsStart_", iso_list[i], "_fit.rds"))
# 
# ## 100 Deaths Start
# temp_overall_100 <- evaluate_country_impact2(original_fit = temp, country_iso = iso_list[i], 
#                                              vaccination_rate = vaccination_rate, 
#                                              bpsv_start_date = bpsv_start_date_100, coverage = 0.6)
# temp_overall_100$start_trigger <- "100_deaths"
# saveRDS(object = temp_overall_100, file = paste0("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual_VaccRate_100DeathsStart/BPSV_counterfactual_VaccRate_100DeathsStart_", iso_list[i], "_fit.rds"))
# 
# ## 1000 Deaths Start
# temp_overall_1000 <- evaluate_country_impact2(original_fit = temp, country_iso = iso_list[i], 
#                                               vaccination_rate = vaccination_rate, 
#                                               bpsv_start_date = bpsv_start_date_1000, coverage = 0.6)
# temp_overall_1000$start_trigger <- "1000_deaths"
# saveRDS(object = temp_overall_1000, file = paste0("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual_VaccRate_1000DeathsStart/BPSV_counterfactual_VaccRate_1000DeathsStart_", iso_list[i], "_fit.rds"))