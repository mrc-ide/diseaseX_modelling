# Load required libraries
library(tidyverse); library(rnaturalearth); library(sf); #library(rgdal); # library(squire.page); 
source("functions/extract_process_generate_fits.R")
source("functions/run_sars_x.R")
source("functions/helper_functions.R")

get_country_draws <- function(country_iso = "IRN") {
  
  # Getting fits for country with excess deaths and extracting the Rt
  excess <- grab_fit(country_iso, TRUE, TRUE)
  out <- squire.page:::generate_draws.rt_optimised(excess)
  daily <- get_deaths_infections_hosps_time(out) %>%
    group_by(replicate) %>%
    arrange(date) %>%
    group_by(date)
  
  # Results generation for each of the 100 draws
  temp_list <- vector(mode = "list", length = length(out$samples))
  for (i in 1:length(out$samples)) {
    
    index <- i
    if (is.null(out$samples[[index]]$prob_severe_multiplier)) {
      out$samples[[index]]$prob_severe_multiplier <- rep(1, length(out$samples[[index]]$tt_prob_hosp_multiplier))
      out$samples[[index]]$tt_prob_severe_multiplier <- out$samples[[index]]$tt_prob_hosp_multiplier
    }
    tt_Rt <- out$samples[[index]]$tt_R0
    tt_R0 <- tt_Rt + 1
    check1 <- squire.page:::run_booster(time_period = max(out$samples[[index]]$tt_R0),
                                        population = squire::get_population(country = out$parameters$country)$n,                                                 
                                        contact_matrix_set = squire::get_mixing_matrix(country = out$parameters$country),                                                   
                                        R0 = out$samples[[index]]$R0,     
                                        tt_R0 = tt_R0, 
                                        hosp_bed_capacity = out$parameters$hosp_bed_capacity,                                     
                                        ICU_bed_capacity = out$parameters$ICU_bed_capacity,                                       
                                        
                                        prob_hosp = out$samples[[index]]$prob_hosp,
                                        prob_severe = out$samples[[index]]$prob_severe,
                                        prob_severe_death_treatment = out$samples[[index]]$prob_severe_death_treatment, 
                                        prob_severe_death_no_treatment = out$samples[[index]]$prob_severe_death_no_treatment, 
                                        prob_non_severe_death_treatment = out$samples[[index]]$prob_non_severe_death_treatment, 
                                        prob_non_severe_death_no_treatment = out$samples[[index]]$prob_non_severe_death_no_treatment, 
                                        dur_R = out$samples[[index]]$dur_R,   
                                        tt_dur_R = out$samples[[index]]$tt_dur_R,
                                        
                                        prob_hosp_multiplier = out$samples[[index]]$prob_hosp_multiplier,
                                        tt_prob_hosp_multiplier = out$samples[[index]]$tt_prob_hosp_multiplier,
                                        prob_severe_multiplier = out$samples[[index]]$prob_severe_multiplier,
                                        tt_prob_severe_multiplier = out$samples[[index]]$tt_prob_severe_multiplier,
                                        dur_V = out$samples[[index]]$dur_V,
                                        tt_dur_V = out$samples[[index]]$tt_dur_V,
                                        vaccine_efficacy_infection = out$samples[[index]]$vaccine_efficacy_infection,
                                        vaccine_efficacy_disease = out$samples[[index]]$vaccine_efficacy_disease,
                                        tt_vaccine_efficacy_disease = out$samples[[index]]$tt_vaccine_efficacy_disease,
                                        tt_vaccine_efficacy_infection = out$samples[[index]]$tt_vaccine_efficacy_infection,
                                        primary_doses = 0,
                                        booster_doses = 0, 
                                        
                                        seeding_cases = NULL,
                                        init = seed_infections(excess$squire_model, excess$parameters$country, out$samples[[index]]$initial_infections))
    check1d <- nimue::format(check1, compartments = "D", summaries = "deaths") %>%
      filter(t > 1, compartment == "deaths") %>%
      mutate(index = i)
    temp_list[[i]] <- check1d
  }
  
  rep_summary <- bind_rows(temp_list) %>%
    group_by(replicate) %>%
    arrange(t) %>%
    group_by(t)  %>%
    summarise(
      across(c(value), ~median(.x, na.rm=TRUE), .names = "no_bpsv_deaths_med"),
      across(c(value), ~quantile(.x, 0.025, na.rm=TRUE), .names = "no_bpsv_deaths_025"),
      across(c(value), ~quantile(.x, 0.975, na.rm=TRUE), .names = "no_bpsv_deaths_975"),
      .groups = "drop") 
  
  return(list(country_iso = country_iso,
              model_fit = rep_summary,
              daily = daily,
              out = list(samples = out$samples, 
                         parameters = out$parameters, 
                         inputs = out$inputs)))
  
}

# Defining function
evaluate_country_impact <- function(original_fit, country_iso = "IRN") {
  
  excess <- grab_fit(country_iso, TRUE, TRUE)
  out <- original_fit$out
  rep_summary <- original_fit$model_fit
  daily <- original_fit$daily
  
  ### Generating the counterfactual impact of a BPSV
  temp <- create_scenarios(R0 = 3, specific_vaccine_start = 500) %>%
    filter(vaccine_scenario == "both_vaccines")
  
  # Generating vaccine coverage matrix
  coverage_spec <- temp$coverage_spec
  coverage_bpsv <- temp$coverage_bpsv
  priority_age_groups <- temp$min_age_group_index_priority:17
  vaccination_age_groups <- temp$min_age_group_index_non_priority:17
  vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(coverage_spec, length(priority_age_groups)), 
                                   rep(0, 17 - length(vaccination_age_groups)), rep(coverage_spec, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)
  
  # Generating vaccine dose series
  detection_time <- 1
  bpsv_start <- 1
  bpsv_protection_delay <- 14
  specific_vaccine_start <- 450
  specific_protection_delay <- 14
  vaccine_doses <- create_vaccination_dose_series(country = out$parameters$country, 
                                                  population_size = sum(squire::get_population(country = out$parameters$country)$n), 
                                                  detection_time = detection_time, 
                                                  vaccine_scenario = "both_vaccines",
                                                  bpsv_start = bpsv_start, 
                                                  bpsv_protection_delay = bpsv_protection_delay,
                                                  specific_vaccine_start = specific_vaccine_start,
                                                  specific_protection_delay = specific_protection_delay,
                                                  vaccination_rate_bpsv = temp$vaccination_rate_bpsv,
                                                  vaccination_rate_spec = temp$vaccination_rate_spec,
                                                  coverage_bpsv = temp$coverage_bpsv,
                                                  coverage_spec = temp$coverage_spec,
                                                  min_age_group_index_priority = temp$min_age_group_index_priority,
                                                  runtime = max(out$samples[[1]]$tt_R0) + 200)
  primary_doses <- vaccine_doses$primary_doses
  second_doses <- vaccine_doses$second_doses
  booster_doses <- vaccine_doses$booster_doses
  
  ## Generating vaccine efficacy matrices
  ve_bpsv <- list(infection = temp$efficacy_infection_bpsv, disease = temp$efficacy_disease_bpsv)
  ve_spec <- list(infection = 0, disease = 0)
  
  ### Infection Efficacy
  ve_i_elderly_bpsv_campaign <- c(ve_bpsv$infection, ve_bpsv$infection, ve_bpsv$infection, 0, ve_spec$infection, ve_spec$infection, 0) 
  ve_i_elderly_spec_campaign <- c(ve_spec$infection, ve_bpsv$infection, ve_bpsv$infection, 0, ve_spec$infection, ve_spec$infection, 0) 
  ve_i_non_elderly <- c(ve_spec$infection, ve_spec$infection, ve_spec$infection, 0, ve_spec$infection, ve_spec$infection, 0)
  vaccine_efficacy_infection_bpsv_campaign <- matrix(c(rep(ve_i_non_elderly, 17 - length(priority_age_groups)), 
                                                       rep(ve_i_elderly_bpsv_campaign, length(priority_age_groups)) ), ncol = 17)
  vaccine_efficacy_infection_spec_campaign <- matrix(c(rep(ve_i_non_elderly, 17 - length(priority_age_groups)), 
                                                       rep(ve_i_elderly_spec_campaign, length(priority_age_groups)) ), ncol = 17)
  
  ### Disease Efficacy
  ve_d_elderly_bpsv_campaign <- c(ve_bpsv$disease, ve_bpsv$disease, ve_bpsv$disease, 0, ve_spec$disease, ve_spec$disease, 0) 
  ve_d_elderly_spec_campaign <- c(ve_spec$disease, ve_bpsv$disease, ve_bpsv$disease, 0, ve_spec$disease, ve_spec$disease, 0) 
  ve_d_non_elderly <- c(ve_spec$disease, ve_spec$disease, ve_spec$disease, 0, ve_spec$disease, ve_spec$disease, 0)
  vaccine_efficacy_disease_bpsv_campaign <- matrix(c(rep(ve_d_non_elderly, 17 - length(priority_age_groups)), 
                                                     rep(ve_d_elderly_bpsv_campaign, length(priority_age_groups)) ), ncol = 17)
  vaccine_efficacy_disease_spec_campaign <- matrix(c(rep(ve_d_non_elderly, 17 - length(priority_age_groups)), 
                                                     rep(ve_d_elderly_spec_campaign, length(priority_age_groups)) ), ncol = 17)
  
  # Duration of vaccination immunity
  dur_V <- matrix(data = c(c(rep(temp$dur_spec, min(priority_age_groups) - 1), rep(temp$dur_bpsv, length(priority_age_groups))),    ## duration of primary series protection secondary dose 1 -> secondary dose 2 compartments (bpsv for elderly, specific vaccine for everyone else in this scenario)
                           c(rep(temp$dur_spec, min(priority_age_groups) - 1), rep(temp$dur_bpsv, length(priority_age_groups))),    ## duration of primary series protection secondary dose 2 -> waned compartments (bpsv for elderly, specific vaccine for everyone else in this scenario)
                           c(rep(temp$dur_spec, min(priority_age_groups) - 1), rep(temp$dur_spec, length(priority_age_groups))),    ## duration of booster protection (specific for elderly, not used for everyone else)
                           c(rep(temp$dur_spec, min(priority_age_groups) - 1), rep(temp$dur_spec, length(priority_age_groups)))),   ## duration of booster protection (specific for elderly, not used for everyone else)
                  nrow = 4, ncol = 17, byrow = TRUE)
  
  ## Boosters are specific vaccine for elderly population, and we're not modelling that here
  vaccine_booster_initial_coverage <- rep(0, 17)
  vaccine_booster_follow_up_coverage <- rep(0, 17)
  
  temp_list_BPSV <- vector(mode = "list", length = 100L)
  for (i in 1:length(out$samples)) {
    index <- i
    tt_Rt <- out$samples[[index]]$tt_R0
    tt_R0 <- tt_Rt + 1
      # excess$squire_mode$parameter_func <- squire.page.sarsX:::nimue_booster_model()$parameter_func
    init <- seed_infections(excess$squire_model, excess$parameters$country, out$samples[[index]]$initial_infections)
    BPSV <- squire.page.sarsX:::run_booster(time_period = max(out$samples[[index]]$tt_R0) + 200,
                                            population = squire::get_population(country = out$parameters$country)$n,                                                 
                                            contact_matrix_set = squire::get_mixing_matrix(country = out$parameters$country),                                                   
                                            R0 = out$samples[[index]]$R0,     
                                            tt_R0 = tt_R0, 
                                            hosp_bed_capacity = out$parameters$hosp_bed_capacity,                                     
                                            ICU_bed_capacity = out$parameters$ICU_bed_capacity,                                       
                                            
                                            prob_hosp = out$samples[[index]]$prob_hosp,
                                            prob_severe = out$samples[[index]]$prob_severe,
                                            prob_severe_death_treatment = out$samples[[index]]$prob_severe_death_treatment, 
                                            prob_severe_death_no_treatment = out$samples[[index]]$prob_severe_death_no_treatment, 
                                            prob_non_severe_death_treatment = out$samples[[index]]$prob_non_severe_death_treatment, 
                                            prob_non_severe_death_no_treatment = out$samples[[index]]$prob_non_severe_death_no_treatment, 
                                            dur_R = out$samples[[index]]$dur_R,   
                                            tt_dur_R = out$samples[[index]]$tt_dur_R,
                                            
                                            prob_hosp_multiplier = out$samples[[index]]$prob_hosp_multiplier,
                                            tt_prob_hosp_multiplier = out$samples[[index]]$tt_prob_hosp_multiplier,
                                            prob_severe_multiplier = out$samples[[index]]$prob_severe_multiplier,
                                            tt_prob_severe_multiplier = out$samples[[index]]$tt_prob_severe_multiplier,
                                            
                                            ## BPSV specific stuff
                                            vaccine_coverage_mat = vaccine_coverage_mat,
                                            primary_doses = primary_doses,  
                                            second_doses = second_doses,
                                            booster_doses = booster_doses,   
                                            vaccine_efficacy_infection = list(vaccine_efficacy_infection_bpsv_campaign, vaccine_efficacy_infection_spec_campaign),
                                            tt_vaccine_efficacy_disease = c(0, temp$specific_vaccine_start - 5),
                                            vaccine_efficacy_disease = list(vaccine_efficacy_disease_bpsv_campaign, vaccine_efficacy_disease_spec_campaign),
                                            tt_vaccine_efficacy_infection = c(0, temp$specific_vaccine_start - 5),
                                            dur_V = dur_V,
                                            vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
                                            vaccine_booster_initial_coverage = vaccine_booster_initial_coverage,
                                            
                                            ## Misc
                                            seeding_cases = NULL,
                                            init = init)
    
    check_BPSV <- nimue::format(BPSV, compartments = "D", summaries = "deaths") %>%
      filter(t > 1, compartment == "deaths") %>%
      mutate(index = i)
    temp_list_BPSV[[i]] <- check_BPSV
  }
  
  rep_bpsv <- bind_rows(temp_list_BPSV) %>%
    group_by(replicate) %>%
    arrange(t) %>%
    group_by(t)  %>%
    summarise(
      across(c(value), ~median(.x, na.rm=TRUE), .names = "bpsv_deaths_med"),
      across(c(value), ~quantile(.x, 0.025, na.rm=TRUE), .names = "bpsv_deaths_025"),
      across(c(value), ~quantile(.x, 0.975, na.rm=TRUE), .names = "bpsv_deaths_975"),
      .groups = "drop")
  
  overall <- rep_bpsv %>%
    left_join(rep_summary, by = "t") %>%
    filter(t < 366)
  
  fit <- daily %>%
    filter(replicate == 1)
  overall$date <- fit$date[1:364] 
  overall <- overall %>%
    pivot_longer(cols = -c(t, date), 
                 names_to = c("scenario", "metric"), 
                 names_pattern = "(.*)_(.*)",
                 values_to = "values") %>%
    pivot_wider(names_from = metric,
                values_from = values)
  
  deaths_data <- out$inputs$data %>%
    mutate(date = date_start + 3, 
           reported_deaths = deaths / 7) %>%
    select(date, reported_deaths)
  
  overall <- overall %>%
    left_join(deaths_data, by = "date")
  
  return(overall) 
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

## Generating BPSV impact (need to go back and figure out when I'm implicitly deploying the BPSV and
##                         make a PDF to check everything looks reasonable)
iso_list <- substr(list.files(path = "outputs/Figure3_SC2_Counterfactual_Impact/raw/"), 5, 7)
deaths_df <- data.frame(iso = rep(NA_character_, length(iso_list)), 
                        empirical_deaths = rep(NA_real_, length(iso_list)), 
                        deaths_BPSV = rep(NA_real_, length(iso_list)))
new_run <- FALSE
for (i in 1:length(iso_list)) {
  if (new_run) {
    temp <- readRDS(paste0("outputs/Figure3_SC2_Counterfactual_Impact/raw/raw_", iso_list[i], "_fit.rds"))
    temp_overall <- evaluate_country_impact(original_fit = temp, country_iso = iso_list[i])
    saveRDS(object = temp_overall, file = paste0("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual/BPSV_counterfactual_", iso_list[i], "_fit.rds"))
  } else {
    temp_overall <- readRDS(paste0("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual/BPSV_counterfactual_", iso_list[i], "_fit.rds"))
  }
  deaths <- temp_overall %>%
    group_by(scenario) %>%
    summarise(med = sum(med)) %>%
    pivot_wider(names_from = scenario, values_from = med)
  deaths_df$iso[i] <- iso_list[i]
  deaths_df$empirical_deaths[i] <- deaths$no_bpsv_deaths
  deaths_df$deaths_BPSV[i] <- deaths$bpsv_deaths
  print(i)
}

data_list <- vector(mode = "list", length = length(iso_list))
for (i in 1:length(iso_list)) {
  bpsv <- readRDS(paste0("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual/BPSV_counterfactual_", iso_list[i], "_fit.rds")) %>%
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
bpsv_ITA <- readRDS("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual/BPSV_counterfactual_ITA_fit.rds") %>%
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
bpsv_IRN <- readRDS("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual/BPSV_counterfactual_IRN_fit.rds") %>%
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

bpsv_COL <- readRDS("outputs/Figure3_SC2_Counterfactual_Impact/BPSV_counterfactual/BPSV_counterfactual_COL_fit.rds") %>%
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
