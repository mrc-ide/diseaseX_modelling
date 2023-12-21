# Load required libraries
library(tidyverse); library(squire.page)

# Sourcing required functions
source("Q_Drive_Copy2/Active_Research_Projects/diseaseX_modelling/functions/extract_process_generate_fits.R")

# Getting fits for Iran with excess deaths and extracting the Rt
excess <- grab_fit("ITA", TRUE, TRUE)
out <- squire.page:::generate_draws.rt_optimised(excess)
daily <- get_deaths_infections_hosps_time(out) %>%
  group_by(replicate) %>%
  arrange(date) %>%
  group_by(date)

# Results generation for each of the 100 draws
temp_list <- vector(mode = "list", length = 100L)
for (i in 1:100) {
  index <- i
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
  print(i)
}

## Checking for random index
index2 <- sample(1:100, 1)
check <- temp_list[[index2]]
fit <- daily %>%
  filter(replicate == index2)
plot(fit$deaths[1:365], type = "l", col = "red")
lines(check$t[1:365], check$value[1:365], type = "l")

## Summarising fitted 
daily_summary <- daily %>%
  select(-infections) %>%
  group_by(replicate) %>%
  arrange(date) %>%
  group_by(date) %>%
  summarise(
    across(c(deaths), ~median(.x, na.rm=TRUE), .names = "{col}_med"),
    across(c(deaths), ~quantile(.x, 0.025, na.rm=TRUE), .names = "{col}_025"),
    across(c(deaths), ~quantile(.x, 0.975, na.rm=TRUE), .names = "{col}_975"),
    .groups = "drop")

rep_summary <- bind_rows(temp_list) %>%
  group_by(replicate) %>%
  arrange(t) %>%
  group_by(t)  %>%
  summarise(
    across(c(value), ~median(.x, na.rm=TRUE), .names = "deaths_med"),
    across(c(value), ~quantile(.x, 0.025, na.rm=TRUE), .names = "deaths_025"),
    across(c(value), ~quantile(.x, 0.975, na.rm=TRUE), .names = "deaths_975"),
    .groups = "drop")

## Comparing outputs 
sum(rep_summary$value_med[1:365]) # my re-running using squire.page - seeding cases not fixed to match how you do it
sum(daily_summary$deaths_med[1:365]) # the original deaths curve (or at least, the median)

## Graphical comparison
plot(daily_summary$deaths_med[1:365], type = "l", col = "red")
lines(rep_summary$t[1:365], rep_summary$value_med[1:365], type = "l")

### Generating the counterfactual impact of a BPSV
source("Q_Drive_Copy2/Active_Research_Projects/diseaseX_modelling/main.R")
source("Q_Drive_Copy2/Active_Research_Projects/diseaseX_modelling/functions/run_sars_x.R")
source("Q_Drive_Copy2/Active_Research_Projects/diseaseX_modelling/functions/helper_functions.R")

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
                                                runtime = max(out$samples[[1]]$tt_R0))
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
for (i in 1:100) {
  index <- i
  tt_Rt <- out$samples[[index]]$tt_R0
  tt_R0 <- tt_Rt + 1
  BPSV <- squire.page.sarsX:::run_booster(time_period = max(out$samples[[index]]$tt_R0),
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
                                          init = seed_infections(excess$squire_model, excess$parameters$country, out$samples[[index]]$initial_infections))
  
  check_BPSV <- nimue::format(BPSV, compartments = "D", summaries = "deaths") %>%
    filter(t > 1, compartment == "deaths") %>%
    mutate(index = i)
  temp_list_BPSV[[i]] <- check_BPSV
  print(i)
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

ggplot(data = overall) +
  geom_point(aes(x = date, y = reported_deaths), col = "grey") +
  geom_ribbon(aes(x = date, ymin = `025`, ymax = `975`, fill = scenario), alpha = 0.2) +
  geom_line(aes(x = date, y = med, col = scenario)) +
  theme_bw() +
  labs(x = "Date", y = "Deaths Per Day")


## Scrap Code
# seed_infections(excess$squire_model, excess$parameters$country, out$samples[[1]]$initial_infections)$E1_0
# x <- out$output[1, 138:154, 1]

plot(rep_summary$value_med[1:365], type = "l")
lines(rep_bpsv$value_025[1:365], col = "blue")
lines(rep_bpsv$value_975[1:365], col = "blue")

# BPSV, fit and my replication for draw 1
check <- temp_list[[1]]
fit <- daily %>%
  filter(replicate == 1)
plot(fit$date[1:365], fit$deaths[1:365], type = "l", col = "red")
lines(fit$date[1:365], check$value[1:365], type = "l")
lines(fit$date[1:365], check_BPSV$value[1:365], col = "purple")

sum(fit$deaths[1:365])
sum(check$value[1:365])
sum(check_BPSV$value[1:365])
