# Load required libraries
library(tidyverse); library(squire.page)

# Sourcing required functions
source("Q_Drive_Copy2/Active_Research_Projects/diseaseX_modelling/functions/extract_process_generate_fits.R")

# Getting fits for Iran
excess <- grab_fit("IRN", TRUE, TRUE)

excess_Rt <- simple_Rt(excess)
excess_Rt_comb <- bind_rows(excess_Rt) %>%
  group_by(date, t) %>%
  summarise(mean_Rt = mean(Rt),
            lower_Rt = quantile(Rt, 0.025),
            upper_Rt = quantile(Rt, 0.975))
plot(excess_Rt_comb$date, excess_Rt_comb$mean_Rt, type = "l")
polygon(c(excess_Rt_comb$date, rev(excess_Rt_comb$date)),
        c(excess_Rt_comb$lower_Rt, rev(excess_Rt_comb$upper_Rt)),
        col = adjustcolor("red", alpha.f = 0.2), border = NA)

# Generating death curves for each of the Rt profiles
out <- squire.page:::generate_draws.rt_optimised(excess)
deaths <- get_deaths_infections_hosps_time(out)

## Generating cumulative and daily deaths
cumulative <- deaths %>%
  group_by(replicate) %>%
  arrange(date) %>%
  mutate(across(c(deaths, infections), ~cumsum(.x))) %>%
  group_by(date) %>%
  summarise(
    across(c(deaths, infections), ~median(.x, na.rm=TRUE), .names = "{col}_med"),
    across(c(deaths, infections), ~quantile(.x, 0.025, na.rm=TRUE), .names = "{col}_025"),
    across(c(deaths, infections), ~quantile(.x, 0.975, na.rm=TRUE), .names = "{col}_975"),
    .groups = "drop")

daily <- deaths %>%
  group_by(replicate) %>%
  arrange(date) %>%
  group_by(date) %>%
  summarise(
    across(c(deaths, infections), ~median(.x, na.rm=TRUE), .names = "{col}_med"),
    across(c(deaths, infections), ~quantile(.x, 0.025, na.rm=TRUE), .names = "{col}_025"),
    across(c(deaths, infections), ~quantile(.x, 0.975, na.rm=TRUE), .names = "{col}_975"),
    .groups = "drop")

plot(daily$date, daily$deaths_med, type = "l")
points(out$inputs$data$date_start + 3, out$inputs$data$deaths / 7)

overall_Rt <- excess_Rt_comb$mean_Rt
overall_tt_Rt <- excess_Rt_comb$t
seeding_cases <- bind_rows(lapply(out$samples, function(x) {
  x$initial_infections
}))
seeding_cases <- median(seeding_cases$initial_infections)
bpsv <- squire.page::run_booster(time_period = max(excess_Rt[[1]]$t),
                                 population = squire::get_population(country = out$parameters$country)$n,                                                 
                                 contact_matrix_set = squire::get_mixing_matrix(country = out$parameters$country),                                                   
                                 R0 = overall_Rt,     
                                 tt_R0 = overall_tt_Rt, 
                                 hosp_bed_capacity =  out$parameters$hosp_bed_capacity,                                     
                                 ICU_bed_capacity = out$parameters$ICU_bed_capacity,                                       
                                 prob_hosp = squire.page:::probs_booster$prob_hosp,
                                 dur_IMild = squire.page:::durs_booster$dur_IMild,
                                 dur_ICase = squire.page:::durs_booster$dur_ICase,
                                 dur_R = squire.page:::durs_booster$dur_R,                                                        
                                 seeding_cases = seeding_cases)


check <- nimue::format(bpsv, compartments = "D", summaries = "deaths") %>%
  filter(t > 1, compartment == "deaths")

check2 <- bpsv$output[, colnames(check2)[grep("E1\\[\\d*,1\\]", colnames(check2))], 1]

check2[1, ]

sum(check$value[1:365])
sum(daily$deaths_med[1:365])


## Running the model with Rt profile AND BPSV
runtime <- max(excess_Rt[[1]]$t)
specific_vaccine_start <- runtime - 50
population <- squire::get_population(country = out$parameters$country)$n                   
mm <- squire::get_mixing_matrix(country = out$parameters$country)
overall_Rt <- excess_Rt_comb$mean_Rt
overall_tt_Rt <- excess_Rt_comb$t
hosp_bed_capacity <- out$parameters$hosp_bed_capacity
ICU_bed_capacity <- out$parameters$ICU_bed_capacity
prob_hosp <- squire.page:::probs_booster$prob_hosp
dur_IMild <- squire.page:::durs_booster$dur_IMild
dur_ICase <- squire.page:::durs_booster$dur_ICase 
dur_R <- squire.page:::durs_booster$dur_R        
seeding_cases <- bind_rows(lapply(out$samples, function(x) {
  x$initial_infections
}))
seeding_cases <- median(seeding_cases$initial_infections)

# overall_Rt <- excess_Rt[[1]]$Rt
# overall_tt_Rt <- excess_Rt[[1]]$t
# seeding_cases <- out$samples[[1]]$initial_infections

vaccine_params <- create_scenarios(R0 = 3, specific_vaccine_start = specific_vaccine_start)
priority_age_groups <- vaccine_params$min_age_group_index_priority[1]:17  
vaccination_age_groups <- vaccine_params$min_age_group_index_non_priority[1]:17
coverage_spec <- 0.01
vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(coverage_spec, length(priority_age_groups)), 
                                 rep(0, 17 - length(vaccination_age_groups)), rep(coverage_spec, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)
vaccine_efficacy_infection <- 0
vaccine_efficacy_disease <- 0
ve_bpsv <- list(infection = vaccine_efficacy_infection, disease = vaccine_efficacy_disease)
ve_spec <- list(infection = 0, disease = 0)
ve_i_elderly_bpsv_campaign <- c(ve_bpsv$infection, ve_bpsv$infection, ve_bpsv$infection, 0, ve_spec$infection, ve_spec$infection, 0) 
ve_i_elderly_spec_campaign <- c(ve_spec$infection, ve_bpsv$infection, ve_bpsv$infection, 0, ve_spec$infection, ve_spec$infection, 0) 
ve_i_non_elderly <- c(ve_spec$infection, ve_spec$infection, ve_spec$infection, 0, ve_spec$infection, ve_spec$infection, 0)
vaccine_efficacy_infection_bpsv_campaign <- matrix(c(rep(ve_i_non_elderly, 17 - length(priority_age_groups)), 
                                                     rep(ve_i_elderly_bpsv_campaign, length(priority_age_groups)) ), ncol = 17)
vaccine_efficacy_infection_spec_campaign <- matrix(c(rep(ve_i_non_elderly, 17 - length(priority_age_groups)), 
                                                     rep(ve_i_elderly_spec_campaign, length(priority_age_groups)) ), ncol = 17)

ve_d_elderly_bpsv_campaign <- c(ve_bpsv$disease, ve_bpsv$disease, ve_bpsv$disease, 0, ve_spec$disease, ve_spec$disease, 0) 
ve_d_elderly_spec_campaign <- c(ve_spec$disease, ve_bpsv$disease, ve_bpsv$disease, 0, ve_spec$disease, ve_spec$disease, 0) 
ve_d_non_elderly <- c(ve_spec$disease, ve_spec$disease, ve_spec$disease, 0, ve_spec$disease, ve_spec$disease, 0)
vaccine_efficacy_disease_bpsv_campaign <- matrix(c(rep(ve_d_non_elderly, 17 - length(priority_age_groups)), 
                                                   rep(ve_d_elderly_bpsv_campaign, length(priority_age_groups)) ), ncol = 17)
vaccine_efficacy_disease_spec_campaign <- matrix(c(rep(ve_d_non_elderly, 17 - length(priority_age_groups)), 
                                                   rep(ve_d_elderly_spec_campaign, length(priority_age_groups)) ), ncol = 17)

rel_infectiousness_vaccinated <- squire.page.sarsX:::probs_booster$rel_infectiousness_vaccinated
rel_infectiousness_vaccinated[rel_infectiousness_vaccinated < 1] <- 1
dur_V <- 3650000
vaccine_booster_follow_up_coverage <- rep(0, 17)
vaccine_booster_initial_coverage <- rep(0, 17)

vaccine_doses <- create_vaccination_dose_series(country = out$parameters$country, 
                                                population_size = sum(population), 
                                                detection_time = 0, 
                                                vaccine_scenario = "both_vaccines",
                                                bpsv_start = 500, 
                                                bpsv_protection_delay = 7,
                                                specific_vaccine_start = specific_vaccine_start,
                                                specific_protection_delay = 10,
                                                vaccination_rate_bpsv = 0.025,
                                                vaccination_rate_spec = 0.025,
                                                coverage_bpsv = 0.8,
                                                coverage_spec = 0.8,
                                                min_age_group_index_priority = vaccine_params$min_age_group_index_priority[1],
                                                runtime = runtime)
primary_doses <- vaccine_doses$primary_doses
second_doses <- vaccine_doses$second_doses
booster_doses <- vaccine_doses$booster_doses

primary_doses[primary_doses > 0] <- 0 
second_doses[second_doses > 0] <- 0 
booster_doses[booster_doses > 0] <- 0 

bpsv <- squire.page.sarsX::run_booster(time_period = runtime,
                                       population = population,                                                 
                                       contact_matrix_set = mm,                                                   
                                       R0 = overall_Rt,     
                                       tt_R0 = overall_tt_Rt, 
                                       hosp_bed_capacity = hosp_bed_capacity,                                     
                                       ICU_bed_capacity = ICU_bed_capacity,                                       
                                       prob_hosp = prob_hosp,
                                       dur_IMild = dur_IMild,
                                       dur_ICase = dur_ICase,
                                       dur_R = dur_R,                                                        
                                       seeding_cases = seeding_cases,
                                       vaccine_coverage_mat = vaccine_coverage_mat,                               
                                       vaccine_efficacy_infection = list(vaccine_efficacy_infection_bpsv_campaign, vaccine_efficacy_infection_spec_campaign),
                                       tt_vaccine_efficacy_infection = c(0, specific_vaccine_start - 5), 
                                       vaccine_efficacy_disease = list(vaccine_efficacy_disease_bpsv_campaign, vaccine_efficacy_disease_spec_campaign),
                                       tt_vaccine_efficacy_disease = c(0, specific_vaccine_start - 5), 
                                       rel_infectiousness_vaccinated = rel_infectiousness_vaccinated, 
                                       dur_V = dur_V,                                              
                                       primary_doses = primary_doses,  
                                       second_doses = second_doses,
                                       booster_doses = booster_doses,                                             
                                       vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
                                       vaccine_booster_initial_coverage = vaccine_booster_initial_coverage)

bpsv <- squire.page::run_booster(time_period = runtime,
                                 population = population,                                                 
                                 contact_matrix_set = mm,                                                   
                                 R0 = overall_Rt,     
                                 tt_R0 = overall_tt_Rt, 
                                 hosp_bed_capacity = hosp_bed_capacity,                                     
                                 ICU_bed_capacity = ICU_bed_capacity,                                       
                                 prob_hosp = prob_hosp,
                                 dur_IMild = dur_IMild,
                                 dur_ICase = dur_ICase,
                                 dur_R = dur_R,                                                        
                                 seeding_cases = seeding_cases)


check <- nimue::format(bpsv, compartments = "D", summaries = "deaths") %>%
  filter(t > 1, compartment == "deaths")

check2 <- bpsv$output[, colnames(check2)[grep("E1\\[\\d*,1\\]", colnames(check2))], 1]

check2[1, ]

sum(check$value[1:365])
sum(daily$deaths_med[1:365])

# deaths_one <- deaths %>%
#   filter(replicate == 1)

plot(2:366, daily$deaths_med[1:365], type = "l", col = "black")
lines(check$t[1:365], check$value[1:365], col = "red")

plot(2:101, daily$deaths_med[1:100], type = "l", col = "black")
lines(check$t[1:100], check$value[1:100], col = "red")


