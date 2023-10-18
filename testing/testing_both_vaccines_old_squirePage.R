# devtools::install_github("mrc-ide/squire.page")

# Load required libraries
source(here::here("main.R"))
set.seed(123)

# Load required functions
source(here::here("functions/run_sars_x.R"))

# Running with old squire.page
test_scenarios <- readRDS("test_scenarios.rds")
test_both <- test_scenarios %>%
  filter(vaccine_scenario == "both_vaccines")
test_both$population_size <- 10^9

standard_pop <- generate_standard_pop(country = test_both$country, population_size = test_both$population_size)
TgVary <- scale_generation_time(target_Tg = test_both$Tg)
TgVary_dur_IMild <- TgVary$dur_IMild
TgVary_dur_ICase <- TgVary$dur_ICase
prob_hosp <- scale_IFR(country = test_both$country, population_size = test_both$population_size, target_IFR = test_both$IFR)

priority_age_groups <- test_both$min_age_group_index_priority:17
vaccination_age_groups <- test_both$min_age_group_index_non_priority:17
vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(test_both$coverage, length(priority_age_groups)), 
                                 rep(0, 17 - length(vaccination_age_groups)),  rep(test_both$coverage, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)


daily_doses <- test_both$vaccination_rate * test_both$population_size / 7
elderly_pop_to_vaccinate <- sum(standard_pop[priority_age_groups]) * test_both$coverage 
time_to_coverage_bpsv <- ceiling(elderly_pop_to_vaccinate/daily_doses)
time_to_coverage_spec <- time_to_coverage_bpsv 

primary_doses <- c(0, daily_doses, 0, daily_doses)
tt_primary_doses <- c(0, 
                      test_both$detection_time + test_both$bpsv_start + test_both$bpsv_protection_delay,
                      test_both$detection_time + test_both$bpsv_start + test_both$bpsv_protection_delay + test_both$time_to_coverage_bpsv,
                      test_both$detection_time + test_both$specific_vaccine_start + test_both$time_to_coverage_spec + test_both$specific_protection_delay)

booster_doses <- c(0, daily_doses)
tt_booster_doses <- c(0,  test_both$detection_time + test_both$specific_vaccine_start + test_both$specific_protection_delay)

ve_bpsv <- list(infection = 0.5, disease = 0.5)
ve_spec <- list(infection = 1, disease = 1)
ve_i_elderly <- c(ve_bpsv$infection, ve_bpsv$infection, ve_bpsv$infection, 0, ve_spec$infection, ve_spec$infection, 0) 
ve_i_non_elderly <- c(ve_spec$infection, ve_spec$infection, ve_spec$infection, 0, ve_spec$infection, ve_spec$infection, 0)
vaccine_efficacy_infection <- matrix(c( rep(ve_i_non_elderly, 17 - length(priority_age_groups)), 
                                        rep(ve_i_elderly, length(priority_age_groups)) ), ncol = 17)

ve_d_elderly <- c(ve_bpsv$disease, ve_bpsv$disease, ve_bpsv$disease, 0, ve_spec$disease, ve_spec$disease, 0) 
ve_d_non_elderly <- c(ve_spec$disease, ve_spec$disease, ve_spec$disease, 0, ve_spec$disease, ve_spec$disease, 0)
vaccine_efficacy_disease <- matrix(c( rep(ve_d_non_elderly, 17 - length(priority_age_groups)), 
                                      rep(ve_d_elderly, length(priority_age_groups)) ), ncol = 17)

vaccine_booster_follow_up_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(0, length(priority_age_groups))) # no booster follow ons
vaccine_booster_initial_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(1, length(priority_age_groups)))

dur_V <- unlist(test_both$dur_V)
dur_R <- test_both$dur_R
tt_Rt <- unlist(test_both$tt_Rt)
Rt <- unlist(test_both$Rt)

mm <- squire::get_mixing_matrix(country = test_both$country)    
rel_infectiousness_vaccinated <- squire.page:::probs_booster$rel_infectiousness_vaccinated
rel_infectiousness_vaccinated[rel_infectiousness_vaccinated < 1] <- 1

runtime <- test_both$runtime
mod_run <- run_booster(time_period = runtime,
                       population = standard_pop,                                                 
                       contact_matrix_set = mm,                                                   
                       R0 = Rt,     
                       tt_R0 = tt_Rt, 
                       hosp_bed_capacity = test_both$hosp_bed_capacity,                                     
                       ICU_bed_capacity = test_both$ICU_bed_capacity,                                       
                       prob_hosp = prob_hosp,
                       dur_IMild = TgVary_dur_IMild,
                       dur_ICase = TgVary_dur_ICase,
                       dur_R = dur_R,                                                        
                       seeding_cases = test_both$seeding_cases,
                       vaccine_coverage_mat = vaccine_coverage_mat,                               
                       vaccine_efficacy_infection = vaccine_efficacy_infection,                   
                       vaccine_efficacy_disease = vaccine_efficacy_disease,     
                       rel_infectiousness_vaccinated = rel_infectiousness_vaccinated, 
                       dur_V = dur_V,             
                       second_dose_delay = 1,
                       protection_delay_rate = 10,
                       primary_doses = primary_doses,  
                       tt_primary_doses = tt_primary_doses,
                       booster_doses = booster_doses,
                       tt_booster_doses = tt_booster_doses,                                             
                       vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
                       vaccine_booster_initial_coverage = vaccine_booster_initial_coverage)

calc_summary_metrics(mod_run)

check_both <- nimue::format(mod_run, compartments = c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_booster_dose"),
                            reduce_age = FALSE) %>%
  filter(t > 1, compartment == "deaths" |  
           compartment == "vaccinated_first_dose" | compartment == "vaccinated_second_dose" | compartment == "vaccinated_booster_dose") %>%
  group_by(replicate, t) 

ggplot() +
  geom_line(data = subset(check_both, age_group %in% c("80+") &
                            compartment == "vaccinated_booster_dose"),
            aes(x = t, y = value, col = compartment)) +
  geom_line(data = subset(check_both, age_group %in% c("80+") &
                            compartment == "vaccinated_first_dose"),
            aes(x = t, y = value, col = compartment)) +
  geom_line(data = subset(check_both, age_group %in% c("80+") &
                            compartment == "vaccinated_second_dose"),
            aes(x = t, y = value, col = compartment)) + ## still need to sort second dose delivery and making sure coveerage is achieved in same time frame as primary doses
  facet_wrap(~age_group) +
  ylim(c(0, 22500))


# ggplot() +
#   geom_line(data = subset(check_both, compartment == "vaccinated_booster_dose"),
#             aes(x = t, y = value, col = compartment)) +
#   geom_line(data = subset(check_both, compartment == "vaccinated_first_dose"),
#             aes(x = t, y = value, col = compartment)) +
#   geom_line(data = subset(check_both, compartment == "vaccinated_second_dose"),
#             aes(x = t, y = value, col = compartment)) + ## still need to sort second dose delivery and making sure coveerage is achieved in same time frame as primary doses
#   facet_wrap(~age_group) 
