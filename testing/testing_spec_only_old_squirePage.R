# devtools::install_github("mrc-ide/squire.page")

# Load required libraries
source(here::here("main.R"))
set.seed(123)

# Load required functions
source(here::here("functions/run_sars_x.R"))

# Running with old squire.page
test_scenarios <- readRDS("test_scenarios.rds")
test_spec <- test_scenarios %>%
  filter(vaccine_scenario == "specific_only")
test_spec$population_size <- 10^6

standard_pop <- generate_standard_pop(country = test_spec$country, population_size = test_spec$population_size)
TgVary <- scale_generation_time(target_Tg = test_spec$Tg)
TgVary_dur_IMild <- TgVary$dur_IMild
TgVary_dur_ICase <- TgVary$dur_ICase
prob_hosp <- scale_IFR(country = test_spec$country, population_size = test_spec$population_size, target_IFR = test_spec$IFR)

priority_age_groups <- test_spec$min_age_group_index_priority:17
vaccination_age_groups <- test_spec$min_age_group_index_non_priority:17
vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(test_spec$coverage, length(priority_age_groups)), 
                                 rep(0, 17 - length(vaccination_age_groups)),  rep(test_spec$coverage, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)

daily_doses <- test_spec$vaccination_rate * test_spec$population_size / 7
primary_doses <- c(0, daily_doses)
tt_primary_doses <- c(0, test_spec$detection_time + test_spec$specific_vaccine_start + test_spec$specific_protection_delay)

booster_doses <- c(0)
tt_booster_doses <- c(0)

ve_i_spec <- c(1, 1, 1, 0, 0, 0, 0)
vaccine_efficacy_infection <- matrix(c(rep(ve_i_spec, 17)), ncol = 17)

ve_d_spec <- c(1, 1, 1, 0, 0, 0, 0)
vaccine_efficacy_disease <- matrix(c(rep(ve_d_spec, 17)), ncol = 17)

vaccine_booster_follow_up_coverage <- rep(0, 17)
vaccine_booster_initial_coverage <- rep(0, 17)

dur_V <- unlist(test_spec$dur_V)
dur_R <- test_spec$dur_R
tt_Rt <- unlist(test_spec$tt_Rt)
Rt <- unlist(test_spec$Rt)

mm <- squire::get_mixing_matrix(country = test_spec$country)    
rel_infectiousness_vaccinated <- squire.page:::probs_booster$rel_infectiousness_vaccinated
rel_infectiousness_vaccinated[rel_infectiousness_vaccinated < 1] <- 1

runtime <- test_spec$runtime
mod_run <- run_booster(time_period = runtime,
                       population = standard_pop,                                                 
                       contact_matrix_set = mm,                                                   
                       R0 = Rt,     
                       tt_R0 = tt_Rt, 
                       hosp_bed_capacity = test_spec$hosp_bed_capacity,                                     
                       ICU_bed_capacity = test_spec$ICU_bed_capacity,                                       
                       prob_hosp = prob_hosp,
                       dur_IMild = TgVary_dur_IMild,
                       dur_ICase = TgVary_dur_ICase,
                       dur_R = dur_R,                                                        
                       seeding_cases = test_spec$seeding_cases,
                       vaccine_coverage_mat = vaccine_coverage_mat,                               
                       vaccine_efficacy_infection = vaccine_efficacy_infection,                   
                       vaccine_efficacy_disease = vaccine_efficacy_disease,     
                       rel_infectiousness_vaccinated = rel_infectiousness_vaccinated, 
                       dur_V = dur_V,             
                       second_dose_delay = 1,
                       protection_delay_rate = 10,
                       primary_doses = primary_doses,  
                       tt_primary_doses = tt_primary_doses,
                       vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
                       vaccine_booster_initial_coverage = vaccine_booster_initial_coverage)

calc_summary_metrics(mod_run)

check_spec <- nimue::format(mod_run, compartments = c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_booster_dose"),
                            reduce_age = FALSE) %>%
  filter(t > 1, compartment == "deaths" |  
           compartment == "vaccinated_first_dose" | compartment == "vaccinated_second_dose" | compartment == "vaccinated_booster_dose") %>%
  group_by(replicate, t) 

ggplot() +
  geom_line(data = subset(check_spec, age_group %in% c("80+") &
                            compartment == "vaccinated_booster_dose"),
            aes(x = t, y = value, col = compartment)) +
  geom_line(data = subset(check_spec, age_group %in% c("80+") &
                            compartment == "vaccinated_first_dose"),
            aes(x = t, y = value, col = compartment)) +
  geom_line(data = subset(check_spec, age_group %in% c("80+") &
                            compartment == "vaccinated_second_dose"),
            aes(x = t, y = value, col = compartment)) + ## still need to sort second dose delivery and making sure coveerage is achieved in same time frame as primary doses
  facet_wrap(~age_group)  +
  ylim(c(0, 22500))

deaths_spec <- check_spec %>%
  filter(compartment == "deaths") %>%
  group_by(age_group) %>%
  summarise(total = sum(value))

deaths_both <- check_both %>%
  filter(compartment == "deaths") %>%
  group_by(age_group) %>%
  summarise(total = sum(value))

sum(deaths_spec$total)
sum(deaths_both$total)

plot(deaths_spec$total)
lines(deaths_both$total)


# ggplot() +
#   geom_line(data = subset(check_both, compartment == "vaccinated_booster_dose"),
#             aes(x = t, y = value, col = compartment)) +
#   geom_line(data = subset(check_both, compartment == "vaccinated_first_dose"),
#             aes(x = t, y = value, col = compartment)) +
#   geom_line(data = subset(check_both, compartment == "vaccinated_second_dose"),
#             aes(x = t, y = value, col = compartment)) + ## still need to sort second dose delivery and making sure coveerage is achieved in same time frame as primary doses
#   facet_wrap(~age_group) 
