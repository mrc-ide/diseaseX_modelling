# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/helper_functions.R"))

## Running the model
run <- run_sars_x(population_size = 1e10, 
                  country = "Argentina",
                  hosp_bed_capacity = 1e10, 
                  ICU_bed_capacity = 1e10,
                  Rt = 2.5, 
                  tt_Rt = 0, 
                  Tg = 6.6, 
                  IFR = 1,
                  vaccine_scenario = "specific_only",       # which scenario to explore
                  detection_time = 14,                      # time at which the pathogen is detected
                  bpsv_start = 14,                          # BPSV distribution start (days after detection time)
                  bpsv_protection_delay = 7,                # time between BPSV dose and protection arising
                  specific_vaccine_start = 300,             # specific vaccine distribution start (days after detection time)
                  specific_protection_delay = 7,            # time between specific dose and protection arising
                  efficacy_infection_bpsv = 0.35,           # vaccine efficacy against infection - BPSV
                  efficacy_disease_bpsv = 0.8,              # vaccine efficacy against disease - BPSV
                  efficacy_infection_spec = 0.55,           # vaccine efficacy against infection - specific vaccine
                  efficacy_disease_spec = 0.9,              # vaccine efficacy against disease - specific vaccine
                  dur_R = 1000 * 365,                       # duration of infection-induced immunity
                  dur_bpsv = 1000 * 365,                    # duration of vaccine-induced immunity for BPSV vaccine
                  dur_spec = 1000 * 365,                    # duration of vaccine-induced immunity for disease-specific vaccines
                  coverage_bpsv = 0.49,                      # proportion of the population vaccinated
                  coverage_spec = 0.5,                     # proportion of the population vaccinated
                  vaccination_rate_bpsv = 0.035,            # bpsv vaccination rate per week as percentage of population (note: percentage of whole population, despite restricted age-groups receiving it)
                  vaccination_rate_spec = 0.015,            # disease-specific vaccination rate per week as percentage of population
                  min_age_group_index_priority = 13,        # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                  min_age_group_index_non_priority = 4,     # index of the youngest age group that *receives* vaccines (4 = 15+)
                  runtime = 600,
                  seeding_cases = 5,
                  output = "full",
                  NPI_int = 0,
                  scenario_index = 0,
                  varied = "") 

check_both <- nimue::format(run$model_output, compartments = c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_booster_dose"),
                            reduce_age = FALSE) %>%
  filter(t > 1, compartment == "deaths" |  
           compartment == "vaccinated_first_dose" | compartment == "vaccinated_second_dose" | compartment == "vaccinated_booster_dose") %>%
  group_by(replicate, t) 

standard_pop <- generate_standard_pop(country = "Argentina", population_size = 1e10)
eighty_plus <- standard_pop[length(standard_pop)]
ggplot() +
  geom_line(data = subset(check_both, age_group %in% c("80+") &
                            compartment == "vaccinated_booster_dose"),
            aes(x = t, y = value / eighty_plus, col = compartment)) +
  geom_line(data = subset(check_both, age_group %in% c("80+") &
                            compartment == "vaccinated_first_dose"),
            aes(x = t, y = value / eighty_plus, col = compartment)) +
  geom_line(data = subset(check_both, age_group %in% c("80+") &
                            compartment == "vaccinated_second_dose"),
            aes(x = t, y = value / eighty_plus, col = compartment)) + ## still need to sort second dose delivery and making sure coveerage is achieved in same time frame as primary doses
  facet_wrap(~age_group) +
  lims(y = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1))


ggplot() +
  geom_line(data = subset(check_both, compartment == "vaccinated_booster_dose"),
            aes(x = t, y = value, col = compartment)) +
  geom_line(data = subset(check_both, compartment == "vaccinated_first_dose"),
            aes(x = t, y = value, col = compartment)) +
  geom_line(data = subset(check_both, compartment == "vaccinated_second_dose"),
            aes(x = t, y = value, col = compartment)) + ## still need to sort second dose delivery and making sure coveerage is achieved in same time frame as primary doses
  facet_wrap(~age_group) 
