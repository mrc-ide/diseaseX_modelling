## testing out the branching process and comparing it to squire.page output


## Running squire.page
population <- squire:::get_population("Argentina")
population_low <- 10^6 * population$n / sum(population$n)
population_high <- 10^9 * population$n / sum(population$n)
mixing_matrix <- squire:::get_mixing_matrix("Argentina")
TgVary <- scale_generation_time(target_Tg = 6.7)
TgVary_dur_IMild <- TgVary$dur_IMild
TgVary_dur_ICase <- TgVary$dur_ICase

mod_run_lowPop <- run_booster(time_period = 365,
                              population = population_low,                                                 
                              contact_matrix_set = mixing_matrix,                                                   
                              R0 = 3.5,     
                              tt_R0 = 0, 
                              hosp_bed_capacity = 10^9,                                     
                              ICU_bed_capacity = 10^9,                                       
                              dur_IMild = TgVary_dur_IMild,
                              dur_ICase = TgVary_dur_ICase,
                              dur_R = 3650000000,                                                        
                              seeding_cases = 1,
                              rel_infectiousness_vaccinated = 1, 
                              dur_V = 3650000000,                                              
                              primary_doses = rep(0, 365),  
                              second_doses = rep(0, 365),
                              booster_doses = rep(0, 365))
mod_run_highPop <- run_booster(time_period = 365,
                               population = population_high,                                                 
                               contact_matrix_set = mixing_matrix,                                                   
                               R0 = 3.5,     
                               tt_R0 = 0, 
                               hosp_bed_capacity = 10^9,                                     
                               ICU_bed_capacity = 10^9,                                       
                               dur_IMild = TgVary_dur_IMild,
                               dur_ICase = TgVary_dur_ICase,
                               dur_R = 3650000000,                                                        
                               seeding_cases = 1,
                               rel_infectiousness_vaccinated = 1, 
                               dur_V = 3650000000,                                              
                               primary_doses = rep(0, 365),  
                               second_doses = rep(0, 365),
                               booster_doses = rep(0, 365))


low_pop_infections <- nimue::format(mod_run_lowPop, compartments = "S", summaries = "hospitalisations") %>%
  filter(compartment != "S") %>%
  select(-replicate) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(pop = "low_pop")

high_pop_infections <- nimue::format(mod_run_highPop, compartments = "S", summaries = "hospitalisations") %>%
  filter(compartment != "S") %>%
  select(-replicate) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(pop = "high_pop")

overall_df <- rbind(low_pop_infections, high_pop_infections)

ggplot(overall_df, aes(x = t, y = value, col = pop)) +
  geom_line() +
  coord_cartesian(xlim = c(0, 60), ylim = c(0, 100))

thresholds <- c(1, 5, 10, 20)
raw_detection_times <- overall_df %>%
  group_by(pop) %>%
  nest() %>%
  mutate(detection_times = map(data, ~ { ## nest creates a column of lists of data
    sapply(thresholds, function(threshold) {
      first(.x$t[.x$value >= threshold])
      # .x$t[which.min(abs(.x$value - threshold))]
    })
  })) %>%
  unnest_longer(col = detection_times, indices_to = "detection_threshold_index") %>%
  mutate(detection_threshold = thresholds[detection_threshold_index]) %>%
  select(-data)

detection_df %>%
  filter(R0 == 3.5, NPI_int == 4, specific_vaccine_start == 100, metric == "Daily Incidence") %>%
  select(R0, detection_time, detection_threshold_hosp) 
  

