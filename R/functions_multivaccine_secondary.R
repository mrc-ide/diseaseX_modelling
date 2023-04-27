
# Specify the country chosen to represent each income group
get_representative_country <- function(income_group){
  case_when(income_group == "HIC" ~ "Canada",
            income_group == "UMIC" ~ "Argentina",
            income_group == "LMIC" ~ "Pakistan",
            income_group == "LIC" ~ "Chad")
  
}



#use an age-distribution defined by probability of death from severe disease if hospitalised
get_ifr_dist <- function(){
  
  prob_severe_death_treatment <- c(0.0038,0.0038,0.0038,0.0038,0.0076,0.0114,0.0152,0.0228,0.0304,0.038,0.076,0.152,0.228,0.304,0.38,0.57,0.95)  
  return(prob_severe_death_treatment)
}

# Run scenarios
run_scenario_2 <- function(target_pop = 1e6,
                           income_group = "HIC",
                           R0 = 3,
                           Rt1 = 0.8,
                           Rt1a = 0.8,
                           Rt1b = 0.8,
                           Rt2 = 3,
                           timing1 = 30,
                           timing1a = 40,
                           timing1b = 70,
                           timing2 = 365,
                           ifr_scaling = 1,
                           coverage = 0, 
                           efficacy_infection_v1 = 0,
                           efficacy_disease_v1 = 0,
                           efficacy_infection_v2 = 0,
                           efficacy_disease_v2 = 0,
                           vaccination_rate = 0.02,
                           duration_R = 1000, # duration of infection-induced immunity
                           duration_V = 5000, # duration of vaccine-induced immunity
                           dur_vacc_delay = 14, # mean duration from vaccination to protection
                           seeding_cases = 5, # define as the number of cases at first sequencing - will need to explore
                           vaccine_1_start = 20, # using stockpiled vaccines - will depend on seeding_cases
                           vaccine_2_start = 100, # 100 day mission for new vaccine
                           lower_priority = 14, #elderly for stockpile of v1
                           lower_vaccine = 4, #15+ for v2
                           two_vaccines = 1, # whether distributing two vaccines or single SARS-3 vaccine
                           runtime = 365*3
){
  
  # Population, mixing and lifespan
  rep_country <- get_representative_country(income_group = income_group)
  pop <- squire::get_population(country = rep_country)$n
  pop_standardise <- target_pop / sum(pop)
  pop <- pop * pop_standardise
  mm <- squire::get_mixing_matrix(country = rep_country)
  life_expectancy <- read.csv("data/life_expectancy_discounted.csv") 
  le_vector <- life_expectancy %>% 
    select(age_group,{{income_group}}) %>%
    rename(life_expect = {{income_group}})
  
  
  
  
  # IFR scaling
  
  prob_severe_death_treatment <- get_ifr_dist()
  
  # scalings as generated in spreadsheet to obtain IFR under uniform attack rate in UK population
  prop_hosp_severe <- case_when(ifr_scaling == 1 ~ 0.5,
                                ifr_scaling == 0.5 ~ 0.1,
                                ifr_scaling == 0.1 ~ 0.1,
                                ifr_scaling == 2 ~ 0.8)
  
  ratio_deaths_severe_to_hosp <- case_when(ifr_scaling == 1 ~ 0.29,
                                           ifr_scaling == 0.5 ~ 0.29,
                                           ifr_scaling == 0.1 ~ 0.06,
                                           ifr_scaling == 2 ~ 0.39)
  
  ratio_deaths_non_severe_to_deaths_severe <- case_when(ifr_scaling == 1 ~ 0.4,
                                                        ifr_scaling == 0.5 ~ 0.215,
                                                        ifr_scaling == 0.1 ~ 0.215,
                                                        ifr_scaling == 2 ~ 0.6)
  
  prob_hosp <- ratio_deaths_severe_to_hosp*prob_severe_death_treatment
  prob_severe <- prop_hosp_severe*prob_hosp
  prob_non_severe_death_treatment <- ratio_deaths_non_severe_to_deaths_severe*prob_severe_death_treatment
  
  
  #vaccine 1 stockpiled and immediately rollout out to elderly populations 
  #vaccine 2 produced after a year with higher efficacy and rolled out first to the elderly then to the rest of the pop
  
  # vaccine efficacy
  ve_v1 <- list(infection = efficacy_infection_v1, disease = efficacy_disease_v1)
  ve_v2 <- list(infection = efficacy_infection_v2, disease = efficacy_disease_v2)
  
  # daily doses of each vaccine - can at a later stage introduce restrictions on stockpile
  v1_daily_doses <- vaccination_rate*sum(pop)/7
  v2_daily_doses <- vaccination_rate*sum(pop)/7
  
  priority_age_groups <- lower_priority:17
  vaccination_age_groups <- lower_vaccine:17
  
  
  if(two_vaccines==1)
  {
    # calculate how long it will take to reach v1 coverage in the target_age_groups
    elderly_pop_to_vaccinate <- sum(pop[priority_age_groups])*coverage
    time_to_coverage_v1 <- elderly_pop_to_vaccinate/v1_daily_doses
    
    #calculate how long it will take to boost the elderly pop (who are vaccinated) 
    time_to_coverage_v2_elderly <- ceiling(elderly_pop_to_vaccinate/v2_daily_doses)
    
    #generate vaccine coverage matrix - rows are priorities for oldest, followed by younger
    vaccine_coverage_mat <- matrix(c( 
      rep(0, 17 - length(priority_age_groups)), rep(coverage, length(priority_age_groups)), 
      rep(0, 17 - length(vaccination_age_groups)), rep(coverage, length(vaccination_age_groups)) ), 
      ncol = 17, byrow = TRUE)
    
    #format for ve in the booster model (first dose, second dose, second dose, waned, booster, booster , waned) 
    
    ve_i_booster_elderly <- c(ve_v1$infection, ve_v1$infection, ve_v1$infection, 0, ve_v2$infection, ve_v2$infection, 0) 
    ve_i_booster_non_elderly <- c(ve_v2$infection, ve_v2$infection, ve_v2$infection, 0, ve_v2$infection, ve_v2$infection, 0)
    
    #setup as a matrix giving VE per age group
    vaccine_efficacy_infection <- matrix(c( rep(ve_i_booster_non_elderly, 17 - length(priority_age_groups)), 
                                            rep(ve_i_booster_elderly, length(priority_age_groups)) ), ncol = 17)
    
    
    
    
    
    #repeat for disease 
    ve_i_booster_elderly <- c(ve_v1$disease, ve_v1$disease, ve_v1$disease, 0, ve_v2$disease, ve_v2$disease, 0) 
    ve_i_booster_non_elderly <- c(ve_v2$disease, ve_v2$disease, ve_v2$disease, 0, ve_v2$disease, ve_v2$disease, 0)
    
    vaccine_efficacy_disease <- matrix(c( rep(ve_i_booster_non_elderly, 17 - length(priority_age_groups)), 
                                          rep(ve_i_booster_elderly, length(priority_age_groups)) ), ncol = 17)
    
    # set up dosing schedule - here the elderly receive v2 ahead of the rest of the population once it becomes available
    
    primary_doses <- c(0, v1_daily_doses, 0, v2_daily_doses) 
    tt_primary_doses <- c(0, vaccine_1_start, vaccine_1_start + time_to_coverage_v1, vaccine_2_start + time_to_coverage_v2_elderly ) 
    booster_doses <- c(0, v2_daily_doses, 0) 
    tt_booster_doses <- c(0, vaccine_2_start, vaccine_2_start + time_to_coverage_v2_elderly )
    
    # lift NPIs once the elderly receive vaccine 1
    if(vaccine_1_start < runtime) timing2 <- vaccine_1_start + time_to_coverage_v2_elderly
    
    
  } else {
    
    #calculate how long it will take to boost the elderly pop (who are vaccinated) 
    elderly_pop_to_vaccinate <- sum(pop[priority_age_groups])*coverage
    time_to_coverage_v2_elderly <- ceiling(elderly_pop_to_vaccinate/v2_daily_doses)
    
    #generate vaccine coverage matrix
    vaccine_coverage_mat <- matrix(c( 
      rep(0, 17 - length(priority_age_groups)), rep(coverage, length(priority_age_groups)), 
      rep(0, 17 - length(vaccination_age_groups)), rep(coverage, length(vaccination_age_groups)) ), 
      ncol = 17, byrow = TRUE)
    
    
    vaccine_coverage_mat <- matrix(c( 
      rep(0, 17 - length(vaccination_age_groups)), rep(coverage, length(vaccination_age_groups)) ), 
      ncol = 17, byrow = TRUE)
    
    #format for ve in the booster model (first dose, second dose, second dose, waned, booster, booster , waned) 
    ve_i_booster <- c(ve_v2$infection, ve_v2$infection, ve_v2$infection, 0, 0, 0, 0)
    vaccine_efficacy_infection <- matrix(c( rep(ve_i_booster, 17 )), ncol = 17)
    
    #for disease
    ve_i_booster <- c(ve_v2$disease, ve_v2$disease, ve_v2$disease, 0, 0, 0, 0)
    vaccine_efficacy_disease <- matrix(c( rep(ve_i_booster, 17 )), ncol = 17)
    
    # set up dosing schedule - here the elderly receive v2 ahead of the rest of the population once it becomes available
    
    primary_doses <- c(0, v2_daily_doses)
    tt_primary_doses <- c(0, vaccine_2_start) 
    booster_doses <- 0
    tt_booster_doses <- 0
    
    # lift NPIs once the elderly receive vaccine 1
    if(vaccine_2_start < runtime) timing2 <- vaccine_2_start + time_to_coverage_v2_elderly
    
  }
  
  
  #since we're ignoring first doses set a one day delay between first dose and second dose (in primary series)
  second_dose_delay <- 1 
  dur_V <- rep(duration_V/2, 4) #spend half the duration in each compartment giving erlang-2 waning
  
  
  # R
  R0 <- c(R0, Rt1, Rt1a, Rt1b, Rt2)
  tt_R0 <- c(0, timing1, timing1a, timing1b, timing2)
  
  r1 <- run_booster( 
    time_period = runtime, 
    population = pop,
    contact_matrix_set = mm,
    R0 = R0, 
    tt_R0 = tt_R0,  # timing of changes in R0 to mimic NPIs
    hosp_bed_capacity = 1000000, # high value for unconstrained
    ICU_bed_capacity = 1000000,
    prob_hosp = prob_hosp,
    prob_severe = prob_severe,
    prob_non_severe_death_no_treatment = prob_non_severe_death_treatment,
    prob_severe_death_treatment = prob_severe_death_treatment,
    dur_R = duration_R,
    seeding_cases = seeding_cases,
    vaccine_coverage_mat = vaccine_coverage_mat, 
    vaccine_efficacy_infection = vaccine_efficacy_infection, 
    vaccine_efficacy_disease = vaccine_efficacy_disease, 
    second_dose_delay = second_dose_delay, 
    dur_V = dur_V, 
    protection_delay_rate = 1/dur_vacc_delay,
    primary_doses = primary_doses, 
    tt_primary_doses = tt_primary_doses, 
    booster_doses = booster_doses, 
    tt_booster_doses = tt_booster_doses )
  
  
  # Create output wrt time and age
  
  o1 <- nimue::format(r1,
                      compartments = c("S", "E", "IMild", "ICase", "IICU", "IHospital", "IRec", "R", "D"),
                      summaries = c("deaths", "infections", "vaccinated_second_dose", "vaccinated_booster_dose", "hospitalisations", "hospital_occupancy", "ICU_occupancy", "hospital_demand_cum"),
                      reduce_age = TRUE) %>%
    filter(t>1) 
  
  if(two_vaccines==1) {
    o1a <- nimue::format(r1, 
                         compartments = c("vaccinated_second_dose", "vaccinated_booster_dose"), 
                         reduce_age = FALSE) %>%
      filter(t>1) %>%
      filter(compartment == "vaccinated_second_dose" | compartment == "vaccinated_booster_dose") %>%
      mutate(
        vaccine = if_else(
          compartment == "vaccinated_second_dose" & age_group %in% sort(unique(age_group))[priority_age_groups],
          "V1",
          "V2"
        )
      ) %>%
      group_by(replicate, t, vaccine) %>%
      summarise(
        value = sum(value),
        .groups = "drop"
      ) %>%
      mutate( compartment = vaccine) %>%
      select(-vaccine)
  } else {
    o1a <- nimue::format(r1, 
                         compartments = c("vaccinated_second_dose", "vaccinated_booster_dose"), 
                         reduce_age = FALSE) %>%
      filter(t>1) %>%
      filter(compartment == "vaccinated_second_dose") %>%
      mutate(vaccine = "V2") %>%
      group_by(replicate, t, vaccine) %>%
      summarise(
        value = sum(value),
        .groups = "drop"
      ) %>%
      mutate( compartment = vaccine) %>%
      select(-vaccine) 
  }
  
  o1 <- rbind(o1,o1a)
  
  o2 <- nimue::format(r1,
                      compartments = c("S", "E", "IMild", "ICase", "IICU", "IHospital", "IRec", "R", "D"),
                      summaries = c("deaths", "infections", "vaccinated_second_dose", "vaccinated_booster_dose", "hospitalisations"),
                      reduce_age = FALSE) %>%
    filter(t>1)
  
  if(two_vaccines==1) {
    o2a <- nimue::format(r1, 
                         compartments = c("vaccinated_second_dose", "vaccinated_booster_dose"), 
                         reduce_age = FALSE) %>%
      filter(t>1) %>%
      filter(compartment == "vaccinated_second_dose" | compartment == "vaccinated_booster_dose") %>%
      mutate(
        vaccine = if_else(
          compartment == "vaccinated_second_dose" & age_group %in% sort(unique(age_group))[priority_age_groups],
          "V1",
          "V2"
        )
      ) %>%
      group_by(replicate, t, vaccine, age_group) %>%
      summarise(
        value = sum(value),
        .groups = "drop"
      ) %>%
      mutate( compartment = vaccine) %>%
      select(-vaccine) 
  } else {
    o2a <- nimue::format(r1, 
                         compartments = c("vaccinated_second_dose", "vaccinated_booster_dose"), 
                         reduce_age = FALSE) %>%
      filter(t>1) %>%
      filter(compartment == "vaccinated_second_dose") %>%
      mutate(vaccine = "V2") %>%
      group_by(replicate, t, vaccine, age_group) %>%
      summarise(
        value = sum(value),
        .groups = "drop"
      ) %>%
      mutate( compartment = vaccine) %>%
      select(-vaccine) 
  }
  o2b <- rbind(o2,o2a)
  le_vector$age_group <- factor(le_vector$age_group)
  o2 <- left_join(o2b, le_vector, by="age_group")
  
  final <- tibble(output = list(o1), output_age = list(o2), timing_2_output = list(timing2))
}


# Format output
format_out <- function(out, scenarios){
  # Combine_inputs and outputs
  out1 <- bind_cols(scenarios, bind_rows(out)) %>%
    mutate(timing2 = unlist(timing_2_output)) %>%
    select( - timing_2_output) %>%
    mutate(contact_reduction = timing2*(1-Rt1/R0)) 
  # from ISARIC review
  
  outcf <- filter(out1, (vaccine_1_start == runtime & vaccine_2_start == runtime)) %>%
    mutate(contact_reduction_cf = timing2*(1-Rt1/R0)) %>%
    mutate(timing2_cf = timing2) %>%
    select(-vaccine_1_start, -vaccine_2_start, -timing2, -contact_reduction) %>%
    rename(output_cf = output,
           output_age_cf = output_age) %>%
    unique()
  # Combine runs and counterfactual and estimate summaries
  # out2 <- filter(out1, vaccine_start < 365)
  summaries <- left_join(out1, outcf,by=c("target_pop","income_group","R0","Rt1","Rt1a", "Rt1b", "Rt2","timing1", "timing1a", "timing1b", "ifr_scaling","coverage", 
                                          "vaccination_rate", "duration_R", "duration_V", "dur_vacc_delay", "seeding_cases",
                                          "efficacy_infection_v1","efficacy_disease_v1", "efficacy_infection_v2", "efficacy_disease_v2",
                                          "lower_priority", "lower_vaccine", "runtime","two_vaccines"))
  summaries <- left_join(summaries, vsl, by=c("income_group"))
  summaries_age <- summarise_age_outputs(summaries, median_hospital_days) 
  final <- select(summaries_age,-contains("output_age"))
}


## reads in data set x 
summarise_age_outputs <- function(x, median_hospital_days) {
  mutate(x, 
         infections = round(map_dbl(output_age, pull_total, outcome = "infections"), 2),
         hospitalisations = round(map_dbl(output_age, pull_total, outcome = "hospitalisations"),2),
         deaths = round(map_dbl(output_age, pull_total, outcome = "deaths"), 2),
         yll = round(map_dbl(output_age, summarise_yll), 2),
         yll_prod= round(map_dbl(output_age, summarise_yll_prod), 2),
         ifr = deaths/infections,
         infections_cf = round(map_dbl(output_age_cf, pull_total, outcome = "infections"), 2),
         hospitalisations_cf = round(map_dbl(output_age_cf, pull_total, outcome = "hospitalisations"), 2),
         deaths_cf = round(map_dbl(output_age_cf, pull_total, outcome = "deaths"), 2),
         yll_cf = round(map_dbl(output_age_cf, summarise_yll), 2),
         yll_prod_cf = round(map_dbl(output_age_cf, summarise_yll_prod), 2),
         ifr_cf = deaths_cf/infections_cf,
         infections_averted = infections_cf - infections,
         hospitalisations_averted = hospitalisations_cf - hospitalisations,
         coi_sum = hospitalisations_averted*median_hospital_days*coi_1,
         deaths_averted = deaths_cf - deaths,
         deaths_averted_prop = deaths_averted /deaths_cf,
         years_life_saved = yll_cf - yll,
         years_prod_saved = yll_prod_cf-yll_prod,
         vsly_sum = years_life_saved * vsly,
         prod_loss =years_prod_saved *gni,
         vsl_sum = deaths_averted*vsl)#deaths_averted*vsly, )#,
  #    vaccine_1 = round(map_dbl(output_age, pull_max, outcome = "V1")), #not sure if this is doing the right thing but not using the output currently
  #    vaccine_2 = round(map_dbl(output_age, pull_max, outcome = "V2"))
}

# Pull sum totals
pull_total <- function(x, outcome){
  filter(x, compartment == outcome) %>%
    pull(value) %>%
    sum()
}

pull_max <- function(x, outcome){
  filter(x, compartment == outcome) %>%
    pull(value) %>%
    max()
}

# Estimate total years of life lost
# now have life expectancy in each age-group

summarise_yll <- function(x, lifespan){
  filter(x, compartment == "deaths") %>%
    mutate(yll = life_expect*value) %>%
    pull(yll) %>%
    sum()
}

summarise_yll_prod <- function(x){
  
  prod_ages <- c("15-20","20-25","25-30","30-34","35-40","40-45","45-50","50-55","55-60","60-65")
  
  working_le <- read.csv ("./data/working_life_expectacy.csv") %>% select (c("age_group",disc_le))
  
  filter(x, compartment == "deaths") %>%
    filter(age_group %in% prod_ages ) %>% left_join(working_le, by="age_group") %>%
    mutate(yll_p = disc_le *value) %>%
    pull(yll_p) %>%
    sum()
}

