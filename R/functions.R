
# Specify the country chosen to represent each income group
get_representative_country <- function(income_group){
  case_when(income_group == "HIC" ~ "Canada",
            income_group == "UMIC" ~ "Argentina",
            income_group == "LMIC" ~ "Pakistan",
            income_group == "LIC" ~ "Chad")
}

# Set vaccine efficacy against infection or disease
set_efficacy <- function(efficacy_infection, efficacy_disease){
  out <- list(
      vaccine_efficacy_infection = rep(efficacy_infection,17),
      vaccine_efficacy_disease = rep(efficacy_disease, 17))
  return(out)
}

#use an age-distribution defined by probability of death from severe disease if hospitalised
get_ifr_dist <- function(){

  prob_severe_death_treatment <- c(0.0038,0.0038,0.0038,0.0038,0.0076,0.0114,0.0152,0.0228,0.0304,0.038,0.076,0.152,0.228,0.304,0.38,0.57,0.95)  
  return(prob_severe_death_treatment)
}



# Run scenarios
run_scenario <- function(target_pop = 1e6,
                         income_group = "HIC",
                         R0 = 3,
                         Rt1 = 0.8,
                         Rt2 = 3,
                         timing1 = 30,
                         timing2 = 365,
                         ifr_scaling = 1,
                         coverage = 0, 
                         vaccine_coverage_mat = "Elderly",   # this is the prioritisation matrix 
                         efficacy_infection = 0,
                         efficacy_disease = 0,
                         vaccine_period = 30,
                         duration_R = 1000, # duration of infection-induced immunity
                         duration_V = 5000, # duration of vaccine-induced immunity
                         dur_vacc_delay = 14, # mean duration from vaccination to protection
                         seeding_cases = 5, # define as the number of cases at first sequencing - will need to explore
                         vaccine_start = 20, # 20 days or 100 days after start of the epidemic - will depend on seeding_cases
                         runtime = 365*3
){
  # R
  R0 <- c(R0, Rt1, Rt2)
  tt_R0 <- c(0, timing1, timing2)
  
  # Population and mixing
  rep_country <- get_representative_country(income_group = income_group)
  pop <- squire::get_population(country = rep_country)$n
  pop_standardise <- target_pop / sum(pop)
  pop <- pop * pop_standardise
  mm <- squire::get_mixing_matrix(country = rep_country)
  
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
  
  
  # vaccine efficacy
  ve <- set_efficacy(efficacy_infection, efficacy_disease)
  
  # vaccine prioritisation strategy
  m1 <- strategy_matrix(vaccine_coverage_mat, max_coverage = coverage)
  max_vaccine = c(0, max(coverage) * target_pop / vaccine_period)
  tt_vaccine = c(0, vaccine_start)
  
  # Run
  r1 <- nimue::run(
    time_period = runtime,
    R0 = R0, 
    tt_R0 = tt_R0,  # timing of changes in R0 to mimic NPIs
    population = pop,
    contact_matrix_set = mm,
    hosp_bed_capacity = 1000000, # high value for unconstrained
    ICU_bed_capacity = 1000000,
    prob_hosp = prob_hosp,
    prob_severe = prob_severe,
    prob_non_severe_death_no_treatment = prob_non_severe_death_treatment,
    prob_severe_death_treatment = prob_severe_death_treatment,
    seeding_cases = seeding_cases,
    seed = 1,
    max_vaccine = max_vaccine,
    tt_vaccine = tt_vaccine,
    dur_V = duration_V,
    dur_vaccine_delay = dur_vacc_delay,
    vaccine_efficacy_infection = ve$vaccine_efficacy_infection,
    vaccine_efficacy_disease = ve$vaccine_efficacy_disease,
    dur_R = duration_R,
    vaccine_coverage_mat = m1
  )
  
  # Create output wrt time and age
  
  o1 <- nimue::format(r1,
                      compartments = c("S", "E", "IMild", "ICase", "IICU", "IHospital", "IRec", "R", "D"),
                      summaries = c("deaths", "infections", "vaccines", "hospitalisations"),
                      reduce_age = TRUE) %>%
    filter(t>1)

  o2 <- nimue::format(r1,
                      compartments = c("S", "E", "IMild", "ICase", "IICU", "IHospital", "IRec", "R", "D"),
                      summaries = c("deaths", "infections", "vaccines", "hospitalisations"),
                      reduce_age = FALSE) %>%
    filter(t>1)
  final <- tibble(output = list(o1), output_age = list(o2))
}

# Format output
format_out <- function(out, scenarios){
  # Combine_inputs and outputs
  out1 <- bind_cols(scenarios, bind_rows(out)) %>%
    mutate(contact_reduction = timing2*(1-Rt1/R0))
    
  outcf <- filter(out1, vaccine_start == 365*2) %>%
      select(-vaccine_start,-timing2) %>%
      rename(output_cf = output,
             output_age_cf = output_age) %>%
      unique()
  # Combine runs and counterfactual and estimate summaries
 # out2 <- filter(out1, vaccine_start < 365)
  summaries <- left_join(out1, outcf,by=c("target_pop","income_group","R0","Rt1","Rt2","timing1","ifr_scaling","efficacy_infection","efficacy_disease"))
  summaries_age <- summarise_age_outputs(summaries) 
  final <- select(summaries_age,-contains("output_age"))
}

## reads in data set x 
summarise_age_outputs <- function(x) {
  mutate(x, 
         infections = round(map_dbl(output_age, pull_total, outcome = "infections"), 2),
         hospitalisations = round(map_dbl(output_age, pull_total, outcome = "hospitalisations"),2),
         deaths = round(map_dbl(output_age, pull_total, outcome = "deaths"), 2),
         yll = round(map_dbl(output_age, summarise_yll), 2),
         ifr = deaths/infections,
         infections_cf = round(map_dbl(output_age_cf, pull_total, outcome = "infections"), 2),
         hospitalisations_cf = round(map_dbl(output_age_cf, pull_total, outcome = "hospitalisations"), 2),
         deaths_cf = round(map_dbl(output_age_cf, pull_total, outcome = "deaths"), 2),
         yll_cf = round(map_dbl(output_age_cf, summarise_yll), 2),
         ifr_cf = deaths_cf/infections_cf,
         infections_averted = infections_cf - infections,
         hospitalisations_averted = hospitalisations_cf - hospitalisations,
         deaths_averted = deaths_cf - deaths,
         deaths_averted_prop = deaths_averted / deaths_cf,
         years_life_saved = yll_cf - yll,
         vaccine_n = round(map_dbl(output_age, pull_total, outcome = "vaccines"))
 )
}

# Pull sum totals
pull_total <- function(x, outcome){
  filter(x, compartment == outcome) %>%
    pull(value) %>%
    sum()
}

# Estimate total years of life lost
summarise_yll <- function(x, lifespan = 86.6){
  filter(x, compartment == "deaths") %>%
    mutate(mid_age = (((as.integer(age_group) - 1) * 5) + 2.5),
           yll = pmax(0, (lifespan - mid_age) * value)) %>%
    pull(yll) %>%
    sum()
}







