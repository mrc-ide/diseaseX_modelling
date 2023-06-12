## Notes - talk to Azra about how to adjust the IFR.


run_sars_x <- function(## Demographic Parameters
                       population_size = 1e6,
                       country = "Argentina",
                       
                       ## Healthcare-Related Parameters
                       hosp_bed_capacity = NULL,
                       ICU_bed_capacity = NULL,
                       
                       ## Epidemiological Parameters
                       Rt = 3,
                       tt_Rt = 0, 
                       Tg = 6.6,
                       IFR = 1,
                       
                       ## Vaccine-Related Parameters
                       vaccine_scenario = "specific_only",       # which scenario to explore
                       bpsv_start = 50,                          # BPSV distribution start
                       specific_vaccine_start = 150,             # specific vaccine distribution start
                       efficacy_infection_bpsv = 0.35,           # vaccine efficacy against infection - BPSV
                       efficacy_disease_bpsv = 0.8,              # vaccine efficacy against disease - BPSV
                       efficacy_infection_spec = 0.55,           # vaccine efficacy against infection - specific vaccine
                       efficacy_disease_spec = 0.9,              # vaccine efficacy against disease - specific vaccine
                       dur_R = 1000 * 365,                       # duration of infection-induced immunity
                       dur_V = rep(1000 * 365, 4),              # duration of vaccine-induced immunity for both vaccines
                       second_dose_delay = 1,                  # controls how many days after "1st dose" people receive second dose; see here: https://github.com/mrc-ide/squire.page/blob/main/inst/odin/nimue_booster.R#L427-L430
                       dur_vacc_delay = 1,                     # mean duration from vaccination to protection
                       coverage = 0.75,                           # proportion of the population vaccinated
                       vaccination_rate = 0.04,                  # vaccination rate per week as percentage of population
                       min_age_group_index_priority = 13,        # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                       min_age_group_index_non_priority = 4,     # index of the youngest age group that *receives* vaccines (4 = 15+)
                       
                       ## Miscellaneous Parameters
                       runtime = 600,
                       seeding_cases = 5) {
  
  ## Country-specific population estimate
  raw_pop <- squire::get_population(country = country)$n                   
  standard_pop <- round(raw_pop * population_size / sum(raw_pop)) 
  
  ## Country-specific healthcare capacity
  if (is.null(hosp_bed_capacity)) {
    hosp_bed_capacity <- population_size
  }
  if (is.null(ICU_bed_capacity)) {
    ICU_bed_capacity <- population_size
  }
  
  ## Country-specific mixing matrix 
  mm <- squire::get_mixing_matrix(country = country)    
  
  ## Adapting Generation Time
  current_Tg <- squire.page:::durs_booster$dur_IMild + squire.page:::durs_booster$dur_ICase
  Tg_ratio <- Tg / current_Tg
  TgVary_dur_IMild <- Tg_ratio * squire.page:::durs_booster$dur_IMild
  TgVary_dur_ICase <- Tg_ratio * squire.page:::durs_booster$dur_ICase
  
  ## Adjusting IFR (NEEDS CHANGING AND CHECKING)
  contact_rates <- apply(mm, 1, sum)
  contact_rates <- c(contact_rates, "80+" = unname(contact_rates[16]))
  pop_contact_weighting <- standard_pop * contact_rates
  non_severe_deaths <- nimue:::probs$prob_hosp * (1 - nimue:::probs$prob_severe) * nimue:::probs$prob_non_severe_death_treatment
  severe_deaths <- nimue:::probs$prob_hosp * nimue:::probs$prob_severe * nimue:::probs$prob_severe_death_treatment
  raw_IFR <- 100 * sum(((non_severe_deaths + severe_deaths) * pop_contact_weighting/sum(pop_contact_weighting))) 
  IFR_scaling_factor <- IFR / raw_IFR
  prob_hosp <- IFR_scaling_factor * nimue:::probs$prob_hosp
  if (sum(prob_hosp > 1) > 0) {
    stop("IFR Scaling Factor too high - results in prob_hosp > 1 for at least one age-group")
  }
  
  ## Setting Up Vaccination Stuff
  if (!vaccine_scenario %in% c("specific_only", "both_vaccines")) {
    stop("parameter vaccine_scenario must be either 'specific_only' or 'both_vaccines'")
  }
  
  if (vaccine_scenario == "specific_only") { ## in this scenario, first vaccine is virus-specific, for all ages
    
    ## Setting vaccination rates and coverages
    daily_doses <- vaccination_rate * sum(standard_pop) / 7         # rate of vaccination with specific vaccine
    priority_age_groups <- min_age_group_index_priority:17          # priority age groups  
    vaccination_age_groups <- min_age_group_index_non_priority:17   # vaccinable age-groups
    vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(coverage, length(priority_age_groups)), 
                                     rep(0, 17 - length(vaccination_age_groups)), rep(coverage, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)

    ## Creating matrices of infection and disease efficacy for each age-group and vaccine received
    ### Infection Efficacy
    ve_i_spec <- c(efficacy_infection_spec, efficacy_infection_spec, efficacy_infection_spec, 0, 0, 0, 0) ## note that booster protection set to 0 as we're assuming virus-spec gives robust protection - no waning, so no booster required
    vaccine_efficacy_infection <- matrix(c(rep(ve_i_spec, 17)), ncol = 17)
    
    ### Disease Efficacy
    ve_d_spec <- c(efficacy_disease_spec, efficacy_disease_spec, efficacy_disease_spec, 0, 0, 0, 0) ## note that booster protection set to 0 as we're assuming virus-spec gives robust protection - no waning, so no booster required
    vaccine_efficacy_disease <- matrix(c(rep(ve_d_spec, 17)), ncol = 17)
    
    # Set up dosing schedule - no boosters here
    primary_doses <- c(0, daily_doses)
    tt_primary_doses <- c(0, specific_vaccine_start) 
    booster_doses <- 0
    tt_booster_doses <- 0
    
    ## Eligibility for Boosters - No One
    vaccine_booster_follow_up_coverage <- rep(0, 17)
    vaccine_booster_initial_coverage <- rep(0, 17)
    
    
  } else if (vaccine_scenario == "both_vaccines") { ## in this scenario, first vaccine is BPSV for 60+, virus-specific for <60, and booster for 60+ is the virus-specific
    
    ## Setting vaccination rates and coverages
    bpsv_daily_doses <- vaccination_rate * sum(standard_pop) / 7    # rate of vaccination with primary series
    spec_daily_doses <- vaccination_rate * sum(standard_pop) / 7    # rate of vaccination with booster series
    priority_age_groups <- min_age_group_index_priority:17          # priority age groups  
    vaccination_age_groups <- min_age_group_index_non_priority:17   # vaccinable age-groups
    vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(coverage, length(priority_age_groups)), 
                                     rep(0, 17 - length(vaccination_age_groups)),  rep(coverage, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)
    # matrix (col = age group, row = prioritisation step) of how to assign vaccines to different age-groups according to priority
    
    ## Calculating time taken to vaccine elderly with 1) BPSV; then 2) virus-specific (given they get priority)
    elderly_pop_to_vaccinate <- sum(standard_pop[priority_age_groups]) * coverage # 60+s receive primary (BNPCV) and booster (diseaseX-specific); under 60s receive just primary (diseaseX-specific)
    time_to_coverage_bpsv <- ceiling(elderly_pop_to_vaccinate / bpsv_daily_doses)     # calculated so that we only vaccinate 60+s with the stockpiled BNPSC
    minimum_spec_development_time_allowed <- (time_to_coverage_bpsv + second_dose_delay + dur_vacc_delay)
    if (minimum_spec_development_time_allowed > (specific_vaccine_start - bpsv_start)) {
      stop("Virus-specific vaccine developed too soon given speed of vaccination campaign")
    }
    time_to_coverage_spec_elderly <- ceiling(elderly_pop_to_vaccinate/spec_daily_doses)
    primary_doses <- c(0, bpsv_daily_doses, 0, spec_daily_doses) 
    tt_primary_doses <- c(0, bpsv_start, ceiling(bpsv_start + time_to_coverage_bpsv), specific_vaccine_start + time_to_coverage_spec_elderly) 
    booster_doses <- c(0, spec_daily_doses, 0) 
    tt_booster_doses <- c(0, specific_vaccine_start, specific_vaccine_start + time_to_coverage_spec_elderly) 
    
    ## Creating matrices of infection and diasease efficacy for each age-group and vaccine received
    ### Infection Efficacy
    ve_bpsv <- list(infection = efficacy_infection_bpsv, disease = efficacy_disease_bpsv)
    ve_spec <- list(infection = efficacy_infection_spec, disease = efficacy_disease_spec)
    ve_i_elderly <- c(ve_bpsv$infection, ve_bpsv$infection, ve_bpsv$infection, 0, ve_spec$infection, ve_spec$infection, 0) 
    ve_i_non_elderly <- c(ve_spec$infection, ve_spec$infection, ve_spec$infection, 0, ve_spec$infection, ve_spec$infection, 0)
    vaccine_efficacy_infection <- matrix(c( rep(ve_i_non_elderly, 17 - length(priority_age_groups)), 
                                            rep(ve_i_elderly, length(priority_age_groups)) ), ncol = 17)
    ### Disease Efficacy
    ve_d_elderly <- c(ve_bpsv$disease, ve_bpsv$disease, ve_bpsv$disease, 0, ve_spec$disease, ve_spec$disease, 0) 
    ve_d_non_elderly <- c(ve_spec$disease, ve_spec$disease, ve_spec$disease, 0, ve_spec$disease, ve_spec$disease, 0)
    vaccine_efficacy_disease <- matrix(c( rep(ve_d_non_elderly, 17 - length(priority_age_groups)), 
                                          rep(ve_d_elderly, length(priority_age_groups)) ), ncol = 17)
    
    ## Eligibility for Boosters (Elderly Only)
    vaccine_booster_follow_up_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(1, length(priority_age_groups)))
    vaccine_booster_initial_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(1, length(priority_age_groups)))
    }
  
  ## Running the Model 
  mod_run <- run_booster(time_period = runtime,
                         population = standard_pop,                                                 
                         contact_matrix_set = mm,                                                   
                         R0 = Rt,     
                         tt_R0 = tt_Rt, 
                         hosp_bed_capacity = hosp_bed_capacity,                                     
                         ICU_bed_capacity = ICU_bed_capacity,                                       
                         prob_hosp = prob_hosp,                                                     
                         dur_R = dur_R,                                                        
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
                         tt_booster_doses = tt_booster_doses,                                       
                         vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
                         vaccine_booster_initial_coverage = vaccine_booster_initial_coverage)       
  
  return(list(
         model_output = mod_run,
         model_arguments = list(population_size = population_size, standard_pop = standard_pop, country = country, hosp_bed_capacity = hosp_bed_capacity, ICU_bed_capacity = ICU_bed_capacity,
                                Rt = Rt, tt_Rt = tt_Rt, Tg = Tg, IFR = IFR, vaccine_scenario = vaccine_scenario, bpsv_start = bpsv_start, specific_vaccine_start = specific_vaccine_start,             
                                efficacy_infection_bpsv = efficacy_infection_bpsv, efficacy_disease_bpsv = efficacy_disease_bpsv, efficacy_infection_spec = efficacy_infection_spec,
                                efficacy_disease_spec = efficacy_disease_spec, dur_R = dur_R, dur_V = dur_V, second_dose_delay = second_dose_delay, dur_vacc_delay = dur_vacc_delay,                     
                                coverage = coverage, vaccination_rate = vaccination_rate, min_age_group_index_priority = min_age_group_index_priority, min_age_group_index_non_priority = min_age_group_index_non_priority,     
                                runtime = runtime, seeding_cases = seeding_cases,
                                tt_primary_doses = tt_primary_doses, tt_booster_doses = tt_booster_doses,
                                vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage, 
                                vaccine_booster_initial_coverage = vaccine_booster_initial_coverage)))
}

y <- run_sars_x(vaccine_scenario = "both_vaccines",
                IFR = 0.01,
                specific_vaccine_start = 250,
                dur_vacc_delay = 1,
                second_dose_delay = 1)  
x <- y$model_output

check <- nimue::format(x, compartments = c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_booster_dose"), 
                       reduce_age = FALSE) %>%
  filter(t > 1,
         compartment == "vaccinated_first_dose" | compartment == "vaccinated_second_dose" | compartment == "vaccinated_booster_dose") %>%
  group_by(replicate, t, age_group) 

pop_df <- data.frame(age_group = sort(unique(check$age_group)),
                     population = y$model_arguments$standard_pop)

check2 <- check %>%
  left_join(pop_df, by = "age_group") %>%
  mutate(prop = value / population) 
ggplot() +
  geom_line(data = check2, aes(x = t, y = prop, col = compartment)) +
  facet_wrap(~age_group) +
  lims(y = c(0, 1), x = c(0, 600)) 


ggplot() +
  geom_line(data = check2[check2$age_group == "65-70", ], aes(x = t, y = prop, col = compartment)) +
  facet_wrap(~age_group) +
  lims(y = c(0, 1), x = c(0, 600)) 
