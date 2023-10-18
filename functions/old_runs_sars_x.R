## Generate vaccination milestone timings based on parameter combinations
### (deprecated because we needed different coverages for bpsv and spec and this wasn't set up to handle that)
### In the old version, for elderly, primary was bpsv, secondary was bpsv and booster was spec.
### In the new version, for elderly, primary is bpsv during bpsv campaign and spec during spec campaign, secondary is always bpsv, and booster is always spec
vaccine_milestone_timings_Old <- function(country, 
                                          population_size, 
                                          detection_time, 
                                          bpsv_start, 
                                          specific_vaccine_start,
                                          vaccination_rate_bpsv,
                                          vaccination_rate_spec,
                                          coverage,
                                          min_age_group_index_priority) {
  
  # Calculating time to coverage from parameter inputs
  standard_pop <- generate_standard_pop(country = country, population_size = population_size)
  priority_age_groups <- min_age_group_index_priority:17            
  elderly_pop_to_vaccinate <- sum(standard_pop[priority_age_groups]) * coverage # 60+s receive primary (BNPCV) and booster (diseaseX-specific); under 60s receive just primary (diseaseX-specific)
  
  daily_doses_bpsv <- vaccination_rate_bpsv * population_size / 7    # rate of vaccination with BPSV
  time_to_coverage_bpsv <- ceiling(elderly_pop_to_vaccinate/daily_doses_bpsv) + 1
  daily_doses_spec <- vaccination_rate_spec * population_size / 7    # rate of vaccination with specific vaccine
  time_to_coverage_spec <-  ceiling(elderly_pop_to_vaccinate/daily_doses_spec) + 1
  
  # Returning dataframe of parameters and key vaccination milestone timings
  temp <- data.frame(country = country,
                     population_size = population_size,
                     detection_time = detection_time,
                     bpsv_start = bpsv_start,
                     specific_vaccine_start = specific_vaccine_start,
                     vaccination_rate_bpsv = vaccination_rate_bpsv, 
                     vaccination_rate_spec = vaccination_rate_spec, 
                     coverage = coverage,
                     min_age_group_index_priority = min_age_group_index_priority,
                     time_to_coverage_bpsv = time_to_coverage_bpsv,
                     time_to_coverage_spec = time_to_coverage_spec)
  return(temp)
}

## Create primary vaccination series
### (deprecated because we needed different coverages for bpsv and spec and this wasn't set up to handle that)
### In the old version, for elderly, primary was bpsv, secondary was bpsv and booster was spec.
### In the new version, for elderly, primary is bpsv during bpsv campaign and spec during spec campaign, secondary is always bpsv, and booster is always spec
create_vaccination_dose_series_Old <- function(country, 
                                               population_size, 
                                               detection_time, 
                                               vaccine_scenario,
                                               bpsv_start, 
                                               bpsv_protection_delay,
                                               specific_vaccine_start,
                                               specific_protection_delay,
                                               vaccination_rate_bpsv,
                                               vaccination_rate_spec,
                                               coverage,
                                               min_age_group_index_priority,
                                               runtime) {
  
  ## Setting Up Vaccination Stuff
  if (!vaccine_scenario %in% c("specific_only", "both_vaccines")) {
    stop("parameter vaccine_scenario must be either 'specific_only' or 'both_vaccines'")
  }
  
  if (vaccine_scenario == "specific_only") {
    
    standard_pop <- generate_standard_pop(country = country, population_size = population_size)
    daily_doses_spec <- vaccination_rate_spec * population_size / 7    # rate of vaccination with primary series
    primary_doses <- 
      c(rep(0, detection_time),
        rep(0, specific_vaccine_start),
        rep(0, specific_protection_delay),        
        rep(daily_doses_spec, runtime - detection_time - specific_vaccine_start - specific_protection_delay))
    
    ## Second Doses
    ## We model full protection as arising the moment you've had the first dose, so second dose is 
    ## redundant given the vaccine-specific protection_delays. We arbitrarily set it to 1 here and therefore stagger
    ## second doses to be a day after primary doses.
    second_doses <- c(0, primary_doses[1:(length(primary_doses) - 1)]) # second dose 1 day after first
    
    ## Booster Doses (only for elderly, this is the disease specific vaccine)
    booster_doses <- rep(0, runtime)
    
  } else {
    
    # Calculating time to coverage from parameter inputs
    standard_pop <- generate_standard_pop(country = country, population_size = population_size)
    priority_age_groups <- min_age_group_index_priority:17            
    elderly_pop_to_vaccinate <- ceiling(sum(standard_pop[priority_age_groups]) * coverage) # 60+s receive primary (BNPCV) and booster (diseaseX-specific); under 60s receive just primary (diseaseX-specific)
    
    daily_doses_bpsv <- vaccination_rate_bpsv * population_size / 7    # rate of vaccination with BPSV
    time_to_coverage_bpsv <- ceiling(elderly_pop_to_vaccinate/daily_doses_bpsv) + 1
    daily_doses_spec <- vaccination_rate_spec * population_size / 7    # rate of vaccination with specific vaccine
    time_to_coverage_spec <-  ceiling(elderly_pop_to_vaccinate/daily_doses_spec) + 1
    
    ## Checking there's no overlap in vaccination campaigns
    ## I think technically we want to avoid BPSV elderly and spec everyone else campaigns (as they're the same series in the model)
    ## so I *think* saying BPSV elderly has to be complete before spec elderly (which is before spec everyone else)
    ## is a conservative estimate
    ## think even more conservative because not counting spec delay in protection development
    minimum_spec_development_time_allowed <- (time_to_coverage_bpsv + bpsv_protection_delay)
    if (minimum_spec_development_time_allowed > (specific_vaccine_start - bpsv_start)) {
      stop("Virus-specific vaccine developed too soon given speed of vaccination campaign")
    }
    
    ## Primary Doses (initial campaign is BPSV for elderly, then specific for everyone who isn't elderly)
    ### Note that we need to check whether another "+1" is needed somewhere to make sure we always vaccinate all the elderly
    primary_doses <- 
      c(rep(0, detection_time),                     ## time between epidemic start and detection
        rep(0, bpsv_start),                         ## time between detection and initiation of BPSV campaign
        rep(0, bpsv_protection_delay),              ## time between initiation of BPSV campaign and people first being protected by that first dose
        rep(daily_doses_bpsv, time_to_coverage_bpsv),    ## protection (if any) emerges in BPSV-vaccinated primary vaccinated folks
        rep(0, specific_vaccine_start - time_to_coverage_bpsv - bpsv_protection_delay - bpsv_start), # specific vaccine campaign starts specific_vaccine_start days after detection
        rep(0, time_to_coverage_spec),              ## no specific vaccine for non-elderly whilst that vaccination campaign is ongoing
        rep(0, specific_protection_delay),          ## time between initiation of specific vaccine campaign and people first being protected by that first dose
        rep(daily_doses_spec, runtime - specific_protection_delay - time_to_coverage_spec - specific_vaccine_start - detection_time)) # specific vaccination of all other ages until end of runtime
    
    ## Second Doses
    ## We model full protection as arising the moment you've had the first dose, so second dose is 
    ## redundant given the vaccine-specific protection_delays. We arbitrarily set it to 1 here and therefore stagger
    ## second doses to be a day after primary doses.
    second_doses <- c(0, primary_doses[1:(length(primary_doses) - 1)]) # second dose 1 day after first
    
    # primary_doses[(detection_time+bpsv_start+bpsv_protection_delay+1):
    #                 (detection_time+bpsv_start+bpsv_protection_delay+50)]
    
    # second_doses[(detection_time+bpsv_start+bpsv_protection_delay+1):
    #                (detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+1)] <- daily_doses
    
    ## extra second doses to get second coverage up to coverage target - why are we not getting there?
    ## unclear to me why, but as long as we're only vaccinating elderly with primary doses
    ## during the initial phase, they're the only folks eligible for second doses during that
    ## initial campaign and so it doesn't make much difference. Just need to watch out
    ## for when detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+4 is very
    ## close to the detection_time + specific_vaccine_start + spec_protection_delay time 
    second_doses[(detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+1)] <- daily_doses_bpsv
    second_doses[(detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+2)] <- daily_doses_bpsv
    second_doses[(detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+3)] <- daily_doses_bpsv
    second_doses[(detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+4)] <- daily_doses_bpsv
    second_doses[(detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+5)] <- daily_doses_bpsv
    # 
    
    ## Booster Doses (only for elderly, this is the disease specific vaccine)
    ## because we now have booster coverage eligibility indicator in, we don't have to worry
    ## about stopping the booster doses - only elderly will ever get it
    booster_doses <- 
      c(rep(0, detection_time),
        rep(0, specific_vaccine_start),
        rep(0, specific_protection_delay),
        rep(daily_doses_spec, runtime - specific_protection_delay - specific_vaccine_start - detection_time))
    
  }
  
  return(list(primary_doses = primary_doses,
              second_doses = second_doses,
              booster_doses = booster_doses))
  
}

## Running the model and summarising outputs
### (deprecated because we needed different coverages for bpsv and spec and this wasn't set up to handle that)
### In the old version, for elderly, primary was bpsv, secondary was bpsv and booster was spec.
### In the new version, for elderly, primary is bpsv during bpsv campaign and spec during spec campaign, secondary is always bpsv, and booster is always spec
run_sars_x_Old <- function(## Demographic Parameters
  population_size = 1e6,
  country = "Argentina",
  
  ## Healthcare-Related Parameters
  hosp_bed_capacity = NULL,
  ICU_bed_capacity = NULL,
  
  ## Epidemiological Parameters
  Rt = 2.5,
  tt_Rt = 0, 
  Tg = 6.6,
  IFR = 0.5,
  
  ## Vaccine-Related Parameters
  vaccine_scenario = "specific_only",       # which scenario to explore
  detection_time = 14,
  bpsv_start = 50,                          # BPSV distribution start (days after detection time)
  bpsv_protection_delay = 7,                # time between BPSV dose and protection arising
  specific_vaccine_start = 250,             # specific vaccine distribution start (days after detection time)
  specific_protection_delay = 7,            # time between specific dose and protection arising
  efficacy_infection_bpsv = 0.35,           # vaccine efficacy against infection - BPSV
  efficacy_disease_bpsv = 0.8,              # vaccine efficacy against disease - BPSV
  efficacy_infection_spec = 0.55,           # vaccine efficacy against infection - specific vaccine
  efficacy_disease_spec = 0.9,              # vaccine efficacy against disease - specific vaccine
  dur_R = 1000 * 365,                       # duration of infection-induced immunity
  dur_bpsv = 1000 * 365,                    # duration of vaccine-induced immunity for BPSV vaccine
  dur_spec = 1000 * 365,                    # duration of vaccine-induced immunity for disease-specific vaccines
  coverage_bpsv = 0.75,                     # proportion of the population vaccinated
  coverage_spec = 0.75,                     # proportion of the population vaccinated
  vaccination_rate_bpsv = 0.01,             # bpsv vaccination rate per week as percentage of population (note: percentage of whole population, despite restricted age-groups receiving it)
  vaccination_rate_spec = 0.01,             # disease-specific vaccination rate per week as percentage of population
  min_age_group_index_priority = 13,        # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
  min_age_group_index_non_priority = 4,     # index of the youngest age group that *receives* vaccines (4 = 15+)
  
  ## Miscellaneous Parameters
  runtime = 600,
  seeding_cases = 5,
  output = "summary",
  NPI_int = 0,
  scenario_index = 0,
  varied = "",
  ...) {
  
  # print(scenario_index)
  
  set.seed(123)
  
  ## Input Checking
  if (!approach %in% c("old", "new")) {
    stop("parameter approach must be either 'old' or 'new'")
  }
  if (!output %in% c("summary", "full")) {
    stop("parameter output must be either 'summary' or 'full'")
  }
  if (is.null(hosp_bed_capacity)) {
    hosp_bed_capacity <- population_size
  }
  if (is.null(ICU_bed_capacity)) {
    ICU_bed_capacity <- population_size
  }
  
  ## Country-specific population estimate
  standard_pop <- generate_standard_pop(country = country, population_size = population_size)
  
  ## Adapting Generation Time
  TgVary <- scale_generation_time(target_Tg = Tg)
  TgVary_dur_IMild <- TgVary$dur_IMild
  TgVary_dur_ICase <- TgVary$dur_ICase
  
  ## Adjusting IFR
  prob_hosp <- scale_IFR(country = country, population_size = population_size, target_IFR = IFR)
  if (sum(prob_hosp > 1) > 0) {
    stop("IFR Scaling Factor too high - results in prob_hosp > 1 for at least one age-group")
  }
  
  ## Setting Up Vaccination Stuff
  if (!vaccine_scenario %in% c("specific_only", "both_vaccines")) {
    stop("parameter vaccine_scenario must be either 'specific_only' or 'both_vaccines'")
  }
  
  if (vaccine_scenario == "specific_only") { ## in this scenario, first vaccine is virus-specific, for all ages
    
    ## Setting vaccination rates and coverages
    priority_age_groups <- min_age_group_index_priority:17          # priority age groups  
    vaccination_age_groups <- min_age_group_index_non_priority:17   # vaccinable age-groups
    vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(coverage_spec, length(priority_age_groups)), 
                                     rep(0, 17 - length(vaccination_age_groups)), rep(coverage_spec, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)
    
    ## Creating matrices of infection and disease efficacy for each age-group and vaccine received
    ### Infection Efficacy
    ve_i_spec <- c(efficacy_infection_spec, efficacy_infection_spec, efficacy_infection_spec, 0, 0, 0, 0) ## note that booster protection set to 0 as we're assuming virus-spec gives robust protection - no waning, so no booster required
    vaccine_efficacy_infection <- matrix(c(rep(ve_i_spec, 17)), ncol = 17)
    
    ### Disease Efficacy
    ve_d_spec <- c(efficacy_disease_spec, efficacy_disease_spec, efficacy_disease_spec, 0, 0, 0, 0) ## note that booster protection set to 0 as we're assuming virus-spec gives robust protection - no waning, so no booster required
    vaccine_efficacy_disease <- matrix(c(rep(ve_d_spec, 17)), ncol = 17)
    
    ### Vaccine protection waning
    dur_V <- matrix(data = c(rep(dur_spec, 17),   ## duration of primary series protection (specific vaccine in this scenario)
                             rep(dur_spec, 17),   ## duration of primary series protection (specific vaccine in this scenario)
                             rep(dur_spec, 17),   ## duration of booster protection (not used here)
                             rep(dur_spec, 17)),  ## duration of booster protection (not used here)
                    nrow = 4, ncol = 17, byrow = TRUE)
    
    # Set up dosing schedule - no boosters here
    vaccine_doses <- create_vaccination_dose_series(country = country, 
                                                    population_size = population_size, 
                                                    detection_time = detection_time, 
                                                    vaccine_scenario = vaccine_scenario,
                                                    bpsv_start = bpsv_start, 
                                                    bpsv_protection_delay = bpsv_protection_delay,
                                                    specific_vaccine_start = specific_vaccine_start,
                                                    specific_protection_delay = specific_protection_delay,
                                                    vaccination_rate_bpsv = vaccination_rate_bpsv,
                                                    vaccination_rate_spec = vaccination_rate_spec,
                                                    coverage = coverage_spec,
                                                    min_age_group_index_priority = min_age_group_index_priority,
                                                    runtime = runtime)
    primary_doses <- vaccine_doses$primary_doses
    second_doses <- vaccine_doses$second_doses
    booster_doses <- vaccine_doses$booster_doses
    
    ## Eligibility for Boosters - No One
    vaccine_booster_follow_up_coverage <- rep(0, 17)
    vaccine_booster_initial_coverage <- rep(0, 17)
    
  } else if (vaccine_scenario == "both_vaccines") { ## in this scenario, first vaccine is BPSV for 60+, virus-specific for <60, and booster for 60+ is the virus-specific
    
    ## Setting vaccination rates and coverages
    priority_age_groups <- min_age_group_index_priority:17          # priority age groups  
    vaccination_age_groups <- min_age_group_index_non_priority:17   # vaccinable age-groups
    vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(coverage_spec, length(priority_age_groups)), 
                                     rep(0, 17 - length(vaccination_age_groups)),  rep(coverage_spec, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)
    # matrix (col = age group, row = prioritisation step) of how to assign vaccines to different age-groups according to priority
    
    ## Calculating vaccine dosing schedules for this scenario
    vaccine_doses <- create_vaccination_dose_series(country = country, 
                                                    population_size = population_size, 
                                                    detection_time = detection_time, 
                                                    vaccine_scenario = vaccine_scenario,
                                                    bpsv_start = bpsv_start, 
                                                    bpsv_protection_delay = bpsv_protection_delay,
                                                    specific_vaccine_start = specific_vaccine_start,
                                                    specific_protection_delay = specific_protection_delay,
                                                    vaccination_rate_bpsv = vaccination_rate_bpsv,
                                                    vaccination_rate_spec = vaccination_rate_spec,
                                                    coverage = coverage_spec,
                                                    min_age_group_index_priority = min_age_group_index_priority,
                                                    runtime = runtime)
    primary_doses <- vaccine_doses$primary_doses
    second_doses <- vaccine_doses$second_doses
    booster_doses <- vaccine_doses$booster_doses
    # note that we have the "+1"s in the underlying function because otherwise we don't *quite* get to vaccinating all of the elderly in
    # the alloted time that Azra was previously calculating. I thought doing ceiling() would be enough and in most cases it is but not all.
    # This way (i.e. by adding +1 to the times) we do! 
    
    ## Creating matrices of infection and disease efficacy for each age-group and vaccine received
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
    
    ### Vaccine protection waning
    dur_V <- matrix(data = c(c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_bpsv, length(priority_age_groups))),    ## duration of primary series protection (bpsv for elderly, specific vaccine for everyone else in this scenario)
                             c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_bpsv, length(priority_age_groups))),    ## duration of primary series protection (bpsv for elderly, specific vaccine for everyone else in this scenario)
                             c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_spec, length(priority_age_groups))),           ## duration of booster protection (specific for elderly, not used for everyone else)
                             c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_spec, length(priority_age_groups)))),         ## duration of booster protection (specific for elderly, not used for everyone else)
                    nrow = 4, ncol = 17, byrow = TRUE)
    
    ## Eligibility for Boosters (Elderly Only)
    vaccine_booster_follow_up_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(0, length(priority_age_groups))) # no booster follow ons
    vaccine_booster_initial_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(1, length(priority_age_groups))) ## note that we set coverage to be 100% as it's already bounded by whatever coverage the 1st/2nd dose coverage got to
    
  } 
  
  ## Ensuring things that could be passed in as lists (when running in parallel) are coerced to right format
  if (is.list(tt_Rt)) {
    tt_Rt <- unlist(tt_Rt)
  }
  if (is.list(Rt)) {
    Rt <- unlist(Rt)
  }
  
  ## Running the Model 
  mm <- squire::get_mixing_matrix(country = country)    
  rel_infectiousness_vaccinated <- squire.page:::probs_booster$rel_infectiousness_vaccinated
  rel_infectiousness_vaccinated[rel_infectiousness_vaccinated < 1] <- 1
  
  if (vaccine_scenario == "specific_only") {
    mod_run <- run_booster(time_period = runtime,
                           population = standard_pop,                                                 
                           contact_matrix_set = mm,                                                   
                           R0 = Rt,     
                           tt_R0 = tt_Rt, 
                           hosp_bed_capacity = hosp_bed_capacity,                                     
                           ICU_bed_capacity = ICU_bed_capacity,                                       
                           prob_hosp = prob_hosp,
                           dur_IMild = TgVary_dur_IMild,
                           dur_ICase = TgVary_dur_ICase,
                           dur_R = dur_R,                                                        
                           seeding_cases = seeding_cases,
                           vaccine_coverage_mat = vaccine_coverage_mat,                               
                           vaccine_efficacy_infection = vaccine_efficacy_infection,                   
                           vaccine_efficacy_disease = vaccine_efficacy_disease,     
                           rel_infectiousness_vaccinated = rel_infectiousness_vaccinated, 
                           dur_V = dur_V,                                              
                           primary_doses = primary_doses,  
                           second_doses = second_doses,
                           booster_doses = booster_doses,                                             
                           vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
                           vaccine_booster_initial_coverage = vaccine_booster_initial_coverage)
  } else if (vaccine_scenario == "both_vaccines") {
    mod_run <- run_booster(time_period = runtime,
                           population = standard_pop,                                                 
                           contact_matrix_set = mm,                                                   
                           R0 = Rt,     
                           tt_R0 = tt_Rt, 
                           hosp_bed_capacity = hosp_bed_capacity,                                     
                           ICU_bed_capacity = ICU_bed_capacity,                                       
                           prob_hosp = prob_hosp,
                           dur_IMild = TgVary_dur_IMild,
                           dur_ICase = TgVary_dur_ICase,
                           dur_R = dur_R,                                                        
                           seeding_cases = seeding_cases,
                           vaccine_coverage_mat = vaccine_coverage_mat,                               
                           vaccine_efficacy_infection = vaccine_efficacy_infection,                   
                           vaccine_efficacy_disease = vaccine_efficacy_disease,     
                           rel_infectiousness_vaccinated = rel_infectiousness_vaccinated, 
                           dur_V = dur_V,                                              
                           primary_doses = primary_doses,  
                           second_doses = second_doses,
                           booster_doses = booster_doses,                                             
                           vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
                           vaccine_booster_initial_coverage = vaccine_booster_initial_coverage)
  }
  
  ## Generating summary metrics for that model run
  summary_metrics <- calc_summary_metrics(mod_run)
  
  if (output == "summary") {
    return(list(
      model_arguments = list(scenario_index = scenario_index, population_size = population_size, standard_pop = standard_pop, country = country, hosp_bed_capacity = hosp_bed_capacity, ICU_bed_capacity = ICU_bed_capacity,
                             Rt = Rt, tt_Rt = tt_Rt, Tg = Tg, IFR = IFR, vaccine_scenario = vaccine_scenario, bpsv_start = bpsv_start, specific_vaccine_start = specific_vaccine_start,             
                             efficacy_infection_bpsv = efficacy_infection_bpsv, efficacy_disease_bpsv = efficacy_disease_bpsv, efficacy_infection_spec = efficacy_infection_spec,
                             efficacy_disease_spec = efficacy_disease_spec, dur_R = dur_R, dur_bpsv = dur_bpsv, dur_spec = dur_spec, bpsv_protection_delay = bpsv_protection_delay, specific_protection_delay = specific_protection_delay,                     
                             coverage_spec = coverage_spec, 
                             coverage_bpsv = if (approach == "old") coverage_spec else coverage_bpsv, 
                             vaccination_rate_bpsv = vaccination_rate_bpsv,
                             vaccination_rate_spec = vaccination_rate_spec,
                             min_age_group_index_priority = min_age_group_index_priority, min_age_group_index_non_priority = min_age_group_index_non_priority,     
                             runtime = runtime, seeding_cases = seeding_cases, detection_time = detection_time,
                             primary_doses = primary_doses, second_doses = second_doses, booster_doses = booster_doses,
                             vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage, 
                             vaccine_booster_initial_coverage = vaccine_booster_initial_coverage,
                             NPI_int = NPI_int, varied = varied),
      summary_metrics = summary_metrics)) 
  } else {
    return(list(
      model_output = mod_run,
      model_arguments = list(scenario_index = scenario_index, population_size = population_size, standard_pop = standard_pop, country = country, hosp_bed_capacity = hosp_bed_capacity, ICU_bed_capacity = ICU_bed_capacity,
                             Rt = Rt, tt_Rt = tt_Rt, Tg = Tg, IFR = IFR, vaccine_scenario = vaccine_scenario, bpsv_start = bpsv_start, specific_vaccine_start = specific_vaccine_start,             
                             efficacy_infection_bpsv = efficacy_infection_bpsv, efficacy_disease_bpsv = efficacy_disease_bpsv, efficacy_infection_spec = efficacy_infection_spec,
                             efficacy_disease_spec = efficacy_disease_spec, dur_R = dur_R, dur_bpsv = dur_bpsv, dur_spec = dur_spec, bpsv_protection_delay = bpsv_protection_delay, specific_protection_delay = specific_protection_delay,               
                             coverage_spec = coverage_spec, 
                             coverage_bpsv = if (approach == "old") coverage_spec else coverage_bpsv, 
                             vaccination_rate_bpsv = vaccination_rate_bpsv,
                             vaccination_rate_spec = vaccination_rate_spec,
                             min_age_group_index_priority = min_age_group_index_priority, min_age_group_index_non_priority = min_age_group_index_non_priority,     
                             runtime = runtime, seeding_cases = seeding_cases, detection_time = detection_time,
                             primary_doses = primary_doses, second_doses = second_doses, booster_doses = booster_doses,
                             vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage, 
                             vaccine_booster_initial_coverage = vaccine_booster_initial_coverage,
                             NPI_int = NPI_int, varied = varied),
      summary_metrics = summary_metrics)) 
  }
}

## Running the model and summarising outputs
# run_sars_x <- function(## Demographic Parameters
#   population_size = 1e6,
#   country = "Argentina",
#   
#   ## Healthcare-Related Parameters
#   hosp_bed_capacity = NULL,
#   ICU_bed_capacity = NULL,
#   
#   ## Epidemiological Parameters
#   Rt = 2.5,
#   tt_Rt = 0, 
#   Tg = 6.6,
#   IFR = 0.5,
#   
#   ## Vaccine-Related Parameters
#   vaccine_scenario = "specific_only",       # which scenario to explore
#   detection_time = 14,
#   bpsv_start = 50,                          # BPSV distribution start (days after detection time)
#   bpsv_protection_delay = 7,                # time between BPSV dose and protection arising
#   specific_vaccine_start = 250,             # specific vaccine distribution start (days after detection time)
#   specific_protection_delay = 7,            # time between specific dose and protection arising
#   efficacy_infection_bpsv = 0.35,           # vaccine efficacy against infection - BPSV
#   efficacy_disease_bpsv = 0.8,              # vaccine efficacy against disease - BPSV
#   efficacy_infection_spec = 0.55,           # vaccine efficacy against infection - specific vaccine
#   efficacy_disease_spec = 0.9,              # vaccine efficacy against disease - specific vaccine
#   dur_R = 1000 * 365,                       # duration of infection-induced immunity
#   dur_bpsv = 1000 * 365,                    # duration of vaccine-induced immunity for BPSV vaccine
#   dur_spec = 1000 * 365,                    # duration of vaccine-induced immunity for disease-specific vaccines
#   coverage_bpsv = 0.75,                     # proportion of the population vaccinated
#   coverage_spec = 0.75,                     # proportion of the population vaccinated
#   vaccination_rate_bpsv = 0.01,             # bpsv vaccination rate per week as percentage of population (note: percentage of whole population, despite restricted age-groups receiving it)
#   vaccination_rate_spec = 0.01,             # disease-specific vaccination rate per week as percentage of population
#   min_age_group_index_priority = 13,        # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
#   min_age_group_index_non_priority = 4,     # index of the youngest age group that *receives* vaccines (4 = 15+)
#   
#   ## Miscellaneous Parameters
#   runtime = 600,
#   approach = "old",
#   seeding_cases = 5,
#   output = "summary",
#   NPI_int = 0,
#   scenario_index = 0,
#   varied = "",
#   ...) {
#   
#   # print(scenario_index)
#   
#   set.seed(123)
#   
#   ## Input Checking
#   if (!approach %in% c("old", "new")) {
#     stop("parameter approach must be either 'old' or 'new'")
#   }
#   if (!output %in% c("summary", "full")) {
#     stop("parameter output must be either 'summary' or 'full'")
#   }
#   if (is.null(hosp_bed_capacity)) {
#     hosp_bed_capacity <- population_size
#   }
#   if (is.null(ICU_bed_capacity)) {
#     ICU_bed_capacity <- population_size
#   }
#   
#   ## Country-specific population estimate
#   standard_pop <- generate_standard_pop(country = country, population_size = population_size)
#   
#   ## Adapting Generation Time
#   TgVary <- scale_generation_time(target_Tg = Tg)
#   TgVary_dur_IMild <- TgVary$dur_IMild
#   TgVary_dur_ICase <- TgVary$dur_ICase
#   
#   ## Adjusting IFR
#   prob_hosp <- scale_IFR(country = country, population_size = population_size, target_IFR = IFR)
#   if (sum(prob_hosp > 1) > 0) {
#     stop("IFR Scaling Factor too high - results in prob_hosp > 1 for at least one age-group")
#   }
#   
#   ## Setting Up Vaccination Stuff
#   if (!vaccine_scenario %in% c("specific_only", "both_vaccines")) {
#     stop("parameter vaccine_scenario must be either 'specific_only' or 'both_vaccines'")
#   }
#   
#   if (vaccine_scenario == "specific_only") { ## in this scenario, first vaccine is virus-specific, for all ages
#     
#     ## Setting vaccination rates and coverages
#     priority_age_groups <- min_age_group_index_priority:17          # priority age groups  
#     vaccination_age_groups <- min_age_group_index_non_priority:17   # vaccinable age-groups
#     vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(coverage_spec, length(priority_age_groups)), 
#                                      rep(0, 17 - length(vaccination_age_groups)), rep(coverage_spec, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)
#     
#     ## Creating matrices of infection and disease efficacy for each age-group and vaccine received
#     ### Infection Efficacy
#     ve_i_spec <- c(efficacy_infection_spec, efficacy_infection_spec, efficacy_infection_spec, 0, 0, 0, 0) ## note that booster protection set to 0 as we're assuming virus-spec gives robust protection - no waning, so no booster required
#     vaccine_efficacy_infection <- matrix(c(rep(ve_i_spec, 17)), ncol = 17)
#     
#     ### Disease Efficacy
#     ve_d_spec <- c(efficacy_disease_spec, efficacy_disease_spec, efficacy_disease_spec, 0, 0, 0, 0) ## note that booster protection set to 0 as we're assuming virus-spec gives robust protection - no waning, so no booster required
#     vaccine_efficacy_disease <- matrix(c(rep(ve_d_spec, 17)), ncol = 17)
#     
#     ### Vaccine protection waning
#     dur_V <- matrix(data = c(rep(dur_spec, 17),   ## duration of primary series protection (specific vaccine in this scenario)
#                              rep(dur_spec, 17),   ## duration of primary series protection (specific vaccine in this scenario)
#                              rep(dur_spec, 17),   ## duration of booster protection (not used here)
#                              rep(dur_spec, 17)),  ## duration of booster protection (not used here)
#                     nrow = 4, ncol = 17, byrow = TRUE)
#     
#     # Set up dosing schedule - no boosters here
#     vaccine_doses <- create_vaccination_dose_series(country = country, 
#                                                     population_size = population_size, 
#                                                     detection_time = detection_time, 
#                                                     vaccine_scenario = vaccine_scenario,
#                                                     bpsv_start = bpsv_start, 
#                                                     bpsv_protection_delay = bpsv_protection_delay,
#                                                     specific_vaccine_start = specific_vaccine_start,
#                                                     specific_protection_delay = specific_protection_delay,
#                                                     vaccination_rate_bpsv = vaccination_rate_bpsv,
#                                                     vaccination_rate_spec = vaccination_rate_spec,
#                                                     coverage = coverage_spec,
#                                                     min_age_group_index_priority = min_age_group_index_priority,
#                                                     runtime = runtime)
#     primary_doses <- vaccine_doses$primary_doses
#     second_doses <- vaccine_doses$second_doses
#     booster_doses <- vaccine_doses$booster_doses
#     
#     ## Eligibility for Boosters - No One
#     vaccine_booster_follow_up_coverage <- rep(0, 17)
#     vaccine_booster_initial_coverage <- rep(0, 17)
#     
#   } else if (vaccine_scenario == "both_vaccines" & approach == "old") { ## in this scenario, first vaccine is BPSV for 60+, virus-specific for <60, and booster for 60+ is the virus-specific
#     
#     ## Setting vaccination rates and coverages
#     priority_age_groups <- min_age_group_index_priority:17          # priority age groups  
#     vaccination_age_groups <- min_age_group_index_non_priority:17   # vaccinable age-groups
#     vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(coverage_spec, length(priority_age_groups)), 
#                                      rep(0, 17 - length(vaccination_age_groups)),  rep(coverage_spec, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)
#     # matrix (col = age group, row = prioritisation step) of how to assign vaccines to different age-groups according to priority
#     
#     ## Calculating vaccine dosing schedules for this scenario
#     vaccine_doses <- create_vaccination_dose_series(country = country, 
#                                                     population_size = population_size, 
#                                                     detection_time = detection_time, 
#                                                     vaccine_scenario = vaccine_scenario,
#                                                     bpsv_start = bpsv_start, 
#                                                     bpsv_protection_delay = bpsv_protection_delay,
#                                                     specific_vaccine_start = specific_vaccine_start,
#                                                     specific_protection_delay = specific_protection_delay,
#                                                     vaccination_rate_bpsv = vaccination_rate_bpsv,
#                                                     vaccination_rate_spec = vaccination_rate_spec,
#                                                     coverage = coverage_spec,
#                                                     min_age_group_index_priority = min_age_group_index_priority,
#                                                     runtime = runtime)
#     primary_doses <- vaccine_doses$primary_doses
#     second_doses <- vaccine_doses$second_doses
#     booster_doses <- vaccine_doses$booster_doses
#     # note that we have the "+1"s in the underlying function because otherwise we don't *quite* get to vaccinating all of the elderly in
#     # the alloted time that Azra was previously calculating. I thought doing ceiling() would be enough and in most cases it is but not all.
#     # This way (i.e. by adding +1 to the times) we do! 
#     
#     ## Creating matrices of infection and disease efficacy for each age-group and vaccine received
#     ### Infection Efficacy
#     ve_bpsv <- list(infection = efficacy_infection_bpsv, disease = efficacy_disease_bpsv)
#     ve_spec <- list(infection = efficacy_infection_spec, disease = efficacy_disease_spec)
#     ve_i_elderly <- c(ve_bpsv$infection, ve_bpsv$infection, ve_bpsv$infection, 0, ve_spec$infection, ve_spec$infection, 0) 
#     ve_i_non_elderly <- c(ve_spec$infection, ve_spec$infection, ve_spec$infection, 0, ve_spec$infection, ve_spec$infection, 0)
#     vaccine_efficacy_infection <- matrix(c( rep(ve_i_non_elderly, 17 - length(priority_age_groups)), 
#                                             rep(ve_i_elderly, length(priority_age_groups)) ), ncol = 17)
#     ### Disease Efficacy
#     ve_d_elderly <- c(ve_bpsv$disease, ve_bpsv$disease, ve_bpsv$disease, 0, ve_spec$disease, ve_spec$disease, 0) 
#     ve_d_non_elderly <- c(ve_spec$disease, ve_spec$disease, ve_spec$disease, 0, ve_spec$disease, ve_spec$disease, 0)
#     vaccine_efficacy_disease <- matrix(c( rep(ve_d_non_elderly, 17 - length(priority_age_groups)), 
#                                           rep(ve_d_elderly, length(priority_age_groups)) ), ncol = 17)
#     
#     ### Vaccine protection waning
#     dur_V <- matrix(data = c(c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_bpsv, length(priority_age_groups))),    ## duration of primary series protection (bpsv for elderly, specific vaccine for everyone else in this scenario)
#                              c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_bpsv, length(priority_age_groups))),    ## duration of primary series protection (bpsv for elderly, specific vaccine for everyone else in this scenario)
#                              c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_spec, length(priority_age_groups))),           ## duration of booster protection (specific for elderly, not used for everyone else)
#                              c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_spec, length(priority_age_groups)))),         ## duration of booster protection (specific for elderly, not used for everyone else)
#                     nrow = 4, ncol = 17, byrow = TRUE)
#     
#     ## Eligibility for Boosters (Elderly Only)
#     vaccine_booster_follow_up_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(0, length(priority_age_groups))) # no booster follow ons
#     vaccine_booster_initial_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(1, length(priority_age_groups))) ## note that we set coverage to be 100% as it's already bounded by whatever coverage the 1st/2nd dose coverage got to
#     
#   } else if (vaccine_scenario == "both_vaccines" & approach == "new") { ## in this scenario, first vaccine is BPSV for 60+, virus-specific for <60, and booster for 60+ is the virus-specific
#     
#     ## Setting vaccination rates and coverages
#     priority_age_groups <- min_age_group_index_priority:17          # priority age groups  
#     vaccination_age_groups <- min_age_group_index_non_priority:17   # vaccinable age-groups
#     vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(coverage_spec, length(priority_age_groups)), 
#                                      rep(0, 17 - length(vaccination_age_groups)),  rep(coverage_spec, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)
#     # matrix (col = age group, row = prioritisation step) of how to assign vaccines to different age-groups according to priority
#     ### New structure is that we compress the duration of the bpsv campaign so that we only achieve bpsv_stockpile_coverage by the end of it - different coverage (applied to specific vaccine) is still applied to vaccine_coverage_matrix
#     ###   we then change the vaccine_efficacy matrix just before specific campaign starts, primary dose efficacy for elderly goes from bpsv -> specific
#     ###   we also turn off second doses so that everyone who gets vaccinated with specific primary doesn't get a second dose
#     ###   we then start the specific vaccination campaign and
#     ###      1) boost all the elderly that got vaccinated with bpsv
#     ###      2) give primary doses (now specific vaccine efficacy profile) to all elderly who *didn't* get the bpsv
#     ###      3) give primary doses to everyone else
#     ###   I think 2) and 3) should be handled automatically within the optimisation steps, but let's see
#     
#     ## Calculating vaccine dosing schedules for this scenario
#     vaccine_doses <- create_vaccination_dose_seriesNew(country = country, 
#                                                        population_size = population_size, 
#                                                        detection_time = detection_time, 
#                                                        vaccine_scenario = vaccine_scenario,
#                                                        bpsv_start = bpsv_start, 
#                                                        bpsv_protection_delay = bpsv_protection_delay,
#                                                        specific_vaccine_start = specific_vaccine_start,
#                                                        specific_protection_delay = specific_protection_delay,
#                                                        vaccination_rate_bpsv = vaccination_rate_bpsv,
#                                                        vaccination_rate_spec = vaccination_rate_spec,
#                                                        coverage_bpsv = coverage_bpsv,
#                                                        coverage_spec = coverage_spec,
#                                                        min_age_group_index_priority = min_age_group_index_priority,
#                                                        runtime = runtime)
#     primary_doses <- vaccine_doses$primary_doses
#     second_doses <- vaccine_doses$second_doses
#     booster_doses <- vaccine_doses$booster_doses
#     # note that we have the "+1"s in the underlying function because otherwise we don't *quite* get to vaccinating all of the elderly in
#     # the alloted time that Azra was previously calculating. I thought doing ceiling() would be enough and in most cases it is but not all.
#     # This way (i.e. by adding +1 to the times) we do! 
#     
#     ## Creating matrices of infection and disease efficacy for each age-group and vaccine received
#     ## Major change compared to previous iterations: 
#     ## During period of bpsv campaign, vaccine efficacy matrix has primary, secondary = bpsv and booster = specific for elderly population
#     ## During period of specific campaign, matrix has primary = specific, secondary = bpsv and booster = specific for elderly population
#     ## Required as some eligible elderly won't have received bpsv during initial campaign due to limited size of bpsv stockpile (coverage bpsv < coverage specific)
#     ## No second doses are distributed during the specific campaign. 
#     
#     ## Ordering in efficacy vector
#     ## 1 = first dose received
#     ## 2 = second dose received 1; 3 = second dose received 2
#     ## 4 = primary/secondary waned
#     ## 5 = booster 1; 6 = booster 2
#     ## 7 = booster waned
#     ve_bpsv <- list(infection = efficacy_infection_bpsv, disease = efficacy_disease_bpsv)
#     ve_spec <- list(infection = efficacy_infection_spec, disease = efficacy_disease_spec)
#     
#     ### Infection Efficacy
#     ve_i_elderly_bpsv_campaign <- c(ve_bpsv$infection, ve_bpsv$infection, ve_bpsv$infection, 0, ve_spec$infection, ve_spec$infection, 0) 
#     ve_i_elderly_spec_campaign <- c(ve_spec$infection, ve_bpsv$infection, ve_bpsv$infection, 0, ve_spec$infection, ve_spec$infection, 0) 
#     ve_i_non_elderly <- c(ve_spec$infection, ve_spec$infection, ve_spec$infection, 0, ve_spec$infection, ve_spec$infection, 0)
#     vaccine_efficacy_infection_bpsv_campaign <- matrix(c(rep(ve_i_non_elderly, 17 - length(priority_age_groups)), 
#                                                          rep(ve_i_elderly_bpsv_campaign, length(priority_age_groups)) ), ncol = 17)
#     vaccine_efficacy_infection_spec_campaign <- matrix(c(rep(ve_i_non_elderly, 17 - length(priority_age_groups)), 
#                                                          rep(ve_i_elderly_spec_campaign, length(priority_age_groups)) ), ncol = 17)
#     
#     ### Disease Efficacy
#     ve_d_elderly_bpsv_campaign <- c(ve_bpsv$disease, ve_bpsv$disease, ve_bpsv$disease, 0, ve_spec$disease, ve_spec$disease, 0) 
#     ve_d_elderly_spec_campaign <- c(ve_spec$disease, ve_bpsv$disease, ve_bpsv$disease, 0, ve_spec$disease, ve_spec$disease, 0) 
#     ve_d_non_elderly <- c(ve_spec$disease, ve_spec$disease, ve_spec$disease, 0, ve_spec$disease, ve_spec$disease, 0)
#     vaccine_efficacy_disease_bpsv_campaign <- matrix(c(rep(ve_d_non_elderly, 17 - length(priority_age_groups)), 
#                                                        rep(ve_d_elderly_bpsv_campaign, length(priority_age_groups)) ), ncol = 17)
#     vaccine_efficacy_disease_spec_campaign <- matrix(c(rep(ve_d_non_elderly, 17 - length(priority_age_groups)), 
#                                                        rep(ve_d_elderly_spec_campaign, length(priority_age_groups)) ), ncol = 17)
#     
#     ### Vaccine protection waning
#     #### NOTE THAT THIS MIGHT NEED CHANGING FOR THE FIGURE WHERE WE ALTER THE WANING RATE OF THE BPSV 
#     dur_V <- matrix(data = c(c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_bpsv, length(priority_age_groups))),    ## duration of primary series protection (bpsv for elderly, specific vaccine for everyone else in this scenario)
#                              c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_bpsv, length(priority_age_groups))),    ## duration of primary series protection (bpsv for elderly, specific vaccine for everyone else in this scenario)
#                              c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_spec, length(priority_age_groups))),           ## duration of booster protection (specific for elderly, not used for everyone else)
#                              c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_spec, length(priority_age_groups)))),         ## duration of booster protection (specific for elderly, not used for everyone else)
#                     nrow = 4, ncol = 17, byrow = TRUE)
#     
#     ## Eligibility for Boosters (Elderly Only)
#     vaccine_booster_follow_up_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(0, length(priority_age_groups))) # no booster follow ons
#     vaccine_booster_initial_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(1, length(priority_age_groups))) ## note that we set coverage to be 100% as it's already bounded by whatever coverage the 1st/2nd dose coverage got to
#   }
#   
#   ## Ensuring things that could be passed in as lists (when running in parallel) are coerced to right format
#   if (is.list(tt_Rt)) {
#     tt_Rt <- unlist(tt_Rt)
#   }
#   if (is.list(Rt)) {
#     Rt <- unlist(Rt)
#   }
#   
#   ## Running the Model 
#   mm <- squire::get_mixing_matrix(country = country)    
#   rel_infectiousness_vaccinated <- squire.page:::probs_booster$rel_infectiousness_vaccinated
#   rel_infectiousness_vaccinated[rel_infectiousness_vaccinated < 1] <- 1
#   
#   if (vaccine_scenario == "specific_only") {
#     mod_run <- run_booster(time_period = runtime,
#                            population = standard_pop,                                                 
#                            contact_matrix_set = mm,                                                   
#                            R0 = Rt,     
#                            tt_R0 = tt_Rt, 
#                            hosp_bed_capacity = hosp_bed_capacity,                                     
#                            ICU_bed_capacity = ICU_bed_capacity,                                       
#                            prob_hosp = prob_hosp,
#                            dur_IMild = TgVary_dur_IMild,
#                            dur_ICase = TgVary_dur_ICase,
#                            dur_R = dur_R,                                                        
#                            seeding_cases = seeding_cases,
#                            vaccine_coverage_mat = vaccine_coverage_mat,                               
#                            vaccine_efficacy_infection = vaccine_efficacy_infection,                   
#                            vaccine_efficacy_disease = vaccine_efficacy_disease,     
#                            rel_infectiousness_vaccinated = rel_infectiousness_vaccinated, 
#                            dur_V = dur_V,                                              
#                            primary_doses = primary_doses,  
#                            second_doses = second_doses,
#                            booster_doses = booster_doses,                                             
#                            vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
#                            vaccine_booster_initial_coverage = vaccine_booster_initial_coverage)
#   } else if (vaccine_scenario == "both_vaccines" & approach == "old") {
#     mod_run <- run_booster(time_period = runtime,
#                            population = standard_pop,                                                 
#                            contact_matrix_set = mm,                                                   
#                            R0 = Rt,     
#                            tt_R0 = tt_Rt, 
#                            hosp_bed_capacity = hosp_bed_capacity,                                     
#                            ICU_bed_capacity = ICU_bed_capacity,                                       
#                            prob_hosp = prob_hosp,
#                            dur_IMild = TgVary_dur_IMild,
#                            dur_ICase = TgVary_dur_ICase,
#                            dur_R = dur_R,                                                        
#                            seeding_cases = seeding_cases,
#                            vaccine_coverage_mat = vaccine_coverage_mat,                               
#                            vaccine_efficacy_infection = vaccine_efficacy_infection,                   
#                            vaccine_efficacy_disease = vaccine_efficacy_disease,     
#                            rel_infectiousness_vaccinated = rel_infectiousness_vaccinated, 
#                            dur_V = dur_V,                                              
#                            primary_doses = primary_doses,  
#                            second_doses = second_doses,
#                            booster_doses = booster_doses,                                             
#                            vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
#                            vaccine_booster_initial_coverage = vaccine_booster_initial_coverage)
#   } else if (vaccine_scenario == "both_vaccines" & approach == "new") {
#     mod_run <- run_booster(time_period = runtime,
#                            population = standard_pop,                                                 
#                            contact_matrix_set = mm,                                                   
#                            R0 = Rt,     
#                            tt_R0 = tt_Rt, 
#                            hosp_bed_capacity = hosp_bed_capacity,                                     
#                            ICU_bed_capacity = ICU_bed_capacity,                                       
#                            prob_hosp = prob_hosp,
#                            dur_IMild = TgVary_dur_IMild,
#                            dur_ICase = TgVary_dur_ICase,
#                            dur_R = dur_R,                                                        
#                            seeding_cases = seeding_cases,
#                            vaccine_coverage_mat = vaccine_coverage_mat,
#                            vaccine_efficacy_infection = list(vaccine_efficacy_infection_bpsv_campaign, vaccine_efficacy_infection_spec_campaign),
#                            tt_vaccine_efficacy_infection = c(0, specific_vaccine_start - 5), # ADD CHECK TO MAKE SURE THIS doesn't bleed into time when bpsv is being given out
#                            vaccine_efficacy_disease = list(vaccine_efficacy_disease_bpsv_campaign, vaccine_efficacy_disease_spec_campaign),
#                            tt_vaccine_efficacy_disease = c(0, specific_vaccine_start - 5), # ADD CHECK TO MAKE SURE THIS doesn't bleed into time when bpsv is being given out
#                            rel_infectiousness_vaccinated = rel_infectiousness_vaccinated, 
#                            dur_V = dur_V,                                              
#                            primary_doses = primary_doses,  
#                            second_doses = second_doses,
#                            booster_doses = booster_doses,                                             
#                            vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
#                            vaccine_booster_initial_coverage = vaccine_booster_initial_coverage)
#   }
#   
#   ## Generating summary metrics for that model run
#   summary_metrics <- calc_summary_metrics(mod_run)
#   
#   if (output == "summary") {
#     return(list(
#       model_arguments = list(scenario_index = scenario_index, population_size = population_size, standard_pop = standard_pop, country = country, hosp_bed_capacity = hosp_bed_capacity, ICU_bed_capacity = ICU_bed_capacity,
#                              Rt = Rt, tt_Rt = tt_Rt, Tg = Tg, IFR = IFR, vaccine_scenario = vaccine_scenario, bpsv_start = bpsv_start, specific_vaccine_start = specific_vaccine_start,             
#                              efficacy_infection_bpsv = efficacy_infection_bpsv, efficacy_disease_bpsv = efficacy_disease_bpsv, efficacy_infection_spec = efficacy_infection_spec,
#                              efficacy_disease_spec = efficacy_disease_spec, dur_R = dur_R, dur_bpsv = dur_bpsv, dur_spec = dur_spec, bpsv_protection_delay = bpsv_protection_delay, specific_protection_delay = specific_protection_delay,                     
#                              coverage_spec = coverage_spec, 
#                              coverage_bpsv = if (approach == "old") coverage_spec else coverage_bpsv, 
#                              vaccination_rate_bpsv = vaccination_rate_bpsv,
#                              vaccination_rate_spec = vaccination_rate_spec,
#                              min_age_group_index_priority = min_age_group_index_priority, min_age_group_index_non_priority = min_age_group_index_non_priority,     
#                              runtime = runtime, seeding_cases = seeding_cases, detection_time = detection_time,
#                              primary_doses = primary_doses, second_doses = second_doses, booster_doses = booster_doses,
#                              vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage, 
#                              vaccine_booster_initial_coverage = vaccine_booster_initial_coverage,
#                              NPI_int = NPI_int, varied = varied),
#       summary_metrics = summary_metrics)) 
#   } else {
#     return(list(
#       model_output = mod_run,
#       model_arguments = list(scenario_index = scenario_index, population_size = population_size, standard_pop = standard_pop, country = country, hosp_bed_capacity = hosp_bed_capacity, ICU_bed_capacity = ICU_bed_capacity,
#                              Rt = Rt, tt_Rt = tt_Rt, Tg = Tg, IFR = IFR, vaccine_scenario = vaccine_scenario, bpsv_start = bpsv_start, specific_vaccine_start = specific_vaccine_start,             
#                              efficacy_infection_bpsv = efficacy_infection_bpsv, efficacy_disease_bpsv = efficacy_disease_bpsv, efficacy_infection_spec = efficacy_infection_spec,
#                              efficacy_disease_spec = efficacy_disease_spec, dur_R = dur_R, dur_bpsv = dur_bpsv, dur_spec = dur_spec, bpsv_protection_delay = bpsv_protection_delay, specific_protection_delay = specific_protection_delay,               
#                              coverage_spec = coverage_spec, 
#                              coverage_bpsv = if (approach == "old") coverage_spec else coverage_bpsv, 
#                              vaccination_rate_bpsv = vaccination_rate_bpsv,
#                              vaccination_rate_spec = vaccination_rate_spec,
#                              min_age_group_index_priority = min_age_group_index_priority, min_age_group_index_non_priority = min_age_group_index_non_priority,     
#                              runtime = runtime, seeding_cases = seeding_cases, detection_time = detection_time,
#                              primary_doses = primary_doses, second_doses = second_doses, booster_doses = booster_doses,
#                              vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage, 
#                              vaccine_booster_initial_coverage = vaccine_booster_initial_coverage,
#                              NPI_int = NPI_int, varied = varied),
#       summary_metrics = summary_metrics)) 
#   }
# }