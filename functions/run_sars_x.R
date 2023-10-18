## Generate vaccination milestone timings based on parameter combinations
### For elderly, primary is bpsv during bpsv campaign and spec during spec campaign, secondary is always bpsv, and booster is always spec
### Note: This supersedes the old version in old_run_sars_x.R, which didn't have the ability to vary the coverage of the bpsv (a proxy for stockpile size)
vaccine_milestone_timings <- function(country, 
                                      population_size, 
                                      detection_time, 
                                      bpsv_start, 
                                      specific_vaccine_start,
                                      vaccination_rate_bpsv,
                                      vaccination_rate_spec,
                                      coverage_bpsv,
                                      coverage_spec,
                                      min_age_group_index_priority) {
  
  # Calculating time to coverage
  ## Note: This is used in the NPI function where we define coverage (for purpose of lifting NPIs) as based on elderly population so this is what is calculated
  standard_pop <- generate_standard_pop(country = country, population_size = population_size)
  priority_age_groups <- min_age_group_index_priority:17       
  
  # Calculating time to BPSV coverage for the elderly population
  daily_doses_bpsv <- vaccination_rate_bpsv * population_size / 7
  elderly_pop_to_vaccinate_bpsv <- ceiling(sum(standard_pop[priority_age_groups]) * coverage_bpsv) 
  time_to_coverage_bpsv <- ceiling(elderly_pop_to_vaccinate_bpsv/daily_doses_bpsv) + 1 ## might need to remove the +1s, TBD
  
  # Calculating time to disease-specific coverage for the elderly population
  # - Note that tis done separately for the elderly population who DID and DID NOT get the BPSV - with daily doses allocated proportionally 
  #   between these groups (i.e. not prioritising a particular elderly sub-group with the disease-specific vaccine)

  ## For the elderly who DID receive the BPSV initially (i.e. up to the limit allowed by the size of the stockpile)
  daily_doses_spec_received_bpsv <- (vaccination_rate_spec * population_size / 7) * coverage_bpsv
  elderly_pop_to_vaccinate_spec_received_bpsv <- ceiling(sum(standard_pop[priority_age_groups]) * coverage_bpsv)
  time_to_coverage_spec_received_bpsv <- ceiling(elderly_pop_to_vaccinate_spec_received_bpsv/daily_doses_spec_received_bpsv) + 1 # +1 added to ensure we reach the coverage targets in non-integer instances
  
  ## For the elderly who DID NOT receive the BPSV initially (despite wanting to i.e. limited by stockpile size, not general willingness)
  ### (coverage_spec - coverage_bpsv)  - proportion who would have got bpsv if stockpile had been big enough (i.e. are willing to be vaccinated)
  daily_doses_spec_no_bpsv <- (vaccination_rate_spec * population_size / 7) * (coverage_spec - coverage_bpsv)
  elderly_pop_to_vaccinate_spec_no_bpsv <- ceiling(sum(standard_pop[priority_age_groups]) * (coverage_spec - coverage_bpsv))
  time_to_coverage_spec_no_bpsv <- ifelse(elderly_pop_to_vaccinate_spec_no_bpsv == 0, 0, ceiling(elderly_pop_to_vaccinate_spec_no_bpsv/daily_doses_spec_no_bpsv) + 1) # +1 added to ensure we reach the coverage targets in non-integer instances
  ## TO DO: FOR REASONS LAID OUT IN CREATE_DOSE_SERIES, NOT SURE THIS IS *QUITE* RIGHT JUST YET
  ##        AS NEED TO ACCOUNT PRIMARY VS BOOSTER DISEASE-SPEC GROUPS FINISHING FIRST (AND HENCE THE EXTRA DOSES BEING AVAILABLE FOR THE OTHER GROUP, THUS SPEEDING UP COMPLETION)
  ## OR MAYBE NONE OF THIS MATTERS BECAUSE IT'S ALL RATIOS AND PROPORTIONAL TO THE POPULATION???

  # Generate and return dataframe of model parameters and associated key vaccination milestone timings
  vaccine_milestone_timings_df <- data.frame(country = country,
                                             population_size = population_size,
                                             detection_time = detection_time,
                                             bpsv_start = bpsv_start,
                                             specific_vaccine_start = specific_vaccine_start,
                                             vaccination_rate_bpsv = vaccination_rate_bpsv, 
                                             vaccination_rate_spec = vaccination_rate_spec, 
                                             coverage_bpsv = coverage_bpsv,
                                             coverage_spec = coverage_spec,
                                             min_age_group_index_priority = min_age_group_index_priority,
                                             time_to_coverage_bpsv = time_to_coverage_bpsv,
                                             time_to_coverage_spec = max(c(time_to_coverage_spec_received_bpsv, time_to_coverage_spec_no_bpsv))) # whichever is later is when all elderly are vaccinated with disease-specific vaccine
  return(vaccine_milestone_timings_df)
}

## Create vaccination dose series
### For elderly, when bpsv is available primary is bpsv during bpsv campaign and spec during spec campaign, secondary is always bpsv, and booster is always spec
### Note: This supersedes the old version in old_run_sars_x.R, which didn't have the ability to vary the coverage of the bpsv (a proxy for stockpile size)
create_vaccination_dose_series <- function(country, 
                                           population_size, 
                                           detection_time, 
                                           vaccine_scenario,
                                           bpsv_start, 
                                           bpsv_protection_delay,
                                           specific_vaccine_start,
                                           specific_protection_delay,
                                           vaccination_rate_bpsv,
                                           vaccination_rate_spec,
                                           coverage_bpsv,
                                           coverage_spec,
                                           min_age_group_index_priority,
                                           runtime) {
  
  # Checking that the vaccine scenario and coverage targets arecorrectly specified
  if (!vaccine_scenario %in% c("specific_only", "both_vaccines")) {
    stop("parameter vaccine_scenario must be either 'specific_only' or 'both_vaccines'")
  }
  if (coverage_spec < coverage_bpsv) {
    step("coverage with the disease-specific vaccine cannot be lower than coverage with the BPSV vaccine - adjust model inputs accordingly")
  }
  
  # If the vaccine scenario only involves the disease-specific vaccine, then primary doses are the disease-specific vaccine for everyone
  if (vaccine_scenario == "specific_only") {
    
    # Primary Doses - Calculating daily number of doses available given a vaccination rate and population size
    ## Disease-specific vaccine becomes available specific_vaccine_start days after pathogen detected - we bake in the specific_protection_delay
    ## directly into the delivery of doses (in contrast to prior approaches where this was in squire.page as a parameter)
    daily_doses_spec <- vaccination_rate_spec * population_size / 7    # rate of vaccination with primary series
    primary_doses <- 
      c(rep(0, detection_time),
        rep(0, specific_vaccine_start),
        rep(0, specific_protection_delay),        
        rep(daily_doses_spec, runtime - detection_time - specific_vaccine_start - specific_protection_delay))
    
    # Second Doses
    ## We model full protection as arising the moment you've had the first dose, so second dose is redundant. 
    ## We arbitrarily set it to 1 here and therefore stagger ## second doses to be a day after primary doses.
    second_doses <- c(0, primary_doses[1:(length(primary_doses) - 1)]) # second dose 1 day after first
    
    ## Booster Doses (none here as the disease-specific vaccine is the only one available; and the primary doses are the disease-specific vaccine in this scenario)
    booster_doses <- rep(0, runtime)
  
  # If the vaccine scenario involves both vaccines, 60+s receive BPSV and disease-specific; under 60s receive just disease-specific
  # - In the model structure, for 60+ age-groups, primary is BPSV initially (during BPSV campaign) and then primary becomes disease-specific (during disease-specific campaign).
  #   Secondary doses are always BPSV and we only have them in BPSV campaign (to get BPSV-vaxxed elderly into second dose compartments), and then when
  #   disease-specific campaign occurs, we distribute it via primary doses (to those who didn't get BPSV-vaxxed) AND booster doses (to those who did).
  #   This is all for 60+; for all other age-groups (who don't get BPSV), primary doses are the disease-specific vaccine.
  } else {
    
    # Calculating daily number of doses available for both vaccines, and the associated time to coverage (of elderly population) from parameter inputs
    standard_pop <- generate_standard_pop(country = country, population_size = population_size)
    priority_age_groups <- min_age_group_index_priority:17       
    
    # Daily number of bpsv doses available and associated time to vaccine elderly population  
    daily_doses_bpsv <- vaccination_rate_bpsv * population_size / 7
    elderly_pop_to_vaccinate_bpsv <- ceiling(sum(standard_pop[priority_age_groups]) * coverage_bpsv)
    time_to_coverage_bpsv <- ceiling(elderly_pop_to_vaccinate_bpsv/daily_doses_bpsv) + 1 # +1 added to ensure we reach the coverage targets in non-integer instances 
    
    # Daily number of disease-specific doses available and associated time to vaccine elderly population
    # - Note that tis done separately for the elderly population who DID and DID NOT get the BPSV - with daily doses allocated proportionally 
    #   between these groups (i.e. not prioritising a particular elderly sub-group with the disease-specific vaccine)
    daily_doses_spec <- vaccination_rate_spec * population_size / 7 ## to general population
    
    ## To those elderly who DID receive the BPSV initially (i.e. up to the limit allowed by the size of the stockpile)
    daily_doses_spec_received_bpsv <- (vaccination_rate_spec * population_size / 7) * coverage_bpsv
    elderly_pop_to_vaccinate_spec_received_bpsv <- ceiling(sum(standard_pop[priority_age_groups]) * coverage_bpsv)
    time_to_coverage_spec_received_bpsv <- ceiling(elderly_pop_to_vaccinate_spec_received_bpsv/daily_doses_spec_received_bpsv) + 1 # +1 added to ensure we reach the coverage targets in non-integer instances

    ## To those elderly who DID NOT receive the BPSV initially (despite wanting to i.e. limited by stockpile size, not general willingness)
    ### coverage_bpsv                    - proportion receiving bpsv 
    ### (1 - coverage_bpsv)              - proportion not receiving bpsv
    ### (coverage_spec - coverage_bpsv)  - proportion who would have got bpsv if stockpile had been big enough (i.e. are willing to be vaccinated)
    daily_doses_spec_no_bpsv <- (vaccination_rate_spec * population_size / 7) * (coverage_spec - coverage_bpsv)
    elderly_pop_to_vaccinate_spec_no_bpsv <- ceiling(sum(standard_pop[priority_age_groups]) * (coverage_spec - coverage_bpsv))
    time_to_coverage_spec_no_bpsv <- ifelse(elderly_pop_to_vaccinate_spec_no_bpsv == 0, 0, ceiling(elderly_pop_to_vaccinate_spec_no_bpsv/daily_doses_spec_no_bpsv) + 1) # +1 added to ensure we reach the coverage targets in non-integer instances
    
    # Checking there isn't any temporal overlap in the BPSV and disease-specific vaccination campaigns
    # - We want to avoid BPSV elderly and disease-specific everyone else campaigns happening concurrently (as they're the same series in the model)
    #   and because we switch what the "primary doses" are (BPSV -> disease-specific) for elderly during the model runs. Making sure the BPSV campaign 
    #   has to be complete before disease-specific campaign begins (which is before the disease-specific campaign for everyone else) I think achieves this.
    minimum_spec_development_time_allowed <- (bpsv_start + time_to_coverage_bpsv + bpsv_protection_delay)
    if (minimum_spec_development_time_allowed > specific_vaccine_start + specific_protection_delay + 7) { # +7 in there just to make it extra conservative and ensure no overlap
      stop("Disease-specific vaccine developed too soon given timing of BPSV campaign and speed of vaccination campaign - adjust model inputs accordingly")
    }
    
    # Primary Doses (These are BPSV for elderly during the BPSV campaign, then disease-specific for everyone during the disease-specific vaccination campaign)
    # - Note that we that for the disease-specific vaccine series, we initially use daily_doses_spec_no_bpsv (to reflect doses being split across this elderly group
    #   and the elderly group receiving the disease-specific vaccine as a booster); and then switch over to daily_doses_spec (to reflect the fact that elderly population
    #   vaccinated with disease specific is completed; and now the full daily_doses_spec is available for just the disease-specific vaccine population receiving it as a primary dose
    #   (i.e. all the younger people who didn't receive a BPSV).
    # - TO DO: Check whether +1 is needed to be added to time_to_coverage_bpsv to sure elderly vaccination targets are achieved
    primary_doses <- 
      c(rep(0, detection_time),                         ## time between epidemic start and detection
        rep(0, bpsv_start),                             ## time between detection and initiation of BPSV campaign
        rep(0, bpsv_protection_delay),                  ## time between initiation of BPSV campaign and people first being protected by that first dose
        rep(daily_doses_bpsv, time_to_coverage_bpsv),   ## protection (if any) emerges in BPSV-vaccinated primary vaccinated folks
        rep(0, specific_vaccine_start - time_to_coverage_bpsv - bpsv_protection_delay - bpsv_start), # specific vaccine becomes available specific_vaccine_start days after detection
        rep(0, specific_protection_delay),                               ## time between initiation of specific vaccine campaign for elderly and them being protected by that vaccine
        rep(daily_doses_spec_no_bpsv, time_to_coverage_spec_no_bpsv),    ## no specific vaccine for non-elderly whilst that vaccination campaign is ongoing
        rep(daily_doses_spec, runtime - time_to_coverage_spec_no_bpsv - specific_protection_delay - specific_vaccine_start - detection_time)) # specific vaccination of all other ages until end of runtime
    ## still not quite sure this is right
    ## time_to_coverage_spec_no_bpsv = 0 won't feature in the vector due to nature of how rep works, so that's fine
    ## but would that not mean the non-elderly disease-specific vaccine starts too soon?
    ## does time_to_coverage_spec_no_bpsv need to be time_to_coverage_spec_received_bpsv or like max(time_to_coverage_spec_no_bpsv, time_to_coverage_spec_no_bpsv) or something?
    ## unclear currently. Need to check.
    ### yeah, think we're still doing more doses than we should be currently
    ### I think it needs to be something like do daily_doses_spec_no_bpsv initially, and then if this finishes early, do daily_doses_spec in the booster elderly population for the remaining
    ### (time_to_coverage_spec_received_bpsv - time_to_coverage_spec_no_bpsv) days (and vice versa). And only then do you start the non-elderly population campaign.
    ### Try implementing this when you're back from your break.
    ## OR MAYBE NONE OF THIS MATTERS BECAUSE IT'S ALL RATIOS AND PROPORTIONAL TO THE POPULATION??? AND SO THEY SHOULD ALWAYS FINISH AT THE SAME TIME ???
    ## YEAH I THINK TIME TO COVERAGE SHOULD BE IDENTICAL, SO time_to_coverage_spec_no_bpsv = time_to_coverage_spec_received_bpsv which means the above is right (I think).
    
    # Second Doses (Only elderly who got BPSV vaccination receive second doses, and they are always the BPSV within the modelling framework)
    # - We model full protection as arising the moment you've had the first dose (in the model) and incorporate protection delays outside the model, as in the above and the below.
    #   The role of the second doses is to transfer BPSV vaccinated elderly into the second dose compartments, to free up the first dose compartments to switch over to being the disease-specific
    #   vaccine for the elderly (required to handle different BPSV and disease-spcific coverage targets; and to make those who did receive the BPSV eligible for a booster, which will be the 
    #   disease-specific vaccine for them). Given this, the second dose timing does not affect model outputs and we therefore arbitrarily set it to 1 here and therefore stagger
    #   second doses to be a day after primary doses.
    second_doses <- 
      c(rep(0, 1),                                  ## arbitrary two day delay between primary and secondary doses
        rep(0, detection_time),                     ## time between epidemic start and detection
        rep(0, bpsv_start),                         ## time between detection and initiation of BPSV campaign
        rep(0, bpsv_protection_delay),              ## time between initiation of BPSV campaign and people first being protected by that first dose
        rep(daily_doses_bpsv, time_to_coverage_bpsv),    ## protection (if any) emerges in BPSV-vaccinated primary vaccinated folks
        rep(0, specific_vaccine_start - time_to_coverage_bpsv - bpsv_protection_delay - bpsv_start), # specific vaccine campaign starts specific_vaccine_start days after detection
        rep(0, specific_protection_delay),          ## time between initiation of specific vaccine campaign and people first being protected by that first dose
        rep(0, time_to_coverage_spec_no_bpsv),    ## no specific vaccine for non-elderly whilst that vaccination campaign is ongoing
        rep(0, runtime - time_to_coverage_spec_no_bpsv - specific_protection_delay - specific_vaccine_start - detection_time)) # specific vaccination of all other ages until end of runtime
    second_doses <- second_doses[1:length(primary_doses)]
    second_doses[(detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+1)] <- daily_doses_bpsv # extra secondary doses to make sure everyone receiving primary goes to secondary
    second_doses[(detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+2)] <- daily_doses_bpsv # we're modelling no additional impact of secondary dose (i.e. you're protected)
    second_doses[(detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+3)] <- daily_doses_bpsv # when you are vaccinated with primary dose. But we do need all the primary dose
                                                                                                                  # elderly folks in the secondary dose compartments when the specific vaccine comes
                                                                                                                  # round, so we can switch the primary to be the specific vaccine (rather than bpsv as it is initially)

    # Booster Doses (only for the elderly population who were vaccinated with the BPSV, this is the disease-specific vaccine)
    # - Because we now have booster coverage eligibility indicator in, we don't have to worry about stopping the booster doses - only elderly will ever get it
    booster_doses <- 
      c(rep(0, detection_time),
        rep(0, specific_vaccine_start),
        rep(0, specific_protection_delay),
        rep(daily_doses_spec_received_bpsv, runtime - specific_protection_delay - specific_vaccine_start - detection_time))
    
  }
  
  return(list(primary_doses = primary_doses, second_doses = second_doses, booster_doses = booster_doses))
}

## Generate default NPI scenarios based on model arguments (vaccination rate, detection time, vaccination start times etc)
### 1 = NPIs reducing Rt < 1 until BPSV campaign is done, followed by minimum mandate until specific campaign is done, then full release
### 2 = NPIs reducing Rt < 1 until BPSV campaign is done, followed by full release
### 3 = NPIs reducing Rt < 1 until spec vaccine campaign is done, followed by full release
### 4 = NPIs minimal mandate until BPSV campaign is done, followed by full release
### 5 = NPIs minimal mandate until spec vaccine campaign is done, followed by full release
### 6 = NPIs reducing Rt < 1 until BPSV campaign is done, followed by gradual release until spec vacc campaign done
### 7 = NPIs minimal mandate until BPSV campaign is done, followed by gradual release until spec vacc campaign done
### 8 = NPIs reducing Rt < 1 until BPSV campaign is done, followed by gradual release from minimal mandate until specific campaign is done
### 9 = No NPIs
default_NPI_scenarios <- function(lockdown_Rt = 0.9,
                                  minimal_mandate_reduction = 0.25,
                                  NPI_scenarios = 1:9,
                                  scenarios) { 
  
  # Unpack the unique parameter combinations in the scenarios dataframe and calculate vaccination milestone timings
  timings <- expand_grid(country = unique(scenarios$country),
                         population_size = unique(scenarios$population_size), 
                         detection_time = unique(scenarios$detection_time),
                         bpsv_start = unique(scenarios$bpsv_start), 
                         specific_vaccine_start = unique(scenarios$specific_vaccine_start),
                         vaccination_rate_bpsv = unique(scenarios$vaccination_rate_bpsv),
                         vaccination_rate_spec = unique(scenarios$vaccination_rate_spec),
                         coverage_bpsv = unique(scenarios$coverage_bpsv),
                         coverage_spec = unique(scenarios$coverage_spec),
                         min_age_group_index_priority = unique(scenarios$min_age_group_index_priority))
  vacc_timings <- pmap(timings, vaccine_milestone_timings)
  vacc_timings <- rbindlist(vacc_timings)
  
  # Checking specific vaccinate starts after bpsv vaccination is completed
  time_diff <- (vacc_timings$bpsv_start + vacc_timings$time_to_coverage_bpsv) - vacc_timings$specific_vaccine_start
  if (sum(time_diff >= 0) > 0) {
    stop("Specific Vaccine start time happens before coverage with BPSV vaccine is finished. Change this.")
  }  
  
  # Generating NPI Scenarios
  NPIs <- expand_grid(vacc_timings, 
                      R0 = unique(scenarios$R0),
                      NPI_int = NPI_scenarios) %>%
    rowwise() %>%
    mutate(temp = ((specific_vaccine_start + time_to_coverage_spec) - (bpsv_start + time_to_coverage_bpsv))) %>%
    mutate(Rt = case_when(NPI_int == 1 ~ list(c(R0, lockdown_Rt, R0 * (1 - minimal_mandate_reduction), R0)),
                          NPI_int == 2 ~ list(c(R0, lockdown_Rt, R0)),
                          NPI_int == 3 ~ list(c(R0, lockdown_Rt, R0)),
                          NPI_int == 4 ~ list(c(R0, R0 * (1 - minimal_mandate_reduction), R0)),
                          NPI_int == 5 ~ list(c(R0, R0 * (1 - minimal_mandate_reduction), R0)),
                          NPI_int == 6 ~ list(c(R0, lockdown_Rt, seq(from = lockdown_Rt, to = R0, length.out = temp), R0)), 
                          NPI_int == 7 ~ list(c(R0, R0 * (1 - minimal_mandate_reduction), seq(from = R0 * (1 - minimal_mandate_reduction), to = R0, length.out = temp), R0)),  
                          NPI_int == 8 ~ list(c(R0, lockdown_Rt, seq(from = R0 * (1 - minimal_mandate_reduction), to = R0, length.out = temp), R0)), 
                          NPI_int == 9 ~ list(R0))) %>%
    mutate(tt_Rt = case_when(NPI_int == 1 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv, detection_time + specific_vaccine_start + time_to_coverage_spec)),
                             NPI_int == 2 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv)),
                             NPI_int == 3 ~ list(c(0, detection_time, detection_time + specific_vaccine_start + time_to_coverage_spec)),
                             NPI_int == 4 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv)),
                             NPI_int == 5 ~ list(c(0, detection_time, detection_time + specific_vaccine_start + time_to_coverage_spec)),
                             NPI_int == 6 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv, seq(from = detection_time + bpsv_start + time_to_coverage_bpsv, to = detection_time + specific_vaccine_start + time_to_coverage_spec - 1))), # possible move the -1 to the first detection_time + bpsv_start + time_to_coverage_bpsv
                             NPI_int == 7 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv, seq(from = detection_time + bpsv_start + time_to_coverage_bpsv, to = detection_time + specific_vaccine_start + time_to_coverage_spec - 1))), # possible move the -1 to the first detection_time + bpsv_start + time_to_coverage_bpsv
                             NPI_int == 8 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv, seq(from = detection_time + bpsv_start + time_to_coverage_bpsv, to = detection_time + specific_vaccine_start + time_to_coverage_spec - 1))), # possible move the -1 to the first detection_time + bpsv_start + time_to_coverage_bpsv
                             NPI_int == 9 ~ list(0))) %>%
    ### check with the lines above whether it's appropriate to have the "-1" at the end for NPIs 6, 7, 8
    select(-temp)
  
  return(NPIs)
}

## Running the model and summarising outputs
### For elderly, primary is bpsv during bpsv campaign and spec during spec campaign, secondary is always bpsv, and booster is always spec
### Note: This supersedes the old version in old_run_sars_x.R, which didn't have the ability to vary the coverage of the bpsv (a proxy for stockpile size)
run_sars_x <- function(## Demographic Parameters
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
  
  ## Sedding the seed
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
    
  }  else if (vaccine_scenario == "both_vaccines" & approach == "new") { ## in this scenario, first vaccine is BPSV for 60+ during bpsv campaign, then virus-specific for 60+ during spec campaign, virus-specific for <60, and booster for 60+ is the virus-specific
    
    ## Setting vaccination rates and coverages
    priority_age_groups <- min_age_group_index_priority:17          # priority age groups  
    vaccination_age_groups <- min_age_group_index_non_priority:17   # vaccinable age-groups
    vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(coverage_spec, length(priority_age_groups)), 
                                     rep(0, 17 - length(vaccination_age_groups)),  rep(coverage_spec, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)
    # matrix (col = age group, row = prioritisation step) of how to assign vaccines to different age-groups according to priority
    ### New structure is that we compress the duration of the bpsv campaign so that we only achieve bpsv_stockpile_coverage by the end of it - different coverage (applied to specific vaccine) is still applied to vaccine_coverage_matrix
    ###   we then change the vaccine_efficacy matrix just before specific campaign starts, primary dose efficacy for elderly goes from bpsv -> specific
    ###   we also turn off second doses so that everyone who gets vaccinated with specific primary doesn't get a second dose
    ###   we then start the specific vaccination campaign and
    ###      1) boost all the elderly that got vaccinated with bpsv
    ###      2) give primary doses (now specific vaccine efficacy profile) to all elderly who *didn't* get the bpsv
    ###      3) give primary doses to everyone else
    ###   I think 2) and 3) should be handled automatically within the optimisation steps, but let's see
    
    ## Calculating vaccine dosing schedules for this scenario
    vaccine_doses <- create_vaccination_dose_seriesNew(country = country, 
                                                       population_size = population_size, 
                                                       detection_time = detection_time, 
                                                       vaccine_scenario = vaccine_scenario,
                                                       bpsv_start = bpsv_start, 
                                                       bpsv_protection_delay = bpsv_protection_delay,
                                                       specific_vaccine_start = specific_vaccine_start,
                                                       specific_protection_delay = specific_protection_delay,
                                                       vaccination_rate_bpsv = vaccination_rate_bpsv,
                                                       vaccination_rate_spec = vaccination_rate_spec,
                                                       coverage_bpsv = coverage_bpsv,
                                                       coverage_spec = coverage_spec,
                                                       min_age_group_index_priority = min_age_group_index_priority,
                                                       runtime = runtime)
    primary_doses <- vaccine_doses$primary_doses
    second_doses <- vaccine_doses$second_doses
    booster_doses <- vaccine_doses$booster_doses
    # note that we have the "+1"s in the underlying function because otherwise we don't *quite* get to vaccinating all of the elderly in
    # the alloted time that Azra was previously calculating. I thought doing ceiling() would be enough and in most cases it is but not all.
    # This way (i.e. by adding +1 to the times) we do! 
    
    ## Creating matrices of infection and disease efficacy for each age-group and vaccine received
    ## Major change compared to previous iterations: 
    ## During period of bpsv campaign, vaccine efficacy matrix has primary, secondary = bpsv and booster = specific for elderly population
    ## During period of specific campaign, matrix has primary = specific, secondary = bpsv and booster = specific for elderly population
    ## Required as some eligible elderly won't have received bpsv during initial campaign due to limited size of bpsv stockpile (coverage bpsv < coverage specific)
    ## No second doses are distributed during the specific campaign. 
    
    ## Ordering in efficacy vector
    ## 1 = first dose received
    ## 2 = second dose received 1; 3 = second dose received 2
    ## 4 = primary/secondary waned
    ## 5 = booster 1; 6 = booster 2
    ## 7 = booster waned
    ve_bpsv <- list(infection = efficacy_infection_bpsv, disease = efficacy_disease_bpsv)
    ve_spec <- list(infection = efficacy_infection_spec, disease = efficacy_disease_spec)
    
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
    
    ### Vaccine protection waning
    #### NOTE THAT THIS MIGHT NEED CHANGING FOR THE FIGURE WHERE WE ALTER THE WANING RATE OF THE BPSV 
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
                           vaccine_efficacy_infection = list(vaccine_efficacy_infection_bpsv_campaign, vaccine_efficacy_infection_spec_campaign),
                           tt_vaccine_efficacy_infection = c(0, specific_vaccine_start - 5), # ADD CHECK TO MAKE SURE THIS doesn't bleed into time when bpsv is being given out
                           vaccine_efficacy_disease = list(vaccine_efficacy_disease_bpsv_campaign, vaccine_efficacy_disease_spec_campaign),
                           tt_vaccine_efficacy_disease = c(0, specific_vaccine_start - 5), # ADD CHECK TO MAKE SURE THIS doesn't bleed into time when bpsv is being given out
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

## Consider changing summary metrics to accommodate not finishing with same R0
## Be wary of when runtime is shorter than some of the numbers automatically spat out when creating tt_Rt 
## Calculate summary metrics (deaths and time under NPIs) for individual model runs
### Note: Breaks when runtime < final tt_Rt value.
### TO DO: Check how this works with secondary country where NPIs might be implemented before the pathogen arrives (i.e. first value in Rt vector < R0)
### TO DO: Check if it works when last value in Rt vector < R0 as well.
calc_summary_metrics <- function(model_output) { # summary metrics = total deaths, time under any NPIs, composite NPI function, 
  
  ## Calculating Deaths
  check <- nimue::format(model_output, compartments = "D", summaries = "deaths") %>%
    filter(t > 1, compartment == "deaths")
  deaths <- sum(check$value)
  
  ## Time Under NPIs
  R0 <- max(model_output$parameters$R0) # changed from model_output$parameters$R0[1] - check whether this is alright
  which_NPI <-  which(model_output$parameters$R0 < R0)
  NPI_times <- model_output$parameters$tt_R0[c(which_NPI, max(which_NPI) + 1)]
  if (length(model_output$parameters$R0) == 1) {
    time_under_NPIs <- 0
  } else {
    time_under_NPIs <- max(NPI_times) - min(NPI_times)
  }

  ## Composite Measure of NPI Stringency and Time
  if (length(model_output$parameters$R0) == 1) {
    composite <- 0
  } else {
    composite <- 0
    rel_stringency <- (1 - model_output$parameters$R0 / R0)[which_NPI]
    for (i in 1:length(rel_stringency)) {
      composite_measure <- (NPI_times[i+1] - NPI_times[i]) * rel_stringency[i]
      composite <- composite + composite_measure
    }
  } 
  
  return(list(deaths = deaths,
              time_under_NPIs = time_under_NPIs,
              composite_NPI = composite))
}
