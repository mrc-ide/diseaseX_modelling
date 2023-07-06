## Identify number of columns with >1 unique values in a dataframe
variable_columns <- function(df) {
  result <- sapply(df, FUN = function(x) {
    length(unique(x)) > 1
  })
  columns <- names(which(result))
  columns <- columns[!(columns %in% c("index", "vaccine_scenario", "deaths", "time_under_NPIs", "composite_NPI"))]
  return(columns)
}

## Generate standardised population size from raw demography
generate_standard_pop <- function(country = "Argentina", population_size = 1e6) {
  raw_pop <- squire::get_population(country = country)$n                   
  standard_pop <- round(raw_pop * population_size / sum(raw_pop)) 
  return(standard_pop)
}

## Generate generation time 
scale_generation_time <- function(target_Tg) {
  current_Tg <- squire.page:::durs_booster$dur_IMild + squire.page:::durs_booster$dur_ICase
  Tg_ratio <- target_Tg / current_Tg
  TgVary_dur_IMild <- Tg_ratio * squire.page:::durs_booster$dur_IMild
  TgVary_dur_ICase <- Tg_ratio * squire.page:::durs_booster$dur_ICase
  return(list(dur_IMild = TgVary_dur_IMild, dur_ICase = TgVary_dur_ICase))
}

## Scale IFR (needs checking/changing)
scale_IFR <- function(country, population_size, target_IFR) {
  
  ## Get mixing matrix for country and use to calculate weighting for infections/deaths
  standard_pop <- generate_standard_pop(country = country, population_size = population_size)
  mm <- squire::get_mixing_matrix(country = country)    
  contact_rates <- apply(mm, 1, sum)
  contact_rates <- c(contact_rates, "80+" = unname(contact_rates[16]))
  pop_contact_weighting <- standard_pop * contact_rates
  
  ## Calculate approx IFR based on demography 
  non_severe_deaths <- nimue:::probs$prob_hosp * (1 - nimue:::probs$prob_severe) * nimue:::probs$prob_non_severe_death_treatment
  severe_deaths <- nimue:::probs$prob_hosp * nimue:::probs$prob_severe * nimue:::probs$prob_severe_death_treatment
  raw_IFR <- 100 * sum(((non_severe_deaths + severe_deaths) * pop_contact_weighting/sum(pop_contact_weighting))) 
  IFR_scaling_factor <- target_IFR / raw_IFR
  prob_hosp <- IFR_scaling_factor * nimue:::probs$prob_hosp
  return(prob_hosp)
}

## Generate vaccination milestone timings based on parameter combinations
vaccine_milestone_timings <- function(country, 
                                      population_size, 
                                      detection_time, 
                                      bpsv_start, 
                                      specific_vaccine_start,
                                      vaccination_rate,
                                      coverage,
                                      min_age_group_index_priority) {
  
  # Calculating time to coverage from parameter inputs
  standard_pop <- generate_standard_pop(country = country, population_size = population_size)
  daily_doses <- vaccination_rate * population_size / 7    # rate of vaccination with primary series
  priority_age_groups <- min_age_group_index_priority:17            
  elderly_pop_to_vaccinate <- sum(standard_pop[priority_age_groups]) * coverage # 60+s receive primary (BNPCV) and booster (diseaseX-specific); under 60s receive just primary (diseaseX-specific)
  time_to_coverage_bpsv <- ceiling(elderly_pop_to_vaccinate/daily_doses) + 1
  time_to_coverage_spec <- time_to_coverage_bpsv
  
  # Returning dataframe of parameters and key vaccination milestone timings
  temp <- data.frame(country = country,
                     population_size = population_size,
                     detection_time = detection_time,
                     bpsv_start = bpsv_start,
                     specific_vaccine_start = specific_vaccine_start,
                     vaccination_rate = vaccination_rate, 
                     coverage = coverage,
                     min_age_group_index_priority = min_age_group_index_priority,
                     time_to_coverage_bpsv = time_to_coverage_bpsv,
                     time_to_coverage_spec = time_to_coverage_spec)
  return(temp)
}

## Create primary vaccination series
# country <- "Argentina"
# population_size <- 10^6
# detection_time <- 14
# vaccine_scenario <- "both_vaccines"
# bpsv_start <- 14
# bpsv_protection_delay <- 7
# specific_vaccine_start <- 200
# specific_protection_delay <- 7
# vaccination_rate <- 0.035
# coverage <- 0.8
# min_age_group_index_priority <- 13
# runtime <- 500
# x <- create_vaccination_dose_series("Argentina", 10^6, 14, "both_vaccines", 14, 7, 200, 7, 0.035, 0.8, 13, 500)
# plot(x$primary_doses, pch = 20)# , xlim = c(220, 260))
# lines(x$second_doses, type = "l")
# lines(x$booster_doses, type = "l", col = "red")
# x <- create_vaccination_dose_series("Argentina", 10^6, 14, "specific_only", 14, 7, 200, 7, 0.035, 0.8, 13, 500)
# plot(x$primary_doses, pch = 20)# , xlim = c(220, 260))
# lines(x$second_doses, type = "l")
# lines(x$booster_doses, type = "l", col = "red")

create_vaccination_dose_series <- function(country, 
                                           population_size, 
                                           detection_time, 
                                           vaccine_scenario,
                                           bpsv_start, 
                                           bpsv_protection_delay,
                                           specific_vaccine_start,
                                           specific_protection_delay,
                                           vaccination_rate,
                                           coverage,
                                           min_age_group_index_priority,
                                           runtime) {
  
  ## Setting Up Vaccination Stuff
  if (!vaccine_scenario %in% c("specific_only", "both_vaccines")) {
    stop("parameter vaccine_scenario must be either 'specific_only' or 'both_vaccines'")
  }
  
  if (vaccine_scenario == "specific_only") {
    
    standard_pop <- generate_standard_pop(country = country, population_size = population_size)
    daily_doses <- vaccination_rate * population_size / 7    # rate of vaccination with primary series
    primary_doses <- 
      c(rep(0, detection_time),
        rep(0, specific_vaccine_start),
        rep(daily_doses, runtime - detection_time - specific_vaccine_start))

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
    daily_doses <- vaccination_rate * population_size / 7    # rate of vaccination with primary series
    priority_age_groups <- min_age_group_index_priority:17            
    elderly_pop_to_vaccinate <- sum(standard_pop[priority_age_groups]) * coverage # 60+s receive primary (BNPCV) and booster (diseaseX-specific); under 60s receive just primary (diseaseX-specific)
    time_to_coverage_bpsv <- ceiling(elderly_pop_to_vaccinate/daily_doses)
    time_to_coverage_spec <- time_to_coverage_bpsv + 1
    
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
        rep(daily_doses, time_to_coverage_bpsv),    ## protection (if any) emerges in BPSV-vaccinated primary vaccinated folks
        rep(0, specific_vaccine_start - time_to_coverage_bpsv - bpsv_protection_delay - bpsv_start), # specific vaccine campaign starts specific_vaccine_start days after detection
        rep(0, time_to_coverage_spec),              ## no specific vaccine for non-elderly whilst that vaccination campaign is ongoing
        rep(0, specific_protection_delay),          ## time between initiation of specific vaccine campaign and people first being protected by that first dose
        rep(daily_doses, runtime - specific_protection_delay - time_to_coverage_spec - specific_vaccine_start - detection_time)) # specific vaccination of all other ages until end of runtime
    
    ## Second Doses
    ## We model full protection as arising the moment you've had the first dose, so second dose is 
    ## redundant given the vaccine-specific protection_delays. We arbitrarily set it to 1 here and therefore stagger
    ## second doses to be a day after primary doses.
    second_doses <- c(0, primary_doses[1:(length(primary_doses) - 1)]) # second dose 1 day after first
    
    # primary_doses[(detection_time+bpsv_start+bpsv_protection_delay+1):
    #                 (detection_time+bpsv_start+bpsv_protection_delay+50)]
    
    # second_doses[(detection_time+bpsv_start+bpsv_protection_delay+1):
    #                (detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+1)] <- daily_doses
    
    ## extra second doses to get up to coverage - why are we not getting there?
    # second_doses[(detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+1)] <- daily_doses
    # second_doses[(detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+2)] <- daily_doses
    # second_doses[(detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+3)] <- daily_doses
    # second_doses[(detection_time+bpsv_start+bpsv_protection_delay+time_to_coverage_bpsv+1+4)] <- daily_doses
    # 
    
    ## Booster Doses (only for elderly, this is the disease specific vaccine)
    ## because we now have booster coverage eligibility indicator in, we don't have to worry
    ## about stopping the booster doses - only elderly will ever get it
    booster_doses <- 
      c(rep(0, detection_time),
        rep(0, specific_vaccine_start),
        rep(0, specific_protection_delay),
        rep(daily_doses, runtime - specific_protection_delay - specific_vaccine_start - detection_time))

  }
  
  return(list(primary_doses = primary_doses,
              second_doses = second_doses,
              booster_doses = booster_doses))

}



## Generates Rt and tt_Rt vectors for - note, this could be replaced with just making a full vector where every timepont has
##                          an associated Rt (instead of just delineating periods as we do currently when dealing with instantaneous changes)
generate_Rt_scenario <- function(Rt, tt_Rt, change_type) {
  
  ## Checking for appropriate Rt change types
  if (any(!(change_type %in% c("gradual", "instant")))) {
    stop("runtime variable can only take the values 'gradual' or 'instant'")
  }
  
  ## Checking whether Rt actually changes or not 
  if (length(Rt) == 1) {
    Rt <- Rt
    tt_Rt <- 0
    return(list(Rt = Rt, tt_Rt = 0))
  }
  
  ## Checking lengths of different model inputs 
  if (length(change_type) != (length(Rt) - 1)) {
    stop("length of change_type must be 1 less than length of Rt vector")
  }
  if (length(change_type) != (length(tt_Rt) - 1)) {
    stop("length of change_type must be 1 less than length of tt_Rt vector")
  }
  if (length(Rt) != length(tt_Rt)) {
    stop("length of Rt vector must be same length as tt_Rt vector")
  }
  
  ## Generating Rt and tt_Rt vectors
  temp_Rt <- Rt[1]
  temp_tt_Rt <- tt_Rt[1]
  Rt_counter <- 1
  tt_Rt_counter <- 1
  for (i in 1:length(change_type)) {
    if (change_type[i] == "instant") {
      Rt_counter <- Rt_counter + 1
      tt_Rt_counter <- tt_Rt_counter + 1
      temp_Rt <- c(temp_Rt, Rt[Rt_counter])
      temp_tt_Rt <- c(temp_tt_Rt, tt_Rt[tt_Rt_counter])
    } else {
      length_output <- length(tt_Rt[tt_Rt_counter]:tt_Rt[tt_Rt_counter + 1])
      temp_Rt <- c(temp_Rt, seq(from = Rt[Rt_counter], to = Rt[Rt_counter + 1], length.out = length_output)[-1])
      temp_tt_Rt <- c(temp_tt_Rt, (tt_Rt[tt_Rt_counter] + 1):(tt_Rt[tt_Rt_counter + 1]))
      Rt_counter <- Rt_counter + 1
      tt_Rt_counter <- tt_Rt_counter + 1  
    }
  }
  return(list(Rt = temp_Rt, tt_Rt = temp_tt_Rt))
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
  timings <- expand_grid(country =  unique(scenarios$country),
                         population_size = unique(scenarios$population_size), 
                         detection_time = unique(scenarios$detection_time),
                         bpsv_start = unique(scenarios$bpsv_start), 
                         specific_vaccine_start = unique(scenarios$specific_vaccine_start),
                         vaccination_rate = unique(scenarios$vaccination_rate),
                         coverage = unique(scenarios$coverage),
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
                             NPI_int == 6 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv, seq(from = detection_time + bpsv_start + time_to_coverage_bpsv, to = detection_time + specific_vaccine_start + time_to_coverage_spec - 1))),
                             NPI_int == 7 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv, seq(from = detection_time + bpsv_start + time_to_coverage_bpsv, to = detection_time + specific_vaccine_start + time_to_coverage_spec - 1))),  
                             NPI_int == 8 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv, seq(from = detection_time + bpsv_start + time_to_coverage_bpsv, to = detection_time + specific_vaccine_start + time_to_coverage_spec - 1))),
                             NPI_int == 9 ~ list(0))) %>%
    ### check with the lines above whether it's appropriate to have the "-1" at the end for NPIs 6, 7, 8
    select(-temp)

  return(NPIs)
}


# Generate combinations of scenarios
create_scenarios <- function(population_size = 1e6,
                             country = "Argentina",
                             hosp_bed_capacity = 1e6,                                         
                             ICU_bed_capacity = 1e6,
                             R0 = c(1.5, 2, 3), 
                             Tg = 7,
                             IFR = c(0.5, 1.5),
                             vaccine_scenario = c("specific_only", "both_vaccines"),
                             detection_time = 14, 
                             bpsv_start = 14,
                             bpsv_protection_delay = 7,
                             specific_vaccine_start = 100,
                             specific_protection_delay = 7,
                             efficacy_infection_bpsv = 0.35,
                             efficacy_disease_bpsv = 0.8, 
                             efficacy_infection_spec = 0.55, 
                             efficacy_disease_spec = 0.9,
                             dur_R = 365000, 
                             dur_V = 365000, 
                             coverage = 0.8,
                             vaccination_rate = 0.035, 
                             min_age_group_index_priority = 13,
                             min_age_group_index_non_priority = 4,
                             runtime = 730,
                             seeding_cases = 2) {
  
  dur_V_temp <- dur_V
  baseline_scenarios <- expand_grid(population_size = population_size,
                                    country = country,
                                    hosp_bed_capacity = hosp_bed_capacity,                                         
                                    ICU_bed_capacity = ICU_bed_capacity,
                                    R0 = R0, 
                                    Tg = Tg,
                                    IFR = IFR,
                                    vaccine_scenario = vaccine_scenario,
                                    detection_time = detection_time, 
                                    bpsv_start = bpsv_start,
                                    bpsv_protection_delay = bpsv_protection_delay, 
                                    specific_vaccine_start = specific_vaccine_start,
                                    specific_protection_delay = specific_protection_delay,
                                    efficacy_infection_bpsv = efficacy_infection_bpsv,
                                    efficacy_disease_bpsv = efficacy_disease_bpsv, 
                                    efficacy_infection_spec = efficacy_infection_spec, 
                                    efficacy_disease_spec = efficacy_disease_spec,
                                    dur_R = dur_R, 
                                    dur_V = 1, 
                                    coverage = coverage,
                                    vaccination_rate = vaccination_rate, 
                                    min_age_group_index_priority = min_age_group_index_priority,
                                    min_age_group_index_non_priority = min_age_group_index_non_priority,
                                    runtime = runtime,
                                    seeding_cases = seeding_cases) 
  baseline_scenarios <- baseline_scenarios %>%
    mutate(dur_V = list(rep(dur_V_temp, 4)))
  
  varying <- variable_columns(baseline_scenarios) 
  
  baseline_scenarios <- baseline_scenarios %>%
    mutate(varied = list(varying))
  
  return(baseline_scenarios)
  
}

## Running the model and summarising outputs
### CHECK - NEED TO ADD IN rel_infectiousness_vaccinated - NOT IN CURRENTLY AND THE DEFAULT IS VERY HIGH (50% REDUCTION)
# population_size = test_both$population_size
# country = test_both$country
# hosp_bed_capacity = test_both$hosp_bed_capacity
# ICU_bed_capacity = test_both$ICU_bed_capacity
# Rt = test_both$Rt
# tt_Rt = test_both$tt_Rt
# Tg = test_both$Tg
# IFR = test_both$IFR
# vaccine_scenario = test_both$vaccine_scenario
# detection_time = test_both$detection_time
# bpsv_start = test_both$bpsv_start
# bpsv_protection_delay = test_both$bpsv_protection_delay
# specific_vaccine_start = test_both$specific_vaccine_start
# specific_protection_delay = test_both$specific_protection_delay
# efficacy_infection_bpsv = 0.1 #test_both$efficacy_infection_bpsv,
# efficacy_disease_bpsv = 0.1
# efficacy_infection_spec = test_both$efficacy_infection_spec
# efficacy_disease_spec = test_both$efficacy_disease_spec
# dur_R = test_both$dur_R
# dur_V = test_both$dur_V
# coverage = test_both$coverage
# vaccination_rate = test_both$vaccination_rate
# min_age_group_index_priority = test_both$min_age_group_index_priority
# min_age_group_index_non_priority = test_both$min_age_group_index_non_priority
# runtime = test_both$runtime
# seeding_cases = test_both$seeding_cases
# output = "full"
# NPI_int = 0
# scenario_index = 0
# varied = ""

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
                       dur_V = rep(1000 * 365, 4),              # duration of vaccine-induced immunity for both vaccines
                       coverage = 0.75,                           # proportion of the population vaccinated
                       vaccination_rate = 0.01,                  # vaccination rate per week as percentage of population
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
  
  set.seed(123)
  
  ## Input Checking
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
    vaccine_doses <- create_vaccination_dose_series(country = country, 
                                                    population_size = population_size, 
                                                    detection_time = detection_time, 
                                                    vaccine_scenario = vaccine_scenario,
                                                    bpsv_start = bpsv_start, 
                                                    bpsv_protection_delay = bpsv_protection_delay,
                                                    specific_vaccine_start = specific_vaccine_start,
                                                    specific_protection_delay = specific_protection_delay,
                                                    vaccination_rate = vaccination_rate,
                                                    coverage = coverage,
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
    vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(coverage, length(priority_age_groups)), 
                                     rep(0, 17 - length(vaccination_age_groups)),  rep(coverage, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)
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
                                                    vaccination_rate = vaccination_rate,
                                                    coverage = coverage,
                                                    min_age_group_index_priority = min_age_group_index_priority,
                                                    runtime = runtime)
    primary_doses <- vaccine_doses$primary_doses
    second_doses <- vaccine_doses$second_doses
    booster_doses <- vaccine_doses$booster_doses
    # note that we have the "+1"s in the underlying function because otherwise we don't *quite* get to vaccinating all of the elderly in
    # the alloted time that Azra was previously calculating. I thought doing ceiling() would be enough and in most cases it is but not all.
    # This way (i.e. by adding +1 to the times) we do! 
    
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
  
  ## Ensuring things that could be passed in as lists (when running in parallel) are coerced to right format
  if (is.list(dur_V)) {
    dur_V <- unlist(dur_V)
    if(length(dur_V) != 4) {
      stop("dur_V must be length 4")
    }
  }
  if (is.list(tt_Rt)) {
    tt_Rt <- unlist(tt_Rt)
  }
  if (is.list(Rt)) {
    Rt <- unlist(Rt)
  }
  if(Rt[1] != Rt[length(Rt)]) { 
    stop("Scenario must start and finish with same R (to reflect reopening)")
  }
  
  ## Running the Model 
  mm <- squire::get_mixing_matrix(country = country)    
  rel_infectiousness_vaccinated <- squire.page:::probs_booster$rel_infectiousness_vaccinated
  rel_infectiousness_vaccinated[rel_infectiousness_vaccinated < 1] <- 1
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
  
  ## Generating summary metrics for that model run
  summary_metrics <- calc_summary_metrics(mod_run)
  
  if (output == "summary") {
    return(list(
      model_arguments = list(scenario_index = scenario_index, population_size = population_size, standard_pop = standard_pop, country = country, hosp_bed_capacity = hosp_bed_capacity, ICU_bed_capacity = ICU_bed_capacity,
                             Rt = Rt, tt_Rt = tt_Rt, Tg = Tg, IFR = IFR, vaccine_scenario = vaccine_scenario, bpsv_start = bpsv_start, specific_vaccine_start = specific_vaccine_start,             
                             efficacy_infection_bpsv = efficacy_infection_bpsv, efficacy_disease_bpsv = efficacy_disease_bpsv, efficacy_infection_spec = efficacy_infection_spec,
                             efficacy_disease_spec = efficacy_disease_spec, dur_R = dur_R, dur_V = dur_V, bpsv_protection_delay = bpsv_protection_delay, specific_protection_delay = specific_protection_delay,                     
                             coverage = coverage, vaccination_rate = vaccination_rate, min_age_group_index_priority = min_age_group_index_priority, min_age_group_index_non_priority = min_age_group_index_non_priority,     
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
                             efficacy_disease_spec = efficacy_disease_spec, dur_R = dur_R, dur_V = dur_V, bpsv_protection_delay = bpsv_protection_delay, specific_protection_delay = specific_protection_delay,               
                             coverage = coverage, vaccination_rate = vaccination_rate, min_age_group_index_priority = min_age_group_index_priority, min_age_group_index_non_priority = min_age_group_index_non_priority,     
                             runtime = runtime, seeding_cases = seeding_cases, detection_time = detection_time,
                             primary_doses = primary_doses, second_doses = second_doses, booster_doses = booster_doses,
                             vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage, 
                             vaccine_booster_initial_coverage = vaccine_booster_initial_coverage,
                             NPI_int = NPI_int, varied = varied),
      summary_metrics = summary_metrics)) 
  }
}

# Summarise and format multiple model simulations with different parameter combinations
format_multirun_output <- function(output_list, parallel = FALSE, cores = NA) {
  
  ## Creating overall dataframe of model outputs
  if (parallel == FALSE) {
    data <- lapply(output_list, function(x) {
      y <- tibble(scenario_index = x$model_arguments$scenario_index, 
                  deaths = x$summary_metrics$deaths, 
                  time_under_NPIs = x$summary_metrics$time_under_NPIs, 
                  composite_NPI = x$summary_metrics$composite_NPI, 
                  country = x$model_arguments$country,
                  population_size = x$model_arguments$population_size,
                  hosp_bed_capacity = x$model_arguments$hosp_bed_capacity,
                  ICU_bed_capacity = x$model_arguments$ICU_bed_capacity,
                  R0 = x$model_arguments$Rt[1],
                  Tg = x$model_arguments$Tg,
                  IFR = x$model_arguments$IFR,
                  vaccine_scenario = x$model_arguments$vaccine_scenario,
                  detection_time = x$model_arguments$detection_time,
                  bpsv_start = ifelse(x$model_arguments$vaccine_scenario == "specific_only", NA, x$model_arguments$bpsv_start),
                  bpsv_protection_delay = ifelse(x$model_arguments$vaccine_scenario == "specific_only", NA, x$model_arguments$bpsv_protection_delay), 
                  specific_vaccine_start = x$model_arguments$specific_vaccine_start,
                  specific_protection_delay = x$model_arguments$specific_protection_delay,
                  efficacy_infection_bpsv = x$model_arguments$efficacy_infection_bpsv,
                  efficacy_disease_bpsv = x$model_arguments$efficacy_disease_bpsv,
                  efficacy_infection_spec = x$model_arguments$efficacy_infection_spec,
                  efficacy_disease_spec = x$model_arguments$efficacy_disease_spec,
                  dur_R = x$model_arguments$dur_R,
                  dur_V = x$model_arguments$dur_V[1],
                  coverage = x$model_arguments$coverage,
                  vaccination_rate = x$model_arguments$vaccination_rate,
                  min_age_group_index_priority = x$model_arguments$min_age_group_index_priority,
                  min_age_group_index_non_priority = x$model_arguments$min_age_group_index_non_priority,
                  runtime = x$model_arguments$runtime,
                  seeding_cases = x$model_arguments$seeding_cases,
                  NPI_int = x$model_arguments$NPI_int,
                  varied = list(x$model_arguments$varied))})
    combined_data <- rbindlist(data)
    
  } else {
    cl <- makeCluster(cores)
    clusterEvalQ(cl, {
      library(data.table)
      library(tibble)
    })
    data <- parLapply(cl, output_list, function(x) {
      y <- tibble(scenario_index = x$model_arguments$scenario_index, 
                  deaths = x$summary_metrics$deaths, 
                  time_under_NPIs = x$summary_metrics$time_under_NPIs, 
                  composite_NPI = x$summary_metrics$composite_NPI, 
                  country = x$model_arguments$country,
                  population_size = x$model_arguments$population_size,
                  hosp_bed_capacity = x$model_arguments$hosp_bed_capacity,
                  ICU_bed_capacity = x$model_arguments$ICU_bed_capacity,
                  R0 = x$model_arguments$Rt[1],
                  Tg = x$model_arguments$Tg,
                  IFR = x$model_arguments$IFR,
                  vaccine_scenario = x$model_arguments$vaccine_scenario,
                  detection_time = x$model_arguments$detection_time,
                  bpsv_start = ifelse(x$model_arguments$vaccine_scenario == "specific_only", NA, x$model_arguments$bpsv_start),
                  bpsv_protection_delay = ifelse(x$model_arguments$vaccine_scenario == "specific_only", NA, x$model_arguments$bpsv_protection_delay), 
                  specific_vaccine_start = x$model_arguments$specific_vaccine_start,
                  specific_protection_delay = x$model_arguments$specific_protection_delay,
                  efficacy_infection_bpsv = x$model_arguments$efficacy_infection_bpsv,
                  efficacy_disease_bpsv = x$model_arguments$efficacy_disease_bpsv,
                  efficacy_infection_spec = x$model_arguments$efficacy_infection_spec,
                  efficacy_disease_spec = x$model_arguments$efficacy_disease_spec,
                  dur_R = x$model_arguments$dur_R,
                  dur_V = x$model_arguments$dur_V[1],
                  coverage = x$model_arguments$coverage,
                  vaccination_rate = x$model_arguments$vaccination_rate,
                  min_age_group_index_priority = x$model_arguments$min_age_group_index_priority,
                  min_age_group_index_non_priority = x$model_arguments$min_age_group_index_non_priority,
                  runtime = x$model_arguments$runtime,
                  seeding_cases = x$model_arguments$seeding_cases,
                  NPI_int = x$model_arguments$NPI_int,
                  varied = list(x$model_arguments$varied))})
    stopCluster(cl) 
    combined_data <- rbindlist(data)
  }
  
  ## Separating out specific only and BPSV/specific scenarios and then left-joining
  both_vax <- combined_data %>%
    filter(vaccine_scenario == "both_vaccines") %>%
    rename(deaths_bpsv = deaths,
           time_under_NPIs_bpsv = time_under_NPIs, ## think we can get rid of this as this isn't specific to both vs vaccine specific scenario
           composite_NPI_bpsv = composite_NPI)     ## think we can get rid of this as this isn't specific to both vs vaccine specific scenario
  
  spec_vax <- combined_data %>%
    filter(vaccine_scenario == "specific_only") %>%
    select(scenario_index, deaths, time_under_NPIs, composite_NPI) %>% 
    rename(deaths_spec = deaths,
           time_under_NPIs_spec = time_under_NPIs,  ## think we can get rid of this as this isn't specific to both vs vaccine specific scenario
           composite_NPI_spec = composite_NPI)      ## think we can get rid of this as this isn't specific to both vs vaccine specific scenario
  
  joined <- both_vax %>%
    left_join(spec_vax, by = "scenario_index") %>% 
    mutate(bpsv_deaths_averted = deaths_spec - deaths_bpsv)
  
  return(joined)
}

### need to come back to this and figure out how this should work when:
#### 1) runtime is < the final tt_Rt value
#### 2) you don't finish up with R0 at the end
## Consider changing summary metrics to accommodate not finishing with same R0
## Be wary of when runtime is shorter than some of the numbers automatically spat out when creating tt_Rt 
calc_summary_metrics <- function(model_output) { # summary metrics = total deaths, time under any NPIs, composite NPI function, 
  
  ## Calculating Deaths
  check <- nimue::format(model_output, compartments = "D", summaries = "deaths") %>%
    filter(t > 1, compartment == "deaths")
  deaths <- sum(check$value)
  
  ## Time Under NPIs
  if (length(model_output$parameters$R0) == 1) {
    time_under_NPIs <- 0
  } else {
    R0 <- model_output$parameters$R0[1]
    which_NPI <-  which(model_output$parameters$R0 < R0)
    NPI_times <- model_output$parameters$tt_R0[c(which_NPI, max(which_NPI) + 1)]
    time_under_NPIs <- max(NPI_times) - min(NPI_times)
  }

  ## Composite of Stringency and Time
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


# generate_Rt_scenario(Rt = c(2), tt_Rt = c(0), change_type = c("instant"))
# generate_Rt_scenario(Rt = c(2, 1), tt_Rt = c(0, 10), change_type = c("gradual"))
# generate_Rt_scenario(Rt = c(2, 1), tt_Rt = c(0, 10), change_type = c("instant"))
# generate_Rt_scenario(Rt = c(2, 1, 1.5), tt_Rt = c(0, 10, 20), change_type = c("gradual", "instant"))
# generate_Rt_scenario(Rt = c(2, 1, 1.5), tt_Rt = c(0, 10, 20), change_type = c("instant", "gradual"))
# generate_Rt_scenario(Rt = c(2, 1, 1.5, 1.2), tt_Rt = c(0, 10, 20, 25), change_type = c("instant", "gradual", "instant"))
# generate_Rt_scenario(Rt = c(2, 1, 1.5, 1.2), tt_Rt = c(0, 10, 20, 25), change_type = c("gradual", "instant", "gradual"))
# generate_Rt_scenario(Rt = c(2, 1, 1.5, 1.2), tt_Rt = c(0, 10, 20, 25), change_type = c("gradual", "gradual", "instant"))
# generate_Rt_scenario(Rt = c(2, 1, 1.5, 1.2), tt_Rt = c(0, 10, 20, 25), change_type = c("instant", "instant", "gradual"))
# default_NPI_scenarios <- function(lockdown_Rt = 0.9,
#                                   minimal_mandate_reduction = 0.25,
#                                   scenarios) {
#   
#   # Unpacking Model Arguments
#   country <- unique(scenarios$country)
#   population_size <- unique(scenarios$population_size)
#   detection_time <- unique(scenarios$detection_time)
#   bpsv_start <- unique(scenarios$bpsv_start)
#   specific_vaccine_start <- unique(scenarios$specific_vaccine_start)
#   vaccination_rate <- unique(scenarios$vaccination_rate)
#   coverage <- unique(scenarios$coverage)
#   R0 <- unique(scenarios$R0)
#   min_age_group_index_priority <- unique(scenarios$min_age_group_index_priority)
#   
#   # Creating Inputs from Model Arguments
#   ## change $ refs in variables that have been unpacked
#   ## Note that this'll probably break if min_age_group_index_priority is varied - come back to change this
#   standard_pop <- generate_standard_pop(country = country, population_size = unique(population_size))
#   daily_doses <- unique(scenarios$vaccination_rate) * unique(scenarios$population_size) / 7    # rate of vaccination with primary series
#   priority_age_groups <- unique(scenarios$min_age_group_index_priority):17            
#   elderly_pop_to_vaccinate <- sum(standard_pop[priority_age_groups]) * coverage # 60+s receive primary (BNPCV) and booster (diseaseX-specific); under 60s receive just primary (diseaseX-specific)
#   time_to_coverage_bpsv <- ceiling(elderly_pop_to_vaccinate/daily_doses) + 1
#   time_to_coverage_spec <- time_to_coverage_bpsv
#   
#   # Checking specific vaccinate starts after bpsv vaccination is completed
#   time_diff <- (bpsv_start + time_to_coverage_bpsv) - specific_vaccine_start
#   if (sum(time_diff >= 0) > 0) {
#     stop("Specific Vaccine start time happens before coverage with BPSV vaccine is finished. Change this.")
#   }
#   
#   # Generating NPI Scenarios
#   NPIs <- expand_grid(country = country,
#                       population_size = population_size, 
#                       coverage = coverage, 
#                       vaccination_rate = vaccination_rate,
#                       min_age_group_index_priority = min_age_group_index_priority,
#                       detection_time = detection_time,
#                       bpsv_start = bpsv_start, 
#                       specific_vaccine_start = specific_vaccine_start,
#                       time_to_coverage_bpsv = time_to_coverage_bpsv,
#                       time_to_coverage_spec = time_to_coverage_spec,
#                       R0 = R0,
#                       NPI_int = 1:9) %>%
#     rowwise() %>%
#     mutate(temp = ((detection_time + specific_vaccine_start + time_to_coverage_spec) - (detection_time + bpsv_start + time_to_coverage_bpsv))) %>%
#     mutate(NPI_scenario = case_when(NPI_int == 1 ~ "LockdownBPSVFinish_minMandateSpecFinish_fullRelease",
#                                     NPI_int == 2 ~ "LockdownBPSVFinish_fullRelease",
#                                     NPI_int == 3 ~ "LockdownSpecFinish_fullRelease",
#                                     NPI_int == 4 ~ "minMandateBPSVFinish_fullRelease",
#                                     NPI_int == 5 ~ "minMandateSpecFinish_fullRelease",
#                                     NPI_int == 6 ~ "LockdownBPSVFinish_gradualReleaseSpecFinish",  
#                                     NPI_int == 7 ~ "minMandateBPSVFinish_gradualReleaseSpecFinish", 
#                                     NPI_int == 8 ~ "LockdownBPSVFinish_minMandategradualReleaseSpecFinish",  
#                                     NPI_int == 9 ~ "Nothing")) %>%
#     mutate(Rt = case_when(NPI_int == 1 ~ list(c(R0, lockdown_Rt, R0 * (1 - minimal_mandate_reduction), R0)),
#                           NPI_int == 2 ~ list(c(R0, lockdown_Rt, R0)),
#                           NPI_int == 3 ~ list(c(R0, lockdown_Rt, R0)),
#                           NPI_int == 4 ~ list(c(R0, R0 * (1 - minimal_mandate_reduction), R0)),
#                           NPI_int == 5 ~ list(c(R0, R0 * (1 - minimal_mandate_reduction), R0)),
#                           NPI_int == 6 ~ list(c(R0, lockdown_Rt, seq(from = lockdown_Rt, to = R0, length.out = temp), R0)), 
#                           NPI_int == 7 ~ list(c(R0, R0 * (1 - minimal_mandate_reduction), seq(from = R0 * (1 - minimal_mandate_reduction), to = R0, length.out = temp), R0)),  
#                           NPI_int == 8 ~ list(c(R0, lockdown_Rt, seq(from = R0 * (1 - minimal_mandate_reduction), to = R0, length.out = temp), R0)), 
#                           NPI_int == 9 ~ list(R0))) %>%
#     mutate(tt_Rt = case_when(NPI_int == 1 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv, detection_time + specific_vaccine_start + time_to_coverage_spec)),
#                              NPI_int == 2 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv)),
#                              NPI_int == 3 ~ list(c(0, detection_time, detection_time + specific_vaccine_start + time_to_coverage_spec)),
#                              NPI_int == 4 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv)),
#                              NPI_int == 5 ~ list(c(0, detection_time, detection_time + specific_vaccine_start + time_to_coverage_spec)),
#                              NPI_int == 6 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv, seq(from = detection_time + bpsv_start + time_to_coverage_bpsv, to = detection_time + specific_vaccine_start + time_to_coverage_spec - 1))),
#                              NPI_int == 7 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv, seq(from = detection_time + bpsv_start + time_to_coverage_bpsv, to = detection_time + specific_vaccine_start + time_to_coverage_spec - 1))),  
#                              NPI_int == 8 ~ list(c(0, detection_time, detection_time + bpsv_start + time_to_coverage_bpsv, seq(from = detection_time + bpsv_start + time_to_coverage_bpsv, to = detection_time + specific_vaccine_start + time_to_coverage_spec - 1))),
#                              NPI_int == 9 ~ list(0))) %>%
#     select(-temp)
#   
#   return(NPIs)
#   
# }

