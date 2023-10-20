## Generates NPI scenarios for secondary country (i.e. country where pathogen doesn't emerge and so there is potentially advance warning relative to importation)
### 4 = NPIs minimal mandate until BPSV campaign is done, followed by full release
### 7 = NPIs minimal mandate until BPSV campaign is done, followed by gradual release until spec vacc campaign done
### 8 = NPIs reducing Rt < 1 until BPSV campaign is done, followed by gradual release from minimal mandate until specific campaign is done
### 9 = No NPIs
default_NPI_scenarios_secondary <- function(lockdown_Rt = 0.9,
                                            minimal_mandate_reduction = 0.25,
                                            NPI_scenarios = 1:9,
                                            scenarios) { 
  
  # Looping through scenarios and setting the detection time. This is either:
  #   The day of "Arrival" for instances where arrival in secondary country is after detection in source country.
  #   Day 0 for instances where detection in origin country precedes arrival of pathogen in secondary country (Day 0 in model simulation = day pathogen is detected in source country)
  #   --> The way this works is we pad the days before detection and importation in the origin country with a bunch of days where Rt = 1 (to emulate no infections arriving).
  #   --> This means we can run the model and get people vaccinated before there's been any meaningful transmission (because small number of starting infection + Rt = 1 = minimal impact on dynamics)
  #   --> and in doing so, effectively simulate the "advanced warning" that the secondary country has in this scenario. 
  for (i in 1:nrow(scenarios)) {
    
    ## Where detection in origin country happens AFTER the arrival of the pathogen in the secondary country
    ## BPSV campaign and specific development i.e. following pathogen detection in source country, starts AFTER pathogen arrival in secondary country
    ### Note - Try to avoid scenarios where days_source_detection_is_ahead_arrival_secondary is exactly 0 as it complicates things quite a lot.
    if (scenarios$days_source_detection_is_ahead_arrival_secondary[i] <= 0) {     
      scenarios$detection_time[i] <- -scenarios$days_source_detection_is_ahead_arrival_secondary[i]
      
    ## Where detection in origin country happens BEFORE the arrival of the pathogen in the secondary country
    ## BPSV campaign and specific development i.e. following pathogen detection in source country, starts BEFORE pathogen arrival in secondary country
    } else if (scenarios$days_source_detection_is_ahead_arrival_secondary[i] > 0) { 
      scenarios$detection_time[i] <- 0 
    }
  }
  
  # Unpack the unique parameter combinations in the scenarios dataframe and calculate vaccination milestone timings
  timings <- expand_grid(country = unique(scenarios$country),
                         population_size = unique(scenarios$population_size), 
                         detection_time = unique(scenarios$detection_time),
                         bpsv_start = unique(scenarios$bpsv_start), 
                         specific_vaccine_start = unique(scenarios$specific_vaccine_start),
                         vaccination_rate_bpsv = unique(scenarios$vaccination_rate_bpsv),
                         vaccination_rate_spec = unique(scenarios$vaccination_rate_spec),
                         coverage = unique(scenarios$coverage),
                         min_age_group_index_priority = unique(scenarios$min_age_group_index_priority),
                         days_source_detection_is_ahead_arrival_secondary = unique(scenarios$days_source_detection_is_ahead_arrival_secondary),
                         detection_time_secondary = unique(scenarios$detection_time_secondary)) %>%
    filter((detection_time == -days_source_detection_is_ahead_arrival_secondary) | (days_source_detection_is_ahead_arrival_secondary > 0 & detection_time == 0))
  vacc_timings <- pmap(timings[-c(ncol(timings) - 1, ncol(timings))], vaccine_milestone_timings)
  vacc_timings <- rbindlist(vacc_timings)
  vacc_timings$days_source_detection_is_ahead_arrival_secondary <- timings$days_source_detection_is_ahead_arrival_secondary
  vacc_timings$detection_time_secondary <- timings$detection_time_secondary
  
  # Checking specific vaccinate starts after bpsv vaccination is completed
  time_diff <- (vacc_timings$bpsv_start + vacc_timings$time_to_coverage_bpsv) - vacc_timings$specific_vaccine_start
  if (sum(time_diff >= 0) > 0) {
    stop("Specific Vaccine start time happens before coverage with BPSV vaccine is finished. Change this.")
  }  
  
  ## Timings when pathogen arrives in secondary country BEFORE detection in the source country - secondary country detection as the trigger for NPI imposition
  ## --> detection_time in source country is trigger for BPSV campaign, but detection_time_secondary is trigger for initial NPI imposition
  ##     We make the assumption that source country ALWAYS detects the pathogen first, and that because of enhanced surveillance following detection, the time until
  ##     the secondary country detects the virus is in calendar days (e.g. like 7 days) rather than dependent on number of hospitalisations.
  ## --> Note that detection_time_secondary here is relative to detection_time (which works because we have constraints that detection_time_secondary has to happen after detection_time)
  ## -> TO DO: What remains somewhat unclear to me what would be the trigger in the secondary country for imposition of NPIs - probably detection in secondary country?
  ##           Need to check this with Azra. 
  ## -> TO DO: Because of the way the triggers work (NPIs dependent on detection_time_secondary but BPSV campaign on detection_time), if detection_time_secondary > (bpsv_start + time_to_coverage_bpsv)
  ##           then campaign will be done by the time detection occurs, and hence NPI scenario 4 becomes identical to NPI scenario 9. I think we can get round this by having
  ##           detection_time_secondary to be just a small amount of time after detection_time. I've put a check in here to make sure of that, but probs we can do better.
  NPIs_arrival_before <- expand_grid(vacc_timings[vacc_timings$days_source_detection_is_ahead_arrival_secondary < 0], 
                                    R0 = unique(scenarios$R0),
                                    NPI_int = NPI_scenarios)
  check <- (NPIs_arrival_before$detection_time + NPIs_arrival_before$detection_time_secondary) < (NPIs_arrival_before$detection_time + NPIs_arrival_before$bpsv_start + NPIs_arrival_before$time_to_coverage_bpsv)
  if (sum(check) > 0) {
    stop("For some scenarios where pathogen arrives in secondary country BEFORE detection, detection_time_secondary is occurring after the completion of the BPSV campaign - change this")
  }
  NPIs_arrival_before <- NPIs_arrival_before %>%
    rowwise() %>%
    mutate(temp = ((specific_vaccine_start + time_to_coverage_spec) - (bpsv_start + time_to_coverage_bpsv))) %>%
    mutate(Rt = case_when(NPI_int == 4 ~ list(c(R0, R0 * (1 - minimal_mandate_reduction), R0)),
                          NPI_int == 7 ~ list(c(R0, R0 * (1 - minimal_mandate_reduction), seq(from = R0 * (1 - minimal_mandate_reduction), to = R0, length.out = temp))),  
                          NPI_int == 8 ~ list(c(R0, lockdown_Rt, seq(from = R0 * (1 - minimal_mandate_reduction), to = R0, length.out = temp))), 
                          NPI_int == 9 ~ list(R0))) %>%
    ## Below are fine for now only because detection_time_secondary < (bpsv_start + time_to_coverage_bpsv) - otherwise would end up with non-monotonically tt_Rt vector
    mutate(tt_Rt = case_when(NPI_int == 4 ~ list(c(0, detection_time + detection_time_secondary, detection_time + bpsv_start + time_to_coverage_bpsv)),
                             NPI_int == 7 ~ list(c(0, detection_time + detection_time_secondary, detection_time + bpsv_start + time_to_coverage_bpsv, seq(from = detection_time + bpsv_start + time_to_coverage_bpsv + 1, to = detection_time + specific_vaccine_start + time_to_coverage_spec - 1))),  
                             NPI_int == 8 ~ list(c(0, detection_time + detection_time_secondary, detection_time + bpsv_start + time_to_coverage_bpsv, seq(from = detection_time + bpsv_start + time_to_coverage_bpsv + 1, to = detection_time + specific_vaccine_start + time_to_coverage_spec - 1))),
                             NPI_int == 9 ~ list(0))) %>%
    select(-temp) ### TO DO: Check with the lines above whether it's appropriate to have the "-1" at the end for NPIs 6, 7, 8
    
  
  ## Timings when pathogen arrives in secondary country AFTER detection in the source - secondary country detection as the trigger for NPI imposition
  ### 4 = NPIs minimal mandate until BPSV campaign is done, followed by full release
  ### 7 = NPIs minimal mandate until BPSV campaign is done, followed by gradual release until spec vacc campaign done
  ### 8 = NPIs reducing Rt < 1 until BPSV campaign is done, followed by gradual release from minimal mandate until specific campaign is done
  ### 9 = No NPIs
  ## Note gradient of the Rt increase in gradual lifting will vary sometimes here depending on comparative timing of different events
  NPIs_arrival_after <- expand_grid(vacc_timings[vacc_timings$days_source_detection_is_ahead_arrival_secondary > 0], 
                                     R0 = unique(scenarios$R0),
                                     NPI_int = NPI_scenarios)
  check <- (NPIs_arrival_before$detection_time + NPIs_arrival_before$detection_time_secondary) > (NPIs_arrival_before$detection_time + NPIs_arrival_before$days_source_detection_is_ahead_arrival_secondary)
  if (sum(check) > 0) {
    stop("For some scenarios where pathogen arrives in secondary country BEFORE detection, detection_time_secondary is occurring after the completion of the BPSV campaign - change this")
  }
  NPIs_arrival_after <- NPIs_arrival_after %>%
    rowwise() %>%
    mutate(days_ago_bpsv_finish = detection_time_secondary - (bpsv_start + time_to_coverage_bpsv)) %>%                 ## remaining days of campaign bpsv after detection - if positive, bpsv campaign finishes before detection in secondary; if negative, bpsv campaign finishes after detection in secondary
    mutate(days_ago_spec_finish = detection_time_secondary - (specific_vaccine_start + time_to_coverage_spec)) %>%     ## remaining days of campaign spec after detection - if positive, spec campaign finishes before detection in secondary; if negative, spec campaign finishes after detection in secondary
    mutate(time_between_bpsv_finish_and_spec_finish = (specific_vaccine_start + time_to_coverage_spec) - (bpsv_start + time_to_coverage_bpsv)) %>%
    mutate(time_between_arrival_and_spec_finish = (specific_vaccine_start + time_to_coverage_spec) - (days_source_detection_is_ahead_arrival_secondary)) %>%
    
    mutate(Rt = case_when(NPI_int == 4 ~ ifelse(days_ago_bpsv_finish < 0 & days_ago_spec_finish < 0, list(c(1, R0, R0 * (1 - minimal_mandate_reduction), R0)), ## corresponding to the detection in secondary occurring before bpsv and spec campaigns have finished
                                                ifelse(days_ago_bpsv_finish >= 0 & days_ago_spec_finish < 0, list(c(1, R0)),                                           ## corresponding to detection in secondary occurring after bpsv campaign finished but before spec campaign has finished
                                                       ifelse(days_ago_bpsv_finish > 0 & days_ago_spec_finish > 0, list(c(1, R0)), NA))),                                    ## corresponding to detection in secondary occurring after bpsv campaign and spec campaign have both finished
                          
                          NPI_int == 7 ~ ifelse(days_ago_bpsv_finish < 0 & days_ago_spec_finish < 0, list(c(1, R0, R0 * (1 - minimal_mandate_reduction),  seq(from = R0 * (1 - minimal_mandate_reduction), to = R0, length.out = time_between_bpsv_finish_and_spec_finish))), ## corresponding to the detection in secondary occurring before bpsv and spec campaigns have finished
                                                ifelse(days_ago_bpsv_finish >= 0 & days_ago_spec_finish < 0, list(c(1, R0, seq(from = R0 * (1 - minimal_mandate_reduction), to = R0, length.out = time_between_arrival_and_spec_finish))),                                        ## corresponding to detection in secondary occurring after bpsv campaign finished but before spec campaign has finished
                                                       ifelse(days_ago_bpsv_finish > 0 & days_ago_spec_finish > 0, list(c(1, R0)), NA))),                                                                                                                    ## corresponding to detection in secondary occurring after bpsv campaign and spec campaign have both finished
                          
                          NPI_int == 8 ~ ifelse(days_ago_bpsv_finish < 0 & days_ago_spec_finish < 0, list(c(1, R0, lockdown_Rt,  seq(from = R0 * (1 - minimal_mandate_reduction), to = R0, length.out = time_between_bpsv_finish_and_spec_finish))), ## corresponding to the detection in secondary occurring before bpsv and spec campaigns have finished
                                                ifelse(days_ago_bpsv_finish >= 0 & days_ago_spec_finish < 0, list(c(1, R0, seq(from = R0 * (1 - minimal_mandate_reduction), to = R0, length.out = time_between_arrival_and_spec_finish))),                                        ## corresponding to detection in secondary occurring after bpsv campaign finished but before spec campaign has finished
                                                       ifelse(days_ago_bpsv_finish > 0 & days_ago_spec_finish > 0, list(c(1, R0)), NA))),
                          
                          NPI_int == 9 ~ list(c(1, R0)))) %>%
    
    
    mutate(tt_Rt = case_when(NPI_int == 4 ~ ifelse(days_ago_bpsv_finish < 0 & days_ago_spec_finish < 0, list(c(0, days_source_detection_is_ahead_arrival_secondary, days_source_detection_is_ahead_arrival_secondary + detection_time_secondary, days_source_detection_is_ahead_arrival_secondary + detection_time_secondary - days_ago_bpsv_finish)),          ## corresponding to the detection in secondary occurring before bpsv and spec campaigns have finished
                                                   ifelse(days_ago_bpsv_finish >= 0 & days_ago_spec_finish < 0, list(c(0, days_source_detection_is_ahead_arrival_secondary)),        ## corresponding to detection in secondary occurring after bpsv campaign finished but before spec campaign has finished
                                                          ifelse(days_ago_bpsv_finish > 0 & days_ago_spec_finish > 0, list(c(0, days_source_detection_is_ahead_arrival_secondary)), NA))), ## corresponding to detection in secondary occurring after bpsv campaign and spec campaign have both finished 
                             
                             NPI_int == 7 ~ ifelse(days_ago_bpsv_finish < 0 & days_ago_spec_finish < 0, list(c(0, days_source_detection_is_ahead_arrival_secondary, days_source_detection_is_ahead_arrival_secondary + detection_time_secondary, days_source_detection_is_ahead_arrival_secondary + detection_time_secondary - days_ago_bpsv_finish, seq(from = days_source_detection_is_ahead_arrival_secondary + detection_time_secondary - days_ago_bpsv_finish + 1, to = days_source_detection_is_ahead_arrival_secondary + detection_time_secondary - days_ago_spec_finish - 1))),  ## does the final term need days_source_detection_is_ahead_arrival_secondary????        ## corresponding to the detection in secondary occurring before bpsv and spec campaigns have finished
                                                   ifelse(days_ago_bpsv_finish >= 0 & days_ago_spec_finish < 0, list(c(0, days_source_detection_is_ahead_arrival_secondary, days_source_detection_is_ahead_arrival_secondary + detection_time_secondary, seq(from = days_source_detection_is_ahead_arrival_secondary + detection_time_secondary + 1, to = days_source_detection_is_ahead_arrival_secondary + detection_time_secondary - days_ago_spec_finish - 1))),        ## corresponding to detection in secondary occurring after bpsv campaign finished but before spec campaign has finished
                                                          ifelse(days_ago_bpsv_finish > 0 & days_ago_spec_finish > 0, list(c(0, days_source_detection_is_ahead_arrival_secondary)), NA))),  ## corresponding to detection in secondary occurring after bpsv campaign and spec campaign have both finished 
                             
                             NPI_int == 8 ~ ifelse(days_ago_bpsv_finish < 0 & days_ago_spec_finish < 0, list(c(0, days_source_detection_is_ahead_arrival_secondary, days_source_detection_is_ahead_arrival_secondary + detection_time_secondary, days_source_detection_is_ahead_arrival_secondary + detection_time_secondary - days_ago_bpsv_finish, seq(from = days_source_detection_is_ahead_arrival_secondary + detection_time_secondary - days_ago_bpsv_finish + 1, to = days_source_detection_is_ahead_arrival_secondary + detection_time_secondary - days_ago_spec_finish - 1))),  ## does the final term need days_source_detection_is_ahead_arrival_secondary????        ## corresponding to the detection in secondary occurring before bpsv and spec campaigns have finished
                                                   ifelse(days_ago_bpsv_finish >= 0 & days_ago_spec_finish < 0, list(c(0, days_source_detection_is_ahead_arrival_secondary, days_source_detection_is_ahead_arrival_secondary + detection_time_secondary, seq(from = days_source_detection_is_ahead_arrival_secondary + detection_time_secondary + 1, to = days_source_detection_is_ahead_arrival_secondary + detection_time_secondary - days_ago_spec_finish - 1))),        ## corresponding to detection in secondary occurring after bpsv campaign finished but before spec campaign has finished
                                                          ifelse(days_ago_bpsv_finish > 0 & days_ago_spec_finish > 0, list(c(0, days_source_detection_is_ahead_arrival_secondary)), NA))),  ## corresponding to detection in secondary occurring after bpsv campaign and spec campaign have both finished 
                             
                             NPI_int == 9 ~ list(c(0, days_source_detection_is_ahead_arrival_secondary)))) %>%
    select(-days_ago_bpsv_finish, -days_ago_spec_finish, -time_between_bpsv_finish_and_spec_finish, -time_between_arrival_and_spec_finish)
  
  # sapply(NPIs_arrival_before$Rt, function(x) length(x)) == sapply(NPIs_arrival_before$tt_Rt, function(x) length(x)) # Checking Rt and tt_Rt are the same length
  # sapply(NPIs_arrival_after$Rt, function(x) length(x)) == sapply(NPIs_arrival_after$tt_Rt, function(x) length(x)) # Checking Rt and tt_Rt are the same length
  overall <- rbind(NPIs_arrival_before, NPIs_arrival_after)
  
  return(overall)
}
