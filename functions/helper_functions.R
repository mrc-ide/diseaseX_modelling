## Define default params
define_default_params <- function() {
  default_params <- list(R0 = 2.5,
                         IFR = 1,
                         population_size = 10^10,
                         hosp_bed_capacity = 10^10,
                         ICU_bed_capacity = 10^10, 
                         Tg = 6.7,                                     
                         detection_time = 1,                           
                         bpsv_start = 7,                               
                         bpsv_protection_delay = 7,                    
                         specific_vaccine_start = 250,    
                         specific_protection_delay = 7,                
                         efficacy_infection_bpsv = 0.35,               
                         efficacy_disease_bpsv = 0.75,                 
                         efficacy_infection_spec = 0.55,               
                         efficacy_disease_spec = 0.9,                  
                         dur_R = 365000000,                            
                         dur_bpsv = 365000000,                         
                         dur_spec = 365000000,                         
                         coverage_bpsv = 0.8,                          
                         coverage_spec = 0.8,                          
                         vaccination_rate_bpsv = 0.035,                
                         vaccination_rate_spec = 0.035,                
                         min_age_group_index_priority = 13,
                         min_age_group_index_non_priority = 4,         
                         seeding_cases = 1,
                         lockdown_Rt = 0.9,
                         minimal_mandate_reduction = 0.25,
                         runtime = 730)
  return(default_params)
}

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
generate_standard_pop <- function(country = "Argentina", population_size = 1e9) {
  raw_pop <- squire::get_population(country = country)$n                   
  standard_pop <- round((population_size / sum(raw_pop)) * raw_pop) 
  return(standard_pop)
}

## Generate generation time
scale_generation_time <- function(target_Tg) {
  current_Tg_mild <- squire.page.sarsX:::durs_booster$dur_E + squire.page.sarsX:::durs_booster$dur_IMild
  current_Tg_case <- squire.page.sarsX:::durs_booster$dur_E + squire.page.sarsX:::durs_booster$dur_ICase
  overall_Tg <- 0.98 * current_Tg_mild + 0.02 * current_Tg_case # (assumed approx 2% IHR - might need changing, though won't make much diff in practice as IMild dominates)
  Tg_ratio <- target_Tg / overall_Tg
  TgVary_dur_IMild <- Tg_ratio * squire.page.sarsX:::durs_booster$dur_IMild
  TgVary_dur_ICase <- Tg_ratio * squire.page.sarsX:::durs_booster$dur_ICase
  return(list(dur_IMild = TgVary_dur_IMild, dur_ICase = TgVary_dur_ICase))
}

## Scale IFR (needs checking/changing - not sure this is quite right)
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

# Generate scenarios with different parameter combinations to run many simulations at once
create_scenarios <- function(country,
                             R0, 
                             IFR,
                             Tg,
                             population_size,
                             hosp_bed_capacity,                                         
                             ICU_bed_capacity,
                             vaccine_scenario = c("specific_only", "both_vaccines"),
                             detection_time, 
                             bpsv_start,
                             bpsv_protection_delay,
                             specific_vaccine_start,
                             specific_protection_delay,
                             efficacy_infection_bpsv,
                             efficacy_disease_bpsv, 
                             efficacy_infection_spec, 
                             efficacy_disease_spec,
                             dur_R, 
                             dur_bpsv, 
                             dur_spec,
                             coverage_bpsv,
                             coverage_spec,
                             vaccination_rate_bpsv, 
                             vaccination_rate_spec,
                             min_age_group_index_priority,
                             min_age_group_index_non_priority,
                             runtime,
                             seeding_cases) {
  
  
  default <- define_default_params()
  
  if (missing(country)) {
    country <- "Argentina"
  }
  if (missing(R0)) {
    R0 <- default$R0
  } 
  if (missing(IFR)) {
    IFR <- default$IFR
  }
  if (missing(Tg)) {
    Tg <- default$Tg
  }
  if (missing(population_size)) {
    population_size <- default$population_size
  }
  if (missing(hosp_bed_capacity)) {
    hosp_bed_capacity <- default$hosp_bed_capacity
  }
  if (missing(ICU_bed_capacity)) {
    ICU_bed_capacity <- default$ICU_bed_capacity
  }
  if (missing(detection_time)) {
    detection_time <- default$detection_time
  }
  if (missing(bpsv_start)) {
    bpsv_start <- default$bpsv_start
  }
  if (missing(bpsv_protection_delay)) {
    bpsv_protection_delay <- default$bpsv_protection_delay
  }
  if (missing(specific_vaccine_start)) {
    specific_vaccine_start <- default$specific_vaccine_start
  }
  if (missing(specific_protection_delay)) {
    specific_protection_delay <- default$specific_protection_delay
  }
  if (missing(efficacy_infection_bpsv)) {
    efficacy_infection_bpsv <- default$efficacy_infection_bpsv
  }
  if (missing(efficacy_disease_bpsv)) {
    efficacy_disease_bpsv <- default$efficacy_disease_bpsv
  }
  if (missing(efficacy_infection_spec)) {
    efficacy_infection_spec <- default$efficacy_infection_spec
  }
  if (missing(efficacy_disease_spec)) {
    efficacy_disease_spec <- default$efficacy_disease_spec
  }
  if (missing(dur_R)) {
    dur_R <- default$dur_R
  }
  if (missing(dur_bpsv)) {
    dur_bpsv <- default$dur_bpsv
  }
  if (missing(dur_spec)) {
    dur_spec <- default$dur_spec
  }
  if (missing(coverage_bpsv)) {
    coverage_bpsv <- default$coverage_bpsv
  }
  if (missing(coverage_spec)) {
    coverage_spec <- default$coverage_spec
  }
  if (missing(vaccination_rate_bpsv)) {
    vaccination_rate_bpsv <- default$vaccination_rate_bpsv
  }
  if (missing(vaccination_rate_spec)) {
    vaccination_rate_spec <- default$vaccination_rate_spec
  }
  if (missing(min_age_group_index_priority)) {
    min_age_group_index_priority <- default$min_age_group_index_priority
  }
  if (missing(min_age_group_index_non_priority)) {
    min_age_group_index_non_priority <- default$min_age_group_index_non_priority
  }
  if (missing(runtime)) {
    runtime <- default$runtime
  }
  if (missing(seeding_cases)) {
    seeding_cases <- default$seeding_cases
  }

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
                                    dur_bpsv = dur_bpsv / 2,  # two waning compartments, so halve the duration to double the rate and keep the total time spent with immunity = dur_bpsv on average
                                    dur_spec = dur_spec / 2,  # two waning compartments, so halve the duration to double the rate and keep the total time spent with immunity = dur_spec on average
                                    coverage_bpsv = coverage_bpsv,
                                    coverage_spec = coverage_spec,
                                    vaccination_rate_bpsv = vaccination_rate_bpsv, 
                                    vaccination_rate_spec = vaccination_rate_spec,
                                    min_age_group_index_priority = min_age_group_index_priority,
                                    min_age_group_index_non_priority = min_age_group_index_non_priority,
                                    runtime = runtime,
                                    seeding_cases = seeding_cases) 
  
  # Identify which columns vary in baseline scenarios (i.e. which parameters you're varying) - helpful reminder for later when processing the outputs
  varying <- variable_columns(baseline_scenarios) 
  baseline_scenarios <- baseline_scenarios %>%
    mutate(varied = list(varying))
  
  return(baseline_scenarios)
  
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
                  dur_bpsv = x$model_arguments$dur_bpsv,
                  dur_spec = x$model_arguments$dur_spec,
                  coverage_bpsv = x$model_arguments$coverage_bpsv,
                  coverage_spec = x$model_arguments$coverage_spec,
                  vaccination_rate_bpsv = x$model_arguments$vaccination_rate_bpsv,
                  vaccination_rate_spec = x$model_arguments$vaccination_rate_spec,
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
                  R0 = max(x$model_arguments$Rt),
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
                  dur_bpsv = x$model_arguments$dur_bpsv,
                  dur_spec = x$model_arguments$dur_spec,
                  coverage_bpsv = x$model_arguments$coverage_bpsv,
                  coverage_spec = x$model_arguments$coverage_spec,
                  vaccination_rate_bpsv = x$model_arguments$vaccination_rate_bpsv,
                  vaccination_rate_spec = x$model_arguments$vaccination_rate_spec,
                  min_age_group_index_priority = x$model_arguments$min_age_group_index_priority,
                  min_age_group_index_non_priority = x$model_arguments$min_age_group_index_non_priority,
                  runtime = x$model_arguments$runtime,
                  seeding_cases = x$model_arguments$seeding_cases,
                  NPI_int = x$model_arguments$NPI_int,
                  varied = list(x$model_arguments$varied))})
    stopCluster(cl) 
    combined_data <- rbindlist(data)
  }
  
  ## Separating out specific only and BPSV/specific scenarios and then left-joining so they're together on a single row
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

pareto_frontier <- function(df) {
  df_sorted <- df %>% arrange(NPI_days)
  n <- nrow(df_sorted)
  is_pareto <- rep(TRUE, n)
  for (i in 2:(n-1)) {
    if ((df_sorted$deaths[i-1] < df_sorted$deaths[i])) {
      is_pareto[i] <- FALSE
    }
  }
  return(df_sorted[is_pareto, ])
}

linear_interpolate <- function(df) {
  new_df <- data.frame()
  for (i in 1:(nrow(df) - 1)) {
    x1 <- df$new_NPI_days[i]
    y1 <- df$deaths[i]
    x2 <- df$new_NPI_days[i + 1]
    y2 <- df$deaths[i + 1]
    
    for (x in floor(x1):ceiling(x2)) {
      if (x >= x1 && x <= x2) {
        y <- y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        new_row <- data.frame(deaths = y, new_NPI_days = x)
        new_df <- rbind(new_df, new_row)
      }
    }
  }
  return(new_df)
}
