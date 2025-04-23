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

convert_array <- function(arr, scenario, pathogen, R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan) {
  
  df <- as.data.frame.table(arr, responseName = "epidemic_size")
  # Rename the Var* columns to something more meaningful
  colnames(df) <- c("iteration_i", "R0_i", "vaccine_efficacy_i", "quarantine_efficacy_i", "epidemic_size")
  
  # Map these integer indices back to their actual values
  df <- df %>%
    mutate(
      iteration           = as.integer(iteration_i),
      scenario            = scenario,
      pathogen            = pathogen,
      R0                  = R0_scan[R0_i],
      vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy_i], 
      vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy_i], 
      quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy_i]
    ) %>%
    # Keep columns in a nice order
    dplyr::select(iteration, scenario, pathogen, R0, vaccine_efficacy_infection, vaccine_efficacy_transmission, quarantine_efficacy, epidemic_size)
  
  
  return(df)
}

convert_array_time_to_n <- function(arr, scenario, pathogen, R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan) {
  
  df <- as.data.frame.table(arr, responseName = "time_to_n")
  # Rename the Var* columns to something more meaningful
  colnames(df) <- c("iteration_i", "R0_i", "vaccine_efficacy_i", "quarantine_efficacy_i", "time_to_n")
  
  # Map these integer indices back to their actual values
  df <- df %>%
    mutate(
      iteration           = as.integer(iteration_i),
      scenario            = scenario,
      pathogen            = pathogen,
      R0                  = R0_scan[R0_i],
      vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy_i], 
      vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy_i], 
      quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy_i]
    ) %>%
    # Keep columns in a nice order
    dplyr::select(iteration, scenario, pathogen, R0, vaccine_efficacy_infection, vaccine_efficacy_transmission, quarantine_efficacy, time_to_n)
  
  return(df)
}

convert_array_Reff <- function(arr, scenario, pathogen, R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan) {
  
  df <- as.data.frame.table(arr, responseName = "time_to_n")
  # Rename the Var* columns to something more meaningful
  colnames(df) <- c("iteration_i", "R0_i", "vaccine_efficacy_i", "quarantine_efficacy_i", "Reff")
  
  # Map these integer indices back to their actual values
  df <- df %>%
    mutate(
      iteration           = as.integer(iteration_i),
      scenario            = scenario,
      pathogen            = pathogen,
      R0                  = R0_scan[R0_i],
      vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy_i], 
      vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy_i], 
      quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy_i]
    ) %>%
    # Keep columns in a nice order
    dplyr::select(iteration, scenario, pathogen, R0, vaccine_efficacy_infection, vaccine_efficacy_transmission, quarantine_efficacy, Reff)
  
  return(df)
}

convert_array_R0 <- function(arr, scenario, pathogen, R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan) {
  
  df <- as.data.frame.table(arr, responseName = "time_to_n")
  # Rename the Var* columns to something more meaningful
  colnames(df) <- c("iteration_i", "R0_i", "vaccine_efficacy_i", "quarantine_efficacy_i", "R0_actual")
  
  # Map these integer indices back to their actual values
  df <- df %>%
    mutate(
      iteration           = as.integer(iteration_i),
      scenario            = scenario,
      pathogen            = pathogen,
      R0                  = R0_scan[R0_i],
      vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy_i], 
      vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy_i], 
      quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy_i]
    ) %>%
    # Keep columns in a nice order
    dplyr::select(iteration, scenario, pathogen, R0, vaccine_efficacy_infection, vaccine_efficacy_transmission, quarantine_efficacy, R0_actual)
  
  return(df)
}

calculate_R0 <- function(tdf) {
  R0 <- mean(tdf$n_offspring, na.rm = TRUE)
  return(R0)
}

calculate_Reff <- function(tdf, type = "not_spatial_vax") {
  if (type == "spatial_vax") {
    Reff <- mean(tdf$n_offspring_new_new, na.rm = TRUE)
  } else {
    Reff <- mean(tdf$n_offspring_post_pruning, na.rm = TRUE)
  }
  return(Reff)
}

make_stacked_plot_timetoN <- function(df, patho, qe,
                                      rel_heights = c(1, 3),
                                      palette      = c("#CA2E6B", "#88C5EE", "#236897", "#13496E", "black")) {
  
  df_sub <- df %>% 
    filter(vaccine_efficacy_infection == 0.35,
           quarantine_efficacy        == qe,
           pathogen                   == patho)
  
  ## bottom panel – containment
  p_cont <- ggplot(df_sub,
                   aes(R0, 100 * proportion_contained, colour = scenario)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_colour_manual(values = palette, 
                        labels = c("No Delay", "2 Days", "1 Week", "2 Weeks", "No Vaccination"),
                        name = "Vaccine\nProtection\nDelay",
                        guide = guide_legend(reverse = TRUE)) +
    labs(x = "R0", y = "% Outbreaks Contained") +
    theme(strip.background = element_rect(fill = "white"),
          plot.margin = margin(t = -10, r = 5, b = 0, l = 5),
          legend.position = "none")
  
  ## top panel – time to N
  p_time <- ggplot(df_sub,
                   aes(R0, time_to_n, colour = scenario)) +
    geom_line() + geom_point() +  theme_bw() +
    scale_colour_manual(values = palette,
                        guide  = "none") +      # legend only once
    labs(x = "", y = "") +
    theme(strip.background = element_rect(fill = "white"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 10, r = 5, b = 0, l = 5),
          legend.position = "none")
  
  plot_grid(p_time, p_cont, nrow = 2, rel_heights = rel_heights)
}

#  
# df <- overall_spatial_df
# qe <- 0
# patho <- "SARS-CoV-1"
# vaccine_efficacy_infection_value <- 0.35
# rel_heights <- c(1, 2)
# palette <- c("#E3AFCB", "#D474A4", "#B52F7B", "#9C105A", "#6B0045", "#474747")

make_stacked_plot_Reff_spatialvax <- function(df, patho, qe, vaccine_efficacy_infection_value, 
                                              rel_heights = c(1, 3),
                                              palette      = c("#E3AFCB", "#D474A4", "#B52F7B", "#9C105A", "#6B0045", "#474747")) {
  
  df_sub <- df %>% 
    filter(vaccine_efficacy_infection == vaccine_efficacy_infection_value,
           quarantine_efficacy        == qe,
           pathogen                   == patho)
  
  # constant to make sure everything's aligned
  x_breaks_raw   <- seq(0.75, 2.50, 0.25)   # tick marks
  x_breaks <- seq(1, 2.5, 0.5)
  bar_width  <- 0.2                    # <- same width you pass to geom_bar
  half_bw    <- bar_width / 1.8
  x_limits   <- range(x_breaks_raw) + c(-half_bw, half_bw)
  pd <- position_dodge(width = bar_width)   # shared dodge
  
  ## bottom panel – containment
  p_cont <- ggplot(df_sub,
                   aes(input_R0, 100 * proportion_contained, colour = factor(surveillance))) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_continuous(breaks = x_breaks,
                       limits = x_limits,
                       expand = c(0, 0)) +
    scale_colour_manual(labels = c(paste0(surveillance_scan[-length(surveillance_scan)], " Hosp."), "No\nVaccine"),
                        values = palette,
                        name = "Surveillance\nThreshold\nTrigger") +
    labs(x = "R0", y = "% Outbreaks\nContained") +
    theme(strip.background = element_rect(fill = "white"),
          plot.margin = margin(t = -10, r = 5, b = 0, l = 5),
          legend.position = "none")
  
  ## top panel – Reff
  df_sub2 <- df_sub %>%
    mutate(plot_bar = ifelse(proportion_contained > 0.25, 0, 1)) %>%
    mutate(Reff_plot = ifelse(plot_bar == 1, Reff_mean, 0))
  p_Reff <- ggplot(df_sub2,
                   aes(input_R0, Reff_mean, fill = factor(surveillance))) +
    geom_bar(stat = "identity",
             position = pd,
             width    = bar_width) +              # <- give bars exact width
    geom_errorbar(aes(ymin = Reff_lower, ymax = Reff_upper),
                  position = pd,
                  width    = bar_width * 0.8,
                  size = 0.35) +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.25) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_x_continuous(breaks = x_breaks,
                       limits = x_limits,
                       expand = c(0, 0)) +
    theme_bw() +
    labs(x = "", y = "REff") +
    theme(strip.background = element_rect(fill = "white"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 10, r = 5, b = 0, l = 5),
          legend.position = "none")
  
  ## bottom panel - time to N
  # p_timetoN <- ggplot(df_sub2, aes(R0, time_to_n_relative, colour = scenario)) + 
  #   geom_line() +
  #   geom_point() +
  #   theme_bw() +
  #   scale_colour_manual(values = rev(palette), 
  #                       labels = c("No Delay", "2 Days", "1 Week", "2 Weeks", "No Vaccination"),
  #                       name = "Vaccine\nProtection\nDelay",
  #                       guide = "none") +
  #   labs(x = "R0", y = "Fold extra time to\nreach epidemic threshold") +
  #   theme(strip.background = element_rect(fill = "white"))
  
  plot_grid(p_Reff, p_cont, nrow = 2, rel_heights = rel_heights,
            align = "v", axis = "l")
}

make_stacked_plot_Reff <- function(df, patho, qe,
                                   rel_heights = c(1, 3),
                                   palette      = c("#CA2E6B", "#88C5EE", "#236897", "#13496E", "black")) {
  
  df_sub <- df %>% 
    filter(vaccine_efficacy_infection == 0.35,
           quarantine_efficacy        == qe,
           pathogen                   == patho)
  
  # constant to make sure everything's aligned
  x_breaks_raw   <- seq(0.75, 2.50, 0.25)   # tick marks
  x_breaks <- seq(1, 2.5, 0.5)
  bar_width  <- 0.2                    # <- same width you pass to geom_bar
  half_bw    <- bar_width / 1.8
  x_limits   <- range(x_breaks_raw) + c(-half_bw, half_bw)
  pd <- position_dodge(width = bar_width)   # shared dodge
  
  ## bottom panel – containment
  p_cont <- ggplot(df_sub,
                   aes(R0, 100 * proportion_contained, colour = scenario)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_x_continuous(breaks = x_breaks,
                       limits = x_limits,
                       expand = c(0, 0)) +
    scale_colour_manual(values = palette, 
                        labels = c("No Delay", "2 Days", "1 Week", "2 Weeks", "No Vaccination"),
                        name = "Vaccine\nProtection\nDelay",
                        guide = guide_legend(reverse = TRUE)) +
    labs(x = "R0", y = "% Outbreaks\nContained") +
    theme(strip.background = element_rect(fill = "white"),
          plot.margin = margin(t = -10, r = 5, b = 0, l = 5),
          legend.position = "none")
  
  ## top panel – Reff
  desired_order <- c("zno_vaccination", "evacc_2weeks_protectDelay", "dvacc_1week_protectDelay", "bvacc_2days_protectDelay",  "avacc_no_delay")
  df_sub2 <- df_sub %>% 
    mutate(scenario = factor(scenario, levels = desired_order))
  
  p_Reff <- ggplot(df_sub2,
                   aes(R0, Reff_mean, fill = scenario)) +
    geom_bar(stat = "identity",
             position = pd,
             width    = bar_width) +              # <- give bars exact width
    geom_errorbar(aes(ymin = Reff_lower, ymax = Reff_upper),
                  position = pd,
                  width    = bar_width * 0.8,
                  size = 0.35) +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.25) +
    scale_fill_manual(values = rev(palette), guide = "none") +
    scale_x_continuous(breaks = x_breaks,
                       limits = x_limits,
                       expand = c(0, 0)) +
    theme_bw() +
    labs(x = "", y = "REff") +
    theme(strip.background = element_rect(fill = "white"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 10, r = 5, b = 0, l = 5),
          legend.position = "none")
  
  ## bottom panel - time to N
  # p_timetoN <- ggplot(df_sub2, aes(R0, time_to_n_relative, colour = scenario)) + 
  #   geom_line() +
  #   geom_point() +
  #   theme_bw() +
  #   scale_colour_manual(values = rev(palette), 
  #                       labels = c("No Delay", "2 Days", "1 Week", "2 Weeks", "No Vaccination"),
  #                       name = "Vaccine\nProtection\nDelay",
  #                       guide = "none") +
  #   labs(x = "R0", y = "Fold extra time to\nreach epidemic threshold") +
  #   theme(strip.background = element_rect(fill = "white"))
  
  plot_grid(p_Reff, p_cont, nrow = 2, rel_heights = rel_heights,
            align = "v", axis = "l")
}
