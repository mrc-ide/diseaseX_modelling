# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))

## Branching-process based calculation of detection times
num_iterations <- 100
detection_hosp <- 1:100 # c(5, 10, 50)
cumulative_window <- 7
delay_hosp <- 10 # check this aligns with squire.page 
R0 <- c(1.5, 2, 2.5, 3, 3.5)
runtime <- c(125, 90, 65, 45, 40)
bp_df <- array(data = NA, dim = c(length(R0), num_iterations, length(detection_hosp), 2))

new_bp <- TRUE
if (new_bp) {
  set.seed(456)
  ## Looping over R0
  for (i in 1:length(R0)) {
    
    ## Looping over number of iterations
    for (j in 1:num_iterations) {
      
      ## Controlling for stochastic fadeout (re-running if it fades out)
      try_again <- 1
      while (try_again == 1) {
        chain_sim_eg <- chain_sim_susc(offspring = "pois", 
                                       mn_offspring = R0[i],
                                       t0 = 0, 
                                       tf = runtime[i],
                                       serial = function(n) {
                                         rgamma(n, shape = 13.4, rate = 2)}, # gamma with mean 6.7 (13.4/2)
                                       pop = 10^6,
                                       initial_immune = 0) 
        if (nrow(chain_sim_eg) > 15) {
          try_again <- 0
        }
      }
      
      ## Calculating daily/rolling 7-day cumulative incidence
      incidence <- chain_sim_eg %>%
        mutate(daily = round(time, digits = 0)) %>%
        group_by(daily) %>%
        summarise(incidence = n()) %>%
        tidyr::complete(daily = min(daily):max(daily), fill = list(incidence = 0)) %>%
        mutate(rolling_incidence = zoo::rollsum(x = incidence, k = 7, na.pad = TRUE, align = "right")) 
      
      ## First time daily hospitalisation incidence goes over detection threshold
      first_incidence_over_threshold <- purrr::map(detection_hosp, ~{
        incidence %>%
          filter(incidence > .x) %>%
          summarise(first_day = min(daily))
      })
      incidence_detection <- unlist(first_incidence_over_threshold) + delay_hosp
      
      ## First time rolling 7-day cumulative hospitalisation incidence goes over detection threshold
      first_cumulative_incidence_over_threshold <- purrr::map(detection_hosp, ~{
        incidence %>%
          filter(rolling_incidence > .x) %>%
          summarise(first_day = min(daily))
      })
      cumulative_incidence_detection <- unlist(first_cumulative_incidence_over_threshold) + delay_hosp
      
      ## Adding these to the array
      bp_df[i, j, , 1] <-  unname(incidence_detection)            # time for daily incidence to eclipse that threshold
      bp_df[i, j, , 2] <-  unname(cumulative_incidence_detection) # time for cumulative 7 day rolling window
      
    }
    print(i)
  }
  
  # Convert to a data frame
  bp_df_long <- reshape2::melt(bp_df)
  bp_df_long$R0 <- R0[(bp_df_long$Var1 - 1) %% length(R0) + 1]
  bp_df_long$iteration <- num_iterations[(bp_df_long$Var2 - 1) %% length(num_iterations) + 1]
  bp_df_long$detection <- detection_hosp[(bp_df_long$Var3 - 1) %% length(detection_hosp) + 1]
  bp_df_long$metric <- c("Daily Incidence", "Cumulative 7-Day Incidence")[(bp_df_long$Var4 - 1) %% 2 + 1]
  bp_df_long <- bp_df_long[, -(1:4)]
  bp_df_long <- bp_df_long %>%
    mutate(detection_threshold = case_when(detection == detection_hosp[1] ~ "low",
                                           detection == detection_hosp[2] ~ "moderate",
                                           detection == detection_hosp[3] ~ "high"))
  saveRDS(object = bp_df_long, file = "outputs/Figure1_bp_detection_times.rds")
} else {
  bp_df_long <- readRDS("outputs/Figure1_bp_detection_times.rds")
}

bp_df_long_mean <- bp_df_long %>%
  filter(!is.infinite(value)) %>%
  group_by(R0, detection, metric, detection_threshold) %>%
  summarise(mean = mean(value))
ggplot(bp_df_long_mean, aes(x = R0, y = mean, col = factor(detection))) + 
  geom_line() +
  facet_grid(.~metric)

# Generate parameter combinations for model running
primary_country_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3, 3.5),                   # Basic reproduction number
                                              IFR = 1,                                       # IFR
                                              population_size = 10^10,                       # population size
                                              Tg = 6.7,                                      # Tg
                                              detection_time = 1,                            # detection time - PLACECHOLDER FOR NOW
                                              bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                              bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                              specific_vaccine_start = c(100, 200, 365),     # specific vaccine distribution start (time after detection time)
                                              specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                              efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                              efficacy_disease_bpsv = 0.75,                  # vaccine efficacy against disease - BPSV
                                              efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                              efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                              dur_R = 365000000,                             # duration of infection-induced immunity
                                              dur_bpsv = 365000000,                          # duration of BPSV vaccine immunity
                                              dur_spec = 365000000,                          # duration of disease-specific vaccine immunity
                                              coverage = 0.8,                                # proportion of the population vaccinated
                                              vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                              min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                              min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)

# Linking these scenarios with the R0-specific 
primary_country_scenarios <- expand_grid(primary_country_scenarios,
                                         detection_threshold = c("low", "moderate", "high"))


# NPI Relevant Parameters
lockdown_Rt <- 0.9                   # Rt achieved under lockdown
minimal_mandate_reduction <- 0.25    # Fold-reduction in R0 achieved under minimal mandate restrictions
NPIs_bpsv_eff <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                       NPI_scenarios = 1:9, scenarios = raw_bpsv_efficacy_scenarios)

# Combining NPI scenarios with parameter combos into one overall dataframe for model running 
bpsv_eff_scenarios <- raw_bpsv_efficacy_scenarios %>%
  full_join(NPIs_bpsv_eff, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                  "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")
