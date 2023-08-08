# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))

# Loading in bp based detection and calculating detection times for the the different R0 values
bp_df_long <- readRDS("outputs/Figure1_bp_detection_times.rds")
prob_hosp <- squire.page:::probs_booster$prob_hosp
arg_pop <- squire::get_population("Argentina")
IHR <- sum(prob_hosp * arg_pop$n / sum(arg_pop$n)) 
num_hosp <- 1:20
detection_hosp <- round(num_hosp / IHR, digits = 0)
num_hosp <- c(1, 5, 10, 20)
infection_thresholds <- detection_hosp[num_hosp]
bp_df_mean_subset <- bp_df_long %>%
  filter(!is.infinite(value),
         detection %in% infection_thresholds) %>%
  group_by(R0, detection, metric) %>%
  summarise(mean = round(mean(value))) 

# Generate parameter combinations for model running

#### Generate initial sets of scenarios (note placeholder for detection time)
raw_rel_start_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3, 3.5),                   # Basic reproduction number
                                            IFR = 1,                                       # IFR
                                            population_size = 10^10,
                                            Tg = 6.7,                                      # Tg
                                            detection_time = 1,                            # PLACEHOLDER
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
                                            seeding_cases = 5,                             # number of initial seeding cases
                                            min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                            min_age_group_index_non_priority = 4,          # index of the youngest age group that *receives* vaccines (4 = 15+)
                                            runtime = 1250)

raw_rel_start_scenarios2 <- expand_grid(raw_rel_start_scenarios, days_ahead_arrival = seq(100, -100, -1))

for (i in 1:nrow(raw_rel_start_scenarios2)) {
  if (raw_rel_start_scenarios2$days_ahead_arrival[i] <= 0) { ## BPSV campaign and specific development starts after pathogen arrival
    
    raw_rel_start_scenarios2$detection_time[i] <- -raw_rel_start_scenarios2$days_ahead_arrival[i]
    raw_rel_start_scenarios2$bpsv_start[i] <- 0
    raw_rel_start_scenarios2$bpsv_protection_delay[i] <- 7
    raw_rel_start_scenarios2$specific_vaccine_start[i] <- raw_rel_start_scenarios2$specific_vaccine_start[i]
    raw_rel_start_scenarios2$specific_protection_delay[i] <- 7
    raw_rel_start_scenarios2$Rt[i] <- list(raw_rel_start_scenarios2$R0[i])
    raw_rel_start_scenarios2$tt_Rt[i] <- list(0)
    
  } else if (raw_rel_start_scenarios2$days_ahead_arrival[i] > 0) { ## BPSV campaign and specific development starts ahead of pathogen arrival
    
    raw_rel_start_scenarios2$detection_time[i] <- 0
    raw_rel_start_scenarios2$bpsv_start[i] <- 0
    raw_rel_start_scenarios2$bpsv_protection_delay[i] <- 7
    raw_rel_start_scenarios2$specific_vaccine_start[i] <- raw_rel_start_scenarios2$specific_vaccine_start[i]
    raw_rel_start_scenarios2$specific_protection_delay[i] <- 7
    raw_rel_start_scenarios2$Rt[i] <- list(c(raw_rel_start_scenarios2$R0[i], 0.9, raw_rel_start_scenarios2$R0[i]))
    raw_rel_start_scenarios2$tt_Rt[i] <- list(c(0, 1, raw_rel_start_scenarios2$days_ahead_arrival[i]))

  } else {
    stop("what???")
  }
}

## Creating overall output and index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vars_for_index <- c(variable_columns(raw_rel_start_scenarios2))
rel_start_scenarios <- raw_rel_start_scenarios2 %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
cores <- parallel::detectCores() - 2
fresh_run <- TRUE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(rel_start_scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/Figure4_PrimarySecondary_TimingComparison.rds")
} else {
  model_outputs <- readRDS("outputs/Figure4_PrimarySecondary_TimingComparison.rds")
}

## Joining back in the detection metrics
rel_timing_df <- rel_start_scenarios %>%
  select(scenario_index, all_of(vars_for_index)) %>%
  filter(vaccine_scenario == "specific_only") %>%
  ungroup() %>%
  select(-vaccine_scenario) %>%
  select(R0, scenario_index, specific_vaccine_start, days_ahead_arrival)

model_outputs2 <- model_outputs %>%
  left_join(rel_timing_df, by = c("R0", "scenario_index", "specific_vaccine_start")) %>%
  filter(days_ahead_arrival >= -30)
population <- unique(model_outputs2$population_size)

ggplot(model_outputs2, aes(x = -days_ahead_arrival, y = deaths_spec, col = factor(R0))) +
  geom_line() +
  facet_wrap(.~specific_vaccine_start) 

ggplot(model_outputs2, aes(x = -days_ahead_arrival, 
                           y = bpsv_deaths_averted * 1000 / population, col = factor(R0))) +
  geom_line() +
  facet_wrap(.~specific_vaccine_start) +
  theme_bw() +
  labs(x = "Days Vaccine Development/Deployment is Ahead of Pathogen Arrival",
       y = "Deaths Averted By BPSV (Per 1,000 Population)")














