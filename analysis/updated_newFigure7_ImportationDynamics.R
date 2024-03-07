# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/helper_functions.R"))
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/runs_sars_x_secondary_country.R"))
source(here::here("functions/branching_process.R"))
source(here::here("functions/branching_process2.R"))

default <- define_default_params()

## Figure for Primary Country and Showing How Surveillance Sensitivity Influences BPSV Impact

## Calculate rough IHR (similar value to Knock et al which had about 2% for London)
prob_hosp <- squire.page.sarsX:::probs_booster$prob_hosp
arg_pop <- squire::get_population("Argentina")
IHR <- sum(prob_hosp * arg_pop$n / sum(arg_pop$n)) 

## Branching-process based calculation of detection times
num_iterations <- 10
num_hosp <- 1:25
detection_hosp <- round(num_hosp / IHR, digits = 0)
cumulative_window <- 7
delay_hosp <- 10 # check this aligns with squire.page 
R0 <- c(1.5, 2.5, 3.5)
runtime <- c(140, 70, 55)
bp_df <- array(data = NA, dim = c(length(R0), num_iterations, max(num_hosp), 2))
dim(bp_df)

new_bp <- FALSE
if (new_bp) {
  set.seed(456)
  ## Looping over R0
  for (i in 1:length(R0)) {
    
    ## Looping over number of iterations
    for (j in 8:num_iterations) {
      
      ## Controlling for stochastic fadeout (re-running if it fades out)
      try_again <- 1
      while (try_again == 1) {
        chain_sim_eg <- chain_sim_susc(offspring = "pois", 
                                       mn_offspring = R0[i],
                                       t0 = 0, 
                                       tf = runtime[i],
                                       serial = function(n) {
                                         rgamma(n, shape = 13.4, rate = 2)}, # gamma with mean 6.7 (13.4/2)
                                       pop = 10^9,
                                       initial_immune = 0) 
        if (nrow(chain_sim_eg) > 100) {
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
  bp_df_long <- bp_df_long 
  saveRDS(object = bp_df_long, file = "outputs/Figure7_SecondaryCountry/Fig7_ImportationDynamics_DetectionTimes.rds")
} else {
  bp_df_long <- readRDS("outputs/Figure7_SecondaryCountry/Fig7_ImportationDynamics_DetectionTimes.rds")
}

# Plotting all the results for the branching process detection times
bp_df_long_mean <- bp_df_long %>%
  filter(!is.infinite(value)) %>%
  group_by(R0, detection, metric) %>%
  summarise(mean = mean(value)) %>%
  mutate(detection_hosp = round(detection * IHR))
bp_detection_time_plot <- ggplot(bp_df_long_mean, aes(x = R0, y = mean, col = factor(detection_hosp))) + 
  geom_line() +
  facet_grid(.~metric)
ggplot(bp_df_long_mean, aes(x = detection_hosp, y = mean, col = factor(R0))) + 
  geom_line() +
  facet_grid(.~metric)

### OPTION 2: ALL IN TERMS OF CALENDAR DAYS, DIFFERENT NPIs CONSIDERED
##### NEED TO TRY OUT ADDING IN R0-DEPENDENT TIME OF DETECTION

#### Generate initial sets of scenarios (note placeholder for detection time)
raw_rel_start_scenarios <- create_scenarios(R0 = R0, specific_vaccine_start = c(100, 250))
sequence <-  seq(100, -100, -5)
raw_rel_start_scenarios2 <- expand_grid(raw_rel_start_scenarios, days_source_detection_is_ahead_arrival_secondary = sequence[-which(sequence == 0)])
raw_rel_start_scenarios2$detection_time_secondary <- 7 # Are we okay for this to be in calendar days and independent of R0? (Potentially yes if we think enhanced surveillance is the play?)

## Generating NPIs based on specific detection times, R0, and other vaccine-related events
NPIs_raw_rel_start <- default_NPI_scenarios_secondary(lockdown_Rt = 0.9, minimal_mandate_reduction = 0.25, 
                                                      NPI_scenarios = c(4, 7, 8), 
                                                      scenarios = raw_rel_start_scenarios2)
rel_start_scenarios <- raw_rel_start_scenarios2 %>%
  full_join(NPIs_raw_rel_start, by = c("R0", "country", "population_size", "bpsv_start", "detection_time_secondary",   # joining by all columns which influence NPI scenario timing
                                       "days_source_detection_is_ahead_arrival_secondary", "specific_vaccine_start", "vaccination_rate_bpsv", "vaccination_rate_spec", 
                                       "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all") %>%
  select(-detection_time.x) %>%
  rename(detection_time = detection_time.y)

## Creating overall output and index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vars_for_index <- c(variable_columns(rel_start_scenarios))
rel_start_scenarios2 <- rel_start_scenarios %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())
# R0 * NPI * spec start * 2 vaccination scenarios * (length(sequence) - 1)
3 * 3 * 2 * 2 * (length(sequence) - 1)

## Running the model and summarising the output
cores <- parallel::detectCores() - 2
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(rel_start_scenarios2, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/Figure7_SecondaryCountry/updated_NEW_Figure7_ImportationDynamics_ScenarioRuns.rds")
} else {
  model_outputs <- readRDS("outputs/Figure7_SecondaryCountry/updated_NEW_Figure7_ImportationDynamics_ScenarioRuns.rds")
}

## Joining back in the detection metrics
rel_timing_df <- rel_start_scenarios2 %>%
  select(scenario_index, days_source_detection_is_ahead_arrival_secondary) %>%
  filter(vaccine_scenario == "specific_only") %>%
  ungroup() %>%
  select(-vaccine_scenario)
model_outputs2 <- model_outputs %>%
  left_join(rel_timing_df, by = c("scenario_index")) 

## Plotting secondary country days ahead advantage vs bpsv lives saved
population_size <- unique(model_outputs2$population_size)
model_outputs2 <- model_outputs2 %>%
  mutate(NPI_scenario = case_when(NPI_int == 4 ~ "bMinimal NPIs", 
                                  NPI_int == 7 ~ "cModerate NPIs",
                                  NPI_int == 8 ~ "dStringent NPIs")) %>%
  mutate(NPI_scenario = factor(NPI_scenario, levels = c("bMinimal NPIs", "cModerate NPIs", "dStringent NPIs")))

text_data <- data.frame(
  x = rep(-100, nlevels(model_outputs2$NPI_scenario)), # x position for the text
  y = rep(4, nlevels(model_outputs2$NPI_scenario)), # y position for the text
  label = rep("Before\nImportation", nlevels(model_outputs2$NPI_scenario)),
  NPI_scenario = levels(model_outputs2$NPI_scenario), # Repeat for each level of NPI_scenario
  specific_vaccine_start = 100
)

model_outputs2 <- model_outputs2 %>%
  mutate(before_after = ifelse(-days_source_detection_is_ahead_arrival_secondary < 0, "Before", "After")) %>%
  filter(NPI_int != 9 & 
           -days_source_detection_is_ahead_arrival_secondary > -60 & 
           -days_source_detection_is_ahead_arrival_secondary < 60)

deaths_averted <- ggplot(model_outputs2) +
  annotate("rect", xmin = -105, xmax = 0, ymin = -Inf, ymax = Inf, fill = "white", alpha = 0.1) +
  annotate("rect", xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.1) +
  geom_line(aes(x = days_source_detection_is_ahead_arrival_secondary, 
                y = bpsv_deaths_averted * 1000 / population_size, 
                col = interaction(factor(R0), before_after))) +
  geom_point(aes(x = days_source_detection_is_ahead_arrival_secondary, 
                 y = bpsv_deaths_averted * 1000 / population_size,
                 col = interaction(factor(R0), before_after))) +
  facet_grid(NPI_scenario ~ specific_vaccine_start,
             labeller = as_labeller(c(`100`='Specific Vaccine In 100 Days', 
                                      `250`='Specific Vaccine In 250 Days',
                                      `bMinimal NPIs`="Minimal NPIs",
                                      `cModerate NPIs`="Moderate NPIs", 
                                      `dStringent NPIs`="Stringent NPIs"))) +
  theme_bw() +
  scale_colour_manual(values = c( "#EBF2D3", "#DBE8B9", "#BED78B", "#E7CAAA", "#E3B273", "#DB9939")) +
  scale_y_continuous(position = "left") +
  scale_x_continuous(breaks = c(-60, -30, 0, 30, 60),
                     labels = c("+60", "+30", "0", "-30", "-60")) +
  labs(x = "Days Pathogen Detection (in Source Country) is Ahead of Importation (to Secondary Country)",
       y = "Deaths Averted By BPSV in Secondary Country (Per 1,000 Population)") +
  geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed") +
  theme(strip.placement = "outside",
        legend.position = "none",
        strip.background = element_rect(fill="white")) +
  coord_cartesian(xlim = c(-60, 60))

# ggsave(filename = "figures/Figure_7_SecondaryCountry_Exploration/NEW_Figure7_Secondary_TimingComparison_DifferentNPIs_CalendarDays.pdf",
#        plot = deaths_averted,
#        width = 9.5,
#        height = 6)
# 
# secondary_legend <- cowplot::plot_grid(deaths_averted + theme(legend.position = "right"))
# ggsave(filename = "figures/Figure_7_SecondaryCountry_Exploration/updated_NEW_Figure7_legend.pdf",
#        plot = secondary_legend,
#        width = 9.5,
#        height = 6)

deaths_averted2 <- deaths_averted + 
  facet_grid(specific_vaccine_start ~ NPI_scenario,
             labeller = as_labeller(c(`100`='Specific Vaccine In 100 Days', 
                                      `250`='Specific Vaccine In 250 Days',
                                      `bMinimal NPIs`="Minimal NPIs",
                                      `cModerate NPIs`="Moderate NPIs", 
                                      `dStringent NPIs`="Stringent NPIs")))
ggsave(filename = "figures/Figure_7_SecondaryCountry_Exploration/updated_NEW_Figure7_legend_altOrientation.pdf",
       plot = deaths_averted2,
       width = 10.5,
       height = 5)
