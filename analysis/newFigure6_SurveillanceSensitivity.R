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

new_bp <- TRUE
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
  saveRDS(object = bp_df_long, file = "outputs/Figure6_Surveillance_Exploration/NEW_Fig6_surveillanceSensitivity.rds")
} else {
  bp_df_long <- readRDS("outputs/Figure6_Surveillance_Exploration/NEW_Fig6_surveillanceSensitivity.rds")
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

## Running
detection_times <- bp_df_long_mean %>%
  filter(metric == "Daily Incidence") %>%
  rename(detection_time = mean) %>%
  select(R0, detection_time, detection_hosp)
detection_times$detection_time  <- round(detection_times$detection_time)

raw_baseline_scenarios <- create_scenarios(R0 = R0, specific_vaccine_start = c(100, 250, 365)) %>%
  select(-detection_time)
baseline_scenarios <- expand_grid(raw_baseline_scenarios,
                                  detection_hosp = unique(detection_times$detection_hosp)) %>%
  left_join(detection_times, by = c("R0", "detection_hosp")) 

NPIs <- default_NPI_scenarios(lockdown_Rt = default$lockdown_Rt, minimal_mandate_reduction = default$minimal_mandate_reduction, 
                              NPI_scenarios = 1:9, scenarios = baseline_scenarios)
scenarios_NPIs <- baseline_scenarios %>%
  full_join(NPIs, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence timing of NPI scenarios
                         "specific_vaccine_start", "vaccination_rate_bpsv", "vaccination_rate_spec",
                         "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all")

detection_times2 <- detection_times %>%
  select(R0, detection_time, detection_hosp)
baseline_scenarios_reduced <- scenarios_NPIs %>%
  semi_join(detection_times2, by = c("R0", "detection_hosp", "detection_time"))
# R0 * NPI * spec start * 2 vaccination scenarios * 25 detection thresholds 3 * 3 * 3 * 2 * 25

## Creating index for output
vars_for_index <- c(variable_columns(baseline_scenarios_reduced), "NPI_int")
scenarios <- baseline_scenarios_reduced %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
fresh_run <- TRUE
if (fresh_run) {
  plan(multisession, workers = 6) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = 2)
  saveRDS(model_outputs, "outputs/Figure6_Surveillance_Exploration/NEW_Figure6_PrimaryCountry_Surveillance.rds")
} else {
  model_outputs <- readRDS("outputs/Figure6_Surveillance_Exploration/NEW_Figure6_PrimaryCountry_Surveillance.rds")
}

scenarios2 <- scenarios %>%
  select(scenario_index, detection_hosp) %>%
  group_by(scenario_index) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  ungroup() %>%
  select(scenario_index, detection_hosp)

model_outputs2 <- model_outputs %>%
  filter(specific_vaccine_start %in% c(100, 250)) %>%
  filter(R0 %in% c(1.5, 2.5, 3.5)) %>%
  filter(NPI_int %in% c(4, 7, 8)) %>%
  left_join(scenarios2, by = "scenario_index")

a <- ggplot(subset(bp_df_long_mean, metric == "Daily Incidence"), 
            aes(x = detection_hosp, y = mean, fill = factor(R0))) + 
  geom_line(aes(col = factor(R0))) +
  geom_point(pch = 21, size = 2) +
  facet_wrap(.~metric) +
  guides(col = guide_legend(title = "R0"),
         fill = "none") +
  theme_bw() +
  labs(x = "Hospitalisation Detection Threshold", y = "Time to Detection (Days)") +
  theme(strip.background = element_rect(fill="#F5F5F5"),
        legend.position = "none")

b <- ggplot(model_outputs2) +
  geom_line(aes(x = detection_hosp, y = 1000 * bpsv_deaths_averted / unique(baseline_scenarios$population_size), 
                col = factor(R0)), size = 1) +
  geom_point(aes(x = detection_hosp, y = 1000 * bpsv_deaths_averted / unique(baseline_scenarios$population_size), 
                fill = factor(R0)), size = 1.5, pch = 21) +
  theme_bw() +
  facet_grid(specific_vaccine_start ~ NPI_int,
             labeller = as_labeller(c(`4`='Minimal NPIs', `7`='Moderate NPIs', `8`='Stringent NPIs', 
                                      `100` = "Specific Vaccine in 100 Days", `250` = "Specific Vaccine in 250 Days"))) +
  labs(x = "Hospitalisation Detection Threshold", y = "Deaths Averted by BPSV\n (Per 1,000 Population)") +
  guides(col = guide_legend(title = "R0"),
         fill = "none") +
  theme(strip.background = element_rect(fill="#F5F5F5"),
        legend.position = "right")

surv_plot <- cowplot::plot_grid(a, b, nrow = 1, rel_widths = c(1, 2),
                                labels = c("A", "B"))
ggsave(filename = "figures/Figure_6_Surveillance/NEW_Fig6_Surveillance_sensitivity.pdf",
       plot = surv_plot,
       width = 11,
       height = 4.8)


### OPTION 2: ALL IN TERMS OF CALENDAR DAYS, DIFFERENT NPIs CONSIDERED
##### NEED TO TRY OUT ADDING IN R0-DEPENDENT TIME OF DETECTION

#### Generate initial sets of scenarios (note placeholder for detection time)
raw_rel_start_scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5),                         # Basic reproduction number
                                            IFR = 1,                                       # IFR
                                            population_size = 10^10,                       # population size
                                            Tg = 6.7,                                      # Tg
                                            detection_time = 1,                            # PLACEHOLDER
                                            bpsv_start = 0,                                # BPSV distribution start (time after detection time)
                                            bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                            specific_vaccine_start = c(100, 250, 365),     # specific vaccine distribution start (time after detection time)
                                            specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                            efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                            efficacy_disease_bpsv = 0.75,                  # vaccine efficacy against disease - BPSV
                                            efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                            efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                            dur_R = 365000000,                             # duration of infection-induced immunity
                                            dur_bpsv = 365000000,                          # duration of BPSV vaccine immunity
                                            dur_spec = 365000000,                          # duration of disease-specific vaccine immunity
                                            coverage_bpsv = 0.8,                           # proportion of the population vaccinated
                                            coverage_spec = 0.8,                           # proportion of the population vaccinated
                                            vaccination_rate_bpsv = 0.035, 
                                            vaccination_rate_spec = 0.035,
                                            seeding_cases = 5,                             # number of initial seeding cases
                                            min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                            min_age_group_index_non_priority = 4,          # index of the youngest age group that *receives* vaccines (4 = 15+)
                                            runtime = 1250)

sequence <-  seq(100, -100, -5)
raw_rel_start_scenarios2 <- expand_grid(raw_rel_start_scenarios, days_source_detection_is_ahead_arrival_secondary = sequence[-which(sequence == 0)])
raw_rel_start_scenarios2$detection_time_secondary <- 7 # Are we okay for this to be in calendar days and independent of R0? (Potentially yes if we think enhanced surveillance is the play?)

## Generating NPIs based on specific detection times, R0, and other vaccine-related events
NPIs_raw_rel_start <- default_NPI_scenarios_secondary(lockdown_Rt = 0.9, minimal_mandate_reduction = 0.25, 
                                                      NPI_scenarios = c(4, 7, 8, 9), 
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

## Running the model and summarising the output
cores <- parallel::detectCores() - 2
fresh_run <- TRUE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(rel_start_scenarios2, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/Figure_7_SecondaryCountry/Figure7_PrimarySecondary_TimingComparison_DifferentNPIs_CalendarDays.rds")
} else {
  model_outputs <- readRDS("outputs/Figure_7_SecondaryCountry/Figure7_PrimarySecondary_TimingComparison_DifferentNPIs_CalendarDays.rds")
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
                                  NPI_int == 8 ~ "dStringent NPIs",
                                  NPI_int == 9 ~ "aNo NPIs")) %>%
  mutate(NPI_scenario = factor(NPI_scenario, levels = c("aNo NPIs", "bMinimal NPIs", "cModerate NPIs", "dStringent NPIs")))

text_data <- data.frame(
  x = rep(-100, nlevels(model_outputs2$NPI_scenario)), # x position for the text
  y = rep(4, nlevels(model_outputs2$NPI_scenario)), # y position for the text
  label = rep("Before\nImportation", nlevels(model_outputs2$NPI_scenario)),
  NPI_scenario = levels(model_outputs2$NPI_scenario), # Repeat for each level of NPI_scenario
  specific_vaccine_start = 100
)

low <- -60
mid <- 10
high <- 60

annotation_data1 <- data.frame(
  NPI_scenario = c("aNo NPIs", "bMinimal NPIs", "cModerate NPIs", "dStringent NPIs"),
  specific_vaccine_start = 200,
  xmin = c(-high - 5, -high - 5), # x min position
  xmax = c(-high + 5, -high + 5), # x max position
  ymin = c(-Inf, -Inf), # y min position
  ymax = c(Inf, Inf) # y max position
)
annotation_data2 <- data.frame(
  NPI_scenario = c("aNo NPIs", "bMinimal NPIs", "cModerate NPIs", "dStringent NPIs"),
  specific_vaccine_start = 200,
  xmin = c(-mid - 5, -mid - 5), # x min position
  xmax = c(-mid + 5, -mid + 5), # x max position
  ymin = c(-Inf, -Inf), # y min position
  ymax = c(Inf, Inf) # y max position
)
annotation_data3 <- data.frame(
  NPI_scenario = c("aNo NPIs", "bMinimal NPIs", "cModerate NPIs", "dStringent NPIs"),
  specific_vaccine_start = 200,
  xmin = c(-low - 5, -low - 5), # x min position
  xmax = c(-low + 5, -low + 5), # x max position
  ymin = c(-Inf, -Inf), # y min position
  ymax = c(Inf, Inf) # y max position
)

deaths_averted <- ggplot(subset(model_outputs2, NPI_int != 9 & 
                                  -days_source_detection_is_ahead_arrival_secondary > -60 & 
                                  -days_source_detection_is_ahead_arrival_secondary < 60)) +
  annotate("rect", xmin = -105, xmax = 0, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.25) +
  annotate("rect", xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, fill = "white", alpha = 0.1) +
  geom_line(aes(x = -days_source_detection_is_ahead_arrival_secondary, 
                y = bpsv_deaths_averted * 1000 / population_size, col = factor(R0))) +
  geom_point(aes(x = -days_source_detection_is_ahead_arrival_secondary, 
                 y = bpsv_deaths_averted * 1000 / population_size, col = factor(R0))) +
  facet_grid(NPI_scenario ~ specific_vaccine_start,
             labeller = as_labeller(c(`100`='Specific Vaccine In 100 Days', 
                                      `220`='Specific Vaccine In 250 Days',
                                      `365`='Specific Vaccine In 365 Days',
                                      `bMinimal NPIs`="Minimal NPIs",
                                      `cModerate NPIs`="Moderate NPIs", 
                                      `dStringent NPIs`="Stringent NPIs",
                                      `aNo NPIs`="No NPIs"))) +
  theme_bw() +
  scale_colour_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
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

ggsave(filename = "figures/Figure_7_SecondaryCountry_Exploration/NEW_Figure7_Secondary_TimingComparison_DifferentNPIs_CalendarDays.pdf",
       plot = deaths_averted,
       width = 9.5,
       height = 6)

secondary_legend <- cowplot::plot_grid(deaths_averted + theme(legend.position = "right"))
ggsave(filename = "figures/Figure_7_SecondaryCountry_Exploration/NEW_Figure7_legend.pdf",
       plot = secondary_legend,
       width = 9.5,
       height = 6)
