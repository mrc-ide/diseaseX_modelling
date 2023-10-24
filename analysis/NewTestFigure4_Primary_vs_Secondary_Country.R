# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/helper_functions.R"))
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/runs_sars_x_secondary_country.R"))

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
                                            specific_vaccine_start = c(100, 200, 365),     # specific vaccine distribution start (time after detection time)
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
cores <- parallel::detectCores() - 3
fresh_run <- TRUE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(rel_start_scenarios2, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/Figure4_PrimarySecondary_TimingComparison_DifferentNPIs_CalendarDays.rds")
} else {
  model_outputs <- readRDS("outputs/Figure4_PrimarySecondary_TimingComparison_DifferentNPIs_CalendarDays.rds")
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

deaths_averted <- ggplot(model_outputs2) +
  annotate("rect", xmin = -105, xmax = 0, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.25) +
  annotate("rect", xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, fill = "white", alpha = 0.1) +
  geom_line(aes(x = -days_source_detection_is_ahead_arrival_secondary, 
                y = bpsv_deaths_averted * 1000 / population_size, col = factor(R0))) +
  geom_point(aes(x = -days_source_detection_is_ahead_arrival_secondary, 
                 y = bpsv_deaths_averted * 1000 / population_size, col = factor(R0))) +
  facet_grid(NPI_scenario ~ specific_vaccine_start,
             labeller = as_labeller(c(`100`='Specific Vaccine In 100 Days', 
                                      `200`='Specific Vaccine In 200 Days',
                                      `365`='Specific Vaccine In 365 Days',
                                      `bMinimal NPIs`="Minimal NPIs",
                                      `cModerate NPIs`="Moderate NPIs", 
                                      `dStringent NPIs`="Stringent NPIs",
                                      `aNo NPIs`="No NPIs")),
             switch = "y") +
  theme_bw() +
  scale_colour_manual(values = c("#B8336A", "#726DA8", "#42A1B6")) +
  scale_y_continuous(position = "right") +
  labs(x = "Time of Pathogen Detection in Source Country\nRelative to Pathogen Arrival in Importing Country",
       y = "Deaths Averted By BPSV (Per 1,000 Population)") +
  geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed") +
  theme(strip.placement = "outside",
        legend.position = "right",
        strip.background = element_rect(fill="white")) +
  coord_cartesian(xlim = c(-100, 100))
