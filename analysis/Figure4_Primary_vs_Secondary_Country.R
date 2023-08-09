# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))

# OPTION 1: ALL IN TERMS OF CALENDAR DAYS, NO NPIs JUST UNMITIGATED

#### Generate initial sets of scenarios (note placeholder for detection time)
raw_rel_start_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3, 3.5),                   # Basic reproduction number
                                            IFR = 1,                                       # IFR
                                            population_size = 10^10,
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
                                            coverage = 0.8,                                # proportion of the population vaccinated
                                            vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                            seeding_cases = 5,                             # number of initial seeding cases
                                            min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                            min_age_group_index_non_priority = 4,          # index of the youngest age group that *receives* vaccines (4 = 15+)
                                            runtime = 1250)

raw_rel_start_scenarios2 <- expand_grid(raw_rel_start_scenarios, days_source_detection_is_ahead_arrival_secondary = seq(100, -100, -10))
raw_rel_start_scenarios2$detection_time_secondary <- 10 # need to change this to make it R0 specific

for (i in 1:nrow(raw_rel_start_scenarios2)) {
  if (raw_rel_start_scenarios2$days_source_detection_is_ahead_arrival_secondary[i] <= 0) { ## BPSV campaign and specific development starts after pathogen arrival
    raw_rel_start_scenarios2$detection_time[i] <- -raw_rel_start_scenarios2$days_source_detection_is_ahead_arrival_secondary[i]
    raw_rel_start_scenarios2$Rt[i] <- list(raw_rel_start_scenarios2$R0[i])
    raw_rel_start_scenarios2$tt_Rt[i] <- list(0)
  } else if (raw_rel_start_scenarios2$days_source_detection_is_ahead_arrival_secondary[i] > 0) { ## BPSV campaign and specific development starts ahead of pathogen arrival
    raw_rel_start_scenarios2$detection_time[i] <- 0
    raw_rel_start_scenarios2$Rt[i] <- list(c(1, raw_rel_start_scenarios2$R0[i]))
    raw_rel_start_scenarios2$tt_Rt[i] <- list(c(0, raw_rel_start_scenarios2$days_source_detection_is_ahead_arrival_secondary[i]))
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
  saveRDS(model_outputs, "outputs/Figure4_PrimarySecondary_TimingComparison_NoNPIs_CalendarDays.rds")
} else {
  model_outputs <- readRDS("outputs/Figure4_PrimarySecondary_TimingComparison_NoNPIs_CalendarDays.rds")
}

## Joining back in the detection metrics
rel_timing_df <- rel_start_scenarios %>%
  select(scenario_index, days_source_detection_is_ahead_arrival_secondary) %>%
  filter(vaccine_scenario == "specific_only") %>%
  ungroup() %>%
  select(-vaccine_scenario)
model_outputs2 <- model_outputs %>%
  left_join(rel_timing_df, by = c("scenario_index")) %>%
  filter(days_source_detection_is_ahead_arrival_secondary > -60)

## Plotting secondary country days ahead advantage vs bpsv lives saved
population_size <- unique(model_outputs2$population_size)
temp_plot <- ggplot(model_outputs2, aes(x = -days_source_detection_is_ahead_arrival_secondary, y = bpsv_deaths_averted * 1000 / population_size, col = factor(R0))) +
  geom_line() +
  facet_wrap(specific_vaccine_start ~ ., 
             labeller = as_labeller(c(`100`='Specific Vaccine In 100 Days', 
                                      `200`='Specific Vaccine In 200 Days',
                                      `365`='Specific Vaccine In 365 Days'))) +
  theme_bw() +
  labs(x = "Days Pathogen Detection in Source Country is Ahead of Pathogen Arrival in Importing Country",
       y = "Deaths Averted By BPSV\n(Per 1,000 Population)") +
  theme(strip.placement = "outside",
        legend.position = "none",
        strip.background = element_rect(fill="white")) 


### Plotting the dynamics in source, early and late secondary countries
population <- squire:::get_population("Argentina")
population <- 10^6 * population$n / sum(population$n)
mixing_matrix <- squire:::get_mixing_matrix("Argentina")
test_run <- run_booster(time_period = 1300,
                        contact_matrix_set = mixing_matrix,
                        population = population,
                        R0 = 1.15,     
                        tt_R0 = 0, 
                        hosp_bed_capacity = 10^9,                                     
                        ICU_bed_capacity = 10^9,                                       
                        dur_R = 365000000000,                                                        
                        seeding_cases = 10,
                        dur_V = 365000000000,                                              
                        primary_doses = rep(0, 1300),  
                        second_doses = rep(0, 1300),
                        booster_doses = rep(0, 1300))


buffer <- 125
infections <- nimue::format(test_run, compartments = "S", summaries = "infections") %>%
  filter(compartment != "S") %>%
  mutate(t = t + buffer) %>%
  select(-replicate) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(R0 = 1.15) %>%
  complete(t = seq(0, max(t), 1),
           fill = list(value = 0.0920248, R0 = 1.15, compartment = "infections"))

secondary_early_start <- 200
infections2 <- infections %>%
  mutate(t = t + secondary_early_start) %>%
  complete(t = seq(0, max(t), 1),
           fill = list(value = NA, R0 = 1.15, compartment = "infections"))

secondary_late_start <- 650
infections3 <- infections %>%
  mutate(t = t + secondary_late_start) %>%
  complete(t = seq(0, max(t), 1),
           fill = list(value = NA, R0 = 1.15, compartment = "infections"))

base_plot <- ggplot(data = infections3, aes(x = t, y = value)) +
  geom_hline(yintercept = 0, color = "black", linetype="solid") + # Adds vertical line at x=0
  geom_line(linewidth = 1, col = "#70C243") +
  geom_line(data = infections2, aes(x = t, y = value), linewidth = 1, col = "#e5a445") +
  geom_line(data = infections, aes(x = t, y = value), linewidth = 1, col = "#ca4d3d") +
  geom_point(x = 0, y = 0, pch = 21, fill = "#ca4d3d", col = "black", size = 2) +
  geom_point(x = secondary_early_start, y = 0, pch = 21, fill = "#e5a445", col = "black", size = 2) +
  geom_point(x = secondary_late_start, y = 0, pch = 21, fill = "#70C243", col = "black", size = 2) +
  
  theme_bw() +
  coord_cartesian(xlim = c(0, 1500), ylim = c(-200, 1000)) +
  
  geom_segment(x = 0, xend = secondary_early_start, 
               y = -100, yend = -100, 
               linewidth = 1, col = "black", linetype = "solid") +
  geom_point(x = 0, y = -100, pch = 21, fill = "#ca4d3d", col = "black", size = 2) +
  geom_point(x = secondary_early_start, y = -100, pch = 21, fill = "#e5a445", col = "black", size = 2) +
  annotate("text", x = secondary_early_start + 25, y = -100, label = "Earlier Importation", color = "black", hjust = 0) +
  
  geom_segment(x = 0, xend = secondary_late_start, 
               y = -200, yend = -200, 
               linewidth = 1, col = "black", linetype = "solid") +
  geom_point(x = 0, y = -200, pch = 21, fill = "#ca4d3d", col = "black", size = 2) +
  geom_point(x = secondary_late_start, y = -200, pch = 21, fill = "#70C243", col = "black", size = 2) +
  annotate("text", x = secondary_late_start + 25, y = -200, label = "Later Importation", color = "black", hjust = 0) +
  
  labs(x = "", y = "Daily Incidence") +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank()) +
  geom_rect(xmin = 1300, xmax = 1325, ymin = 950, ymax = 1000, fill = "#ca4d3d", col = "#ca4d3d") +
  geom_rect(xmin = 1300, xmax = 1325, ymin = 875, ymax = 925, fill = "#e5a445", col = "#e5a445") +
  geom_rect(xmin = 1300, xmax = 1325, ymin = 800, ymax = 850, fill = "#70C243", col = "#70C243") +
  
  annotate("text", x = 1340, y = 975, label = "Source Country", color = "black", hjust = 0) +
  annotate("text", x = 1340, y = 900, label = "Early Importer", color = "black", hjust = 0) +
  annotate("text", x = 1340, y = 825, label = "Late Importer", color = "black", hjust = 0)


image_path <- magick::image_read_pdf(path = "test_primarySecondary_figure.pdf")
base_ft_inset <- ggdraw() +
  draw_image(image_path, 
             x = 0.0, y = 0.85, hjust = 0, vjust = 1, scale = 2, width = 0.25, height = 0.25) +
  draw_plot(base_plot) 
# 10.91 * 3.45 --> 11 * 3.5

cowplot::plot_grid(base_ft_inset, temp_plot, axis = "hv", align = "tblr",
                   nrow = 2, rel_heights = c(1, 2), labels = c("A", "B"))

### OPTION 2: ALL IN TERMS OF CALENDAR DAYS, DIFFERENT NPIs CONSIDERED

#### Generate initial sets of scenarios (note placeholder for detection time)
raw_rel_start_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3, 3.5),                   # Basic reproduction number
                                            IFR = 1,                                       # IFR
                                            population_size = 10^10,
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
                                            coverage = 0.8,                                # proportion of the population vaccinated
                                            vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                            seeding_cases = 5,                             # number of initial seeding cases
                                            min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                            min_age_group_index_non_priority = 4,          # index of the youngest age group that *receives* vaccines (4 = 15+)
                                            runtime = 1250)

raw_rel_start_scenarios2 <- expand_grid(raw_rel_start_scenarios, days_source_detection_is_ahead_arrival_secondary = seq(100, -100, -2))
raw_rel_start_scenarios2$detection_time_secondary <- 10 # need to change this to make it R0 specific - also be wary when this is longer than time_to_bpsv_coverage - need to change NPIs_arrival_after subset in the NPIs generator

## Generating NPIs based on specific detection times, R0, and other vaccine-related events
NPIs_raw_rel_start <- default_NPI_scenarios_secondary(lockdown_Rt = 0.9, minimal_mandate_reduction = 0.25, 
                                                      NPI_scenarios = c(4, 7, 8, 9), scenarios = raw_rel_start_scenarios2)
rel_start_scenarios <- raw_rel_start_scenarios2 %>%
  full_join(NPIs_raw_rel_start, by = c("R0", "country", "population_size", "bpsv_start", "detection_time_secondary",   # joining by all columns which influence NPI scenario timing
                                       "days_source_detection_is_ahead_arrival_secondary", "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all") %>%
  select(-detection_time.x) %>%
  rename(detection_time = detection_time.y)

## happens with specific 100, 200, 365 and R0 1.5 and days_source_detection_is_ahead_arrival_secondary of 16
## happens with specific 200, 365 (maybe 100, tbd) R0 2 and days_source_detection_is_ahead_arrival_secondary of 16
##days_source_detection_is_ahead_arrival_secondary of 16 appears to be the common denominator

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
  left_join(rel_timing_df, by = c("scenario_index")) # %>%
  #filter(days_source_detection_is_ahead_arrival_secondary > -60)

## Plotting secondary country days ahead advantage vs bpsv lives saved
population_size <- unique(model_outputs2$population_size)
ggplot(model_outputs2, aes(x = days_source_detection_is_ahead_arrival_secondary, y = bpsv_deaths_averted * 1000 / population_size, col = factor(R0))) +
  geom_line() +
  facet_grid(NPI_int ~ specific_vaccine_start) +
  theme_bw() +
  labs(x = "Days Pathogen Detection in Source Country is Ahead of Pathogen Arrival in Importing Country",
       y = "Deaths Averted By BPSV\n(Per 1,000 Population)") +
  theme(strip.placement = "outside",
        legend.position = "none",
        strip.background = element_rect(fill="white")) 




