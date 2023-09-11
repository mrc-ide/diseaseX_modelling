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
  summarise(mean = mean(value)) 

# NPI Relevant Parameters
lockdown_Rt <- 0.9                   # Rt achieved under lockdown
minimal_mandate_reduction <- 0.25    # Fold-reduction in R0 achieved under minimal mandate restrictions

# Calculating vaccination rates for use in sensitivity analyses
days_to_bpsv_coverage <- seq(20, 120, 20)
population_size <- 10^10
min_age_group_index_priority <- 13
standard_pop <- generate_standard_pop(country = "Argentina", population_size = population_size)
coverage <- 0.8
elderly_pop_to_vaccinate <- sum(standard_pop[min_age_group_index_priority:17]) * coverage 
bpsv_vaccination_rate <- elderly_pop_to_vaccinate / (days_to_bpsv_coverage * population_size / 7)
bpsv_coverage_time_df <- data.frame(days_to_bpsv_coverage  = days_to_bpsv_coverage, vaccination_rate_bpsv = bpsv_vaccination_rate)

## Generating the scenarios
raw_vacc_rate_scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5),                         # Basic reproduction number
                                            IFR = 1,                                       # IFR
                                            population_size = 10^10,                       # population size
                                            Tg = 7,                                        # Tg
                                            detection_time = 1,                            # PLACEHOLDER FOR NOW
                                            bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                            bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                            specific_vaccine_start = c(150, 220, 365),     # specific vaccine distribution start (time after detection time)
                                            specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                            efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                            efficacy_disease_bpsv = 0.75,                  # vaccine efficacy against disease - BPSV
                                            efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                            efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                            dur_R = 365000000,                             # duration of infection-induced immunity
                                            dur_bpsv = 365000000,                          # duration of BPSV vaccine immunity
                                            dur_spec = 365000000,                          # duration of disease-specific vaccine immunity
                                            coverage = coverage,                           # proportion of the population vaccinated
                                            vaccination_rate = bpsv_vaccination_rate,      # vaccination rate per week as percentage of population
                                            vaccination_rate_spec = 0.035,                 # disease-specific vaccination rate per week as percentage of population
                                            min_age_group_index_priority = min_age_group_index_priority,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                            min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)

## Join detection time dataframe (note currently have all combos of R0 and detection time and we want specific pairings)
raw_vacc_rate_scenarios2 <- expand_grid(raw_vacc_rate_scenarios,
                                        detection_threshold = unique(bp_df_mean_subset$detection)) %>%
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean) 

## Generating NPIs based on specific detection times, R0, and other vaccine-related events
NPIs_vacc_rate <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, 
                                        minimal_mandate_reduction = minimal_mandate_reduction, 
                                        NPI_scenarios = c(4, 7, 8, 9), 
                                        scenarios = raw_vacc_rate_scenarios2)
vacc_rate_scenarios <- raw_vacc_rate_scenarios2 %>%
  full_join(NPIs_vacc_rate, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",         # joining by all columns which influence NPI scenario timing
                                   "specific_vaccine_start", "vaccination_rate_bpsv", "vaccination_rate_spec", 
                                   "coverage", "min_age_group_index_priority"), multiple = "all")

## Filtering the above to only select R0 and detection time pairs that actually occurred (the expand grid call above generated all combos)
R0_detection_time_pairs <- bp_df_mean_subset %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time, metric)
vacc_rate_scenarios2 <- vacc_rate_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "immunity_duration") 

## Creating overall output and index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vars_for_index <- c(variable_columns(vacc_rate_scenarios), "NPI_int")
vacc_rate_scenarios3 <- vacc_rate_scenarios2 %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
cores <- parallel::detectCores() - 2
fresh_run <- TRUE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(vacc_rate_scenarios3, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/vaccination_rate_exploration_scenarios.rds")
} else {
  model_outputs <- readRDS("outputs/vaccination_rate_exploration_scenarios.rds")
}

## Joining back in the detection metrics
detection_df <- vacc_rate_scenarios3 %>%
  select(scenario_index, specific_vaccine_start, all_of(vars_for_index)) %>%
  filter(vaccine_scenario == "specific_only") %>%
  ungroup() %>%
  select(-vaccine_scenario) %>%
  select(scenario_index, specific_vaccine_start, NPI_int, detection_time, detection_threshold, metric)

model_outputs2 <- model_outputs %>%
  left_join(detection_df, by = c("scenario_index", "specific_vaccine_start", "detection_time", "NPI_int")) %>%
  mutate(detection_threshold_hosp = round(detection_threshold * IHR)) %>%
  mutate(detection_timing = case_when(detection_threshold_hosp == 1 ~ "Early",
                                      detection_threshold_hosp == 5 ~ "Intermediate",
                                      detection_threshold_hosp == 10 ~ "Late",
                                      detection_threshold_hosp == 20 ~ "Very Late"))

specific_vaccine_start_fixed <- 365
vacc_rate_plotting <- model_outputs2 %>%
  filter(metric == "Daily Incidence") %>%
  group_by(R0, specific_vaccine_start, NPI_int, detection_time, detection_threshold_hosp, vaccination_rate_bpsv) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted) * 1000 / population_size,
            max_deaths_averted = max(bpsv_deaths_averted) * 1000 / population_size,
            central_deaths_averted = bpsv_deaths_averted * 1000 / population_size,
            perc_deaths_averted = 100 * bpsv_deaths_averted / deaths_spec,
            total_deaths_spec = deaths_spec * 1000 / population_size,
            total_deaths_bpsv = deaths_bpsv * 1000 / population_size,
            time_under_NPIs_bpsv = time_under_NPIs_bpsv,
            composite_NPI_bpsv = composite_NPI_bpsv) %>%
  left_join(bpsv_coverage_time_df, by = "vaccination_rate_bpsv") 

ggplot(subset(vacc_rate_plotting, detection_threshold_hosp == 5)) +
  geom_line(aes(x = days_to_bpsv_coverage, y = central_deaths_averted, col = factor(specific_vaccine_start)), size = 1) +
  facet_grid(R0~NPI_int)


##### Delay to Access Figure #####
owid_vacc_data <- read.csv("data/owid-covid-data.csv") %>%
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>%
  filter(!is.na(total_vaccinations_per_hundred)) %>%
  filter(date > as.Date("01/11/2020", format = "%d/%m/%Y"))

time_to_one <- owid_vacc_data %>%
  filter(continent != "") %>%
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>%
  group_by(continent, location) %>%
  filter(!is.na(total_vaccinations_per_hundred) & total_vaccinations_per_hundred >= 1) %>%
  summarise(time = min(date)) %>%
  ungroup() %>%
  mutate(delay = as.numeric(time - min(time)))

time_to_one$continent <- factor(time_to_one$continent, levels = rev(levels(factor(time_to_one$continent))))
access_delay_empirical_boxplot <- ggplot(time_to_one, aes(x = continent, y = delay)) + 
  geom_boxplot(aes(col = continent), position = position_dodge(0.8), outlier.shape = NA, linewidth = 1.5) +
  geom_jitter(aes(fill = continent), position = position_jitterdodge(1.75), size = 2, pch = 21) +
  theme_bw() +
  labs(x = "", y = "Delay to 1% Population Vaccinated (Days)") +
  coord_flip() +
  lims(y = c(0, 300)) +
  theme(legend.position = "none")

## Generating the scenarios
raw_vacc_delay_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3, 3.5),          # Basic reproduction number
                                             IFR = 1,                                       # IFR
                                             population_size = 10^10,                       # population size
                                             Tg = 7,                                        # Tg
                                             detection_time = 1,                            # PLACEHOLDER FOR NOW
                                             bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                             bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                             specific_vaccine_start = 100 + seq(0, 300, 25),# specific vaccine distribution start (time after detection time)
                                             specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                             efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                             efficacy_disease_bpsv = 0.75,                  # vaccine efficacy against disease - BPSV
                                             efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                             efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                             dur_R = 365000000,                             # duration of infection-induced immunity
                                             dur_bpsv = 365000000,                          # duration of BPSV vaccine immunity
                                             dur_spec = 365000000,                          # duration of disease-specific vaccine immunity
                                             coverage = coverage,                           # proportion of the population vaccinated
                                             vaccination_rate_bpsv = 0.035, 
                                             vaccination_rate_spec = 0.035,                 # disease-specific vaccination rate per week as percentage of population
                                             min_age_group_index_priority = min_age_group_index_priority,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                             min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)

## Join detection time dataframe (note currently have all combos of R0 and detection time and we want specific pairings)
raw_vacc_delay_scenarios2 <- expand_grid(raw_vacc_delay_scenarios,
                                         detection_threshold = unique(bp_df_mean_subset$detection)) %>%
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean) 

## Generating NPIs based on specific detection times, R0, and other vaccine-related events
NPIs_vacc_delay <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, 
                                         minimal_mandate_reduction = minimal_mandate_reduction, 
                                         NPI_scenarios = c(4, 7, 8), 
                                         scenarios = raw_vacc_delay_scenarios2) 

vacc_delay_scenarios <- raw_vacc_delay_scenarios2 %>%
  full_join(NPIs_vacc_delay, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",         # joining by all columns which influence NPI scenario timing
                                   "specific_vaccine_start", "vaccination_rate_bpsv", "vaccination_rate_spec", 
                                   "coverage", "min_age_group_index_priority"), multiple = "all")

## Filtering the above to only select R0 and detection time pairs that actually occurred (the expand grid call above generated all combos)
R0_detection_time_pairs <- bp_df_mean_subset %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time, metric)
vacc_delay_scenarios2 <- vacc_delay_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "immunity_duration") 

## Creating overall output and index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vars_for_index <- c(variable_columns(vacc_delay_scenarios2))
vacc_delay_scenarios3 <- vacc_delay_scenarios2 %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
cores <- parallel::detectCores() - 2
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(vacc_delay_scenarios3, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/vaccination_access_delay_exploration_scenarios.rds")
} else {
  model_outputs <- readRDS("outputs/vaccination_access_delay_exploration_scenarios.rds")
}

## Joining back in the detection metrics
detection_df <- vacc_delay_scenarios3 %>%
  select(scenario_index, specific_vaccine_start, all_of(vars_for_index)) %>%
  filter(vaccine_scenario == "specific_only") %>%
  ungroup() %>%
  select(-vaccine_scenario) %>%
  select(scenario_index, specific_vaccine_start, NPI_int, detection_time, detection_threshold, metric)

model_outputs2 <- model_outputs %>%
  left_join(detection_df, by = c("scenario_index", "specific_vaccine_start", "detection_time", "NPI_int")) %>%
  mutate(detection_threshold_hosp = round(detection_threshold * IHR)) %>%
  mutate(detection_timing = case_when(detection_threshold_hosp == 1 ~ "Early",
                                      detection_threshold_hosp == 5 ~ "Intermediate",
                                      detection_threshold_hosp == 10 ~ "Late",
                                      detection_threshold_hosp == 20 ~ "Very Late")) %>%
  filter(metric == "Daily Incidence") 

time_to_one_continent_df <- time_to_one %>%
  group_by(continent) %>%
  summarise(median_delay = median(delay))

labeller_lookup <- c(`4` = "Minimal NPIs", `7` = "Intermediate NPIs", `8` = "Stringent NPIs")

x <- ggplot(data = subset(model_outputs2, detection_timing == "Intermediate" & NPI_int %in% c(4, 7, 8))) +
  geom_line(aes(x = specific_vaccine_start - min(model_outputs2$specific_vaccine_start), 
                y = bpsv_deaths_averted * 1000 / unique(model_outputs2$population_size), col = factor(R0))) +
  scale_colour_manual(values = c("#F2E9E4", 
                                 "#C9ADA7", 
                                 "#9A8C98",
                                 "#4A4E69",
                                 "#22223B")) +
  new_scale("colour") +
  geom_vline(data = time_to_one_continent_df, aes(xintercept = median_delay, col = continent)) +
  labs(y = "Deaths Averted by Stockpiled BPSV", x = "Delay to Access Disease-Specific Vaccine (Days)") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(NPI_int ~ ., nrow = 3, labeller = labeller(NPI_int = labeller_lookup), strip.position = "right")

delay_access_plot <- cowplot::plot_grid(access_delay_empirical_boxplot, x, nrow = 2, align = "v", axis = "lr", rel_heights = c(1, 2.5),
                                        labels = c("A", "B"))
