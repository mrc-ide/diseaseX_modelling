# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))

# Generate parameter combinations for model running
raw_primary_country_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3, 3.5),                   # Basic reproduction number
                                                  IFR = 1,                                       # IFR
                                                  population_size = 10^10,                       # population size
                                                  Tg = 6.7,                                      # Tg
                                                  detection_time = seq(10, 150, 10),             # detection time - PLACECHOLDER FOR NOW
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
                                                  min_age_group_index_non_priority = 4,           # index of the youngest age group that *receives* vaccines (4 = 15+)
                                                  seeding_cases = 1,
                                                  runtime = 1460)         

# NPI Relevant Parameters
lockdown_Rt <- 0.9                   # Rt achieved under lockdown
minimal_mandate_reduction <- 0.25    # Fold-reduction in R0 achieved under minimal mandate restrictions
NPIs_primary_country <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                              NPI_scenarios = 1:9, scenarios = raw_primary_country_scenarios)

# Combining NPI scenarios with parameter combos into one overall dataframe for model running 
all_combos_primary_country_scenarios <- raw_primary_country_scenarios %>%
  full_join(NPIs_primary_country, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                         "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")


## Creating index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vars_for_index <- c(variable_columns(all_combos_primary_country_scenarios), "NPI_int")
all_combos_primary_country_scenarios <- all_combos_primary_country_scenarios %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
cores <- parallel::detectCores() - 2
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(all_combos_primary_country_scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/test_Figure1_primaryCountry_outputs.rds")
} else {
  model_outputs <- readRDS("outputs/test_Figure1_primaryCountry_outputs.rds")
}


## Selecting which NPIs to include
colour_func <- scales::hue_pal()(max(model_outputs$NPI_int))
population_size <- unique(model_outputs$population_size)
runtime <- unique(model_outputs$runtime)
NPI_to_include <- 1:9

## NPI Plot
initial_NPI_df <- NPIs_primary_country %>%
  filter(R0 == 2.5, specific_vaccine_start == 200, NPI_int %in% NPI_to_include) 

table(initial_NPI_df$detection_time)
NPI_df <- initial_NPI_df %>%
  filter(detection_time == 30) %>%
  select(R0, detection_time, bpsv_start, specific_vaccine_start, time_to_coverage_bpsv, time_to_coverage_spec, NPI_int, Rt, tt_Rt) %>%
  rowwise() %>%
  mutate(scenario_info = list(tibble(Rt = Rt, tt_Rt = tt_Rt))) %>%
  select(-Rt, -tt_Rt) %>%
  unnest(cols = c(scenario_info)) %>%
  mutate(scenario = paste0("Scenario ", NPI_int)) %>%
  group_by(scenario) %>%
  mutate(next_time = lead(tt_Rt),
         next_value = lead(Rt)) %>%
  mutate(next_time = ifelse(is.na(next_time), runtime, next_time),
         next_value = ifelse(is.na(next_value), R0, next_value))
overplot_factor <- 1

NPI_plot <- ggplot(NPI_df, aes(x = tt_Rt - overplot_factor, colour = scenario)) +
  geom_hline(aes(yintercept = 1), linewidth = 0.2) +
  geom_hline(aes(yintercept = lockdown_Rt), linetype = "dashed", linewidth = 0.2) +
  geom_hline(aes(yintercept = R0 * (1 - minimal_mandate_reduction)), linetype = "dashed", linewidth = 0.2) +
  geom_vline(aes(xintercept = 0), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time + bpsv_start + time_to_coverage_bpsv), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time + specific_vaccine_start + time_to_coverage_spec), linewidth = 0.2) +
  geom_segment(aes(xend = next_time + overplot_factor, y = Rt, yend = Rt), size = 1) +
  geom_segment(aes(x = next_time, xend = next_time, y = Rt, yend = next_value), size = 1) +
  theme_bw() +
  scale_colour_manual(values = colour_func) +
  scale_x_continuous(breaks = c(0, unique(NPI_df$detection_time),
                                unique(NPI_df$detection_time) + unique(NPI_df$bpsv_start) + unique(NPI_df$time_to_coverage_bpsv),
                                unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start) + unique(NPI_df$time_to_coverage_spec)),
                     labels = c("", "", "BPSV\nFinish", "Spec\nFinish")) +
  scale_y_continuous(breaks = c(0, 1, unique(NPI_df$R0)),
                     labels = c("", "1", "R0")) +
  facet_wrap(scenario~., nrow = 9) +
  coord_cartesian(xlim = c(0, unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start) + unique(NPI_df$time_to_coverage_spec) + 10),
                  ylim = c(0.5, unique(NPI_df$R0) + 0.5)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill="white"))

centralSpec_deathsAverted <- model_outputs %>%
  filter(specific_vaccine_start == 100,
         R0 %in% c(1.5, 2, 2.5, 3, 3.5))

x <- ggplot(data = centralSpec_deathsAverted) +
  geom_bar(aes(x = factor(detection_time), y = deaths_spec, fill = factor(NPI_int)), 
           stat = "identity", position = "dodge") +
  facet_grid(NPI_int~R0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

y <- ggplot(data = centralSpec_deathsAverted) +
  geom_bar(aes(x = factor(detection_time), y = deaths_bpsv, fill = factor(NPI_int)), 
           stat = "identity", position = "dodge") +
  facet_grid(NPI_int~R0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

z <- ggplot(data = centralSpec_deathsAverted) +
  geom_bar(aes(x = factor(detection_time), y = bpsv_deaths_averted, fill = factor(NPI_int)), 
           stat = "identity", position = "dodge") +
  facet_grid(NPI_int~R0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# check why there's an increase in deaths going from 10 -> 150 detection time for scenario 9
## because specific vaccine is being developed at day 200 after detection
## early detection for R0 1.5 means decent chance that you implement vaccine before
## epidemic properly takes off :)
a <- cowplot::plot_grid(NPI_plot, x, nrow = 1, rel_widths = c(1, 5))

b <- cowplot::plot_grid(NPI_plot, y, nrow = 1, rel_widths = c(1, 5))

c <- cowplot::plot_grid(NPI_plot, z, nrow = 1, rel_widths = c(1, 5))

## R0 vs Detection Time Plot (discrete detection time)
R0_detectionFactor_plot <- ggplot(data = subset(centralSpec_deathsAverted, metric == "Daily Incidence")) +
  geom_bar(aes(x = factor(R0), y = bpsv_deaths_averted * 1000 / population_size, 
               colour = factor(NPI_int), fill = factor(NPI_int)), 
           stat = "identity", position = "dodge") +
  facet_wrap(.~detection_timing, strip.position = "top",
             labeller = as_labeller(c(Early = "Early\nDetection", 
                                      Intermediate = "Intermediate\nDetection",
                                      Late = "Late\nDetection",
                                      `Very Late` = "Very Late\nDetection"))) +
  labs(x = "Reproduction Number (R0)", y = "Deaths Averted by BPSV (Per 1000)") + 
  scale_fill_manual(values = NPI_colours) +
  scale_colour_manual(values = NPI_colours) +
  theme_classic() +
  theme(strip.placement = "outside",
        legend.position = "none")

complete_R0_detectionFactor_plot <- cowplot::plot_grid(NPI_plot, R0_detectionFactor_plot, 
                                                       nrow = 1, rel_widths = c(1, 4), labels = c("C", "D"))