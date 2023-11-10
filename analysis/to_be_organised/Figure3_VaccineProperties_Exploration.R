##### THINGS TO CHECK
### 1) ref_infectious_vaccinated was previously set to 50% - have defaulted back to 100% for now. What should we do?
### 2) What to do about variable seeding cases?
### 3) What to do about healthcare capacity? 
### 4) Bug - have to have second dose coverage last a few days longer than primary dose coverage to get
###    up to same overall level of coverage.
### 5) Bug - dynamics of first dose delivery vs booster dose delivery are different (especially so for small population sizes)
###    result is that both vs spec-only scenarios are non-identical (as in former, booster = spec for elderly and in latter first = spec for elderly)
###    ==> this is fixed when you make population sizes really big (10^8 minimum - probably best to work in 10^9 and scale down to per 1000)
### 6) Bug - time_to_coverage_bpsv <- ceiling(elderly_pop_to_vaccinate/daily_doses) + 1 - have to do "+1" to make sure all elderly get
###    vaxxed - does result in a little overshoot into younger age-groups, but doesn't influence deaths basically at all
### 7) Concern: The way dose_pops/vaccination_cov is calculated includes small but non-zero ineligible groups of people like those in hospital.
###    Probably it shouldn't.
### 8) IFR calculation under the hood on my end not quite right, and need to sort (chat to Azra about this)

### Collecting intuition on bits of behaviour
# 1) Only people who can get 2nd dose can be boostered. But in both-vaccines scenario, some of the BPSV'd elderly will die
#    leaving <coverage % of them left for receipt of the specific vaccine.
# 2) Some of the folks 

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

# Generate parameter combinations for model running

### BPSV Efficacy Against Disease

#### Generate initial sets of scenarios (note placeholder for detection time)
raw_bpsv_efficacy_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3, 3.5),                   # Basic reproduction number
                                                IFR = 1,                                       # IFR
                                                population_size = 10^10,
                                                Tg = 5.5,                                      # Tg
                                                detection_time = 1,                            # PLACEHOLDER FOR NOW
                                                bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                                bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                                specific_vaccine_start = c(100, 200, 365),     # specific vaccine distribution start (time after detection time)
                                                specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                                efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                                efficacy_disease_bpsv = seq(0.05, 1, 0.05),    # vaccine efficacy against disease - BPSV
                                                efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                                efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                                dur_R = 365000000,                             # duration of infection-induced immunity
                                                dur_bpsv = 365000000,                          # duration of BPSV vaccine immunity
                                                dur_spec = 365000000,                          # duration of disease-specific vaccine immunity
                                                coverage = 0.8,                                # proportion of the population vaccinated
                                                vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                                min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                                min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)

## Join detection time dataframe (note currently have all combos of R0 and detection time and we want specific pairings)
raw_bpsv_efficacy_scenarios2 <- expand_grid(raw_bpsv_efficacy_scenarios,
                                            detection_threshold = unique(bp_df_mean_subset$detection)) %>%
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean) 

## Generating NPIs based on specific detection times, R0, and other vaccine-related events
NPIs_bpsv_eff <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                       NPI_scenarios = c(4, 7, 8), scenarios = raw_bpsv_efficacy_scenarios2)
bpsv_eff_scenarios <- raw_bpsv_efficacy_scenarios2 %>%
  full_join(NPIs_bpsv_eff, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                  "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")

## Filtering the above to only select R0 and detection time pairs that actually occurred (the expand grid call above generated all combos)
R0_detection_time_pairs <- bp_df_mean_subset %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time, metric)
final_bpsv_eff_scenarios <- bpsv_eff_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "disease_efficacy")
# R0 * spec vax * vaccine scenarios * detection * detection metric * efficacy * NPIs
5 * 3 * 2 * 3 * 2 * 20 * 3

### BPSV Delay Between Vaccination Receipt & Protection
raw_bpsv_delay_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3, 3.5),                   # Basic reproduction number
                                             IFR = 1,                                       # IFR
                                             population_size = 10^10,
                                             Tg = 5.5,                                      # Tg
                                             detection_time = 1,                            # PLACEHOLDER FOR NOW
                                             bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                             bpsv_protection_delay = seq(2, 42, 2),         # delay between receipt of BPSV dose and protection
                                             specific_vaccine_start = c(100, 200, 365),     # specific vaccine distribution start (time after detection time)
                                             specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                             efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                             efficacy_disease_bpsv = 0.75,                  # vaccine efficacy against disease - BPSV
                                             efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                             efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                             dur_R = 365000000,                                # duration of infection-induced immunity
                                             dur_bpsv = 365000000,                          # duration of BPSV vaccine immunity
                                             dur_spec = 365000000,                          # duration of disease-specific vaccine immunity
                                             coverage = 0.8,                                # proportion of the population vaccinated
                                             vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                             min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                             min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)

## Join detection time dataframe (note currently have all combos of R0 and detection time and we want specific pairings)
raw_bpsv_delay_scenarios2 <- expand_grid(raw_bpsv_delay_scenarios,
                                         detection_threshold = unique(bp_df_mean_subset$detection)) %>%
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean) 

## Generating NPIs based on specific detection times, R0, and other vaccine-related events
NPIs_bpsv_delay <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                         NPI_scenarios = c(4, 7, 8), scenarios = raw_bpsv_delay_scenarios2)
bpsv_delay_scenarios <- raw_bpsv_delay_scenarios2 %>%
  full_join(NPIs_bpsv_delay, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                  "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")

## Filtering the above to only select R0 and detection time pairs that actually occurred (the expand grid call above generated all combos)
R0_detection_time_pairs <- bp_df_mean_subset %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time, metric)
final_bpsv_delay_scenarios <- bpsv_delay_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "protection_delay")
# R0 * spec vax * vaccine scenarios * detection * detection metric * efficacy * NPIs
5 * 3 * 2 * 3 * 2 * 21 * 3

### BPSV Duration of Immunity
raw_dur_bpsv_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3, 3.5),                   # Basic reproduction number
                                           IFR = 1,                                       # IFR
                                           population_size = 10^10,
                                           Tg = 5.5,                                      # Tg
                                           detection_time = 1,                            # PLACEHOLDER FOR NOW
                                           bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                           bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                           specific_vaccine_start = c(100, 200, 365),     # specific vaccine distribution start (time after detection time)
                                           specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                           efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                           efficacy_disease_bpsv = 0.75,                  # vaccine efficacy against disease - BPSV
                                           efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                           efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                           dur_R = 365000000,                             # duration of infection-induced immunity
                                           dur_bpsv = seq(20, 365, 15),                   # duration of BPSV vaccine immunity
                                           dur_spec = 365000000,                          # duration of disease-specific vaccine immunity
                                           coverage = 0.8,                                # proportion of the population vaccinated
                                           vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                           min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                           min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)

## Join detection time dataframe (note currently have all combos of R0 and detection time and we want specific pairings)
raw_dur_bpsv_scenarios2 <- expand_grid(raw_dur_bpsv_scenarios,
                                       detection_threshold = unique(bp_df_mean_subset$detection)) %>%
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean) 

## Generating NPIs based on specific detection times, R0, and other vaccine-related events
NPIs_bpsv_dur <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                         NPI_scenarios = c(4, 7, 8), scenarios = raw_dur_bpsv_scenarios2)
bpsv_dur_scenarios <- raw_dur_bpsv_scenarios2 %>%
  full_join(NPIs_bpsv_dur, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                  "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")

## Filtering the above to only select R0 and detection time pairs that actually occurred (the expand grid call above generated all combos)
R0_detection_time_pairs <- bp_df_mean_subset %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time, metric)
final_bpsv_dur_scenarios <- bpsv_dur_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "immunity_duration") 
# R0 * spec vax * vaccine scenarios * detection * detection metric * efficacy * NPIs
5 * 3 * 2 * 3 * 2 * 24 * 3

## Creating overall output and index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vaccine_property_scenarios <- rbind(final_bpsv_eff_scenarios, final_bpsv_delay_scenarios, final_bpsv_dur_scenarios)
vars_for_index <- c(variable_columns(vaccine_property_scenarios), "NPI_int")
vaccine_property_scenarios <- vaccine_property_scenarios %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
cores <- parallel::detectCores() - 2
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(vaccine_property_scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/vaccine_properties_exploration_scenarios.rds")
} else {
  model_outputs <- readRDS("outputs/vaccine_properties_exploration_scenarios.rds")
}

## Joining back in the detection metrics
detection_df <- vaccine_property_scenarios %>%
  select(scenario_index, all_of(vars_for_index)) %>%
  filter(vaccine_scenario == "specific_only") %>%
  ungroup() %>%
  select(-vaccine_scenario) %>%
  select(R0, scenario_index, specific_vaccine_start, NPI_int, detection_time, detection_threshold, metric)

model_outputs2 <- model_outputs %>%
  left_join(detection_df, by = c("R0", "scenario_index", "specific_vaccine_start", "detection_time", "NPI_int")) %>%
  mutate(detection_threshold_hosp = round(detection_threshold * IHR)) %>%
  mutate(detection_timing = case_when(detection_threshold_hosp == 1 ~ "Early",
                                      detection_threshold_hosp == 5 ~ "Intermediate",
                                      detection_threshold_hosp == 10 ~ "Late",
                                      detection_threshold_hosp == 20 ~ "Very Late"))

## Plotting the output
colour_func <- scales::hue_pal()(max(model_outputs$NPI_int))
NPI_colours <- c("#C64191", "#F0803C", "#0D84A9")
population_size <- unique(model_outputs$population_size)
runtime <- unique(model_outputs$runtime)
NPI_to_include <- c(4, 7, 8) # c(2, 4, 5, 7, 8)

NPI_df <- NPIs_bpsv_eff %>%
  filter(R0 == 2.5, detection_time == 40, specific_vaccine_start == 200, NPI_int %in% NPI_to_include) %>%
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
  scale_colour_manual(values = NPI_colours) +
  scale_x_continuous(breaks = c(0, unique(NPI_df$detection_time),
                                unique(NPI_df$detection_time) + unique(NPI_df$bpsv_start) + unique(NPI_df$time_to_coverage_bpsv),
                                unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start) + unique(NPI_df$time_to_coverage_spec)),
                     labels = c("", "", "BPSV\nFinish", "Spec\nFinish")) +
  scale_y_continuous(breaks = c(0, 1, unique(NPI_df$R0)),
                     labels = c("", "1", "R0")) +
  facet_wrap(scenario~., nrow = 3,
             labeller = as_labeller(c(`Scenario 4`='Minimal', `Scenario 7`='Moderate', `Scenario 8`='Stringent'))) +
  coord_cartesian(xlim = c(0, unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start) + unique(NPI_df$time_to_coverage_spec) + 10),
                  ylim = c(0.5, unique(NPI_df$R0) + 0.5)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill="#F5F5F5"))

## Disease Efficacy plot
IFR_fixed <- 1
specific_vaccine_start_fixed <- 200
disease_efficacy_plotting <- model_outputs2 %>%
  filter(IFR == IFR_fixed, NPI_int %in% NPI_to_include, specific_vaccine_start == specific_vaccine_start_fixed, 
         metric == "Daily Incidence") %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "specific_vaccine_start", "efficacy_disease_bpsv")))) %>%
  group_by(R0, specific_vaccine_start, efficacy_disease_bpsv, NPI_int, detection_threshold_hosp) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted) * 1000 / population_size,
            max_deaths_averted = max(bpsv_deaths_averted) * 1000 / population_size,
            central_deaths_averted = bpsv_deaths_averted * 1000 / population_size,
            perc_deaths_averted = 100 * bpsv_deaths_averted / deaths_spec,
            total_deaths_spec = deaths_spec * 1000 / population_size,
            total_deaths_bpsv = deaths_bpsv * 1000 / population_size,
            time_under_NPIs_bpsv = time_under_NPIs_bpsv,
            composite_NPI_bpsv = composite_NPI_bpsv)

ribbon_plotting_eff <- disease_efficacy_plotting %>%
  filter(R0 != 2 & R0 != 3 & detection_threshold_hosp == 5) %>%
  group_by(efficacy_disease_bpsv, R0, detection_threshold_hosp) %>%
  summarise(lower = ifelse(min(central_deaths_averted) - 0.2 < 0, 0, min(central_deaths_averted) - 0.2),
            upper = max(central_deaths_averted) + 0.2)

disease_efficacy_plot <- ggplot(subset(disease_efficacy_plotting, R0 != 2 & R0 != 3 & detection_threshold_hosp == 5)) +
  geom_ribbon(data = ribbon_plotting_eff, aes(x = 100 * efficacy_disease_bpsv, ymin = lower, ymax = upper, group = R0), 
              alpha = 0.1, colour = "black", linetype = "dashed") +
  geom_line(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted, 
                col = interaction(factor(R0), factor(NPI_int))), size = 1) +
  geom_jitter(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted,
                  fill = interaction(factor(R0), factor(NPI_int))),
              size = 3, pch = 21, width = 0, height = 0) +
  theme_bw() +
  scale_colour_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                      n_colours = 3))))  +
  scale_fill_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                    n_colours = 3)),
                               rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                    n_colours = 3)),
                               rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                    n_colours = 3)))) +
  scale_x_continuous(breaks = c(25, 50, 75, 100), labels = paste0(c(25, 50, 75, 100), "%")) +
  annotate("text", x = 102, y = 0.65, label = expression(paste("R"[0], " 1.5")), color = "black", size = 5, hjust = 0) +
  annotate("text", x = 102, y = 5.5, label = expression(paste("R"[0], " 2.5")), color = "black", size = 5, hjust = 0) +
  annotate("text", x = 102, y = 6.75, label = expression(paste("R"[0], " 3.5")), color = "black", size = 5, hjust = 0) +
  coord_cartesian(ylim = c(0, 9),
                  xlim = c(min(disease_efficacy_plotting$efficacy_disease_bpsv) * 100, 110)) +
  labs(x = "BPSV Disease Efficacy", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none")

disease_efficacy_plot2 <- disease_efficacy_plot + 
  annotation_custom(
    ggplotGrob(NPI_plot),
    xmin = 0, xmax = 35, ymin = 3.8, ymax = 9.25)

eff_plot <- ggplot(subset(disease_efficacy_plotting, R0 != 2 & R0 != 3 & detection_threshold_hosp == 5)) +
  geom_ribbon(data = subset(ribbon_plotting_eff, R0 == 2.5 & detection_threshold_hosp == 5),
              aes(x = 100 * efficacy_disease_bpsv,
                  ymin = lower, ymax = upper, group = R0),
              alpha = 0.1, col = "black",linetype = "dashed") +
  geom_ribbon(data = subset(ribbon_plotting_eff, R0 == 1.5 & detection_threshold_hosp == 5),
              aes(x = 100 * efficacy_disease_bpsv,
                  ymin = lower, ymax = upper, group = R0),
              alpha = 0.05) +
              # alpha = 0.03, col = "black",linetype = "dashed") +
  geom_ribbon(data = subset(ribbon_plotting_eff, R0 == 3.5 & detection_threshold_hosp == 5),
              aes(x = 100 * efficacy_disease_bpsv,
                  ymin = lower, ymax = upper, group = R0),
              alpha = 0.05) +
              # alpha = 0.05, col = "black",linetype = "dashed") +
  geom_line(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted, 
                col = interaction(factor(R0), factor(NPI_int))), size = 1) +
  geom_jitter(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted,
                  fill = interaction(factor(R0), factor(NPI_int))),
              size = 3, pch = 21, width = 0, height = 0) +
  theme_bw() +
  scale_alpha_manual(values = c(0.1, 0.5, 0.1)) +
  scale_colour_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                      n_colours = 3))))  +
  scale_fill_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                    n_colours = 3)),
                               rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                    n_colours = 3)),
                               rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                    n_colours = 3)))) +
  scale_x_continuous(breaks = c(25, 50, 75, 100), labels = paste0(c(25, 50, 75, 100), "%")) +
  annotate("text", x = 102, y = 0.65, label = expression(paste("R"[0], " 1.5")), color = "black", size = 5, hjust = 0) +
  annotate("text", x = 102, y = 5.5, label = expression(paste("R"[0], " 2.5")), color = "black", size = 5, hjust = 0) +
  annotate("text", x = 102, y = 6.75, label = expression(paste("R"[0], " 3.5")), color = "black", size = 5, hjust = 0) +
  coord_cartesian(ylim = c(0, 9),
                  xlim = c(min(disease_efficacy_plotting$efficacy_disease_bpsv) * 100, 110)) +
  labs(x = "BPSV Disease Efficacy", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none")
eff_plot2 <- eff_plot + 
  annotation_custom(
    ggplotGrob(NPI_plot),
    xmin = 0, xmax = 35, ymin = 3.8, ymax = 9.25)

## Duration of protection plot
dur_protect_plotting <- model_outputs2 %>%
  filter(IFR == IFR_fixed, 
         NPI_int %in% NPI_to_include, 
         specific_vaccine_start == specific_vaccine_start_fixed,
         metric == "Daily Incidence") %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "specific_vaccine_start", "dur_bpsv")))) %>%
  group_by(R0, specific_vaccine_start, dur_bpsv , NPI_int, detection_threshold_hosp) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted) * 1000 / population_size,
            max_deaths_averted = max(bpsv_deaths_averted) * 1000 / population_size,
            central_deaths_averted = bpsv_deaths_averted * 1000 / population_size,
            perc_deaths_averted = 100 * bpsv_deaths_averted / deaths_spec,
            total_deaths_spec = deaths_spec * 1000 / population_size,
            total_deaths_bpsv = deaths_bpsv * 1000 / population_size,
            time_under_NPIs_bpsv = time_under_NPIs_bpsv,
            composite_NPI_bpsv = composite_NPI_bpsv)

dur_protect_plot <- ggplot(subset(dur_protect_plotting, R0 == 2.5 & detection_threshold_hosp == 5)) +
  geom_line(aes(x = dur_bpsv, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
  geom_point(aes(x = dur_bpsv, y = central_deaths_averted, fill = factor(NPI_int)), 
             size = 2, pch = 21, col = "black") +
  scale_colour_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                      n_colours = 3))[2],
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3))[2],
                                 rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                      n_colours = 3))[2]))  +
  scale_fill_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                    n_colours = 3))[2],
                               rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                    n_colours = 3))[2],
                               rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                    n_colours = 3))[2]))  +
  theme_bw() +
  lims(y = c(0, max(subset(dur_protect_plotting, R0 == 2.5)$central_deaths_averted))) +
  labs(x = "BPSV Immunity Duration (Days)", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none")

ribbon_plotting_dur <- dur_protect_plotting %>%
  filter(R0 != 2 & R0 != 3 & detection_threshold_hosp == 5) %>%
  group_by(dur_bpsv, R0) %>%
  summarise(lower = ifelse(min(central_deaths_averted) - 0.05 < 0, 0, min(central_deaths_averted) - 0.05),
            upper = max(central_deaths_averted) + 0.05)

dur_protect_plot2 <- ggplot(subset(dur_protect_plotting, R0 != 2 & R0 != 3 & detection_threshold_hosp == 5)) +
  geom_ribbon(data = ribbon_plotting_dur, aes(x = dur_bpsv, ymin = lower, ymax = upper, group = R0), 
              alpha = 0.1, colour = "black", linetype = "dashed") +
  geom_line(aes(x = dur_bpsv, y = central_deaths_averted, 
                col = interaction(factor(R0), factor(NPI_int))), size = 1) +
  geom_point(aes(x = dur_bpsv, y = central_deaths_averted, 
                 fill = interaction(factor(R0), factor(NPI_int))), 
             size = 2, pch = 21, col = "black") +
  scale_colour_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                      n_colours = 3))))  +
  scale_fill_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                    n_colours = 3)),
                               rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                    n_colours = 3)),
                               rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                    n_colours = 3))))  +
  theme_bw() +
  labs(x = "BPSV Immunity Duration (Days)", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none")

## Delay to protection plot
delay_protect_plotting <- model_outputs2 %>%
  filter(IFR == IFR_fixed, 
         NPI_int %in% NPI_to_include,
         specific_vaccine_start == specific_vaccine_start_fixed,
         metric == "Daily Incidence") %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "specific_vaccine_start", "bpsv_protection_delay")))) %>%
  group_by(R0, specific_vaccine_start, bpsv_protection_delay, NPI_int, detection_threshold_hosp) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted) * 1000 / population_size,
            max_deaths_averted = max(bpsv_deaths_averted) * 1000 / population_size,
            central_deaths_averted = bpsv_deaths_averted * 1000 / population_size,
            perc_deaths_averted = 100 * bpsv_deaths_averted / deaths_spec,
            total_deaths_spec = deaths_spec * 1000 / population_size,
            total_deaths_bpsv = deaths_bpsv * 1000 / population_size,
            time_under_NPIs_bpsv = time_under_NPIs_bpsv,
            composite_NPI_bpsv = composite_NPI_bpsv)

delay_protect_plot <- ggplot(subset(delay_protect_plotting, R0 %in% c(1.5, 2.5, 3.5) & NPI_int == 7 &
                                      detection_threshold_hosp %in% c(1, 10))) +
  geom_line(aes(x = bpsv_protection_delay, y = central_deaths_averted, col = factor(R0)), size = 1) +
  geom_point(aes(x = bpsv_protection_delay, y = central_deaths_averted, fill = factor(R0)), 
             size = 2, pch = 21, col = "black") +
  scale_colour_manual(values = c(rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3))[1],
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3))[2],
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3))[3]))  +
  scale_fill_manual(values = c(rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3))[1],
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3))[2],
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3))[3]))  +
  facet_grid(detection_threshold_hosp ~ .,
             labeller = as_labeller(c(`1` = "Early Detect", `10` = "Late Detect"))) +
  theme_bw() +
  lims(y = c(0, max(subset(delay_protect_plotting, R0 == 3.5)$central_deaths_averted))) +
  labs(x = "Protection Delay (Days)", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none")

ribbon_plotting_delay <- delay_protect_plotting %>%
  filter(R0 != 1.5 & R0 != 3) %>%
  group_by(bpsv_protection_delay, R0) %>%
  summarise(lower = ifelse(min(central_deaths_averted) - 0.2 < 0, 0, min(central_deaths_averted) - 0.2),
            upper = max(central_deaths_averted) + 0.2)

delay_protect_plot2 <- ggplot(subset(delay_protect_plotting, R0 != 1.5 & R0 != 3)) +
  geom_ribbon(data = ribbon_plotting_delay, aes(x = bpsv_protection_delay, ymin = lower, ymax = upper, group = R0), 
              alpha = 0.1, colour = "black", linetype = "dashed") +
  geom_line(aes(x = bpsv_protection_delay, y = central_deaths_averted, 
                col = interaction(factor(R0), factor(NPI_int))), size = 1) +
  geom_point(aes(x = bpsv_protection_delay, y = central_deaths_averted, 
                 fill = interaction(factor(R0), factor(NPI_int))), size = 2, pch = 21, col = "black") +
  scale_colour_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                      n_colours = 3))))  +
  scale_fill_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                      n_colours = 3))))  +
  theme_bw() +
  lims(y = c(0, max(subset(delay_protect_plotting, R0 != 1.5)$central_deaths_averted) + 0.2)) +
  labs(x = "Protection Delay (Days)", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none")

## Altogether Plotting
temp2 <- cowplot::plot_grid(dur_protect_plot, delay_protect_plot, nrow = 2, ncol = 1, labels = c("B", "C"))
overall2 <- cowplot::plot_grid(disease_efficacy_plot2, temp2, ncol = 2, rel_widths = c(1.5, 1), labels = c("A", ""))
ggsave(filename = "figures/Figure3_VaccineProperties_Exploration.pdf",
       plot = overall2,
       height = 7.75,
       width = 9.92)

#####
# 
# delay_protect_plot <- ggplot(subset(delay_protect_plotting, R0 != 1.5 & NPI_int == 7)) +
#   geom_line(aes(x = bpsv_protection_delay, y = central_deaths_averted, col = factor(R0)), size = 1) +
#   scale_colour_manual(values = rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
#                                                     n_colours = 3, view_palette = TRUE))) +
#   theme_bw() +
#   labs(x = "Protection Delay (Days)", y = "Deaths Averted By BPSV Per 1000") +
#   guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
#   theme(legend.position = "none")

# ribbon_plotting <- disease_efficacy_plotting %>%
#   filter(R0 != 1.5) %>%
#   group_by(efficacy_disease_bpsv, R0) %>%
#   summarise(lower = ifelse(min(central_deaths_averted) - 0.2 < 0, 0, min(central_deaths_averted) - 0.2),
#             upper = max(central_deaths_averted) + 0.2)
# 
# disease_efficacy_plot <- ggplot(subset(disease_efficacy_plotting, R0 != 1.5)) +
#   geom_ribbon(data = ribbon_plotting, aes(x = 100 * efficacy_disease_bpsv, ymin = lower, ymax = upper, group = R0), 
#               alpha = 0.1, colour = "black", linetype = "dashed") +
#   geom_line(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted, 
#                 col = interaction(factor(R0), factor(NPI_int))), size = 1) +
#   geom_jitter(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted,
#                   fill = interaction(factor(R0), factor(NPI_int))),
#               size = 3, pch = 21, width = 0, height = 0) +
#   theme_bw() +
#   scale_colour_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
#                                                       n_colours = 3, view_palette = TRUE)),
#                                  rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
#                                                       n_colours = 3, view_palette = TRUE)),
#                                  rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
#                                                       n_colours = 3, view_palette = TRUE))))  +
#   scale_fill_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
#                                                     n_colours = 3, view_palette = TRUE)),
#                                rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
#                                                     n_colours = 3, view_palette = TRUE)),
#                                rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
#                                                     n_colours = 3, view_palette = TRUE)))) +
#   scale_x_continuous(breaks = c(25, 50, 75, 100), labels = paste0(c(25, 50, 75, 100), "%")) +
#   annotate("text", x = 102, y = 2, label = "R0=2.0", color = "black", size = 5, hjust = 0) +
#   annotate("text", x = 102, y = 5.5, label = "R0=2.5", color = "black", size = 5, hjust = 0) +
#   annotate("text", x = 102, y = 6.75, label = "R0=3.5", color = "black", size = 5, hjust = 0) +
#   coord_cartesian(ylim = c(0, 8),
#                   xlim = c(min(disease_efficacy_plotting$efficacy_disease_bpsv) * 100, 110)) +
#   labs(x = "BPSV Disease Efficacy", y = "Deaths Averted By BPSV Per 1000") +
#   guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
#   theme(legend.position = "none")
# 
# disease_efficacy_plot2 <- disease_efficacy_plot + 
#   annotation_custom(
#     ggplotGrob(NPI_plot),
#     xmin = 0, xmax = 25, ymin = 3.3, ymax = 8.25)
# ggsave(filename = "figures/disease_efficacy_plot.pdf",
#        plot = disease_efficacy_plot2,
#        height = 6.5,
#        width = 6.5)
# 
# ## Vaccine Protection Delay Plot
# 
# 
# 
# 
# ggplot(subset(disease_efficacy_plotting, specific_vaccine_start != 100 & R0 != 1.5)) +
#   # geom_point(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted, fill = R0),
#   #            size = 3, pch = 21) +
#   geom_jitter(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted, fill = R0),
#               size = 3, pch = 21, width = 5) +
#   theme_bw() +
#   scale_x_continuous(breaks = c(0, 25, 50, 75, 100), labels = paste0(c(0, 25, 50, 75, 100), "%")) +
#   scale_fill_viridis_c(name = "Vaccine Efficacy", option = "magma") +
#   labs(x = "BPSV Disease Efficacy", y = "Deaths Averted By BPSV Per 1000") +
#   guides(colour = guide_legend("NPI\nScenario")) +
#   theme(legend.position = "bottom")
# 
# # ggplot(disease_efficacy_plotting) +
# #   geom_point(aes(x = composite_NPI_bpsv, y = central_deaths_averted, col = 100 * efficacy_disease_bpsv), size = 1) +
# #   theme_bw() +
# #   labs(x = "BPSV Disease Efficacy", y = "Deaths Averted By BPSV Per 1000") +
# #   guides(colour = guide_legend("NPI\nScenario")) +
# #   theme(legend.position = "none")
# 
# # ggplot(disease_efficacy_plotting) +
# #   geom_point(aes(x = composite_NPI_bpsv, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
# #   facet_grid(R0 ~ specific_vaccine_start, scales = "free_y") +
# #   scale_colour_manual(values = NPI_colours) +
# #   scale_fill_manual(values = NPI_colours) +
# #   scale_x_continuous(breaks = c(0, 25, 50, 70, 100), labels = paste0(c(0, 25, 50, 75, 100), "%")) +
# #   theme_bw() +
# #   labs(x = "BPSV Disease Efficacy", y = "Deaths Averted By BPSV Per 1000") +
# #   guides(colour = guide_legend("NPI\nScenario")) +
# #   theme(legend.position = "none")
# 
# x <- ggplot(disease_efficacy_plotting) +
#   geom_line(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
#   facet_grid(R0 ~ specific_vaccine_start, scales = "free_y") +
#   scale_colour_manual(values = NPI_colours) +
#   scale_fill_manual(values = NPI_colours) +
#   scale_x_continuous(breaks = c(0, 25, 50, 75, 100), labels = paste0(c(0, 25, 50, 75, 100), "%")) +
#   theme_bw() +
#   labs(x = "BPSV Disease Efficacy", y = "Deaths Averted By BPSV Per 1000") +
#   guides(colour = guide_legend("NPI\nScenario")) +
#   theme(legend.position = "right")
# 
# test <- disease_efficacy_plotting %>%
#   filter(round(efficacy_disease_bpsv, 2) == 0.75)
# 
# y <- ggplot(disease_efficacy_plotting) +
#   geom_line(aes(x = 100 * efficacy_disease_bpsv, y = perc_deaths_averted, col = factor(NPI_int)), size = 1) +
#   facet_grid(R0 ~ specific_vaccine_start, scales = "free_y") +
#   scale_colour_manual(values = NPI_colours) +
#   scale_fill_manual(values = NPI_colours) +
#   scale_x_continuous(breaks = c(0, 25, 50, 75, 100), labels = paste0(c(0, 25, 50, 75, 100), "%")) +
#   theme_bw() +
#   labs(x = "BPSV Disease Efficacy", y = "% of Deaths Averted By BPSV") +
#   guides(colour = guide_legend("NPI\nScenario")) +
#   theme(legend.position = "none")
# 
# z <- ggplot(subset(disease_efficacy_plotting, efficacy_disease_bpsv == 0.5)) +
#   geom_bar(aes(x = factor(NPI_int), y = total_deaths_spec, fill = factor(NPI_int)), stat = "identity", size = 1) +
#   facet_grid(R0 ~ specific_vaccine_start) +
#   scale_colour_manual(values = NPI_colours) +
#   scale_fill_manual(values = NPI_colours) +
#   theme_bw() +
#   labs(x = "BPSV Disease Efficacy", y = "Total Deaths Baseline Scenario") +
#   guides(colour = guide_legend("NPI\nScenario")) +
#   theme(legend.position = "none")
# 
# cowplot::plot_grid(NPI_plot, x, nrow = 1, rel_widths = c(1, 3))
# 
# ## Delay plot
# delay_plotting <- model_outputs %>%
#   filter(map_lgl(varied, ~ setequal(., c("R0", "IFR", "bpsv_protection_delay")))) %>%
#   group_by(bpsv_protection_delay, NPI_int) %>%
#   summarise(min_deaths_averted = min(bpsv_deaths_averted),
#             max_deaths_averted = max(bpsv_deaths_averted),
#             central_deaths_averted = bpsv_deaths_averted[R0 == 2.5 & IFR == 1])
# delay_plot <- ggplot(delay_plotting) +
#   geom_line(aes(x = bpsv_protection_delay, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
#   geom_ribbon(aes(x = bpsv_protection_delay, ymin = min_deaths_averted, ymax = max_deaths_averted,
#                   fill = factor(NPI_int)), alpha = 0.1) +
#   geom_line(aes(x = bpsv_protection_delay, y = min_deaths_averted, col = factor(NPI_int)), size = 0.1) +
#   geom_line(aes(x = bpsv_protection_delay, y = max_deaths_averted, col = factor(NPI_int)), size = 0.1) + 
#   scale_colour_manual(values = NPI_colours) +
#   scale_fill_manual(values = NPI_colours) +
#   theme_bw() +
#   labs(x = "Delay Between Vaccination & Protection", y = "Deaths Averted By BPSV")
# 
# ## Efficacy x delay heatmap
# efficacy_delay_joint_plotting <-  model_outputs %>%
#   filter(map_lgl(varied, ~ setequal(., c("R0", "IFR", "efficacy_disease_bpsv", "dur_vacc_delay")))) %>%
#   filter(NPI_int == 2, R0 == 2.5, IFR == 1) %>%
#   select(efficacy_disease_bpsv, dur_vacc_delay, bpsv_deaths_averted)
# 
# ggplot(data = efficacy_delay_joint_plotting, aes(x = 100 * efficacy_disease_bpsv, y = dur_vacc_delay, fill = bpsv_deaths_averted)) +
#   geom_raster() +
#   scale_fill_gradient2(low = "white", high = "#7A6F9B", mid = "#8B85C1", 
#                        midpoint = 2500, space = "Lab", 
#                        name = "Deaths\nAverted") +
#   scale_x_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"),
#                      expand = c(0, 0)) +
#   scale_y_continuous(breaks = seq(0, 14, 1), expand = c(0, 0)) +
#   theme_bw() + 
#   theme(panel.grid = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   labs(x = "BPSV Disease Efficacy", y = "Delay from Vaccination to Protection")
# 
# 
# 
# ## Plotting the NPI Scenarios
# 
# # 
# # # Plotting model outputted deaths and time under NPIs
# # model_outputs <- model_outputs %>%
# #   mutate(specific_vaccine_start = factor(specific_vaccine_start),
# #          detection_time = factor(detection_time),
# #          NPI_int = factor(NPI_int))
# # 
# # absolute_deaths_plot <- ggplot() +
# #   geom_segment(data = model_outputs, aes(x = composite_NPI_spec, xend = composite_NPI_spec, 
# #                                          y = deaths_spec, yend = deaths_bpsv + 75, group = factor(NPI_int)),
# #                arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
# #   geom_point(data = model_outputs, aes(x = composite_NPI_spec, y = deaths_spec, fill = factor(NPI_int)), 
# #              shape = 4, colour = "black", size = 2, pch = 21) +
# #   geom_point(data = model_outputs, aes(x = composite_NPI_spec, y = deaths_bpsv, fill = factor(NPI_int)), 
# #              colour = "black", size = 4, pch = 21) +
# #   theme_bw() +
# #   labs(x = "NPI Days (Composite Duration+Stringency", y = "Disease Deaths") +
# #   guides(fill = guide_legend(title = "Scenario"))
# # 
# # deaths_averted_plot <- ggplot() +
# #   geom_bar(data = model_outputs, aes(x = factor(NPI_int), y = bpsv_deaths_averted, fill = factor(NPI_int)), stat = "identity") +
# #   labs(x = "", y = "Deaths Averted") +
# #   theme_bw() +
# #   theme(legend.position = "none",
# #         axis.title.x = element_blank(),
# #         plot.background = element_rect(colour = "black"))
# # 
# # inset_prop_y <- 0.4
# # min_deaths <- min(model_outputs$deaths_bpsv)
# # max_deaths <- max(model_outputs$deaths_spec)
# # inset_ymin <- min_deaths + (1 - inset_prop_y) * (max_deaths - min_deaths)
# # inset_ymax <- max_deaths
# # 
# # inset_prop_x <- 0.5
# # min_NPI <- min(model_outputs$composite_NPI_spec)
# # max_NPI <- max(model_outputs$composite_NPI_spec)
# # inset_xmin <- min_NPI + (1 - inset_prop_x) * (max_NPI - min_NPI)
# # inset_xmax <- max_NPI
# # 
# # overall_deaths_plot <- absolute_deaths_plot + 
# #   annotation_custom(
# #     ggplotGrob(deaths_averted_plot), 
# #     xmin = inset_xmin, xmax = inset_xmax, ymin = inset_ymin, ymax = inset_ymax)
# # cowplot::plot_grid(NPI_plot, overall_deaths_plot)
