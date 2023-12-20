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

# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/helper_functions.R"))

# Loading in bp based detection and calculating detection times for the the different R0 values
bp_df_long <- readRDS("outputs/Figure1_bp_detection_times.rds")
prob_hosp <- squire.page.sarsX:::probs_booster$prob_hosp
arg_pop <- squire::get_population("Argentina")
IHR <- sum(prob_hosp * arg_pop$n / sum(arg_pop$n)) 
num_hosp <- 1:20
detection_hosp <- round(num_hosp / IHR, digits = 0)
num_hosp <- c(1, 5, 10)
infection_thresholds <- detection_hosp[num_hosp]
bp_df_mean_subset <- bp_df_long %>%
  filter(!is.infinite(value),
         detection %in% infection_thresholds) %>%
  group_by(R0, detection, metric) %>%
  summarise(mean = mean(value)) 

## Generating default parameters
default <- define_default_params()

# Generate parameter combinations for model running

### BPSV Efficacy Against Disease

#### Generate initial sets of scenarios (note placeholder for detection time)
#### STILL NEED TO WORK ON THIS!!!!! #####
raw_bpsv_efficacy_scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5),                         # Basic reproduction number
                                                efficacy_disease_bpsv = seq(0.05, 1, 0.05))    # vaccine efficacy against disease - BPSV
length(c(1.5, 2.5, 3.5)) * length(seq(0.05, 1, 0.05)) * 2                                            


raw_bpsv_efficacy_scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5),                         # Basic reproduction number
                                                IFR = default$IFR,                                       # IFR
                                                population_size = default$population_size,
                                                hosp_bed_capacity = default$hosp_bed_capacity,
                                                ICU_bed_capacity = default$ICU_bed_capacity,
                                                Tg = default$Tg,                                      # Tg
                                                detection_time = 1,                            # PLACEHOLDER FOR NOW
                                                bpsv_start = default$bpsv_start,                                # BPSV distribution start (time after detection time)
                                                bpsv_protection_delay = default$bpsv_protection_delay,                     # delay between receipt of BPSV dose and protection
                                                specific_vaccine_start = default$specific_vaccine_start,                  # specific vaccine distribution start (time after detection time)
                                                specific_protection_delay = default$specific_protection_delay,                 # delay between receipt of specific dose and protection
                                                efficacy_infection_bpsv = default$efficacy_infection_bpsv,                # vaccine efficacy against infection - BPSV
                                                efficacy_disease_bpsv = seq(0.05, 1, 0.05),    # vaccine efficacy against disease - BPSV
                                                efficacy_infection_spec = default$efficacy_infection_spec,                # vaccine efficacy against infection - specific vaccine
                                                efficacy_disease_spec = default$efficacy_disease_spec,                   # vaccine efficacy against disease - specific vaccine
                                                dur_R = default$dur_R,                             # duration of infection-induced immunity
                                                dur_bpsv = default$dur_bpsv,                          # duration of BPSV vaccine immunity
                                                dur_spec = default$dur_spec,                          # duration of disease-specific vaccine immunity
                                                coverage_bpsv = 0.8,                           # proportion of the population vaccinated
                                                coverage_spec = 0.8,                           # proportion of the population vaccinated
                                                vaccination_rate_bpsv = 0.035,                 # vaccination rate per week as percentage of population
                                                vaccination_rate_spec = 0.035,                 # vaccination rate per week as percentage of population
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
                                  "specific_vaccine_start", "vaccination_rate_bpsv", "vaccination_rate_spec",
                                  "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all")

## Filtering the above to only select R0 and detection time pairs that actually occurred (the expand grid call above generated all combos)
R0_detection_time_pairs <- bp_df_mean_subset %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time, metric)
final_bpsv_eff_scenarios <- bpsv_eff_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "disease_efficacy")
# R0 * vaccine scenarios * detection * detection metric * efficacy * NPIs
3 * 2 * 3 * 2 * 20 * 3

### BPSV Duration of Immunity
raw_dur_bpsv_scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5),                         # Basic reproduction number
                                           IFR = 1,                                       # IFR
                                           population_size = 10^10,
                                           hosp_bed_capacity = 10^10,
                                           ICU_bed_capacity = 10^10,
                                           Tg = 6.7,                                      # Tg
                                           detection_time = 1,                            # PLACEHOLDER FOR NOW
                                           bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                           bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                           specific_vaccine_start = 220,                  # specific vaccine distribution start (time after detection time)
                                           specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                           efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                           efficacy_disease_bpsv = 0.75,                  # vaccine efficacy against disease - BPSV
                                           efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                           efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                           dur_R = 365000000,                             # duration of infection-induced immunity
                                           dur_bpsv = seq(20, 365, 15),                   # duration of BPSV vaccine immunity
                                           dur_spec = 365000000,                          # duration of disease-specific vaccine immunity
                                           coverage_bpsv = 0.8,                           # proportion of the population vaccinated
                                           coverage_spec = 0.8,                           # proportion of the population vaccinated
                                           vaccination_rate_bpsv = 0.035,                 # vaccination rate per week as percentage of population
                                           vaccination_rate_spec = 0.035,                 # vaccination rate per week as percentage of population
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
                                  "specific_vaccine_start", "vaccination_rate_bpsv", "vaccination_rate_spec",
                                  "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all")

## Filtering the above to only select R0 and detection time pairs that actually occurred (the expand grid call above generated all combos)
R0_detection_time_pairs <- bp_df_mean_subset %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time, metric)
final_bpsv_dur_scenarios <- bpsv_dur_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "immunity_duration") 
# R0 * vaccine scenarios * detection * detection metric * efficacy * NPIs
3 * 2 * 3 * 2 * 24 * 3

### Specific Vaccination Development Time
raw_spec_vax_scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5),                         # Basic reproduction number
                                           IFR = 1,                                       # IFR
                                           population_size = 10^10,
                                           hosp_bed_capacity = 10^10,
                                           ICU_bed_capacity = 10^10,
                                           Tg = 6.7,                                      # Tg
                                           detection_time = 1,                            # PLACEHOLDER FOR NOW
                                           bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                           bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                           specific_vaccine_start = seq(65, 365, 15),     # specific vaccine distribution start (time after detection time)
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
                                           vaccination_rate_bpsv = 0.035,                 # vaccination rate per week as percentage of population
                                           vaccination_rate_spec = 0.035,                 # vaccination rate per week as percentage of population
                                           min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                           min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)

## Join detection time dataframe (note currently have all combos of R0 and detection time and we want specific pairings)
raw_spec_vax_scenarios2 <- expand_grid(raw_spec_vax_scenarios,
                                       detection_threshold = unique(bp_df_mean_subset$detection)) %>%
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean) 

## Generating NPIs based on specific detection times, R0, and other vaccine-related events
NPIs_spec_vax <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                       NPI_scenarios = c(4, 7, 8), scenarios = raw_spec_vax_scenarios2)
spec_vax_scenarios <- raw_spec_vax_scenarios2 %>%
  full_join(NPIs_spec_vax, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                  "specific_vaccine_start", "vaccination_rate_bpsv", "vaccination_rate_spec",
                                  "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all")

## Filtering the above to only select R0 and detection time pairs that actually occurred (the expand grid call above generated all combos)
R0_detection_time_pairs <- bp_df_mean_subset %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time, metric)
final_spec_vax_scenarios <- spec_vax_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "spec_development_time") 
# R0 * vaccine scenarios * detection * detection metric * efficacy * NPIs
3 * 2 * 3 * 2 * 19 * 3

## Creating overall output and index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vaccine_property_scenarios <- rbind(final_bpsv_eff_scenarios, final_bpsv_dur_scenarios, final_spec_vax_scenarios)
vars_for_index <- c(variable_columns(vaccine_property_scenarios), "NPI_int")
vaccine_property_scenarios <- vaccine_property_scenarios %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
cores <- parallel::detectCores() - 2
fresh_run <- TRUE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(vaccine_property_scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/Figure3_VaccineProperties/vaccine_properties_exploration_scenarios.rds")
} else {
  model_outputs <- readRDS("outputs/Figure3_VaccineProperties/vaccine_properties_exploration_scenarios.rds")
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
  filter(R0 == 2.5, detection_time == 40, specific_vaccine_start == 220, NPI_int %in% NPI_to_include) %>%
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
specific_vaccine_start_fixed <- 220
disease_efficacy_plotting <- model_outputs2 %>%
  filter(IFR == IFR_fixed, NPI_int %in% NPI_to_include, specific_vaccine_start == specific_vaccine_start_fixed, 
         metric == "Daily Incidence") %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "efficacy_disease_bpsv")))) %>%
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

NPI_1_colours <- rev(generate_palette(NPI_colours[1], modification = "go_lighter", n_colours = 3))[c(1, 3, 1)]
NPI_2_colours <- rev(generate_palette(NPI_colours[2], modification = "go_lighter", n_colours = 3))[c(1, 3, 1)]
NPI_3_colours <- rev(generate_palette(NPI_colours[3], modification = "go_lighter", n_colours = 3))[c(1, 3, 1)]

disease_efficacy_plot <- ggplot(subset(disease_efficacy_plotting, R0 != 2 & R0 != 3 & detection_threshold_hosp == 5)) +
  geom_ribbon(data = ribbon_plotting_eff, aes(x = 100 * efficacy_disease_bpsv, ymin = lower, ymax = upper, group = R0), 
              alpha = 0.1, colour = "black", linetype = "dashed") +
  geom_line(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted, 
                col = interaction(factor(R0), factor(NPI_int))), size = 1) +
  geom_jitter(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted,
                  fill = interaction(factor(R0), factor(NPI_int))),
              size = 3, pch = 21, width = 0, height = 0) +
  theme_bw() +
  scale_colour_manual(values = c(NPI_1_colours, NPI_2_colours, NPI_3_colours))  +
  scale_fill_manual(values = c(NPI_1_colours, NPI_2_colours, NPI_3_colours)) + 
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

## Duration of protection plot
dur_protect_plotting <- model_outputs2 %>%
  filter(IFR == IFR_fixed, 
         NPI_int %in% NPI_to_include, 
         specific_vaccine_start == specific_vaccine_start_fixed,
         metric == "Daily Incidence") %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "dur_bpsv")))) %>%
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
  scale_colour_manual(values = NPI_colours)  +
  scale_fill_manual(values = NPI_colours)  +
  theme_bw() +
  lims(y = c(0, max(subset(dur_protect_plotting, R0 == 2.5)$central_deaths_averted))) +
  labs(x = "BPSV Immunity Duration (Days)", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none")

## Specific Vaccine Development Time Plotting
spec_dev_plotting <- model_outputs2 %>%
  filter(IFR == IFR_fixed, 
         NPI_int %in% NPI_to_include, 
         metric == "Daily Incidence") %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "specific_vaccine_start")))) %>%
  group_by(R0, specific_vaccine_start, dur_bpsv , NPI_int, detection_threshold_hosp) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted) * 1000 / population_size,
            max_deaths_averted = max(bpsv_deaths_averted) * 1000 / population_size,
            central_deaths_averted = bpsv_deaths_averted * 1000 / population_size,
            perc_deaths_averted = 100 * bpsv_deaths_averted / deaths_spec,
            total_deaths_spec = deaths_spec * 1000 / population_size,
            total_deaths_bpsv = deaths_bpsv * 1000 / population_size,
            time_under_NPIs_bpsv = time_under_NPIs_bpsv,
            composite_NPI_bpsv = composite_NPI_bpsv)

spec_dev_plot <- ggplot(subset(spec_dev_plotting, R0 == 2.5 & detection_threshold_hosp == 5)) +
  geom_line(aes(x = specific_vaccine_start, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
  geom_point(aes(x = specific_vaccine_start, y = central_deaths_averted, fill = factor(NPI_int)), 
             size = 2, pch = 21, col = "black") + 
  scale_colour_manual(values = NPI_colours)  +
  scale_fill_manual(values = NPI_colours)  +
  theme_bw() +
  lims(y = c(0, max(subset(spec_dev_plotting, R0 == 2.5)$central_deaths_averted))) +
  labs(x = "Time to Specific Vaccine Development (Days)", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none")


## Altogether Plotting
temp2 <- cowplot::plot_grid(dur_protect_plot, spec_dev_plot, nrow = 2, ncol = 1, labels = c("B", "C"))
overall2 <- cowplot::plot_grid(disease_efficacy_plot2, temp2, ncol = 2, rel_widths = c(1.5, 1), labels = c("A", ""))
ggsave(filename = "figures/Figure_3_VaccineProperties/Figure3_VaccineProperties_Exploration.pdf",
       plot = overall2,
       height = 7.75,
       width = 9.92)
