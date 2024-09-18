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
bp_df_long <- readRDS("outputs/Figure1_branchingProcess_Containment/Figure1_bp_detection_times.rds")
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
  summarise(mean = mean(value)) %>%
  filter(metric == "Daily Incidence")
R0_detection_time_pairs <- bp_df_mean_subset %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time, metric)

# Generating default parameters
default <- define_default_params()

# Generate parameter combinations for model running

## BPSV Efficacy Against Disease
### 1) Generate initial sets of scenarios (note placeholder for detection time)
raw_bpsv_efficacy_scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5), efficacy_disease_bpsv = seq(0.05, 1, 0.05)) 

### 2) Join detection time dataframe (note currently have all combos of R0 and detection time and we want specific pairings)
raw_bpsv_efficacy_scenarios2 <- expand_grid(raw_bpsv_efficacy_scenarios, detection_threshold = unique(bp_df_mean_subset$detection)) %>% 
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean) 

### 3) Generating NPIs based on specific detection times, R0, and other vaccine-related events
NPIs_bpsv_eff <- default_NPI_scenarios(lockdown_Rt = default$lockdown_Rt, minimal_mandate_reduction = default$minimal_mandate_reduction,
                                       NPI_scenarios = c(4, 7, 8), scenarios = raw_bpsv_efficacy_scenarios2)

### 4) Joining by all columns which influence NPI scenario timing
bpsv_eff_scenarios <- raw_bpsv_efficacy_scenarios2 %>%
  full_join(NPIs_bpsv_eff, by = c("R0", "country", "population_size", "detection_time", "bpsv_start", "specific_vaccine_start",
                                  "vaccination_rate_bpsv", "vaccination_rate_spec", "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all")

### 5) Filtering the above to only select R0 and detection time pairs that actually occurred (the expand grid call above generated all combos)
final_bpsv_eff_scenarios <- bpsv_eff_scenarios %>% 
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "disease_efficacy")
# R0 * vaccine scenarios * detection threshold * efficacy * NPIs = 3 * 2 * 3 * 1 * 20 * 3

## BPSV Duration of Immunity
raw_dur_bpsv_scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5), dur_bpsv = seq(20, 365, 15))

raw_dur_bpsv_scenarios2 <- expand_grid(raw_dur_bpsv_scenarios, detection_threshold = unique(bp_df_mean_subset$detection)) %>%
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean) 

NPIs_bpsv_dur <- default_NPI_scenarios(lockdown_Rt = default$lockdown_Rt, minimal_mandate_reduction = default$minimal_mandate_reduction, 
                                       NPI_scenarios = c(4, 7, 8), scenarios = raw_dur_bpsv_scenarios2)

bpsv_dur_scenarios <- raw_dur_bpsv_scenarios2 %>%
  full_join(NPIs_bpsv_dur, by = c("R0", "country", "population_size", "detection_time", "bpsv_start", "specific_vaccine_start", 
                                  "vaccination_rate_bpsv", "vaccination_rate_spec", "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all")

final_bpsv_dur_scenarios <- bpsv_dur_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "immunity_duration") 
# R0 * vaccine scenarios * detection threshold * efficacy * NPIs = 3 * 2 * 3 * 24 * 3

## BPSV Efficacy Against Infection
raw_bpsv_efficacy_infection_scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5), efficacy_infection_bpsv = seq(0.05, 1, 0.05)) 

raw_bpsv_efficacy_infection_scenarios2 <- expand_grid(raw_bpsv_efficacy_infection_scenarios, detection_threshold = unique(bp_df_mean_subset$detection)) %>% 
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean) 

NPIs_bpsv_eff_inf <- default_NPI_scenarios(lockdown_Rt = default$lockdown_Rt, minimal_mandate_reduction = default$minimal_mandate_reduction,
                                           NPI_scenarios = c(4, 7, 8), scenarios = raw_bpsv_efficacy_infection_scenarios2)

bpsv_eff_inf_scenarios <- raw_bpsv_efficacy_infection_scenarios2 %>%
  full_join(NPIs_bpsv_eff_inf, by = c("R0", "country", "population_size", "detection_time", "bpsv_start", "specific_vaccine_start",
                                      "vaccination_rate_bpsv", "vaccination_rate_spec", "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all")

final_bpsv_eff_inf_scenarios <- bpsv_eff_inf_scenarios %>% 
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "infection_efficacy")
# R0 * vaccine scenarios * detection threshold * efficacy * NPIs = 3 * 2 * 3 * 1 * 20 * 3

## Figure 5 - Varying Disease-Specific Vaccine Development Time (and associated access)
raw_vacc_delay_scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5), specific_vaccine_start = 100 + seq(0, 720, 5), runtime = 1000)

raw_vacc_delay_scenarios2 <- expand_grid(raw_vacc_delay_scenarios, detection_threshold = unique(bp_df_mean_subset$detection)) %>%
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean)

NPIs_vacc_delay <- default_NPI_scenarios(lockdown_Rt = default$lockdown_Rt, minimal_mandate_reduction = default$minimal_mandate_reduction,
                                         NPI_scenarios = c(4, 7, 8), scenarios = raw_vacc_delay_scenarios2)

vacc_delay_scenarios <- raw_vacc_delay_scenarios2 %>%
  full_join(NPIs_vacc_delay, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",
                                    "specific_vaccine_start", "vaccination_rate_bpsv", "vaccination_rate_spec",
                                    "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all")

final_vacc_delay_scenarios <- vacc_delay_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "specific_development_time")

## Creating overall output and index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vaccine_property_scenarios <- rbind(final_bpsv_eff_scenarios, final_bpsv_dur_scenarios, 
                                    final_bpsv_eff_inf_scenarios, final_vacc_delay_scenarios)

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
  saveRDS(model_outputs, "outputs/Figure4_BPSVProperties/updated_NEW_BPSV_properties_exploration_scenarios.rds")
} else {
  model_outputs <- readRDS("outputs/Figure4_BPSVProperties/updated_NEW_BPSV_properties_exploration_scenarios.rds")
}

## Joining back in the detection metrics
detection_df <- vaccine_property_scenarios %>%
  select(scenario_index, specific_vaccine_start, all_of(vars_for_index)) %>%
  filter(vaccine_scenario == "specific_only") %>%
  ungroup() %>%
  select(-vaccine_scenario) %>%
  select(R0, scenario_index, specific_vaccine_start, NPI_int, detection_time, detection_threshold)
model_outputs2 <- model_outputs %>%
  left_join(detection_df, by = c("R0", "scenario_index", "specific_vaccine_start", "detection_time", "NPI_int")) %>%
  mutate(detection_threshold_hosp = round(detection_threshold * IHR)) %>%
  mutate(detection_timing = case_when(detection_threshold_hosp == 1 ~ "Early",
                                      detection_threshold_hosp == 5 ~ "Intermediate",
                                      detection_threshold_hosp == 10 ~ "Late",
                                      detection_threshold_hosp == 20 ~ "Very Late"))

## Fig 5A Inset - NPI Scenarios
colour_func <- scales::hue_pal()(max(model_outputs$NPI_int))
NPI_colours <- c("#C64191", "#F0803C", "#0D84A9")
population_size <- unique(model_outputs$population_size)
runtime <- unique(model_outputs$runtime)
NPI_to_include <- c(4, 7, 8) # c(2, 4, 5, 7, 8)

NPI_df <- NPIs_bpsv_eff %>%
  filter(R0 == 2.5, detection_time == 40, specific_vaccine_start == 250, NPI_int %in% NPI_to_include) %>%
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
         next_value = ifelse(is.na(next_value), R0, next_value)) %>%
  mutate(specific_vaccine_start = 150)
overplot_factor <- 1

NPI_plot <- ggplot(NPI_df, aes(x = tt_Rt - overplot_factor, colour = scenario)) +
  geom_hline(aes(yintercept = 1), linewidth = 0.2) +
  geom_hline(aes(yintercept = default$lockdown_Rt), linetype = "dashed", linewidth = 0.2) +
  geom_hline(aes(yintercept = R0 * (1 - default$minimal_mandate_reduction)), linetype = "dashed", linewidth = 0.2) +
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

## Figure 5A Main - BPSV Disease Efficacy
disease_efficacy_plotting <- model_outputs2 %>%
  filter(IFR == 1, NPI_int %in% NPI_to_include, specific_vaccine_start == 250) %>%
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
  filter(R0 != 2.5 & detection_threshold_hosp == 5) %>%
  group_by(efficacy_disease_bpsv, R0, detection_threshold_hosp) %>%
  summarise(lower = ifelse(min(central_deaths_averted) - 0.2 < 0, 0, min(central_deaths_averted) - 0.2),
            upper = max(central_deaths_averted) + 0.2)

NPI_1_colours <- rev(generate_palette(NPI_colours[1], modification = "go_lighter", n_colours = 3))[c(1, 3, 1)]
NPI_2_colours <- rev(generate_palette(NPI_colours[2], modification = "go_lighter", n_colours = 3))[c(1, 3, 1)]
NPI_3_colours <- rev(generate_palette(NPI_colours[3], modification = "go_lighter", n_colours = 3))[c(1, 3, 1)]

disease_efficacy_plot_supp <- ggplot(subset(disease_efficacy_plotting, R0 != 2.5 & detection_threshold_hosp == 5)) +
  # geom_ribbon(data = ribbon_plotting_eff, aes(x = 100 * efficacy_disease_bpsv, ymin = lower, ymax = upper, group = R0), 
  #             alpha = 0.1, colour = "black", linetype = "dashed") +
  geom_line(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted, 
                col = interaction(factor(R0), factor(NPI_int))), size = 1) +
  geom_jitter(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted,
                  fill = interaction(factor(R0), factor(NPI_int))),
              size = 2, pch = 21, width = 0, height = 0) +
  theme_bw() +
  scale_colour_manual(values = c(NPI_1_colours[1:2], NPI_2_colours[1:2], NPI_3_colours[1:2]))  +
  scale_fill_manual(values = c(NPI_1_colours[1:2], NPI_2_colours[1:2], NPI_3_colours[1:2])) + 
  scale_x_continuous(breaks = c(25, 50, 75, 100), labels = paste0(c(25, 50, 75, 100), "%")) +
  # annotate("text", x = 102, y = 0.65, label = expression(paste("R"[0], " 1.5")), color = "black", size = 4, hjust = 0) +
  # annotate("text", x = 102, y = 6.75, label = expression(paste("R"[0], " 3.5")), color = "black", size = 4, hjust = 0) +
  coord_cartesian(ylim = c(0, max(subset(disease_efficacy_plotting, R0 == 3.5)$central_deaths_averted)),
                  xlim = c(min(disease_efficacy_plotting$efficacy_disease_bpsv) * 100, 110)) +
  labs(x = "BPSV Disease Efficacy", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  facet_wrap(R0 ~ ., nrow = 2,
             labeller = as_labeller(c(`1.5`='R0=1.5', `3.5`='R0=3.5'))) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"))
disease_efficacy_plot2 <- disease_efficacy_plot + 
  annotation_custom(
    ggplotGrob(NPI_plot),
    xmin = 0, xmax = 35, ymin = 3.8, ymax = 9.25)


## Figure 5B - Specific Vaccine Development Time
spec_dev_plotting <- model_outputs2 %>%
  filter(specific_vaccine_start <= 365) %>%
  filter(specific_vaccine_start %% 10 != 0) %>%
  group_by(R0, specific_vaccine_start, NPI_int, detection_threshold_hosp) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted) * 1000 / population_size,
            max_deaths_averted = max(bpsv_deaths_averted) * 1000 / population_size,
            central_deaths_averted = bpsv_deaths_averted * 1000 / population_size,
            perc_deaths_averted = 100 * bpsv_deaths_averted / deaths_spec,
            total_deaths_spec = deaths_spec * 1000 / population_size,
            total_deaths_bpsv = deaths_bpsv * 1000 / population_size,
            time_under_NPIs_bpsv = time_under_NPIs_bpsv,
            composite_NPI_bpsv = composite_NPI_bpsv)

NPI_1_colours <- rev(generate_palette(NPI_colours[1], modification = "go_lighter", n_colours = 3))[c(1, 3, 1)]
NPI_2_colours <- rev(generate_palette(NPI_colours[2], modification = "go_lighter", n_colours = 3))[c(1, 3, 1)]
NPI_3_colours <- rev(generate_palette(NPI_colours[3], modification = "go_lighter", n_colours = 3))[c(1, 3, 1)]

ribbon_plotting_eff <- spec_dev_plotting %>%
  group_by(specific_vaccine_start, R0) %>%
  summarise(lower = ifelse(min(central_deaths_averted) - 0.1 < 0, 0, min(central_deaths_averted) - 0.1),
            upper = max(central_deaths_averted) + 0.1)

NPI_colours_new <- rep(NPI_colours, each = 3)

spec_dev_plot <- ggplot(data = subset(spec_dev_plotting, detection_threshold_hosp == 5 & R0 != 2.5), aes(group = interaction(R0, NPI_int))) +
  # geom_ribbon(data = ribbon_plotting_eff, aes(x = specific_vaccine_start, ymin = lower, ymax = upper, group = R0),
  #             alpha = 0.1, colour = "black", linetype = "dashed") +
  geom_line(aes(x = specific_vaccine_start, y = central_deaths_averted,
                col = interaction(factor(R0), factor(NPI_int))), size = 1) +
  geom_point(aes(x = specific_vaccine_start, y = central_deaths_averted, fill = interaction(factor(R0), factor(NPI_int))),
             size = 2, pch = 21, col = "black") +
  scale_colour_manual(values = c(NPI_1_colours[1:2], NPI_2_colours[1:2], NPI_3_colours[1:2]))  +
  scale_fill_manual(values = c(NPI_1_colours[1:2], NPI_2_colours[1:2], NPI_3_colours[1:2])) + 
  theme_bw() +
  # annotate("text", x = 375, y = 1.3, label = expression(paste("R"[0], " 1.5")), color = "black", size = 4, hjust = 0) +
  # annotate("text", x = 375, y = 4.3, label = expression(paste("R"[0], " 2.5")), color = "black", size = 4, hjust = 0) +
  # annotate("text", x = 375, y = 5.5, label = expression(paste("R"[0], " 3.5")), color = "black", size = 4, hjust = 0) +
  lims(y = c(0, max(subset(spec_dev_plotting, R0 == 3.5)$central_deaths_averted) + 0.2),
       x = c(100, 410)) +
  labs(x = "Time to Specific Vaccine Development (Days)", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none") +
  facet_wrap(R0 ~ ., nrow = 2) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"))


## Figure 5C - Duration of Immune Protection
dur_protect_plotting <- model_outputs2 %>%
  filter(IFR == 1, 
         NPI_int %in% NPI_to_include, 
         specific_vaccine_start == 250) %>%
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

dur_protect_plot_supp <- ggplot(subset(dur_protect_plotting, R0 != 2.5 & detection_threshold_hosp == 5)) +
  geom_line(aes(x = dur_bpsv, y = central_deaths_averted, col = interaction(factor(R0), factor(NPI_int))), size = 1) +
  geom_point(aes(x = dur_bpsv, y = central_deaths_averted, fill = interaction(factor(R0), factor(NPI_int))), 
             size = 2, pch = 21, col = "black") +
  scale_colour_manual(values = c(NPI_1_colours[1:2], NPI_2_colours[1:2], NPI_3_colours[1:2]))  +
  scale_fill_manual(values = c(NPI_1_colours[1:2], NPI_2_colours[1:2], NPI_3_colours[1:2])) + 
  theme_bw() +
  lims(y = c(0, max(subset(dur_protect_plotting, R0 == 3.5)$central_deaths_averted))) +
  labs(x = "BPSV Immunity Duration (Days)", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  facet_wrap(R0 ~ ., nrow = 2,
             labeller = as_labeller(c(`1.5`='R0=1.5', `3.5`='R0=3.5'))) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"))

## Figure 5D - BPSV Infection Efficacy plot
bpsv_inf_efficacy_plotting <- model_outputs2 %>%
  filter(IFR == 1, 
         NPI_int %in% NPI_to_include, 
         specific_vaccine_start == 250) %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "efficacy_infection_bpsv")))) %>%
  group_by(R0, specific_vaccine_start, efficacy_infection_bpsv , NPI_int, detection_threshold_hosp) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted) * 1000 / population_size,
            max_deaths_averted = max(bpsv_deaths_averted) * 1000 / population_size,
            central_deaths_averted = bpsv_deaths_averted * 1000 / population_size,
            perc_deaths_averted = 100 * bpsv_deaths_averted / deaths_spec,
            total_deaths_spec = deaths_spec * 1000 / population_size,
            total_deaths_bpsv = deaths_bpsv * 1000 / population_size,
            time_under_NPIs_bpsv = time_under_NPIs_bpsv,
            composite_NPI_bpsv = composite_NPI_bpsv)

bpsv_inf_efficacy_plot <- ggplot(subset(bpsv_inf_efficacy_plotting, R0 == 2.5 & detection_threshold_hosp == 5)) +
  geom_line(aes(x = 100 * efficacy_infection_bpsv, y = central_deaths_averted, col = interaction(factor(R0), factor(NPI_int))), size = 1) +
  geom_point(aes(x = 100 * efficacy_infection_bpsv, y = central_deaths_averted, fill = interaction(factor(R0), factor(NPI_int))), 
             size = 2, pch = 21, col = "black") +
  scale_colour_manual(values = c(NPI_1_colours[1:2], NPI_2_colours[1:2], NPI_3_colours[1:2]))  +
  scale_fill_manual(values = c(NPI_1_colours[1:2], NPI_2_colours[1:2], NPI_3_colours[1:2])) + 
  theme_bw() +
  lims(y = c(0, max(subset(bpsv_inf_efficacy_plotting, R0 == 2.5)$central_deaths_averted))) +
  labs(x = "BPSV Infection Efficacy", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none")

bpsv_inf_efficacy_plot_supp <- ggplot(subset(bpsv_inf_efficacy_plotting, R0 != 2.5 & detection_threshold_hosp == 5)) +
  geom_line(aes(x = 100 * efficacy_infection_bpsv, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
  geom_point(aes(x = 100 * efficacy_infection_bpsv, y = central_deaths_averted, fill = factor(NPI_int)), 
             size = 2, pch = 21, col = "black") +
  scale_colour_manual(values = NPI_colours)  +
  scale_fill_manual(values = NPI_colours)  +
  theme_bw() +
  lims(y = c(0, max(subset(bpsv_inf_efficacy_plotting, R0 == 3.5)$central_deaths_averted))) +
  labs(x = "BPSV Infection Efficacy", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  facet_wrap(R0 ~ ., nrow = 2,
             labeller = as_labeller(c(`1.5`='R0=1.5', `3.5`='R0=3.5'))) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"))

supp_fig3_first_half <- cowplot::plot_grid(disease_efficacy_plot_supp, dur_protect_plot_supp, bpsv_inf_efficacy_plot_supp,
                                           nrow = 1, labels = c("A", "B", "C"))
saveRDS(supp_fig3_first_half, "outputs/Figure4_BPSVProperties/SuppFig3_altR0_vaccine_properties_figure.rds")
ggsave(filename = "figures/Figure_4_VaccineProperties/SuppFig3_altR0_vaccine_properties_firsthalf.pdf",
       plot = supp_fig3_first_half,
       height = 8,
       width = 4.65)

## Altogether Plotting
temp <- cowplot::plot_grid(bpsv_inf_efficacy_plot, dur_protect_plot, nrow = 2, labels = c("C", "D"))
cowplot::plot_grid(disease_efficacy_plot2, spec_dev_plot, temp, nrow = 1, 
                   labels = c("A", "B", NULL), rel_widths = c(2, 2, 1.5))                                                                                          

temp <- cowplot::plot_grid(bpsv_inf_efficacy_plot, dur_protect_plot, nrow = 2, labels = c("B", "C"))
overall2 <- cowplot::plot_grid(disease_efficacy_plot2, temp, spec_dev_plot, nrow = 1, 
                   labels = c("A", "B", "D"), rel_widths = c(2, 1.5, 2))                                                                                          
ggsave(filename = "figures/Figure_4_VaccineProperties/updated_NEW_Figure4_VaccineProperties_Exploration.pdf",
       plot = overall2,
       height = 7,
       width = 12.5)

## Alternative Figure
disease_efficacy_plot <- ggplot(subset(disease_efficacy_plotting, R0 == 2.5 & detection_threshold_hosp == 5)) +
  geom_line(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted, 
                col = factor(NPI_int)), size = 1) +
  geom_jitter(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted,
                  fill = factor(NPI_int)),
              size = 2, pch = 21, width = 0, height = 0) +
  theme_bw() + 
  scale_colour_manual(values = NPI_colours)  +
  scale_fill_manual(values = NPI_colours) + 
  scale_x_continuous(breaks = c(25, 50, 75, 100), labels = paste0(c(25, 50, 75, 100), "%")) +
  annotate("text", x = 102, y = 5.5, label = expression(paste("R"[0], " 2.5")), color = "black", size = 4, hjust = 0) +
  coord_cartesian(ylim = c(0, 10.5),
                  xlim = c(min(disease_efficacy_plotting$efficacy_disease_bpsv) * 100, 110)) +
  labs(x = "BPSV Disease Efficacy", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none")
disease_efficacy_plot2 <- disease_efficacy_plot + 
  annotation_custom(
    ggplotGrob(NPI_plot),
    xmin = 0, xmax = 50, ymin = 3.8, ymax = 10.75)

fig_4_first_half <- cowplot::plot_grid(disease_efficacy_plot2, temp, 
                                       labels = c("A", NA), rel_widths = c(1.2, 1))
saveRDS(fig_4_first_half,
        file = "figures/Figure_4_VaccineProperties/test_new_Fig4_first_half.rds")

# ## BPSV Coverage plot
# bpsv_coverage_plotting <- model_outputs2 %>%
#   filter(IFR == IFR_fixed, 
#          NPI_int %in% NPI_to_include, 
#          specific_vaccine_start == specific_vaccine_start_fixed) %>%
#   filter(map_lgl(varied, ~ setequal(., c("R0", "coverage_bpsv")))) %>%
#   group_by(R0, specific_vaccine_start, coverage_bpsv , NPI_int, detection_threshold_hosp) %>%
#   summarise(min_deaths_averted = min(bpsv_deaths_averted) * 1000 / population_size,
#             max_deaths_averted = max(bpsv_deaths_averted) * 1000 / population_size,
#             central_deaths_averted = bpsv_deaths_averted * 1000 / population_size,
#             perc_deaths_averted = 100 * bpsv_deaths_averted / deaths_spec,
#             total_deaths_spec = deaths_spec * 1000 / population_size,
#             total_deaths_bpsv = deaths_bpsv * 1000 / population_size,
#             time_under_NPIs_bpsv = time_under_NPIs_bpsv,
#             composite_NPI_bpsv = composite_NPI_bpsv)
# 
# bpsv_coverage_plot <- ggplot(subset(bpsv_coverage_plotting, R0 == 2.5 & detection_threshold_hosp == 5)) +
#   geom_line(aes(x = 100 * coverage_bpsv, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
#   geom_point(aes(x = 100 * coverage_bpsv, y = central_deaths_averted, fill = factor(NPI_int)), 
#              size = 2, pch = 21, col = "black") +
#   scale_colour_manual(values = NPI_colours)  +
#   scale_fill_manual(values = NPI_colours)  +
#   theme_bw() +
#   lims(y = c(0, max(subset(bpsv_coverage_plotting, R0 == 2.5)$central_deaths_averted))) +
#   labs(x = "% Coverage 65+ Population With BPSV", y = "Deaths Averted By BPSV Per 1000") +
#   guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
#   theme(legend.position = "none")
# 
# ## BPSV Vaccination Rate plot
# bpsv_rate_plotting <- model_outputs2 %>%
#   filter(IFR == IFR_fixed, 
#          NPI_int %in% NPI_to_include, 
#          specific_vaccine_start == specific_vaccine_start_fixed) %>%
#   filter(map_lgl(varied, ~ setequal(., c("R0", "vaccination_rate_bpsv")))) %>%
#   group_by(R0, specific_vaccine_start, vaccination_rate_bpsv , NPI_int, detection_threshold_hosp) %>%
#   summarise(min_deaths_averted = min(bpsv_deaths_averted) * 1000 / population_size,
#             max_deaths_averted = max(bpsv_deaths_averted) * 1000 / population_size,
#             central_deaths_averted = bpsv_deaths_averted * 1000 / population_size,
#             perc_deaths_averted = 100 * bpsv_deaths_averted / deaths_spec,
#             total_deaths_spec = deaths_spec * 1000 / population_size,
#             total_deaths_bpsv = deaths_bpsv * 1000 / population_size,
#             time_under_NPIs_bpsv = time_under_NPIs_bpsv,
#             composite_NPI_bpsv = composite_NPI_bpsv) %>%
#   left_join(bpsv_coverage_time_df, by = "vaccination_rate_bpsv")
# 
# bpsv_rate_plot <- ggplot(subset(bpsv_rate_plotting, R0 == 2.5 & detection_threshold_hosp == 5)) +
#   geom_line(aes(x = days_to_bpsv_coverage, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
#   geom_point(aes(x = days_to_bpsv_coverage, y = central_deaths_averted, fill = factor(NPI_int)), 
#              size = 2, pch = 21, col = "black") +
#   scale_colour_manual(values = NPI_colours)  +
#   scale_fill_manual(values = NPI_colours)  +
#   theme_bw() +
#   lims(y = c(0, max(subset(bpsv_rate_plotting, R0 == 2.5)$central_deaths_averted))) +
#   labs(x = "Days to Complete BPSV Campaign", y = "Deaths Averted By BPSV Per 1000") +
#   guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
#   theme(legend.position = "none")