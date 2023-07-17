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

# NPI Relevant Parameters
lockdown_Rt <- 0.9                   # Rt achieved under lockdown
minimal_mandate_reduction <- 0.25    # Fold-reduction in R0 achieved under minimal mandate restrictions

# Generate parameter combinations for model running

### BPSV Efficacy Against Disease
raw_bpsv_efficacy_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3.5),                      # Basic reproduction number
                                                IFR = c(0.5, 1, 1.5),                          # IFR
                                                population_size = 10^10,
                                                Tg = 5.5,                                      # Tg
                                                detection_time = 14,                           # detection time
                                                bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                                bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                                specific_vaccine_start = c(100, 200, 365, 500),# specific vaccine distribution start (time after detection time)
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
raw_bpsv_efficacy_scenarios$main_varied <- "disease_efficacy"
NPIs_bpsv_eff <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                       NPI_scenarios = 1:9, scenarios = raw_bpsv_efficacy_scenarios)
bpsv_eff_scenarios <- raw_bpsv_efficacy_scenarios %>%
  full_join(NPIs_bpsv_eff, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                  "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")

### BPSV Delay Between Vaccination Receipt & Protection
raw_bpsv_delay_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3.5),                      # Basic reproduction number
                                             IFR = c(0.5, 1, 1.5),                          # IFR
                                             population_size = 10^10,
                                             Tg = 5.5,                                      # Tg
                                             detection_time = 14,                           # detection time
                                             bpsv_start = 14,                               # BPSV distribution start (time after detection time)
                                             bpsv_protection_delay = seq(2, 21, 1),         # delay between receipt of BPSV dose and protection
                                             specific_vaccine_start = c(100, 200, 365, 500),# specific vaccine distribution start (time after detection time)
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
raw_bpsv_delay_scenarios$main_varied <- "protection_delay"
NPIs_bpsv_delay <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                         NPI_scenarios = 1:9, scenarios = raw_bpsv_delay_scenarios)
bpsv_delay_scenarios <- raw_bpsv_delay_scenarios %>%
  full_join(NPIs_bpsv_delay, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                    "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")

### BPSV Duration of Immunity
raw_dur_bpsv_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3.5),                      # Basic reproduction number
                                           IFR = c(0.5, 1, 1.5),                          # IFR
                                           population_size = 10^10,
                                           Tg = 5.5,                                      # Tg
                                           detection_time = 14,                           # detection time
                                           bpsv_start = 14,                               # BPSV distribution start (time after detection time)
                                           bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                           specific_vaccine_start = c(100, 200, 365, 500),# specific vaccine distribution start (time after detection time)
                                           specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                           efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                           efficacy_disease_bpsv = 0.75,    # vaccine efficacy against disease - BPSV
                                           efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                           efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                           dur_R = 365000000,                             # duration of infection-induced immunity
                                           dur_bpsv = seq(20, 365, 15),                   # duration of BPSV vaccine immunity
                                           dur_spec = 365000000,                          # duration of disease-specific vaccine immunity
                                           coverage = 0.8,                                # proportion of the population vaccinated
                                           vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                           min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                           min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)
raw_dur_bpsv_scenarios$main_varied <- "immunity_duration"
NPIs_bpsv_dur <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                        NPI_scenarios = 1:9, scenarios = raw_dur_bpsv_scenarios)
bpsv_dur_scenarios <- raw_dur_bpsv_scenarios %>%
  full_join(NPIs_bpsv_dur, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                  "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")

vaccine_property_scenarios <- rbind(bpsv_eff_scenarios, bpsv_delay_scenarios, bpsv_dur_scenarios)

## Creating index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
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

## Plotting the output
colour_func <- scales::hue_pal()(9)
NPI_colours <- colour_func[c(2, 4, 6)]
population_size <- unique(model_outputs$population_size)
runtime <- unique(model_outputs$runtime)

NPI_df <- NPIs_bpsv_eff %>%
  filter(R0 == 2.5) %>%
  filter(specific_vaccine_start == 200) %>%
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
  geom_vline(aes(xintercept = detection_time), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time + bpsv_start + time_to_coverage_bpsv), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time + specific_vaccine_start + time_to_coverage_spec), linewidth = 0.2) +
  geom_segment(aes(xend = next_time + overplot_factor, y = Rt, yend = Rt), size = 1) +
  geom_segment(aes(x = next_time, xend = next_time, y = Rt, yend = next_value), size = 1) +
  theme_bw() +
  scale_colour_manual(values = NPI_colours) +
  facet_wrap(scenario~., ncol = 1) +
  labs(x = "Time Since Spillover (Days)") +
  coord_cartesian(xlim = c(0, unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start) + unique(NPI_df$time_to_coverage_spec) + 10),
                  ylim = c(0.5, unique(NPI_df$R0) + 0.5)) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="#F5F5F5"))

## Efficacy plot
disease_efficacy_plotting <- model_outputs %>%
  filter(IFR == 1) %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "IFR", "specific_vaccine_start", "efficacy_disease_bpsv")))) %>%
  group_by(R0, specific_vaccine_start, efficacy_disease_bpsv, NPI_int) %>%
  filter(NPI_int != 9) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted) * 1000 / population_size,
            max_deaths_averted = max(bpsv_deaths_averted) * 1000 / population_size,
            central_deaths_averted = bpsv_deaths_averted * 1000 / population_size,
            perc_deaths_averted = 100 * bpsv_deaths_averted / deaths_spec,
            total_deaths_spec = deaths_spec * 1000 / population_size,
            total_deaths_bpsv = deaths_bpsv * 1000 / population_size)

x <- ggplot(disease_efficacy_plotting) +
  geom_line(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
  facet_grid(R0 ~ specific_vaccine_start, scales = "free_y") +
  scale_colour_manual(values = NPI_colours) +
  scale_fill_manual(values = NPI_colours) +
  scale_x_continuous(breaks = c(0, 25, 50, 70, 100), labels = paste0(c(0, 25, 50, 75, 100), "%")) +
  theme_bw() +
  labs(x = "BPSV Disease Efficacy", y = "Deaths Averted By BPSV Per 1000") +
  guides(colour = guide_legend("NPI\nScenario")) +
  theme(legend.position = "none")

y <- ggplot(disease_efficacy_plotting) +
  geom_line(aes(x = 100 * efficacy_disease_bpsv, y = perc_deaths_averted, col = factor(NPI_int)), size = 1) +
  facet_grid(R0 ~ specific_vaccine_start, scales = "free_y") +
  scale_colour_manual(values = NPI_colours) +
  scale_fill_manual(values = NPI_colours) +
  scale_x_continuous(breaks = c(0, 25, 50, 70, 100), labels = paste0(c(0, 25, 50, 75, 100), "%")) +
  theme_bw() +
  labs(x = "BPSV Disease Efficacy", y = "% of Deaths Averted By BPSV") +
  guides(colour = guide_legend("NPI\nScenario")) +
  theme(legend.position = "none")

z <- ggplot(subset(disease_efficacy_plotting, efficacy_disease_bpsv == 0.5)) +
  geom_bar(aes(x = factor(NPI_int), y = total_deaths_spec, fill = factor(NPI_int)), stat = "identity", size = 1) +
  facet_grid(R0 ~ specific_vaccine_start) +
  scale_colour_manual(values = NPI_colours) +
  scale_fill_manual(values = NPI_colours) +
  theme_bw() +
  labs(x = "BPSV Disease Efficacy", y = "Total Deaths Baseline Scenario") +
  guides(colour = guide_legend("NPI\nScenario")) +
  theme(legend.position = "none")

cowplot::plot_grid(NPI_plot, x, nrow = 1, rel_widths = c(1, 3))

## Delay plot
delay_plotting <- model_outputs %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "IFR", "bpsv_protection_delay")))) %>%
  group_by(bpsv_protection_delay, NPI_int) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted),
            max_deaths_averted = max(bpsv_deaths_averted),
            central_deaths_averted = bpsv_deaths_averted[R0 == 2.5 & IFR == 1])
delay_plot <- ggplot(delay_plotting) +
  geom_line(aes(x = bpsv_protection_delay, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
  geom_ribbon(aes(x = bpsv_protection_delay, ymin = min_deaths_averted, ymax = max_deaths_averted,
                  fill = factor(NPI_int)), alpha = 0.1) +
  geom_line(aes(x = bpsv_protection_delay, y = min_deaths_averted, col = factor(NPI_int)), size = 0.1) +
  geom_line(aes(x = bpsv_protection_delay, y = max_deaths_averted, col = factor(NPI_int)), size = 0.1) + 
  scale_colour_manual(values = NPI_colours) +
  scale_fill_manual(values = NPI_colours) +
  theme_bw() +
  labs(x = "Delay Between Vaccination & Protection", y = "Deaths Averted By BPSV")

## Efficacy x delay heatmap
efficacy_delay_joint_plotting <-  model_outputs %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "IFR", "efficacy_disease_bpsv", "dur_vacc_delay")))) %>%
  filter(NPI_int == 2, R0 == 2.5, IFR == 1) %>%
  select(efficacy_disease_bpsv, dur_vacc_delay, bpsv_deaths_averted)

ggplot(data = efficacy_delay_joint_plotting, aes(x = 100 * efficacy_disease_bpsv, y = dur_vacc_delay, fill = bpsv_deaths_averted)) +
  geom_raster() +
  scale_fill_gradient2(low = "white", high = "#7A6F9B", mid = "#8B85C1", 
                       midpoint = 2500, space = "Lab", 
                       name = "Deaths\nAverted") +
  scale_x_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 14, 1), expand = c(0, 0)) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(x = "BPSV Disease Efficacy", y = "Delay from Vaccination to Protection")



## Plotting the NPI Scenarios

# 
# # Plotting model outputted deaths and time under NPIs
# model_outputs <- model_outputs %>%
#   mutate(specific_vaccine_start = factor(specific_vaccine_start),
#          detection_time = factor(detection_time),
#          NPI_int = factor(NPI_int))
# 
# absolute_deaths_plot <- ggplot() +
#   geom_segment(data = model_outputs, aes(x = composite_NPI_spec, xend = composite_NPI_spec, 
#                                          y = deaths_spec, yend = deaths_bpsv + 75, group = factor(NPI_int)),
#                arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
#   geom_point(data = model_outputs, aes(x = composite_NPI_spec, y = deaths_spec, fill = factor(NPI_int)), 
#              shape = 4, colour = "black", size = 2, pch = 21) +
#   geom_point(data = model_outputs, aes(x = composite_NPI_spec, y = deaths_bpsv, fill = factor(NPI_int)), 
#              colour = "black", size = 4, pch = 21) +
#   theme_bw() +
#   labs(x = "NPI Days (Composite Duration+Stringency", y = "Disease Deaths") +
#   guides(fill = guide_legend(title = "Scenario"))
# 
# deaths_averted_plot <- ggplot() +
#   geom_bar(data = model_outputs, aes(x = factor(NPI_int), y = bpsv_deaths_averted, fill = factor(NPI_int)), stat = "identity") +
#   labs(x = "", y = "Deaths Averted") +
#   theme_bw() +
#   theme(legend.position = "none",
#         axis.title.x = element_blank(),
#         plot.background = element_rect(colour = "black"))
# 
# inset_prop_y <- 0.4
# min_deaths <- min(model_outputs$deaths_bpsv)
# max_deaths <- max(model_outputs$deaths_spec)
# inset_ymin <- min_deaths + (1 - inset_prop_y) * (max_deaths - min_deaths)
# inset_ymax <- max_deaths
# 
# inset_prop_x <- 0.5
# min_NPI <- min(model_outputs$composite_NPI_spec)
# max_NPI <- max(model_outputs$composite_NPI_spec)
# inset_xmin <- min_NPI + (1 - inset_prop_x) * (max_NPI - min_NPI)
# inset_xmax <- max_NPI
# 
# overall_deaths_plot <- absolute_deaths_plot + 
#   annotation_custom(
#     ggplotGrob(deaths_averted_plot), 
#     xmin = inset_xmin, xmax = inset_xmax, ymin = inset_ymin, ymax = inset_ymax)
# cowplot::plot_grid(NPI_plot, overall_deaths_plot)
