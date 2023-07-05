# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))

# Generate parameter combinations for model running (note Rt and tt_Rt has a placeholder)
raw_bpsv_scenarios <- create_scenarios(R0 = c(2.25, 2.5, 2.75),                       # Basic reproduction number
                                       IFR = c(0.75, 1, 1.25),                        # IFR
                                       Tg = 7,                                        # Tg
                                       detection_time = 14,                           # detection time
                                       bpsv_start = 14,                               # BPSV distribution start (time after detection time)
                                       specific_vaccine_start = 100,                  # specific vaccine distribution start (time after detection time)
                                       efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                       efficacy_disease_bpsv = seq(0.01, 1, 0.02),    # vaccine efficacy against disease - BPSV
                                       efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                       efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                       dur_R = 365000,                                # duration of infection-induced immunity
                                       dur_V = 365000,                                # duration of vaccine-induced immunity for both vaccines
                                       second_dose_delay = 7,                         # controls how many days after "1st dose" people receive second dose; see here: https://github.com/mrc-ide/squire.page/blob/main/inst/odin/nimue_booster.R#L427-L430
                                       dur_vacc_delay = 7,                            # mean duration from vaccination to protection
                                       coverage = 0.8,                                # proportion of the population vaccinated
                                       vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                       min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                       min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)

raw_delay_scenarios <- create_scenarios(R0 = c(2.25, 2.5, 2.75),                       # Basic reproduction number
                                        IFR = c(0.75, 1, 1.25),                        # IFR
                                        Tg = 7,                                        # Tg
                                        detection_time = 14,                           # detection time
                                        bpsv_start = 14,                               # BPSV distribution start (time after detection time)
                                        specific_vaccine_start = 100,                  # specific vaccine distribution start (time after detection time)
                                        efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                        efficacy_disease_bpsv = 0.75,                  # vaccine efficacy against disease - BPSV
                                        efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                        efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                        dur_R = 365000,                                # duration of infection-induced immunity
                                        dur_V = 365000,                                # duration of vaccine-induced immunity for both vaccines
                                        second_dose_delay = 7,                         # controls how many days after "1st dose" people receive second dose; see here: https://github.com/mrc-ide/squire.page/blob/main/inst/odin/nimue_booster.R#L427-L430
                                        dur_vacc_delay = 1:14,                         # mean duration from vaccination to protection
                                        coverage = 0.8,                                # proportion of the population vaccinated
                                        vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                        min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                        min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)

raw_delay_bpsv_scenarios <- create_scenarios(R0 = c(2.25, 2.5, 2.75),                       # Basic reproduction number
                                             IFR = c(0.75, 1, 1.25),                        # IFR
                                             Tg = 7,                                        # Tg
                                             detection_time = 14,                           # detection time
                                             bpsv_start = 14,                               # BPSV distribution start (time after detection time)
                                             specific_vaccine_start = 100,                  # specific vaccine distribution start (time after detection time)
                                             efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                             efficacy_disease_bpsv = seq(0.05, 1, 0.05),    # vaccine efficacy against disease - BPSV
                                             efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                             efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                             dur_R = 365000,                                # duration of infection-induced immunity
                                             dur_V = 365000,                                # duration of vaccine-induced immunity for both vaccines
                                             second_dose_delay = 1:14,                      # controls how many days after "1st dose" people receive second dose; see here: https://github.com/mrc-ide/squire.page/blob/main/inst/odin/nimue_booster.R#L427-L430
                                             dur_vacc_delay = 7,                            # mean duration from vaccination to protection (applied to primary, secondary and boosters)
                                             coverage = 0.8,                                # proportion of the population vaccinated
                                             vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                             min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                             min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)
 

## Generating defuault NPI scenarios (i.e. Rt and tt_Rt for model parameter combinations) and joining to parameter combos
lockdown_Rt <- 0.9
minimal_mandate_reduction <- 0.25

NPIs_bpsv_eff <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, NPI_scenarios = c(2, 4, 6), scenarios = raw_bpsv_scenarios)
bpsv_scenarios <- raw_bpsv_scenarios %>%
  full_join(NPIs_bpsv_eff, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                  "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")

NPIs_delay <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, NPI_scenarios = c(2, 4, 6), scenarios = raw_delay_scenarios)
delay_scenarios <- raw_delay_scenarios %>%
  full_join(NPIs_delay, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                               "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")

NPIs_bpsv_eff_delay <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, NPI_scenarios = c(2, 4, 6), scenarios = raw_delay_bpsv_scenarios)
delay_bpsv_scenarios <- raw_delay_bpsv_scenarios %>%
  full_join(NPIs_bpsv_eff_delay, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                        "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")

all_scenarios <- rbind(bpsv_scenarios, delay_scenarios, delay_bpsv_scenarios)

## Creating index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vars_for_index <- c(variable_columns(all_scenarios), "NPI_int")
all_scenarios <- all_scenarios %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = 60) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(all_scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = 50)
  saveRDS(model_outputs, "outputs/vaccine_properties_exploration_scenarios.rds")
} else {
  model_outputs <- readRDS("outputs/vaccine_properties_exploration_scenarios.rds")
}

## Plotting the output
colour_func <- scales::hue_pal()(9)
NPI_colours <- colour_func[c(2, 4, 6)]

unique(model_outputs$varied)


## Efficacy plot
disease_efficacy_plotting <- model_outputs %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "IFR", "efficacy_disease_bpsv")))) %>%
  group_by(efficacy_disease_bpsv, NPI_int) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted),
            max_deaths_averted = max(bpsv_deaths_averted),
            central_deaths_averted = bpsv_deaths_averted[R0 == 2.5 & IFR == 1])
disease_efficacy_plot <- ggplot(disease_efficacy_plotting) +
  geom_line(aes(x = 100 * efficacy_disease_bpsv, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
  geom_ribbon(aes(x = 100 * efficacy_disease_bpsv, ymin = min_deaths_averted, ymax = max_deaths_averted,
                  fill = factor(NPI_int)), alpha = 0.1) +
  geom_line(aes(x = 100 * efficacy_disease_bpsv, y = min_deaths_averted, col = factor(NPI_int)), size = 0.1) +
  geom_line(aes(x = 100 * efficacy_disease_bpsv, y = max_deaths_averted, col = factor(NPI_int)), size = 0.1) + 
  scale_colour_manual(values = NPI_colours) +
  scale_fill_manual(values = NPI_colours) +
  scale_x_continuous(breaks = c(0, 25, 50, 70, 100), labels = paste0(c(0, 25, 50, 75, 100), "%")) +
  theme_bw() +
  labs(x = "BPSV Disease Efficacy", y = "Deaths Averted By BPSV")


## Delay plot
delay_plotting <- model_outputs %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "IFR", "dur_vacc_delay")))) %>%
  group_by(dur_vacc_delay, NPI_int) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted),
            max_deaths_averted = max(bpsv_deaths_averted),
            central_deaths_averted = bpsv_deaths_averted[R0 == 2.5 & IFR == 1])
delay_plot <- ggplot(delay_plotting) +
  geom_line(aes(x = dur_vacc_delay, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
  geom_ribbon(aes(x = dur_vacc_delay, ymin = min_deaths_averted, ymax = max_deaths_averted,
                  fill = factor(NPI_int)), alpha = 0.1) +
  geom_line(aes(x = dur_vacc_delay, y = min_deaths_averted, col = factor(NPI_int)), size = 0.1) +
  geom_line(aes(x = dur_vacc_delay, y = max_deaths_averted, col = factor(NPI_int)), size = 0.1) + 
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
# NPI_df <- NPIs %>%
#   select(R0, detection_time, bpsv_start, specific_vaccine_start, time_to_coverage_bpsv, time_to_coverage_spec, NPI_int, Rt, tt_Rt) %>%
#   rowwise() %>%
#   mutate(scenario_info = list(tibble(Rt = Rt, tt_Rt = tt_Rt))) %>%
#   select(-Rt, -tt_Rt) %>%
#   unnest(cols = c(scenario_info)) %>%
#   mutate(scenario = paste0("Scenario ", NPI_int)) %>%
#   group_by(scenario) %>%
#   mutate(next_time = lead(tt_Rt),
#          next_value = lead(Rt)) %>%
#   mutate(next_time = ifelse(is.na(next_time), unique(baseline_scenarios$runtime), next_time),
#          next_value = ifelse(is.na(next_value), R0, next_value))
# overplot_factor <- 1
# 
# NPI_plot <- ggplot(NPI_df, aes(x = tt_Rt - overplot_factor, colour = scenario)) +
#   geom_hline(aes(yintercept = 1), linewidth = 0.2) +
#   geom_hline(aes(yintercept = lockdown_Rt), linetype = "dashed", linewidth = 0.2) +
#   geom_hline(aes(yintercept = R0 * (1 - minimal_mandate_reduction)), linetype = "dashed", linewidth = 0.2) +
#   geom_vline(aes(xintercept = detection_time), linewidth = 0.2) +
#   geom_vline(aes(xintercept = detection_time + bpsv_start + time_to_coverage_bpsv), linewidth = 0.2) +
#   geom_vline(aes(xintercept = detection_time + specific_vaccine_start + time_to_coverage_spec), linewidth = 0.2) +
#   geom_segment(aes(xend = next_time + overplot_factor, y = Rt, yend = Rt), size = 1) +
#   geom_segment(aes(x = next_time, xend = next_time, y = Rt, yend = next_value), size = 1) +
#   theme_bw() +
#   facet_wrap(~scenario) +
#   labs(x = "Time Since Spillover (Days)") +
#   coord_cartesian(xlim = c(0, unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start) + unique(NPI_df$time_to_coverage_spec) + 10), 
#                   ylim = c(0.5, unique(NPI_df$R0) + 0.5)) +
#   theme(legend.position = "none",
#         strip.background = element_rect(fill="#F5F5F5"))
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
