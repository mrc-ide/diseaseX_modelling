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

## Generating defuault NPI scenarios (i.e. Rt and tt_Rt for model parameter combinations) and joining to parameter combos
lockdown_Rt <- 0.9
minimal_mandate_reduction <- 0.25
NPIs <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, NPI_scenarios = c(2, 4, 6), scenarios = raw_bpsv_scenarios)
bpsv_scenarios <- raw_bpsv_scenarios %>%
  full_join(NPIs, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenarios
                         "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"),
            multiple = "all")

## Creating index for output
vars_for_index <- c(variable_columns(bpsv_scenarios), "NPI_int")
bpsv_scenarios <- bpsv_scenarios %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = 50) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(bpsv_scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = 50)
  saveRDS(model_outputs, "outputs/univariate_bpsv_efficacy.rds")
} else {
  model_outputs <- readRDS("outputs/univariate_bpsv_efficacy.rds")
}

## Plotting the output
colour_func <- scales::hue_pal()(9)
NPI_colours <- colour_func[c(2, 4, 6)]

bpsv_plotting <- model_outputs %>%
  group_by(efficacy_disease_bpsv, NPI_int) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted),
            max_deaths_averted = max(bpsv_deaths_averted),
            central_deaths_averted = bpsv_deaths_averted[R0 == 2.5 & IFR == 1])
colnames(bpsv_plotting)

ggplot(bpsv_plotting) +
  geom_line(aes(x = efficacy_disease_bpsv, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
  geom_ribbon(aes(x = efficacy_disease_bpsv, ymin = min_deaths_averted, ymax = max_deaths_averted,
                  fill = factor(NPI_int)), alpha = 0.1) +
  geom_line(aes(x = efficacy_disease_bpsv, y = min_deaths_averted, col = factor(NPI_int)), size = 0.1) +
  geom_line(aes(x = efficacy_disease_bpsv, y = max_deaths_averted, col = factor(NPI_int)), size = 0.1) + 
  scale_colour_manual(values = NPI_colours) +
  scale_fill_manual(values = NPI_colours)

ggplot(subset(model_outputs, IFR == 1 & NPI_int == 2), aes(x = efficacy_disease_bpsv, y = bpsv_deaths_averted,
                          col = R0)) +
  geom_point()





## Plotting the NPI Scenarios
NPI_df <- NPIs %>%
  select(R0, detection_time, bpsv_start, specific_vaccine_start, time_to_coverage_bpsv, time_to_coverage_spec, NPI_int, Rt, tt_Rt) %>%
  rowwise() %>%
  mutate(scenario_info = list(tibble(Rt = Rt, tt_Rt = tt_Rt))) %>%
  select(-Rt, -tt_Rt) %>%
  unnest(cols = c(scenario_info)) %>%
  mutate(scenario = paste0("Scenario ", NPI_int)) %>%
  group_by(scenario) %>%
  mutate(next_time = lead(tt_Rt),
         next_value = lead(Rt)) %>%
  mutate(next_time = ifelse(is.na(next_time), unique(baseline_scenarios$runtime), next_time),
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
  facet_wrap(~scenario) +
  labs(x = "Time Since Spillover (Days)") +
  coord_cartesian(xlim = c(0, unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start) + unique(NPI_df$time_to_coverage_spec) + 10), 
                  ylim = c(0.5, unique(NPI_df$R0) + 0.5)) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="#F5F5F5"))

# Plotting model outputted deaths and time under NPIs
model_outputs <- model_outputs %>%
  mutate(specific_vaccine_start = factor(specific_vaccine_start),
         detection_time = factor(detection_time),
         NPI_int = factor(NPI_int))

absolute_deaths_plot <- ggplot() +
  geom_segment(data = model_outputs, aes(x = composite_NPI_spec, xend = composite_NPI_spec, 
                                         y = deaths_spec, yend = deaths_bpsv + 75, group = factor(NPI_int)),
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_point(data = model_outputs, aes(x = composite_NPI_spec, y = deaths_spec, fill = factor(NPI_int)), 
             shape = 4, colour = "black", size = 2, pch = 21) +
  geom_point(data = model_outputs, aes(x = composite_NPI_spec, y = deaths_bpsv, fill = factor(NPI_int)), 
             colour = "black", size = 4, pch = 21) +
  theme_bw() +
  labs(x = "NPI Days (Composite Duration+Stringency", y = "Disease Deaths") +
  guides(fill = guide_legend(title = "Scenario"))

deaths_averted_plot <- ggplot() +
  geom_bar(data = model_outputs, aes(x = factor(NPI_int), y = bpsv_deaths_averted, fill = factor(NPI_int)), stat = "identity") +
  labs(x = "", y = "Deaths Averted") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.background = element_rect(colour = "black"))

inset_prop_y <- 0.4
min_deaths <- min(model_outputs$deaths_bpsv)
max_deaths <- max(model_outputs$deaths_spec)
inset_ymin <- min_deaths + (1 - inset_prop_y) * (max_deaths - min_deaths)
inset_ymax <- max_deaths

inset_prop_x <- 0.5
min_NPI <- min(model_outputs$composite_NPI_spec)
max_NPI <- max(model_outputs$composite_NPI_spec)
inset_xmin <- min_NPI + (1 - inset_prop_x) * (max_NPI - min_NPI)
inset_xmax <- max_NPI

overall_deaths_plot <- absolute_deaths_plot + 
  annotation_custom(
    ggplotGrob(deaths_averted_plot), 
    xmin = inset_xmin, xmax = inset_xmax, ymin = inset_ymin, ymax = inset_ymax)
cowplot::plot_grid(NPI_plot, overall_deaths_plot)


# example_NPI <- subset(NPI_df, scenario == "Scenario 1")
# example_NPI$next_time[length(example_NPI$next_time)] <- 150
# a <- ggplot(example_NPI) +
#   geom_segment(aes(x = tt_Rt - overplot_factor, xend = next_time + overplot_factor, y = Rt, yend = Rt), size = 2) +
#   geom_segment(aes(x = next_time, xend = next_time, y = Rt, yend = next_value), size = 2) +
#   geom_hline(aes(yintercept = 1)) +
#   geom_hline(aes(yintercept = lockdown_Rt), linetype = "dashed") +
#   geom_hline(aes(yintercept = temp_R0 * (1 - minimal_mandate_reduction)), linetype = "dashed") +
#   geom_vline(aes(xintercept = temp_detection_time)) +
#   geom_vline(aes(xintercept = temp_detection_time + bpsv_start + time_to_coverage_bpsv)) +
#   geom_vline(aes(xintercept = temp_detection_time + temp_specific_vaccine_start + time_to_coverage_spec)) +
#   theme_bw() +
#   geom_segment(aes(x = temp_detection_time, y = temp_R0 + 1, xend = temp_detection_time, yend = temp_R0 + 0.7), arrow = arrow(length = unit(0.3, "cm"))) +
#   annotate("text", x = temp_detection_time, y = temp_R0 + 1.2, label = "Detection") +
#   geom_segment(aes(x = temp_detection_time + bpsv_start + time_to_coverage_bpsv, y = temp_R0 + 1, xend = temp_detection_time + bpsv_start + time_to_coverage_bpsv, yend = temp_R0 + 0.7), arrow = arrow(length = unit(0.3, "cm"))) +
#   annotate("text", x = temp_detection_time + bpsv_start + time_to_coverage_bpsv, y = temp_R0 + 1.2, label = "BPSV Vacc.\nCompleted", hjust = 0.5) +
#   geom_segment(aes(x = temp_detection_time + temp_specific_vaccine_start + time_to_coverage_spec, y = temp_R0 + 1, xend = temp_detection_time + temp_specific_vaccine_start + time_to_coverage_spec, yend = temp_R0 + 0.7), arrow = arrow(length = unit(0.3, "cm"))) +
#   annotate("text", x = temp_detection_time + temp_specific_vaccine_start + time_to_coverage_spec, y = temp_R0 + 1.2, label = "Spec. Vacc.\nCompleted", hjust = 0.5) +
#   geom_segment(aes(x = -20, y = temp_R0 * (1 - minimal_mandate_reduction), xend = -10, yend = temp_R0 * (1 - minimal_mandate_reduction)), arrow = arrow(length = unit(0.3, "cm"))) +
#   annotate("text", x = -20, y = temp_R0 * (1 - minimal_mandate_reduction), label = "Min.\nMandate", hjust = 1) +
#   geom_segment(aes(x = -20, y = lockdown_Rt, xend = -10, yend = 0.9), arrow = arrow(length = unit(0.3, "cm"))) +
#   annotate("text", x = -20, y = lockdown_Rt, label = "Lockdown", hjust = 1) +
#   geom_segment(aes(x = -20, y = temp_R0, xend = -10, yend = temp_R0), arrow = arrow(length = unit(0.3, "cm"))) +
#   annotate("text", x = -20, y = temp_R0, label = "R0", hjust = 1) +
#   theme(plot.margin = margin(2.5, 1, 2.5, 2.5, "cm")) +
#   coord_cartesian(clip = 'off', xlim = c(0, 150), ylim = c(0.5, temp_R0 + 0.5)) +
#   scale_y_continuous(position = "right") +
#   labs(x = "Time (Days)")