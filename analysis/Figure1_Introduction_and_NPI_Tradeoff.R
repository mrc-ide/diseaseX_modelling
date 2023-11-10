# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))

# Getting the detection time
R0_subset <- c(1.5, 2.5, 3.5)
detection_theshold_hosp <- 5
bp_df_long <- readRDS("outputs/Figure1_bp_detection_times.rds")
detection_threshold_inf <- unique(bp_df_long$detection)[detection_theshold_hosp]
bp_subset <- bp_df_long %>%
  filter(R0 %in% R0_subset, detection == detection_threshold_inf, metric == "Daily Incidence") %>%
  select(-iteration) %>%
  filter(!is.infinite(value)) %>%
  group_by(R0) %>%
  summarise(time_to_detection = round(mean(value)))
time_to_detection <- round(bp_subset$time_to_detection)

# Generate parameter combinations for model running (note Rt and tt_Rt has a placeholder)
raw_baseline_scenarios <- create_scenarios(R0 = R0_subset,                                # Basic reproduction number
                                           IFR = 1,                                       # IFR
                                           population_size = 10^10,                       # population size
                                           Tg = 6.7,                                      # Tg
                                           detection_time = 1,                            # detection time PLACEHOLDER FOR NOW
                                           bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                           bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                           specific_vaccine_start = c(100, 220, 365),     # specific vaccine distribution start (time after detection time)
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
                                           seeding_cases = 1)

# Linking these scenarios with the R0-specific detection times for a given threshold
baseline_scenarios <- expand_grid(raw_baseline_scenarios,
                                  time_to_detection = time_to_detection) %>%
  left_join(bp_subset, by = c("R0", "time_to_detection")) %>%
  mutate(detection_time = time_to_detection)

## Generating defuault NPI scenarios (i.e. Rt and tt_Rt for model parameter combinations) and joining to parameter combos
lockdown_Rt <- 0.9
minimal_mandate_reduction <- 0.25
NPIs <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                              NPI_scenarios = 1:9, scenarios = baseline_scenarios)
scenarios_NPIs <- baseline_scenarios %>%
  full_join(NPIs, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenarios
                         "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"),
            multiple = "all")

# Filtering the above to only select R0 and detection time pairs that actually occurred (function above produces all pairwise combos of them)
R0_detection_time_pairs <- bp_subset %>%
  mutate(detection_time = round(time_to_detection, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time)
baseline_scenarios_reduced <- scenarios_NPIs %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time"))

## Creating index for output
vars_for_index <- c(variable_columns(baseline_scenarios_reduced), "NPI_int")
scenarios <- baseline_scenarios_reduced %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
fresh_run <- TRUE
if (fresh_run) {
  plan(multisession, workers = 4) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = 2)
  saveRDS(model_outputs, "outputs/Figure2_NPI_Exploration_Outputs.rds")
} else {
  model_outputs <- readRDS("outputs/Figure2_NPI_Exploration_Outputs.rds")
}

## Plotting the NPI Scenarios
NPI_df <- NPIs %>%
  filter(specific_vaccine_start == 220, R0 == 2.5, detection_time == 62) %>%
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

NPI_composite_df <- data.frame(NPI_int = 1:9, 
                               composite = model_outputs2$composite_NPI_bpsv[model_outputs2$specific_vaccine_start == 100 & model_outputs2$R0 == 2.5],
                               scenario = unique(NPI_df$scenario))
NPI_df$NPI_int <- factor(NPI_df$NPI_int, levels = NPI_composite_df$NPI_int[order(NPI_composite_df$composite)])
NPI_df$scenario <- factor(NPI_df$scenario, levels = NPI_composite_df$scenario[order(NPI_composite_df$composite)])
colour_func <- scales::hue_pal()(max(model_outputs$NPI_int))
colour_func2 <- colour_func[order(NPI_composite_df$composite)]
NPI_df$scenario2 <- paste0("Scenario ", as.integer(NPI_df$scenario))

NPI_plot <- ggplot(NPI_df, aes(x = tt_Rt - overplot_factor, colour = scenario)) +
  geom_hline(aes(yintercept = 1), linewidth = 0.2) +
  geom_hline(aes(yintercept = lockdown_Rt), linetype = "dashed", linewidth = 0.2) +
  geom_hline(aes(yintercept = R0 * (1 - minimal_mandate_reduction)), linetype = "dashed", linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time + bpsv_start + time_to_coverage_bpsv), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time + specific_vaccine_start), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time + specific_vaccine_start + time_to_coverage_spec), linewidth = 0.2) +
  geom_segment(aes(xend = next_time + overplot_factor, y = Rt, yend = Rt), size = 1) +
  geom_segment(aes(x = next_time, xend = next_time, y = Rt, yend = next_value), size = 1) +
  theme_bw() +
  scale_colour_manual(values = colour_func2) +
  facet_wrap(~scenario2, nrow = 5) +
  coord_cartesian(xlim = c(0, unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start) + unique(NPI_df$time_to_coverage_spec) + 50), 
                  ylim = c(0.5, unique(NPI_df$R0) + 0.5)) +
  scale_x_continuous(breaks = c(unique(NPI_df$detection_time),
                                unique(NPI_df$detection_time) + unique(NPI_df$bpsv_start) + unique(NPI_df$time_to_coverage_bpsv),
                                unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start),
                                unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start) + unique(NPI_df$time_to_coverage_spec)),
                     labels = c("Detection", "BPSV Finish", "Spec Start", "Spec Finish")) +
  scale_y_continuous(breaks = c(1, 2.5), labels = c(1, "R0")) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="#F5F5F5"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.05),
        axis.title.x = element_blank())

# Plotting model outputted deaths and time under NPIs
model_outputs2 <- model_outputs %>%
  mutate(specific_vaccine_start = factor(specific_vaccine_start),
         detection_time = factor(detection_time),
         NPI_int = factor(NPI_int))

absolute_deaths_plot <- ggplot() +
  geom_segment(data = subset(model_outputs2, R0 == 2.5 & specific_vaccine_start %in% c(100, 220)), 
               aes(x = composite_NPI_spec, xend = composite_NPI_spec, 
                   y = deaths_spec * 1000 / population_size, yend = deaths_bpsv * 1000 / population_size + 0.2, group = factor(NPI_int)),
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_point(data = subset(model_outputs2, R0 == 2.5 & specific_vaccine_start %in% c(100, 220)),
             aes(x = composite_NPI_spec, y = deaths_spec * 1000 / population_size, fill = factor(NPI_int)), 
             shape = 4, colour = "black", size = 2, pch = 21) +
  geom_point(data = subset(model_outputs2, R0 == 2.5 & specific_vaccine_start %in% c(100, 220)),
             aes(x = composite_NPI_spec, y = deaths_bpsv * 1000 / population_size, fill = factor(NPI_int)), 
             colour = "black", size = 4, pch = 21) +
  facet_wrap(specific_vaccine_start ~ ., scales = "free_x",
             nrow = 2, strip.position = "right",
             labeller = as_labeller(c(`100`='Specific Vaccine In 100 Days', `220`='Specific Vaccine In 220 Days'))) +
  scale_fill_manual(values = colour_func) +
  theme_bw() +
  lims(y = c(0, max(model_outputs2$deaths_spec[model_outputs2$R0 == 2.5] * 1000 / population_size))) +
  labs(x = "NPI Days (Composite Duration & Stringency)", y = "Disease Deaths Per 1000 Population") +
  theme(strip.placement = "outside",
        legend.position = "none",
        strip.background = element_rect(fill="white")) +
  guides(fill = guide_legend(title = "Scenario"))

model_outputs2$NPI_int2 <- factor(model_outputs2$NPI_int, 
                                  levels = NPI_composite_df$NPI_int[order(NPI_composite_df$composite)])
deaths_averted_100_plot <- ggplot() +
  geom_bar(data = subset(model_outputs2, R0 == 2.5 & specific_vaccine_start %in% c(100)), 
           aes(x = factor(NPI_int2), y = bpsv_deaths_averted * 1000 / population_size, fill = factor(NPI_int)), stat = "identity") +
  labs(x = "BPSV Deaths Averted\nPer 1000 Pop", y = "Deaths Averted") +
  scale_fill_manual(values = colour_func) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x =  element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.background = element_rect(colour = "black"))

NPI_composite_df2 <- data.frame(NPI_int = 1:9, 
                                composite = model_outputs2$composite_NPI_bpsv[model_outputs2$specific_vaccine_start == 200 & model_outputs2$R0 == 2.5],
                                scenario = unique(NPI_df$scenario))
model_outputs2$NPI_int3 <- factor(model_outputs2$NPI_int, levels = NPI_composite_df2$NPI_int[order(NPI_composite_df2$composite)])
deaths_averted_200_plot <- ggplot() +
  geom_bar(data = subset(model_outputs2, R0 == 2.5 & specific_vaccine_start %in% c(200)), 
           aes(x = factor(NPI_int3), y = bpsv_deaths_averted * 1000 / population_size, fill = factor(NPI_int)), stat = "identity") +
  labs(x = "BPSV Deaths Averted\nPer 1000 Pop", y = "Deaths Averted") +
  scale_fill_manual(values = colour_func) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x =  element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.background = element_rect(colour = "black"))

deaths_averted_plot <- absolute_deaths_plot +
  gg_inset(ggplot2::ggplotGrob(deaths_averted_100_plot), data = data.frame(specific_vaccine_start = 100),
           xmin = ggplot_build(absolute_deaths_plot)$layout$panel_params[[1]]$x.range[1] +
             ((ggplot_build(absolute_deaths_plot)$layout$panel_params[[1]]$x.range[2] -
                 ggplot_build(absolute_deaths_plot)$layout$panel_params[[1]]$x.range[1]) / 1.5), xmax = Inf, 
           ymin = ggplot_build(absolute_deaths_plot)$layout$panel_params[[1]]$y.range[2] / 3, ymax = Inf) +
  gg_inset(ggplot2::ggplotGrob(deaths_averted_200_plot), data = data.frame(specific_vaccine_start = 200),
           xmin = ggplot_build(absolute_deaths_plot)$layout$panel_params[[2]]$x.range[1] +
             ((ggplot_build(absolute_deaths_plot)$layout$panel_params[[2]]$x.range[2] -
                 ggplot_build(absolute_deaths_plot)$layout$panel_params[[2]]$x.range[1]) / 1.5), xmax = Inf, 
           ymin = ggplot_build(absolute_deaths_plot)$layout$panel_params[[2]]$y.range[2] / 3, ymax = Inf)

Figure2 <- cowplot::plot_grid(NPI_plot, deaths_averted_plot,
                              nrow = 1, rel_widths = c(1, 1.4),
                              labels = c("B", "C"))
ggsave(filename = "figures/Figure2BC_NPI_Plot.pdf", plot = Figure2, width = 9.25, height = 7.5)
