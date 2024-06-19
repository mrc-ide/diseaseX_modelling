# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/helper_functions.R"))

## Generating default parameters
default <- define_default_params()

# Getting the detection time
R0_subset <- c(1.5, 2.5, 3.5)
detection_theshold_hosp <- 5
bp_df_long <- readRDS("outputs/Figure1_branchingProcess_Containment/Figure1_bp_detection_times.rds")
detection_threshold_inf <- unique(bp_df_long$detection)[detection_theshold_hosp]
bp_subset <- bp_df_long %>%
  filter(R0 %in% R0_subset, detection == detection_threshold_inf, metric == "Daily Incidence") %>%
  select(-iteration) %>%
  filter(!is.infinite(value)) %>%
  group_by(R0) %>%
  summarise(time_to_detection = round(mean(value)))
time_to_detection <- round(bp_subset$time_to_detection)

# Generate parameter combinations for model running
raw_baseline_scenarios <- create_scenarios(R0 = R0_subset,                                                # Basic reproduction number
                                           specific_vaccine_start = c(100, 250, 365))                     # Specific vaccine distribution start (time after detection time)

# Linking these scenarios with the R0-specific detection times for a given threshold
### Note - this generates all R0 and detection time combinations when there should be a unique pair of R0 & detection time - we sort this below
baseline_scenarios <- expand_grid(raw_baseline_scenarios,
                                  time_to_detection = time_to_detection) %>%
  left_join(bp_subset, by = c("R0", "time_to_detection")) %>%
  mutate(detection_time = time_to_detection)

## Generating default NPI scenarios (i.e. Rt and tt_Rt for model parameter combinations) and joining to parameter combos
NPIs <- default_NPI_scenarios(lockdown_Rt = default$lockdown_Rt, minimal_mandate_reduction = default$minimal_mandate_reduction, 
                              NPI_scenarios = 1:9, scenarios = baseline_scenarios)
scenarios_NPIs <- baseline_scenarios %>%
  full_join(NPIs, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence timing of NPI scenarios
                         "specific_vaccine_start", "vaccination_rate_bpsv", "vaccination_rate_spec",
                         "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all")

# Filtering the above to only select R0 and detection time pairs that actually occurred (function above produces all pairwise combos of them)
R0_detection_time_pairs <- bp_subset %>%
  mutate(detection_time = round(time_to_detection, digits = 0))
baseline_scenarios_reduced <- scenarios_NPIs %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time"))
# R0 * NPI * spec start * 2 vaccination scenarios
3 * 9 * 3 * 2 

## Creating index for output
vars_for_index <- c(variable_columns(baseline_scenarios_reduced), "NPI_int")
scenarios <- baseline_scenarios_reduced %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = 4) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = 2)
  saveRDS(model_outputs, "outputs/Figure2_NPI_Exploration/NEW_Figure2_NPI_Exploration_Outputs.rds")
} else {
  model_outputs <- readRDS("outputs/Figure2_NPI_Exploration/NEW_Figure2_NPI_Exploration_Outputs.rds")
}

## Plotting the NPI Scenarios
NPI_df <- NPIs %>%
  filter(specific_vaccine_start == 250, R0 == 2.5, detection_time == 62) %>%
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
                               composite = model_outputs$composite_NPI_bpsv[model_outputs$specific_vaccine_start == 100 & model_outputs$R0 == 2.5],
                               scenario = unique(NPI_df$scenario))
NPI_df$NPI_int <- factor(NPI_df$NPI_int, levels = NPI_composite_df$NPI_int[order(NPI_composite_df$composite)])
NPI_df$scenario <- factor(NPI_df$scenario, levels = NPI_composite_df$scenario[order(NPI_composite_df$composite)])
colour_func <- scales::hue_pal()(max(model_outputs$NPI_int))
colour_func2 <- colour_func[order(NPI_composite_df$composite)]
NPI_df$scenario2 <- paste0("Scenario ", as.integer(NPI_df$scenario))

NPI_plot <- ggplot(NPI_df, aes(x = tt_Rt - overplot_factor, colour = scenario)) +
  geom_hline(aes(yintercept = 1), linewidth = 0.2) +
  geom_hline(aes(yintercept = default$lockdown_Rt), linetype = "dashed", linewidth = 0.2) +
  geom_hline(aes(yintercept = R0 * (1 - default$minimal_mandate_reduction)), linetype = "dashed", linewidth = 0.2) +
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
spec_dev_scenarios <- c(100, 250)
model_outputs2 <- model_outputs %>%
  mutate(specific_vaccine_start = factor(specific_vaccine_start),
         detection_time = factor(detection_time),
         NPI_int = factor(NPI_int))
absolute_deaths_plot <- ggplot() +
  geom_segment(data = subset(model_outputs2, R0 == 2.5 & specific_vaccine_start %in% spec_dev_scenarios), 
               aes(x = composite_NPI_spec, xend = composite_NPI_spec, 
                   y = deaths_spec * 1000 / default$population_size, yend = deaths_bpsv * 1000 / default$population_size + 0.2, group = factor(NPI_int)),
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_point(data = subset(model_outputs2, R0 == 2.5 & specific_vaccine_start %in% spec_dev_scenarios),
             aes(x = composite_NPI_spec, y = deaths_spec * 1000 / default$population_size, fill = factor(NPI_int)), 
             shape = 4, colour = "black", size = 2, pch = 21) +
  geom_point(data = subset(model_outputs2, R0 == 2.5 & specific_vaccine_start %in% spec_dev_scenarios),
             aes(x = composite_NPI_spec, y = deaths_bpsv * 1000 / default$population_size, fill = factor(NPI_int)), 
             colour = "black", size = 4, pch = 21) +
  facet_wrap(specific_vaccine_start ~ ., scales = "free_x",
             nrow = 2, strip.position = "right",
             labeller = as_labeller(c(`100`='Specific Vaccine In 100 Days', `250`='Specific Vaccine In 250 Days'))) +
  scale_fill_manual(values = colour_func) +
  theme_bw() +
  lims(y = c(0, max(model_outputs2$deaths_spec[model_outputs2$R0 == 2.5] * 1000 / default$population_size))) +
  labs(x = "NPI Days (Composite Duration & Stringency)", y = "Disease Deaths Per 1000 Population") +
  theme(strip.placement = "outside",
        legend.position = "none",
        strip.background = element_rect(fill="white")) +
  guides(fill = guide_legend(title = "Scenario"))

model_outputs2$NPI_int2 <- factor(model_outputs2$NPI_int, levels = NPI_composite_df$NPI_int[order(NPI_composite_df$composite)])
deaths_averted_100_plot <- ggplot() +
  geom_bar(data = subset(model_outputs2, R0 == 2.5 & specific_vaccine_start %in% c(100)), 
           aes(x = factor(NPI_int2), y = bpsv_deaths_averted * 1000 / default$population_size, fill = factor(NPI_int)), stat = "identity") +
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
                                composite = model_outputs2$composite_NPI_bpsv[model_outputs2$specific_vaccine_start == 250 & model_outputs2$R0 == 2.5],
                                scenario = unique(NPI_df$scenario))
model_outputs2$NPI_int3 <- factor(model_outputs2$NPI_int, levels = NPI_composite_df2$NPI_int[order(NPI_composite_df2$composite)])
deaths_averted_200_plot <- ggplot() +
  geom_bar(data = subset(model_outputs2, R0 == 2.5 & specific_vaccine_start %in% c(250)), 
           aes(x = factor(NPI_int3), y = bpsv_deaths_averted * 1000 / default$population_size, fill = factor(NPI_int)), stat = "identity") +
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
  gg_inset(ggplot2::ggplotGrob(deaths_averted_200_plot), data = data.frame(specific_vaccine_start = 250),
           xmin = ggplot_build(absolute_deaths_plot)$layout$panel_params[[2]]$x.range[1] +
             ((ggplot_build(absolute_deaths_plot)$layout$panel_params[[2]]$x.range[2] -
                 ggplot_build(absolute_deaths_plot)$layout$panel_params[[2]]$x.range[1]) / 1.5), xmax = Inf, 
           ymin = ggplot_build(absolute_deaths_plot)$layout$panel_params[[2]]$y.range[2] / 3, ymax = Inf)

Figure2 <- cowplot::plot_grid(NPI_plot, deaths_averted_plot,
                              nrow = 1, rel_widths = c(1, 1.4),
                              labels = c("B", "C"))
ggsave(filename = "figures/Figure_2_FrameworkIntro/NEW_Figure2BC_NPI_Plot.pdf", plot = Figure2, width = 9.25, height = 7.5)


test <- subset(model_outputs2, R0 == 2.5 & specific_vaccine_start %in% spec_dev_scenarios)
x <- data.frame(specific_vaccine_start = test$specific_vaccine_start,
                deaths = test$deaths_spec, 
                NPI_int = test$NPI_int,
                NPI_days = test$composite_NPI_spec)
frontiers <- x %>% 
  group_by(specific_vaccine_start) %>% 
  group_modify(~pareto_frontier(.x))
frontiers$new_NPI_days <- floor(frontiers$NPI_days)
interpolated_frontiers <- frontiers %>% 
  group_by(specific_vaccine_start) %>% 
  group_modify(~linear_interpolate(.x))
test$bpsv_deaths_averted[test$bpsv_deaths_averted < 0] <- 0
test_250 <- test %>%
  filter(specific_vaccine_start == 250)
interpolated_frontiers_250 <- interpolated_frontiers %>%
  filter(specific_vaccine_start == 250)
new_NPI_days_250 <- sapply(test_250$deaths_bpsv, function(x) {
  temp_df <- interpolated_frontiers_250[which.min(abs(x - interpolated_frontiers_250$deaths)), ]
  return(temp_df$new_NPI_days)
})
test_100 <- test %>%
  filter(specific_vaccine_start == 100)
interpolated_frontiers_100 <- interpolated_frontiers %>%
  filter(specific_vaccine_start == 100)
new_NPI_days_100 <- sapply(test_100$deaths_bpsv, function(x) {
  temp_df <- interpolated_frontiers_100[which.min(abs(x - interpolated_frontiers_100$deaths)), ]
  return(temp_df$new_NPI_days)
})

new_NPI_days_250 - test_250$composite_NPI_bpsv
ordering_250 <- NPI_composite_df$NPI_int[order(NPI_composite_df$composite)]
plot((new_NPI_days_250 - test_250$composite_NPI_bpsv)[ordering_250])

new_NPI_days_100 - test_100$composite_NPI_bpsv
ordering_100 <- NPI_composite_df$NPI_int[order(NPI_composite_df$composite)]
plot((new_NPI_days_100 - test_100$composite_NPI_bpsv)[ordering_100])

new <- data.frame(specific_vaccine_start = test$specific_vaccine_start,
           NPI_days_averted = c((new_NPI_days_100 - test_100$composite_NPI_bpsv)[ordering_100],
                                (new_NPI_days_250 - test_250$composite_NPI_bpsv)[ordering_250]),
           NPI_int = test$NPI_int)
new$NPI_days_averted[new$NPI_days_averted < 0] <- 0.1
new$NPI_int2 <- factor(new$NPI_int, levels = NPI_composite_df$NPI_int[order(NPI_composite_df$composite)])
NPI_days_averted <- ggplot(new) +
  geom_bar(aes(x = NPI_int, y = NPI_days_averted, fill = NPI_int), stat = "identity") +
  facet_grid(specific_vaccine_start ~ ., scales = "free_y",
             labeller = as_labeller(c(`100`='Specific Vaccine In 100 Days', `250`='Specific Vaccine In 250 Days'))) +
  scale_fill_manual(values = colour_func2) +
  theme_bw() +
  labs(y = "NPI Days Averted by BPSV") +
  theme(legend.position = "none",
        axis.title.x =  element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="white"))


x <- cowplot::plot_grid(deaths_averted_plot, NPI_days_averted,
                   nrow = 1, rel_widths = c(1.4, 2.4/3.75),
                   labels = c("C", "D"),
                   align = "h", axis = "tb")
Figure2_new <- cowplot::plot_grid(NPI_plot, x,
                        nrow = 1, rel_widths = c(1, 1.4 + 2.4/3.75),
                        labels = c("B", NA))
ggsave(filename = "figures/Figure_2_FrameworkIntro/NEW_Figure2BCD_NPI_Plot_DaysAverted.pdf", 
       plot = Figure2_new, width = 9.25 + (9.25 * 1/3.75), height = 7.5)

## Supplementary Figure
absolute_deaths_plot_supp <- ggplot() +
  geom_segment(data = subset(model_outputs2, R0 != 2.5), 
               aes(x = composite_NPI_spec, xend = composite_NPI_spec, 
                   y = deaths_spec * 1000 / default$population_size, yend = deaths_bpsv * 1000 / default$population_size + 0.2, group = factor(NPI_int)),
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_point(data = subset(model_outputs2, R0 != 2.5),
             aes(x = composite_NPI_spec, y = deaths_spec * 1000 / default$population_size, fill = factor(NPI_int)), 
             shape = 4, colour = "black", size = 2, pch = 21) +
  geom_point(data = subset(model_outputs2, R0 != 2.5),
             aes(x = composite_NPI_spec, y = deaths_bpsv * 1000 / default$population_size, fill = factor(NPI_int)), 
             colour = "black", size = 4, pch = 21) +
  facet_grid(specific_vaccine_start ~ R0, scales = "free_x",
             labeller = as_labeller(c(`100`='VSV In 100 Days', 
                                      `250`='VSV In 250 Days',
                                      `365`='VSV In 365 Days',
                                      `1.5`='R0=1.5',
                                      `3.5`='R0=3.5'))) +
  scale_fill_manual(values = colour_func) +
  theme_bw() +
  lims(y = c(0, max(model_outputs2$deaths_spec[model_outputs2$R0 == 3.5] * 1000 / default$population_size))) +
  labs(x = "NPI Days (Composite Duration & Stringency)", y = "Disease Deaths Per 1000 Population") +
  theme(strip.placement = "outside",
        legend.position = "none",
        strip.background = element_rect(fill="white")) +
  guides(fill = guide_legend(title = "Scenario"))

deaths_averted_supp_plot <- ggplot() +
  geom_bar(data = subset(model_outputs2, R0 != 2.5), 
           aes(x = factor(NPI_int2), y = bpsv_deaths_averted * 1000 / default$population_size, fill = factor(NPI_int)), stat = "identity") +
  labs(x = "NPI Scenario", y = "BPSV Deaths Averted\nPer 1000 Pop") +
  scale_fill_manual(values = colour_func) +
  facet_grid(specific_vaccine_start ~ R0, scales = "free_x",
             labeller = as_labeller(c(`100`='VSV In 100 Days', 
                                      `250`='VSV In 250 Days',
                                      `365`='VSV In 365 Days',
                                      `1.5`='R0=1.5',
                                      `3.5`='R0=3.5'))) +
  scale_x_discrete(labels = 1:9) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill="white"))


NPI_plot_supp <- NPI_plot + 
  facet_wrap(~scenario2, nrow = 1) +
  theme(axis.text.x = element_text(size = 6))

deaths_supp_plot <- cowplot::plot_grid(absolute_deaths_plot_supp, deaths_averted_supp_plot,
                                       nrow = 1, labels = c("B", "C"))
supp_deaths_plot_final <- cowplot::plot_grid(NPI_plot_supp, deaths_supp_plot, nrow = 2, rel_heights = c(1, 3),
                                             labels = c("A", NA))
ggsave(filename = "figures/Figure_2_FrameworkIntro/Supp_Figure1_NPI_Plot_Sensitivity.pdf",
       plot = supp_deaths_plot_final,
       height = 8,
       width = 13)
