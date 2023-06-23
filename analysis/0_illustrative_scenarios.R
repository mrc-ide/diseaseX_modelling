# Notes:
## 1) Double check IFR calc (fine to just alter prob_hosp?)
## 2) Are we fine assuming unlimited healthcare capacity?

# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))

# Generate parameter combinations for model running (note Rt and tt_Rt has a placeholder)
baseline_scenarios <- create_scenarios(#### Demographic Parameters
                                       population_size = 1e6,
                                       country = "Argentina",
                                        
                                       #### Healthcare Parameters
                                       hosp_bed_capacity = 100000000,                                          
                                       ICU_bed_capacity = 100000000, 
                                       
                                       #### Epidemiological Parameters
                                       R0 = c(1.5, 2, 3),
                                       Tg = 7,
                                       IFR = c(0.5, 1.5),
                                       number_NPI_scenarios = 9,
                                       
                                       ## Vaccine-Related Parameters
                                       detection_time = c(14, 28),                    # detection time
                                       bpsv_start = 14,                               # BPSV distribution start (time after detection time)
                                       specific_vaccine_start = c(100, 200, 365),     # specific vaccine distribution start (time after detection time)
                                       efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                       efficacy_disease_bpsv = 0.8,                   # vaccine efficacy against disease - BPSV
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

## Generating NPI scenarios (i.e. Rt and tt_Rt for model parameter combinations)
NPIs <- default_NPI_scenarios(lockdown_Rt = 0.9, minimal_mandate_reduction = 0.25, scenarios = baseline_scenarios)

# Identifying varying NPI-relevant parameters and linking NPI scenarios df and parameters dfs together by these variables
all_varying <- names(which(apply(baseline_scenarios, 2, function(x) {
  length(unique(x))
}) > 1))
NPI_relevant_varying <- all_varying[all_varying %in% c("R0", "detection_time", "bpsv_start", "specific_vaccine_start", 
                                                       "coverage", "vaccination_rate", "min_age_group_index_priority", "NPI_int")]
scenarios <- baseline_scenarios %>%
  left_join(NPIs, by = NPI_relevant_varying)
scenarios$index <- 1:nrow(scenarios)

## Running 
tic()
plan(multisession, workers = 6) # multicore does nothing on windows as multicore isn't supported
system.time({out <- future_pmap(scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
toc()
saveRDS(out, "outputs/test_raw_outputs.rds")

### need to add NPI scenario in here
cl <- makeCluster(6)
clusterEvalQ(cl, {
  library(data.table)
})
tic()
data <- parLapply(cl, out, function(x) {
  y <- data.frame(index = x$model_arguments$index, 
                  x$summary_metrics, 
                  country = x$model_arguments$country,
                  population_size = x$model_arguments$population_size,
                  hosp_bed_capacity = x$model_arguments$hosp_bed_capacity,
                  ICU_bed_capacity = x$model_arguments$ICU_bed_capacity,
                  R0 = x$model_arguments$Rt[1],
                  Tg = x$model_arguments$Tg,
                  IFR = x$model_arguments$IFR,
                  vaccine_scenario = x$model_arguments$vaccine_scenario,
                  detection_time = x$model_arguments$detection_time,
                  bpsv_start = ifelse(x$model_arguments$vaccine_scenario == "specific_only", NA, x$model_arguments$bpsv_start),
                  specific_vaccine_start = x$model_arguments$specific_vaccine_start,
                  efficacy_infection_bpsv = x$model_arguments$efficacy_infection_bpsv,
                  efficacy_disease_bpsv = x$model_arguments$efficacy_disease_bpsv,
                  efficacy_infection_spec = x$model_arguments$efficacy_infection_spec,
                  efficacy_disease_spec = x$model_arguments$efficacy_disease_spec,
                  dur_R = x$model_arguments$dur_R,
                  dur_V = x$model_arguments$dur_V[1],
                  second_dose_delay = x$model_arguments$second_dose_delay,
                  dur_vacc_delay = x$model_arguments$dur_vacc_delay,
                  coverage = x$model_arguments$coverage,
                  vaccination_rate = x$model_arguments$vaccination_rate,
                  min_age_group_index_priority = x$model_arguments$min_age_group_index_priority,
                  min_age_group_index_non_priority = x$model_arguments$min_age_group_index_non_priority,
                  runtime = x$model_arguments$runtime,
                  seeding_cases = x$model_arguments$seeding_cases,
                  NPI_scenario = x$model_arguments$NPI_scenario,
                  NPI_int = x$model_arguments$NPI_int) # newest runs this'll work fine but for current run should be NPI_scenario_int 
})
toc()
stopCluster(cl) 

combined_data <- rbindlist(data)
saveRDS(combined_data, "outputs/bpsv_efficacy_sensitivity.rds")

combined_data <- readRDS("outputs/bpsv_efficacy_sensitivity.rds")

two_vax <- combined_data %>%
  filter(vaccine_scenario == "both_vaccines") %>%
  rename(deaths_bpsv = deaths,
         time_under_NPIs_bpsv = time_under_NPIs,
         composite_NPI_bpsv = composite_NPI)

one_vax <- combined_data %>%
  filter(vaccine_scenario == "specific_only") %>%
  select(R0, Tg, IFR, NPI_scenario, detection_time, specific_vaccine_start, deaths, time_under_NPIs, composite_NPI, efficacy_disease_bpsv) %>% ## IMPORTANT - YOU NEED TO INCLUDE THE VARIABLE YOU'RE DOING THE SENSITIVITY ANALYSIS OVER
  rename(deaths_spec = deaths,
         time_under_NPIs_spec = time_under_NPIs,
         composite_NPI_spec = composite_NPI)

# create vector of all variables that are varying, so that you don't have to worry about figuring out which ones each time
# if you miss one out, this can fuck up silently and you'll be comparing the wrong pairs of scenarios
joined <- two_vax %>%
  left_join(one_vax, by = c("R0", "Tg", "IFR", "NPI_scenario", "detection_time", "specific_vaccine_start", "efficacy_disease_bpsv")) %>% ## IMPORTANT - YOU NEED TO INCLUDE THE VARIABLE YOU'RE DOING THE SENSITIVITY ANALYSIS OVER
  mutate(deaths_saved = deaths_spec - deaths_bpsv)

plot(joined$deaths_bpsv, joined$deaths_spec, xlim = c(0, 17500), ylim = c(0, 17500))
lines(1:17500, 1:17500)
plot(joined$efficacy_disease_bpsv, joined$deaths_saved)

##### create some holder variables here for specific_vaccine_start, R0 etc that will
##### enable you to filter both the results plots and the NPI scenario plots by the same 
##### things and ensure what you present in the nice figure is underpinned by the same set of runs
test_plot <- joined %>%
  filter(R0 == 2,
         Tg == 7,
         IFR == 1.5) %>%
  mutate(specific_vaccine_start = factor(specific_vaccine_start),
         detection_time = factor(detection_time),
         NPI_int = factor(NPI_int))

ggplot(test_plot, aes(x = efficacy_disease_bpsv, y = deaths_saved, colour = specific_vaccine_start)) +
  geom_path() +
  facet_grid(paste0("Detection Time: ", detection_time) ~ paste0("Scenario ", NPI_int))

test_plot2 <- joined %>%
  mutate(pathogen = ifelse(R0 == 3 & Tg == 7 & IFR == 0.5, "SARS-CoV-2-Like", 
                           ifelse(R0 == 1.5 & Tg == 14 & IFR == 1.5, "SARS-CoV-1-Like", "not_relevant"))) %>%
  filter(specific_vaccine_start %in% c(100, 200, 365),
         pathogen != "not_relevant", 
         detection_time == 14) %>%
  group_by(efficacy_disease_bpsv, pathogen, NPI_int) %>%
  mutate(min_ribbon = min(deaths_saved),
         max_ribbon = max(deaths_saved))

ggplot() +
  geom_line(data = test_plot2, aes(x = efficacy_disease_bpsv, y = deaths_saved, 
                  group = interaction(NPI_int, pathogen, factor(specific_vaccine_start)), 
                  colour = factor(specific_vaccine_start))) +
  geom_ribbon(data = test_plot2,
              aes(x = efficacy_disease_bpsv,
                  group = interaction(NPI_int, pathogen, factor(specific_vaccine_start)),
                  ymin = min_ribbon,
                  ymax = max_ribbon), alpha = 0.2) +
  facet_grid(pathogen ~ paste0("Scenario ", NPI_int))

ggplot() +
  geom_line(data = subset(test_plot2, NPI_int == 1 ), 
            aes(x = efficacy_disease_bpsv, y = deaths_saved, 
                group = interaction(NPI_int, pathogen, factor(specific_vaccine_start)), 
                colour = factor(specific_vaccine_start))) # +
  # geom_ribbon(data = subset(test_plot2, NPI_int == 1),
  #             aes(x = efficacy_disease_bpsv,
  #                 group = interaction(NPI_int, pathogen, factor(specific_vaccine_start)),
  #                 ymin = min_ribbon,
  #                 ymax = max_ribbon), alpha = 0.2)

temp_plot <- ggplot() +
  geom_segment(data = subset(test_plot2, pathogen == "SARS-CoV-2-Like" & specific_vaccine_start == 100),
               aes(x = composite_NPI_spec, xend = composite_NPI_spec, y = deaths_spec, yend = deaths_bpsv + 75,
                   group = factor(NPI_int)),
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_point(data = subset(test_plot2, pathogen == "SARS-CoV-2-Like" & specific_vaccine_start == 100), 
             aes(x = composite_NPI_spec, y = deaths_spec, fill = factor(NPI_int)), shape = 4, colour = "black", size = 2, pch = 21) +
  geom_point(data = subset(test_plot2, pathogen == "SARS-CoV-2-Like" & specific_vaccine_start == 100), 
             aes(x = composite_NPI_spec, y = deaths_bpsv, fill = factor(NPI_int)), colour = "black", size = 4, pch = 21) +
  theme_bw() +
  labs(x = "NPI Days (Composite Duration+Stringency", 
       y = "Disease Deaths") +
  guides(fill = guide_legend(title = "Scenario"))
x <- ggplot() +
  geom_bar(data = subset(test_plot2, pathogen == "SARS-CoV-2-Like" & specific_vaccine_start == 200), 
           aes(x = factor(NPI_int), y = deaths_saved, fill = factor(NPI_int)), stat = "identity") +
  labs(x = "", y = "Deaths Averted") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.background = element_rect(colour = "black"))
y <- temp_plot + 
  annotation_custom(
    ggplotGrob(x), 
    xmin = 50, xmax = 90, ymin = 3700, ymax = 5000)
cowplot::plot_grid(b, y)

## Plotting out the NPI scenarios
temp_R0 <- R0[2]
temp_detection_time <- detection_time[1]
temp_specific_vaccine_start <- specific_vaccine_start[1]
NPI_df <- NPIs %>%
  filter(R0 == temp_R0,
         detection_time == temp_detection_time,
         specific_vaccine_start == temp_specific_vaccine_start) %>%
  select(NPI_int, NPI_scenario, Rt, tt_Rt) %>%
  rowwise() %>%
  mutate(scenario_info = list(tibble(Rt = Rt, tt_Rt = tt_Rt))) %>%
  select(-Rt, -tt_Rt) %>%
  unnest(cols = c(scenario_info))
NPI_df$scenario <- paste0("Scenario ", NPI_df$NPI_int)
NPI_df <- NPI_df %>%
  group_by(scenario) %>%
  mutate(next_time = lead(tt_Rt),
         next_value = lead(Rt)) %>%
  mutate(next_time = ifelse(is.na(next_time), unique(baseline_scenarios$runtime), next_time),
         next_value = ifelse(is.na(next_value), temp_R0, next_value))
overplot_factor <- 1
example_NPI <- subset(NPI_df, scenario == "Scenario 1")
example_NPI$next_time[length(example_NPI$next_time)] <- 150
a <- ggplot(example_NPI) +
  geom_segment(aes(x = tt_Rt - overplot_factor, xend = next_time + overplot_factor, y = Rt, yend = Rt), size = 2) +
  geom_segment(aes(x = next_time, xend = next_time, y = Rt, yend = next_value), size = 2) +
  geom_hline(aes(yintercept = 1)) +
  geom_hline(aes(yintercept = lockdown_Rt), linetype = "dashed") +
  geom_hline(aes(yintercept = temp_R0 * (1 - minimal_mandate_reduction)), linetype = "dashed") +
  geom_vline(aes(xintercept = temp_detection_time)) +
  geom_vline(aes(xintercept = temp_detection_time + bpsv_start + time_to_coverage_bpsv)) +
  geom_vline(aes(xintercept = temp_detection_time + temp_specific_vaccine_start + time_to_coverage_spec)) +
  theme_bw() +
  geom_segment(aes(x = temp_detection_time, y = temp_R0 + 1, xend = temp_detection_time, yend = temp_R0 + 0.7), arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = temp_detection_time, y = temp_R0 + 1.2, label = "Detection") +
  geom_segment(aes(x = temp_detection_time + bpsv_start + time_to_coverage_bpsv, y = temp_R0 + 1, xend = temp_detection_time + bpsv_start + time_to_coverage_bpsv, yend = temp_R0 + 0.7), arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = temp_detection_time + bpsv_start + time_to_coverage_bpsv, y = temp_R0 + 1.2, label = "BPSV Vacc.\nCompleted", hjust = 0.5) +
  geom_segment(aes(x = temp_detection_time + temp_specific_vaccine_start + time_to_coverage_spec, y = temp_R0 + 1, xend = temp_detection_time + temp_specific_vaccine_start + time_to_coverage_spec, yend = temp_R0 + 0.7), arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = temp_detection_time + temp_specific_vaccine_start + time_to_coverage_spec, y = temp_R0 + 1.2, label = "Spec. Vacc.\nCompleted", hjust = 0.5) +
  geom_segment(aes(x = -20, y = temp_R0 * (1 - minimal_mandate_reduction), xend = -10, yend = temp_R0 * (1 - minimal_mandate_reduction)), arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = -20, y = temp_R0 * (1 - minimal_mandate_reduction), label = "Min.\nMandate", hjust = 1) +
  geom_segment(aes(x = -20, y = lockdown_Rt, xend = -10, yend = 0.9), arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = -20, y = lockdown_Rt, label = "Lockdown", hjust = 1) +
  geom_segment(aes(x = -20, y = temp_R0, xend = -10, yend = temp_R0), arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = -20, y = temp_R0, label = "R0", hjust = 1) +
  theme(plot.margin = margin(2.5, 1, 2.5, 2.5, "cm")) +
  coord_cartesian(clip = 'off', xlim = c(0, 150), ylim = c(0.5, temp_R0 + 0.5)) +
  scale_y_continuous(position = "right") +
  labs(x = "Time (Days)")

b <- ggplot(NPI_df, aes(x = tt_Rt - overplot_factor, colour = scenario)) +
  geom_hline(aes(yintercept = 1), linewidth = 0.2) +
  geom_hline(aes(yintercept = lockdown_Rt), linetype = "dashed", linewidth = 0.2) +
  geom_hline(aes(yintercept = temp_R0 * (1 - minimal_mandate_reduction)), linetype = "dashed", linewidth = 0.2) +
  geom_vline(aes(xintercept = temp_detection_time), linewidth = 0.2) +
  geom_vline(aes(xintercept = temp_detection_time + bpsv_start + time_to_coverage_bpsv), linewidth = 0.2) +
  geom_vline(aes(xintercept = temp_detection_time + temp_specific_vaccine_start + time_to_coverage_spec), linewidth = 0.2) +
  geom_segment(aes(xend = next_time + overplot_factor, y = Rt, yend = Rt), size = 1) +
  geom_segment(aes(x = next_time, xend = next_time, y = Rt, yend = next_value), size = 1) +
  theme_bw() +
  facet_wrap(~scenario) +
  coord_cartesian(xlim = c(0, 150), ylim = c(0.5, temp_R0 + 0.5)) +
  theme(legend.position = "none",
        strip.background =element_rect(fill="#F5F5F5"))

cowplot::plot_grid(a, b, rel_widths = c(1, 1.2)) # 11.3 x 5.7 dimensions is good
