### CHECK - rel_infectiousness_vaccinated - ASSUMPTIONS AROUND THIS??? DEFAULT IS VERY HIGH (50% REDUCTION)
### CHECK - UNLIMITED HEALTHCARE_CAPACITY???
### CHECK - SEEDING CASES???
### CHECK - WEIRD SECOND DOSES COVERAGE BEING LOWER THAN PRIMARY FOR ELDERLY IN BOTH VACCINES MODEL
###         AND REQUIRING AN EXTRA DAY OR TWO TO GET IN THE RIGHT BALLPARK

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
                                       bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                       specific_vaccine_start = 100,                  # specific vaccine distribution start (time after detection time)
                                       specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                       efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                       efficacy_disease_bpsv = seq(0.01, 1, 0.02),    # vaccine efficacy against disease - BPSV
                                       efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                       efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                       dur_R = 365000,                                # duration of infection-induced immunity
                                       dur_V = 365000,                                # duration of vaccine-induced immunity for both vaccines
                                       coverage = 0.8,                                # proportion of the population vaccinated
                                       vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                       min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                       min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)

raw_delay_scenarios <- create_scenarios(R0 = c(2.25, 2.5, 2.75),                       # Basic reproduction number
                                        IFR = c(0.75, 1, 1.25),                        # IFR
                                        Tg = 7,                                        # Tg
                                        detection_time = 14,                           # detection time
                                        bpsv_start = 14,                               # BPSV distribution start (time after detection time)
                                        bpsv_protection_delay = 1:14,                  # delay between receipt of BPSV dose and protection
                                        specific_vaccine_start = 100,                  # specific vaccine distribution start (time after detection time)
                                        specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                        efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                        efficacy_disease_bpsv = 0.75,                  # vaccine efficacy against disease - BPSV
                                        efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                        efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                        dur_R = 365000,                                # duration of infection-induced immunity
                                        dur_V = 365000,                                # duration of vaccine-induced immunity for both vaccines
                                        coverage = 0.8,                                # proportion of the population vaccinated
                                        vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                        min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                        min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)

raw_delay_bpsv_scenarios <- create_scenarios(R0 = c(2.25, 2.5, 2.75),                       # Basic reproduction number
                                             IFR = c(0.75, 1, 1.25),                        # IFR
                                             Tg = 7,                                        # Tg
                                             detection_time = 14,                           # detection time
                                             bpsv_start = 14,                               # BPSV distribution start (time after detection time)
                                             bpsv_protection_delay = 1:14,                  # delay between receipt of BPSV dose and protection
                                             specific_vaccine_start = 100,                  # specific vaccine distribution start (time after detection time)
                                             specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                             efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                             efficacy_disease_bpsv = seq(0.05, 1, 0.05),    # vaccine efficacy against disease - BPSV
                                             efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                             efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                             dur_R = 365000,                                # duration of infection-induced immunity
                                             dur_V = 365000,                                # duration of vaccine-induced immunity for both vaccines
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

## Testing and checking
test_scenarios <- all_scenarios %>%
  filter(scenario_index == 1)
test_spec <- test_scenarios %>%
  filter(vaccine_scenario == "specific_only")
test_both <- test_scenarios %>%
  filter(vaccine_scenario == "both_vaccines")

specific_only <- run_sars_x(population_size = test_spec$population_size,
                            country = test_spec$country,
                            hosp_bed_capacity = test_spec$hosp_bed_capacity,
                            ICU_bed_capacity = test_spec$ICU_bed_capacity,
                            Rt = test_spec$Rt,
                            tt_Rt = test_spec$tt_Rt,
                            Tg = test_spec$Tg,
                            IFR = test_spec$IFR,
                            vaccine_scenario = test_spec$vaccine_scenario,
                            detection_time = test_spec$detection_time,
                            bpsv_start = test_spec$bpsv_start,
                            bpsv_protection_delay = test_spec$bpsv_protection_delay,
                            specific_vaccine_start = test_spec$specific_vaccine_start,
                            specific_protection_delay = test_spec$specific_protection_delay,
                            efficacy_infection_bpsv = test_spec$efficacy_infection_bpsv,
                            efficacy_disease_bpsv = test_spec$efficacy_disease_bpsv,
                            efficacy_infection_spec = test_spec$efficacy_infection_spec,
                            efficacy_disease_spec = test_spec$efficacy_disease_spec,
                            dur_R = test_spec$dur_R,
                            dur_V = test_spec$dur_V,
                            coverage = test_spec$coverage,
                            vaccination_rate = test_spec$vaccination_rate,
                            min_age_group_index_priority = test_spec$min_age_group_index_priority,
                            min_age_group_index_non_priority = test_spec$min_age_group_index_non_priority,
                            runtime = test_spec$runtime,
                            seeding_cases = test_spec$seeding_cases,
                            output = "full",
                            NPI_int = 0,
                            scenario_index = 0,
                            varied = "")

both_vaccines <- run_sars_x(population_size = test_both$population_size,
                            country = test_both$country,
                            hosp_bed_capacity = test_both$hosp_bed_capacity,
                            ICU_bed_capacity = test_both$ICU_bed_capacity,
                            Rt = test_both$Rt,
                            tt_Rt = test_both$tt_Rt,
                            Tg = test_both$Tg,
                            IFR = test_both$IFR,
                            vaccine_scenario = test_both$vaccine_scenario,
                            detection_time = 7, test_both$detection_time,
                            bpsv_start = test_both$bpsv_start,
                            bpsv_protection_delay = test_both$bpsv_protection_delay,
                            specific_vaccine_start = test_both$specific_vaccine_start,
                            specific_protection_delay = test_both$specific_protection_delay,
                            efficacy_infection_bpsv = 0.01, #test_both$efficacy_infection_bpsv,
                            efficacy_disease_bpsv = 0.01, #test_both$efficacy_disease_bpsv,
                            efficacy_infection_spec = test_both$efficacy_infection_spec,
                            efficacy_disease_spec = test_both$efficacy_disease_spec,
                            dur_R = test_both$dur_R,
                            dur_V = test_both$dur_V,
                            coverage = test_both$coverage,
                            vaccination_rate = test_both$vaccination_rate,
                            min_age_group_index_priority = test_both$min_age_group_index_priority,
                            min_age_group_index_non_priority = test_both$min_age_group_index_non_priority,
                            runtime = test_both$runtime,
                            seeding_cases = test_both$seeding_cases,
                            output = "full",
                            NPI_int = 0,
                            scenario_index = 0,
                            varied = "")

# plot(specific_only$model_arguments$primary_doses, type = "l")
# lines(specific_only$model_arguments$second_doses, col = "red")
# lines(specific_only$model_arguments$booster_doses, col = "blue")
# 
# plot(both_vaccines$model_arguments$primary_doses, type = "l")
# lines(both_vaccines$model_arguments$second_doses, col = "red")
# lines(both_vaccines$model_arguments$booster_doses, col = "blue")

specific_only$summary_metrics$deaths / both_vaccines$summary_metrics$deaths ## this is still not perfect, but it's within a tolerable range (~1% difference - think that's due to combo of ODE solver and something else going on I'm currently missing)

sum(both_vaccines$model_arguments$primary_doses)
sum(both_vaccines$model_arguments$second_doses)

check <- nimue::format(both_vaccines$model_output, compartments = c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_booster_dose"),
                       reduce_age = FALSE) %>%
  filter(t > 1,
         compartment == "vaccinated_first_dose" | compartment == "vaccinated_second_dose" | compartment == "vaccinated_booster_dose") %>%
  group_by(replicate, t)
pop_df <- data.frame(age_group = sort(unique(check$age_group)), population = both_vaccines$model_output$parameters$population)
check <- check %>%
  left_join(pop_df, by = "age_group") %>%
  mutate(prop = value / population)
ggplot(data = subset(check, age_group %in% c("50-55", "80+"))) +
  geom_line(aes(x = t, y = prop, col = compartment)) +
  facet_wrap(~age_group)

ggplot(data = check) + #subset(check, age_group %in% c("50-55", "60-65"))) +
  geom_line(aes(x = t, y = prop, col = compartment)) +
  coord_cartesian(ylim=c(0, 0.02))+
  facet_wrap(~age_group) 

ggplot(data = subset(check, age_group %in% c("50-55", "80+"))) +
  geom_line(aes(x = t, y = value, col = compartment)) +
  facet_wrap(~age_group)

# specific_only$model_output$output
# 
# check <- nimue::format(specific_only$model_output, compartments = c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_booster_dose"),
#                        reduce_age = FALSE) %>%
#   filter(t > 1,
#          compartment == "vaccinated_first_dose" | compartment == "vaccinated_second_dose" | compartment == "vaccinated_booster_dose") %>%
#   group_by(replicate, t)
# pop_df <- data.frame(age_group = sort(unique(check$age_group)), population = specific_only$model_output$parameters$population)
# check <- check %>%
#   left_join(pop_df, by = "age_group") %>%
#   mutate(prop = value / population)
# ggplot() +
#   geom_line(data = check, aes(x = t, y = prop, col = compartment)) +
#   facet_wrap(~age_group)



## Running the model and summarising the output
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = 3) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(all_scenarios[1:20, ], run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = 3)
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
