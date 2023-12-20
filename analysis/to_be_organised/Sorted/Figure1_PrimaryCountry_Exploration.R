# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))

## Calculate rough IHR (similar value to Knock et al which had about 2% for London)
prob_hosp <- squire.page.sarsX:::probs_booster$prob_hosp
arg_pop <- squire::get_population("Argentina")
IHR <- sum(prob_hosp * arg_pop$n / sum(arg_pop$n)) 

## Branching-process based calculation of detection times
num_iterations <- 50
num_hosp <- 1:20
detection_hosp <- round(num_hosp / IHR, digits = 0)
cumulative_window <- 7
delay_hosp <- 10 # check this aligns with squire.page 
R0 <- c(1.5, 2, 2.5, 3, 3.5)
runtime <- c(140, 95, 70, 60, 55)

new_bp <- FALSE
if (new_bp) {
  set.seed(456)
  ## Looping over R0
  for (i in 1:length(R0)) {
    
    ## Looping over number of iterations
    for (j in 1:num_iterations) {
      
      ## Controlling for stochastic fadeout (re-running if it fades out)
      try_again <- 1
      while (try_again == 1) {
        chain_sim_eg <- chain_sim_susc(offspring = "pois", 
                                       mn_offspring = R0[i],
                                       t0 = 0, 
                                       tf = runtime[i],
                                       serial = function(n) {
                                         rgamma(n, shape = 13.4, rate = 2)}, # gamma with mean 6.7 (13.4/2)
                                       pop = 10^9,
                                       initial_immune = 0) 
        if (nrow(chain_sim_eg) > 100) {
          try_again <- 0
        }
      }
      
      ## Calculating daily/rolling 7-day cumulative incidence
      incidence <- chain_sim_eg %>%
        mutate(daily = round(time, digits = 0)) %>%
        group_by(daily) %>%
        summarise(incidence = n()) %>%
        tidyr::complete(daily = min(daily):max(daily), fill = list(incidence = 0)) %>%
        mutate(rolling_incidence = zoo::rollsum(x = incidence, k = 7, na.pad = TRUE, align = "right")) 
      
      ## First time daily hospitalisation incidence goes over detection threshold
      first_incidence_over_threshold <- purrr::map(detection_hosp, ~{
        incidence %>%
          filter(incidence > .x) %>%
          summarise(first_day = min(daily))
      })
      incidence_detection <- unlist(first_incidence_over_threshold) + delay_hosp
      
      ## First time rolling 7-day cumulative hospitalisation incidence goes over detection threshold
      first_cumulative_incidence_over_threshold <- purrr::map(detection_hosp, ~{
        incidence %>%
          filter(rolling_incidence > .x) %>%
          summarise(first_day = min(daily))
      })
      cumulative_incidence_detection <- unlist(first_cumulative_incidence_over_threshold) + delay_hosp
      
      ## Adding these to the array
      bp_df[i, j, , 1] <-  unname(incidence_detection)            # time for daily incidence to eclipse that threshold
      bp_df[i, j, , 2] <-  unname(cumulative_incidence_detection) # time for cumulative 7 day rolling window
      
    }
    print(i)
  }
  
  # Convert to a data frame
  bp_df_long <- reshape2::melt(bp_df)
  bp_df_long$R0 <- R0[(bp_df_long$Var1 - 1) %% length(R0) + 1]
  bp_df_long$iteration <- num_iterations[(bp_df_long$Var2 - 1) %% length(num_iterations) + 1]
  bp_df_long$detection <- detection_hosp[(bp_df_long$Var3 - 1) %% length(detection_hosp) + 1]
  bp_df_long$metric <- c("Daily Incidence", "Cumulative 7-Day Incidence")[(bp_df_long$Var4 - 1) %% 2 + 1]
  bp_df_long <- bp_df_long[, -(1:4)]
  bp_df_long <- bp_df_long 
  saveRDS(object = bp_df_long, file = "outputs/Figure1_bp_detection_times.rds")
} else {
  bp_df_long <- readRDS("outputs/Figure1_bp_detection_times.rds")
}

# Plotting all the results for the branching process detection times
bp_df_long_mean <- bp_df_long %>%
  filter(!is.infinite(value)) %>%
  group_by(R0, detection, metric) %>%
  summarise(mean = mean(value)) 
bp_detection_time_plot <- ggplot(bp_df_long_mean, aes(x = R0, y = mean, col = factor(detection))) + 
  geom_line() +
  facet_grid(.~metric)

## Extracting subset of times to run with
num_hosp <- c(1, 5, 10, 20)
infection_thresholds <- detection_hosp[num_hosp]
bp_df_mean_subset <- bp_df_long_mean %>%
  filter(detection %in% infection_thresholds)

# Generate parameter combinations for model running
raw_primary_country_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3, 3.5),                   # Basic reproduction number
                                                  IFR = 1,                                       # IFR
                                                  population_size = 10^10,                       # population size
                                                  Tg = 6.7,                                      # Tg
                                                  detection_time = 1,                            # detection time - PLACECHOLDER FOR NOW
                                                  bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                                  bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                                  specific_vaccine_start = c(100, 200, 365),     # specific vaccine distribution start (time after detection time)
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

# Linking these scenarios with the R0-specific 
primary_country_scenarios <- expand_grid(raw_primary_country_scenarios,
                                         detection_threshold = unique(bp_df_mean_subset$detection)) %>%
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean) 
# 5 R0 * 3 specific development times * 2 vaccine scenarios * 2 detection scenarios * 4 detection thresholds (* n NPIs)

# NPI Relevant Parameters
lockdown_Rt <- 0.9                   # Rt achieved under lockdown
minimal_mandate_reduction <- 0.25    # Fold-reduction in R0 achieved under minimal mandate restrictions
NPIs_primary_country <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                              NPI_scenarios = c(4, 7, 8), scenarios = primary_country_scenarios)

# Combining NPI scenarios with parameter combos into one overall dataframe for model running 
all_combos_primary_country_scenarios <- primary_country_scenarios %>%
  full_join(NPIs_primary_country, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                         "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")


# Filtering the above to only select R0 and detection time pairs that actually occurred (function above produces all pairwise combos of them)
R0_detection_time_pairs <- bp_df_mean_subset %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time)

final_primary_country_scenarios <- all_combos_primary_country_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time")) %>%
  relocate(c("detection_threshold", "metric"), .after = last_col()) %>%
  rename(detection_threshold_inf = detection_threshold) %>%
  mutate(detection_threshold_hosp = round(detection_threshold_inf * IHR))

## Creating index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vars_for_index <- c(variable_columns(final_primary_country_scenarios), "NPI_int")
final_primary_country_scenarios <- final_primary_country_scenarios %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
cores <- parallel::detectCores() - 2
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(final_primary_country_scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/Figure1_primaryCountry_outputs.rds")
} else {
  model_outputs <- readRDS("outputs/Figure1_primaryCountry_outputs.rds")
}

## Joining back in the detection metrics
detection_df <- final_primary_country_scenarios %>%
  select(scenario_index, all_of(vars_for_index)) %>%
  filter(vaccine_scenario == "specific_only") %>%
  ungroup() %>%
  select(-vaccine_scenario) %>%
  select(R0, scenario_index, specific_vaccine_start, NPI_int, detection_time, detection_threshold_inf, detection_threshold_hosp, metric)

model_outputs2 <- model_outputs %>%
  left_join(detection_df, by = c("R0", "scenario_index", "specific_vaccine_start", "detection_time", "NPI_int")) %>%
  mutate(detection_timing = case_when(detection_threshold_hosp == 1 ~ "Early",
                                      detection_threshold_hosp == 5 ~ "Intermediate",
                                      detection_threshold_hosp == 10 ~ "Late",
                                      detection_threshold_hosp == 20 ~ "Very Late"))

## Selecting which NPIs to include
colour_func <- scales::hue_pal()(max(model_outputs$NPI_int))
NPI_colours <- c("#C64191", "#F0803C", "#0D84A9")
population_size <- unique(model_outputs$population_size)
runtime <- unique(model_outputs$runtime)
NPI_to_include <- c(4, 7, 8) # c(2, 4, 5, 7, 8)

## NPI Plot
initial_NPI_df <- NPIs_primary_country %>%
  filter(R0 == 2.5, specific_vaccine_start == 200, NPI_int %in% NPI_to_include) 

table(initial_NPI_df$detection_time)
NPI_df <- initial_NPI_df %>%
  filter(detection_time == 35) %>%
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
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill="white"))

centralSpec_deathsAverted <- model_outputs2 %>%
  filter(specific_vaccine_start %in% c(100, 200),
         detection_threshold_hosp %in% c(1, 5, 10),
         R0 %in% c(1.5, 2, 2.5, 3, 3.5),
         NPI_int %in% NPI_to_include)

# ggplot(data = subset(centralSpec_deathsAverted, metric == "Cumulative 7-Day Incidence")) +
#   geom_bar(aes(x = factor(detection_timing), y = bpsv_deaths_averted, 
#                fill = interaction(factor(NPI_int), factor(specific_vaccine_start))), 
#            stat = "identity", position = "dodge") +
#   facet_grid(NPI_int~R0)

## R0 vs Detection Time Plot (discrete detection time)
my_labeller <- labeller(
  R0 = function(x) paste0("R0 = ", x),
  specific_vaccine_start = function(x) ifelse(x == 100, "Specific Vaccine\n@ 100 Days", "Specific Vaccine\n@ 200 Days"))

R0_NPI_surv_specStart_plot <- ggplot(data = subset(centralSpec_deathsAverted, metric == "Daily Incidence" &
                                                   R0 %in% c(1.5, 2.5, 3))) +
  geom_bar(aes(x = factor(detection_timing), y = bpsv_deaths_averted * 1000 / population_size, 
               colour = factor(NPI_int), fill = factor(NPI_int)), 
           stat = "identity", position = "dodge") +
  facet_grid(specific_vaccine_start~R0, labeller = my_labeller) +
  labs(x = "Detection Timing", y = "Deaths Averted by BPSV (Per 1000)") + 
  scale_fill_manual(values = NPI_colours) +
  scale_colour_manual(values = NPI_colours) +
  theme_bw() +
  theme(strip.placement = "outside",
        legend.position = "none",
        strip.background = element_rect(fill="white"))

complete_R0_detectionFactor_plot <- cowplot::plot_grid(NPI_plot, R0_NPI_surv_specStart_plot, 
                                                       nrow = 1, rel_widths = c(1, 4), labels = c("C", "D"),
                                                       align = "h", axis = "tb")

## Plotting Epidemic Curves and Timings
population <- squire:::get_population("Argentina")
population <- 10^6 * population$n / sum(population$n)
mixing_matrix <- squire:::get_mixing_matrix("Argentina")

low_R0 <- run_booster(time_period = 365,
                      contact_matrix_set = mixing_matrix,
                      population = population,
                      R0 = 1.25,     
                      tt_R0 = 0, 
                      hosp_bed_capacity = 10^9,                                     
                      ICU_bed_capacity = 10^9,                                       
                      dur_R = 365000000000,                                                        
                      seeding_cases = 50,
                      dur_V = 365000000000,                                              
                      primary_doses = rep(0, 365),  
                      second_doses = rep(0, 365),
                      booster_doses = rep(0, 365))

medium_R0 <- run_booster(time_period = 365,
                         contact_matrix_set = mixing_matrix,
                         population = population,
                         R0 = 1.1,     
                         tt_R0 = 0, 
                         hosp_bed_capacity = 10^9,                                     
                         ICU_bed_capacity = 10^9,                                       
                         dur_R = 365000000000,                                                        
                         seeding_cases = 100,
                         dur_V = 365000000000,                                              
                         primary_doses = rep(0, 365),  
                         second_doses = rep(0, 365),
                         booster_doses = rep(0, 365))

high_R0 <- run_booster(time_period = 365,
                       contact_matrix_set = mixing_matrix,
                       population = population,
                       R0 = 1.75,     
                       tt_R0 = 0, 
                       hosp_bed_capacity = 10^9,                                     
                       ICU_bed_capacity = 10^9,                                       
                       dur_R = 365000000000,                                                        
                       seeding_cases = 100,
                       dur_V = 365000000000,                                              
                       primary_doses = rep(0, 365),  
                       second_doses = rep(0, 365),
                       booster_doses = rep(0, 365))

low_R0_infections <- nimue::format(low_R0, compartments = "S", summaries = "hospitalisations") %>%
  filter(compartment != "S") %>%
  select(-replicate) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(R0 = 1.25)

medium_R0_infections <- nimue::format(medium_R0, compartments = "S", summaries = "hospitalisations") %>%
  filter(compartment != "S") %>%
  select(-replicate) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(R0 = 2)

high_R0_infections <- nimue::format(high_R0, compartments = "S", summaries = "hospitalisations") %>%
  filter(compartment != "S") %>%
  select(-replicate) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(R0 = 3)

thresholds <- c(2, 5, 10, 50)
raw_detection_times <- low_R0_infections %>%
  group_by(R0) %>%
  nest() %>%
  mutate(detection_times = map(data, ~ { ## nest creates a column of lists of data
    sapply(thresholds, function(threshold) {
      first(.x$t[.x$value >= threshold])
      # .x$t[which.min(abs(.x$value - threshold))]
    })
  })) %>%
  unnest_longer(col = detection_times, indices_to = "detection_threshold_index") %>%
  mutate(detection_threshold = thresholds[detection_threshold_index]) %>%
  select(-data)

deploy_time <- 18
bpsv_campaign_dur <- 28
spec_develop_time <- 80

relevant_time <- subset(raw_detection_times, R0 == 1.25 & detection_threshold == 2)$detection_times
time_at_detection <- low_R0_infections$t[low_R0_infections$t == relevant_time]
hosp_at_detection <- low_R0_infections$value[time_at_detection]
time_at_BPSV_start <- low_R0_infections$t[low_R0_infections$t == relevant_time + deploy_time]
hosp_at_BPSV_start <- low_R0_infections$value[low_R0_infections$t == relevant_time + deploy_time]
time_at_BPSV_finish <- low_R0_infections$t[low_R0_infections$t == relevant_time + deploy_time + bpsv_campaign_dur]
hosp_at_BPSV_finish <- low_R0_infections$value[low_R0_infections$t == relevant_time + deploy_time + bpsv_campaign_dur]
time_at_spec_start <- low_R0_infections$t[low_R0_infections$t == relevant_time + spec_develop_time]
hosp_at_spec_start <- low_R0_infections$value[low_R0_infections$t == relevant_time + spec_develop_time]

singleR0_example_plot_withExtra <- ggplot(low_R0_infections, aes(x = t, y = log(value + 1), col = factor(R0))) +
  scale_colour_manual(values = c("#DCDBA8")) +
  # geom_line(data = medium_R0_infections, aes(x = t, y = log(value + 1)), col = "#C9DBBA", linewidth = 3.5, alpha = 0.25) +
  # geom_line(data = high_R0_infections, aes(x = t, y = log(value + 1)), col = "#FAA381", linewidth = 3.5, alpha = 0.25) +
  theme_bw() +
  labs(x = "Time Since Spillover", y = "Daily Incidence") +
  coord_cartesian(xlim = c(0, 300)) +
  geom_line(linewidth = 3.5) +
  
  ## Pathogen Detection Event Indicator Segments
  geom_segment(x = -10, xend = time_at_detection - 5, 
               y = log(hosp_at_detection + 1), yend = log(hosp_at_detection + 1), 
               linewidth = 0.25, col = "black", linetype = "solid",
               arrow = arrow(length = unit(0.015, "npc"), type = "closed")) +
  annotate("text", x = -8, y = log(hosp_at_detection + 1) + 0.145, label = "Pathogen Detection", color = "black", hjust = 0, size = 3.25) +
  geom_point(x = time_at_detection, y = log(hosp_at_detection + 1), pch = 21, fill = "#DCDBA8", col = "black", size = 2) +
  
  ## BPSV Campaign Start Event Indicator Segments
  geom_segment(x = -10, xend = time_at_BPSV_start - 5, 
               y = log(hosp_at_BPSV_start + 1), yend = log(hosp_at_BPSV_start + 1), 
               linewidth = 0.25, col = "black", linetype = "solid",
               arrow = arrow(length = unit(0.015, "npc"), type = "closed")) +
  annotate("text", x = -8, y = log(hosp_at_BPSV_start + 1) + 0.145, label = "BPSV Campaign Start", color = "black", hjust = 0, size = 3.25) +
  geom_point(x = time_at_BPSV_start, y = log(hosp_at_BPSV_start + 1), pch = 21, fill = "#DCDBA8", col = "black", size = 2) +
  
  ## BPSV Campaign Finish Event Indicator Segments
  geom_segment(x = -10, xend = time_at_BPSV_finish - 5, 
               y = log(hosp_at_BPSV_finish + 1), yend = log(hosp_at_BPSV_finish + 1), 
               linewidth = 0.25, col = "black", linetype = "solid",
               arrow = arrow(length = unit(0.015, "npc"), type = "closed")) +
  annotate("text", x = -8, y = log(hosp_at_BPSV_finish + 1) + 0.145, label = "BPSV Campaign Finish", color = "black", hjust = 0, size = 3.25) +
  geom_point(x = time_at_BPSV_finish, y = log(hosp_at_BPSV_finish + 1), pch = 21, fill = "#DCDBA8", col = "black", size = 2) +
  
  ## Specific Vaccine Campaign Start Indicator Segments
  geom_segment(x = -10, xend = time_at_spec_start - 7, 
               y = log(hosp_at_spec_start), yend = log(hosp_at_spec_start), 
               linewidth = 0.25, col = "black", linetype = "solid",
               arrow = arrow(length = unit(0.015, "npc"), type = "closed")) +
  annotate("text", x = -8, y = log(hosp_at_spec_start) + 0.3, label = "Specific Vaccine\nDevelopment Time", color = "black", hjust = 0, size = 3.25) +
  geom_point(x = time_at_spec_start - 2, y = log(hosp_at_spec_start), pch = 21, fill = "#DCDBA8", col = "black", size = 2) +
  
  ## BPSV "Time to Deploy" Period Segments
  geom_segment(x = time_at_detection, xend = time_at_detection + deploy_time, 
               y = log(hosp_at_detection + 0.75), yend = log(hosp_at_detection + 0.75), 
               linewidth = 0.25, col = "black", linetype = "dashed",
               arrow = arrow(length = unit(0.015, "npc"), type = "closed", ends = "both")) +
  geom_segment(x = time_at_detection + deploy_time/4, xend = time_at_detection + 3*deploy_time/4, 
               y = log(hosp_at_detection + 0.5), yend = log(hosp_at_detection + 0.5), 
               linewidth = 0.25, col = "black", linetype = "solid") +
  geom_segment(x = time_at_detection + 7, xend = time_at_detection + 7, 
               y = 0.5, yend = log(hosp_at_detection + 0.5), 
               linewidth = 0.25, col = "black", linetype = "solid") +
  geom_segment(x = time_at_detection + 7, xend = time_at_detection + 42.5, 
               y = 0.5, yend = 0.5, 
               linewidth = 0.25, col = "black", linetype = "solid") +
  annotate("label", x = time_at_detection + 42.5, y = 0.5, label = "Time to Deploy BPSV", color = "black", hjust = 0, size = 3.25) +
  
  ## BPSV "Vaccination Campaign" Period Segments
  geom_segment(x = time_at_detection + deploy_time, xend = time_at_detection + deploy_time + bpsv_campaign_dur, 
               y = log(hosp_at_BPSV_start) + 0.8, yend = log(hosp_at_BPSV_start) + 0.8, 
               linewidth = 0.25, col = "black", linetype = "dashed",
               arrow = arrow(length = unit(0.015, "npc"), type = "closed", ends = "both")) +
  geom_segment(x = time_at_detection + deploy_time + bpsv_campaign_dur/4, xend = time_at_detection + deploy_time + bpsv_campaign_dur - bpsv_campaign_dur/4, 
               y = log(hosp_at_BPSV_start) + 0.7, yend = log(hosp_at_BPSV_start) + 0.7, 
               linewidth = 0.25, col = "black", linetype = "solid") +
  geom_segment(x = time_at_detection + deploy_time + bpsv_campaign_dur/4 + bpsv_campaign_dur/4, xend = time_at_detection + deploy_time + bpsv_campaign_dur/4 + bpsv_campaign_dur/4, 
               y = log(hosp_at_BPSV_start) + 0.2, yend = log(hosp_at_BPSV_start) + 0.7, 
               linewidth = 0.25, col = "black", linetype = "solid") +
  geom_segment(x = time_at_detection + deploy_time + bpsv_campaign_dur/4 + bpsv_campaign_dur/4, xend = time_at_detection + 42.5, 
               y = log(hosp_at_BPSV_start) + 0.2, yend = log(hosp_at_BPSV_start) + 0.2, 
               linewidth = 0.25, col = "black", linetype = "solid") +
  annotate("label", x = time_at_detection + 43.5, y = log(hosp_at_BPSV_start) + 0.2, label = "Time to Complete\nBPSV Vaccination", color = "black", hjust = 0, size = 3.25) +
  
  ## Specific Vaccine Development Segments
  geom_segment(x = time_at_detection, xend = time_at_detection + spec_develop_time - 2, 
               y = log(hosp_at_spec_start) - 0.15, yend = log(hosp_at_spec_start) - 0.15, 
               linewidth = 0.25, col = "black", linetype = "dashed",
               arrow = arrow(length = unit(0.015, "npc"), type = "closed", ends = "both")) +
  geom_segment(x = time_at_detection + spec_develop_time/2 - spec_develop_time/10, xend = time_at_detection + spec_develop_time/2 + spec_develop_time/10, 
               y = log(hosp_at_spec_start) - 0.25, yend = log(hosp_at_spec_start) - 0.25, 
               linewidth = 0.25, col = "black", linetype = "solid") +
  geom_segment(x = time_at_detection + spec_develop_time/2, xend = time_at_detection + spec_develop_time/2, 
               y = 2.5, yend = log(hosp_at_spec_start) - 0.25, 
               linewidth = 0.25, col = "black", linetype = "solid") +
  geom_segment(x = time_at_detection + spec_develop_time/2, xend = time_at_detection + 92.5, 
               y = 2.5, yend = 2.5, 
               linewidth = 0.25, col = "black", linetype = "solid") +
  annotate("label", x = time_at_detection + 92.5, y = 2.4, label = "Time to Develop\nSpecific Vaccine", color = "black", hjust = 0, size = 3.25) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") #+
  # geom_rect(xmin = 210, xmax = 305, ymin = 4.4, ymax = 5.4, fill = "white", col = "black", linewidth = 0.25) +
  # geom_rect(xmin = 215, xmax = 225, ymin = 4.5, ymax = 4.7, fill = "#C9DBBA", col = "#C9DBBA") +
  # annotate("text", x = 230, y = 4.6, label = "R0=1.5", color = "black", hjust = 0) +
  # geom_rect(xmin = 215, xmax = 225, ymin = 4.8, ymax = 5, fill = "#DCDBA8", col = "#DCDBA8") +
  # annotate("text", x = 230, y = 4.9, label = "R0=2.5", color = "black", hjust = 0) +
  # geom_rect(xmin = 215, xmax = 225, ymin = 5.1, ymax = 5.3, fill = "#FAA381", col = "#FAA381") +
  # annotate("text", x = 230, y = 5.2, label = "R0=3.5", color = "black", hjust = 0) 

## Generating Heatmap of Specific Vaccine Development Time vs Detection Time
num_hosp <- 1:20
infection_thresholds <- detection_hosp[num_hosp]
bp_df_mean_subset <- bp_df_long_mean %>%
  filter(detection %in% infection_thresholds)

# Generate parameter combinations for model running
raw_specStart_detectTime_scenarios <- create_scenarios(R0 = c(2.5),                                  # Basic reproduction number
                                                       IFR = 1,                                       # IFR
                                                       population_size = 10^10,                       # population size
                                                       Tg = 6.7,                                      # Tg
                                                       detection_time = 1,                            # detection time - PLACECHOLDER FOR NOW
                                                       bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                                       bpsv_protection_delay = 7,                     # delay between receipt of BPSV dose and protection
                                                       specific_vaccine_start = seq(100, 200, 5),     # specific vaccine distribution start (time after detection time)
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

# Linking these scenarios with the R0-specific 
specStart_detectTime_scenarios <- expand_grid(raw_specStart_detectTime_scenarios,
                                         detection_threshold = unique(bp_df_mean_subset$detection)) %>%
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean) 

# NPI Relevant Parameters
lockdown_Rt <- 0.9                   # Rt achieved under lockdown
minimal_mandate_reduction <- 0.25    # Fold-reduction in R0 achieved under minimal mandate restrictions
NPIs_specStart_detectTime <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                                   NPI_scenarios = 8, scenarios = specStart_detectTime_scenarios)

# Combining NPI scenarios with parameter combos into one overall dataframe for model running 
all_combos_specStart_detectTime_scenarios <- specStart_detectTime_scenarios %>%
  full_join(NPIs_specStart_detectTime, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                              "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")


# Filtering the above to only select R0 and detection time pairs that actually occurred (function above produces all pairwise combos of them)
R0_detection_time_pairs <- bp_df_mean_subset %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time)

final_specStart_detectTime_scenarios <- all_combos_specStart_detectTime_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time")) %>%
  relocate(c("detection_threshold", "metric"), .after = last_col()) %>%
  rename(detection_threshold_inf = detection_threshold) %>%
  mutate(detection_threshold_hosp = round(detection_threshold_inf * IHR))

## Creating index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vars_for_index <- c(variable_columns(final_specStart_detectTime_scenarios), "NPI_int")
final_specStart_detectTime_scenarios <- final_specStart_detectTime_scenarios %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
cores <- parallel::detectCores() - 2
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(final_specStart_detectTime_scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  specStart_detectTime_model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(specStart_detectTime_model_outputs, "outputs/Figure1_specStart_detectTime_sens.rds")
} else {
  specStart_detectTime_model_outputs <- readRDS("outputs/Figure1_specStart_detectTime_sens.rds")
}

## Joining back in the detection metrics
detection_df <- final_specStart_detectTime_scenarios %>%
  select(scenario_index, all_of(vars_for_index), R0) %>%
  filter(vaccine_scenario == "specific_only") %>%
  ungroup() %>%
  select(-vaccine_scenario) %>%
  select(R0, scenario_index, specific_vaccine_start, NPI_int, detection_time, detection_threshold_inf, detection_threshold_hosp, metric)

model_outputs2 <- specStart_detectTime_model_outputs %>%
  left_join(detection_df, by = c("R0", "scenario_index", "specific_vaccine_start", "detection_time", "NPI_int")) %>%
  mutate(detection_timing = case_when(detection_threshold_hosp == 1 ~ "Early",
                                      detection_threshold_hosp == 5 ~ "Intermediate",
                                      detection_threshold_hosp == 10 ~ "Late",
                                      detection_threshold_hosp == 20 ~ "Very Late"))

specStart_detectTime_heatmap_plot <- ggplot(subset(model_outputs2, NPI_int == 8 & R0 == 2.5 & specific_vaccine_start <= 200), 
                            aes(x = specific_vaccine_start, y = detection_threshold_hosp, fill = bpsv_deaths_averted  * 1000 / population_size)) + 
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "#7C9EB2", high = "#37617A", midpoint = 3) +
  theme_bw() +
  labs(x = "Time to Develop Specific Vaccine (Days)", y = "Detection Threshold\n(# Hospitalizations)") +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = c(0.77, 0.75),
        legend.key.size = unit(0.02, 'npc'),
        legend.text = element_text(size = 8)) +
  scale_y_continuous(breaks = 1:20) +
  scale_x_continuous(breaks = seq(100, 200, 10)) +
  coord_cartesian(x = c(102, 198),y = c(1.2, 19.5)) +
  guides(fill = guide_legend(title = "Deaths Averted\n(Per 1000)"))

top_row <- cowplot::plot_grid(singleR0_example_plot_withExtra, specStart_detectTime_heatmap_plot, nrow = 1, rel_widths = 2, labels = c("A", "B"), align = "h", axis = "tb")
final_Figure1 <- cowplot::plot_grid(top_row, complete_R0_detectionFactor_plot, nrow = 2, rel_heights = c(1, 1.5))
ggsave(plot = final_Figure1,
       filename = "figures/Figure1_PrimaryCountry_SpecStart_vs_DetectionTime_Draft.pdf", width = 7.5, height = 8)
