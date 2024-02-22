# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/helper_functions.R"))
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/runs_sars_x_secondary_country.R"))
source(here::here("functions/branching_process.R"))
source(here::here("functions/branching_process2.R"))

default <- define_default_params()

## Figure for Primary Country and Showing How Surveillance Sensitivity Influences BPSV Impact

## Calculate rough IHR (similar value to Knock et al which had about 2% for London)
prob_hosp <- squire.page.sarsX:::probs_booster$prob_hosp
arg_pop <- squire::get_population("Argentina")
IHR <- sum(prob_hosp * arg_pop$n / sum(arg_pop$n)) 

## Branching-process based calculation of detection times
num_iterations <- 25
num_hosp <- c(1, 2, 3, 4, 5, seq(10, 50, 10)) #1:20
detection_hosp <- round(num_hosp / IHR, digits = 0)
cumulative_window <- 7
delay_hosp <- 9.1  
R0 <- 2.5 # c(1.5, 2.5, 3.5)
runtime <- 135 # c(140, 70, 55)
bp_df <- matrix(data = NA, nrow = num_iterations, ncol = length(num_hosp))

new_bp <- FALSE
if (new_bp) {
  
  set.seed(456)
  
  ## Looping over R0
  for (i in 1:num_iterations) {

    ## Controlling for stochastic fadeout (re-running if it fades out)
    try_again <- 1
    while (try_again == 1) {
      chain_sim_eg <- chain_sim_susc_ring_vacc2(offspring = "pois",
                                                mn_offspring = R0,
                                                generation_time = function(n) {
                                                  rgamma(n, shape = 13.4, rate = 2)}, # gamma with mean 6.7 (13.4/2),
                                                t0 = 0, tf = 250, 
                                                pop = 10^9, check_final_size = 25000,
                                                initial_immune = 0,
                                                seeding_cases = 1, 
                                                prop_asymptomatic = 0,
                                                infection_to_onset = function(n) {
                                                  rgamma(n, shape = 0.0001, rate = 1)},
                                                vaccine_start = 1000, vaccine_coverage = 0,
                                                vaccine_efficacy_infection = 0,
                                                vaccine_efficacy_transmission = 0,
                                                vaccine_logistical_delay = 0,
                                                vaccine_protection_delay = 1000)
      if (sum(!is.na(chain_sim_eg$time_infection)) > 100) {
        try_again <- 0
      }
    }
      
    ## Calculating daily/rolling 7-day cumulative incidence
    ## Note the decrease in rolling incidence because check_final_size caps it
    ##  which means some folks who have been infected, but who we haven't generated offspring for
    ##  are present in the dataset
    incidence <- chain_sim_eg %>%
      filter(!is.na(time_infection)) %>%
      mutate(daily = round(time_infection, digits = 0)) %>%
      group_by(daily) %>%
      summarise(incidence = n()) %>%
      tidyr::complete(daily = min(daily):max(daily), fill = list(incidence = 0)) %>%
      mutate(rolling_incidence = zoo::rollsum(x = incidence, k = 7, na.pad = TRUE, align = "right"))
    
    ## First time rolling 7-day cumulative hospitalisation incidence goes over detection threshold
    # first_cumulative_incidence_over_threshold <- purrr::map(detection_hosp, ~{
    #   incidence %>%
    #     filter(rolling_incidence > .x) %>%
    #     summarise(first_day = min(daily))
    # })
    # cumulative_incidence_detection <- unlist(first_cumulative_incidence_over_threshold) + delay_hosp
    
    first_incidence_over_threshold <- purrr::map(detection_hosp, ~{
      incidence %>%
        filter(incidence > .x) %>%
        summarise(first_day = min(daily))
    })
    incidence_detection <- unlist(first_incidence_over_threshold) + delay_hosp
    
    ## Adding these to the array
    bp_df[i, ] <-  unname(incidence_detection) # time for cumulative 7 day rolling window
    
    print(i)
  }

  # Convert to a data frame
  bp_df_long <- reshape2::melt(bp_df) # R0 x iterations x num_hosp
  bp_df_long$iteration <- bp_df_long$Var1
  bp_df_long$detection <- detection_hosp[(bp_df_long$Var2 - 1) %% length(detection_hosp) + 1]
  bp_df_long$metric <- "Daily Incidence"  # "Cumulative 7-Day Incidence"
  bp_df_long <- bp_df_long[, -(1:2)]
  saveRDS(object = bp_df_long, file = "outputs/Figure6_Surveillance_Exploration/NEW_Fig6_surveillanceSensitivity.rds")
} else {
  bp_df_long <- readRDS("outputs/Figure6_Surveillance_Exploration/NEW_Fig6_surveillanceSensitivity.rds")
}

# Plotting all the results for the branching process detection times
bp_df_long_mean <- bp_df_long %>%
  filter(!is.infinite(value)) %>%
  group_by(detection, metric) %>%
  summarise(mean = mean(value)) %>%
  mutate(detection_hosp = round(detection * IHR))

inset_detection_plot <- ggplot(subset(bp_df_long_mean, metric == "Daily Incidence")) + 
  geom_line(aes(x = mean, y = detection_hosp, col = factor(R0))) +
  geom_point(data = subset(bp_df_long_mean, metric == "Daily Incidence" & detection_hosp == 3),
             aes(x = mean, y = detection_hosp, fill = factor(R0)),
             pch = 21, size = 4) +
  guides(col = guide_legend(title = "R0"),
         fill = "none") +
  scale_x_continuous(position = "top") +
  theme_classic() +
  labs(y = "Hosp. Detection Threshold", x = "Time to Detection (Days)") +
  theme(strip.background = element_rect(fill="#F5F5F5"),
        legend.position = "none")

## Running
bp_df_long_mean$R0 <- R0
detection_times <- bp_df_long_mean %>%
  filter(metric == "Daily Incidence") %>%
  rename(detection_time = mean) %>%
  select(R0, detection_time, detection_hosp)
detection_times$detection_time  <- round(detection_times$detection_time)

raw_baseline_scenarios <- create_scenarios(R0 = R0, specific_vaccine_start = c(100, 250)) %>%
  select(-detection_time)
baseline_scenarios <- expand_grid(raw_baseline_scenarios,
                                  detection_hosp = unique(detection_times$detection_hosp)) %>%
  left_join(detection_times, by = c("R0", "detection_hosp")) 

NPIs <- default_NPI_scenarios(lockdown_Rt = default$lockdown_Rt, minimal_mandate_reduction = default$minimal_mandate_reduction, 
                              NPI_scenarios = c(4, 7, 8), scenarios = baseline_scenarios)
scenarios_NPIs <- baseline_scenarios %>%
  full_join(NPIs, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence timing of NPI scenarios
                         "specific_vaccine_start", "vaccination_rate_bpsv", "vaccination_rate_spec",
                         "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all")

detection_times2 <- detection_times %>%
  select(R0, detection_time, detection_hosp)
baseline_scenarios_reduced <- scenarios_NPIs %>%
  semi_join(detection_times2, by = c("R0", "detection_hosp", "detection_time"))
# R0 * NPI * spec start * 2 vaccination scenarios * 25 detection thresholds 3 * 3 * 3 * 2 * 25
3 * 1 * 2 * 2 * 25

## Creating index for output
vars_for_index <- c(variable_columns(baseline_scenarios_reduced), "NPI_int")
scenarios <- baseline_scenarios_reduced %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = 6) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = 2)
  saveRDS(model_outputs, "outputs/Figure6_Surveillance_Exploration/NEW_Figure6_PrimaryCountry_Surveillance.rds")
} else {
  model_outputs <- readRDS("outputs/Figure6_Surveillance_Exploration/NEW_Figure6_PrimaryCountry_Surveillance.rds")
}

## Plotting the sensitivity analyses for detection threshold
scenarios2 <- scenarios %>%
  select(scenario_index, detection_hosp) %>%
  group_by(scenario_index) %>%
  slice_head(n = 1) %>%
  as.data.frame() %>%
  ungroup() %>%
  select(scenario_index, detection_hosp)

model_outputs2 <- model_outputs %>%
  filter(specific_vaccine_start %in% c(100, 250)) %>%
  filter(R0 %in% c(2.5)) %>%
  filter(NPI_int %in% c(4, 7, 8)) %>%
  left_join(scenarios2, by = "scenario_index")

model_outputs2$hosp_check <- ifelse(model_outputs2$detection_hosp == 3, "Yes", "No")
detection_sensitivity_plot <- ggplot(subset(model_outputs2, NPI_int == 7)) +
  geom_bar(aes(x = factor(detection_hosp), y = 1000 * bpsv_deaths_averted / unique(baseline_scenarios$population_size), 
               fill = interaction(factor(specific_vaccine_start), factor(hosp_check))), 
           stat = "identity", size = 1) +
  theme_bw() +
  facet_wrap(specific_vaccine_start ~ ., nrow = 2,
             labeller = as_labeller(c(`100` = "Specific Vaccine in 100 Days", 
                                      `250` = "Specific Vaccine in 250 Days"))) +
  labs(x = "Hospitalisation Detection Threshold", y = "Deaths Averted by BPSV\n (Per 1,000 Population)") +
  scale_fill_manual(values = c("grey", "grey", "#4F96F9", "#D387AB")) +
  guides(col = guide_legend(title = "R0"),
         fill = "none") +
  scale_y_continuous(position = "right") +
  theme(strip.background = element_rect(fill="#F5F5F5"),
        legend.position = "right")

detection_sensitivity_plot100 <- ggplot(subset(model_outputs2, NPI_int == 7 & specific_vaccine_start == 100)) +
  geom_bar(aes(x = factor(detection_hosp), y = 1000 * bpsv_deaths_averted / unique(baseline_scenarios$population_size), 
               fill = interaction(factor(specific_vaccine_start), factor(hosp_check))), 
           stat = "identity", size = 1) +
  theme_bw() +
  facet_wrap(specific_vaccine_start ~ ., nrow = 2,
             labeller = as_labeller(c(`100` = "Specific Vaccine in 100 Days", 
                                      `250` = "Specific Vaccine in 250 Days"))) +
  labs(x = "Hospitalisation Detection Threshold", y = "Deaths Averted by BPSV\n (Per 1,000 Population)") +
  scale_fill_manual(values = c("grey", "#4F96F9")) +
  guides(col = guide_legend(title = "R0"),
         fill = "none") +
  scale_y_continuous(position = "right") +
  theme(strip.background = element_rect(fill="#F5F5F5"),
        legend.position = "right")

detection_sensitivity_plot250 <- ggplot(subset(model_outputs2, NPI_int == 7 & specific_vaccine_start == 250)) +
  geom_bar(aes(x = factor(detection_hosp), y = 1000 * bpsv_deaths_averted / unique(baseline_scenarios$population_size), 
               fill = interaction(factor(specific_vaccine_start), factor(hosp_check))), 
           stat = "identity", size = 1) +
  theme_bw() +
  facet_wrap(specific_vaccine_start ~ ., nrow = 2,
             labeller = as_labeller(c(`100` = "Specific Vaccine in 100 Days", 
                                      `250` = "Specific Vaccine in 250 Days"))) +
  labs(x = "Hospitalisation Detection Threshold", y = "Deaths Averted by BPSV\n (Per 1,000 Population)") +
  scale_fill_manual(values = c("grey", "#D387AB")) +
  guides(col = guide_legend(title = "R0"),
         fill = "none") +
  scale_y_continuous(position = "right") +
  theme(strip.background = element_rect(fill="#F5F5F5"),
        legend.position = "right")


### Running individual simulations
R0_check <- 2.5
detection_hosp_check <- 3
cutoff_time <- 300
NPI_int_check <- 7

### No vaccines
nothing <- baseline_scenarios_reduced %>%
  filter(R0 == R0_check, specific_vaccine_start == 250, NPI_int == NPI_int_check, detection_hosp == detection_hosp_check, vaccine_scenario == "specific_only")
nothing$detection_hosp == 5

no_vaccine <- run_sars_x(population_size = nothing$population_size,
                         country = nothing$country,
                         hosp_bed_capacity = nothing$hosp_bed_capacity,
                         ICU_bed_capacity = nothing$ICU_bed_capacity,
                         Rt = nothing$Rt, 
                         tt_Rt = nothing$tt_Rt,
                         Tg = nothing$Tg,
                         IFR = nothing$IFR,
                         vaccine_scenario = "specific_only", 
                         detection_time = 5,                      
                         bpsv_start = 100000,                          
                         bpsv_protection_delay = 10000,
                         specific_protection_delay = 10, 
                         specific_vaccine_start = 600,                      
                         efficacy_infection_bpsv = 0,           
                         efficacy_disease_bpsv = 0,              
                         efficacy_infection_spec = 0,           
                         efficacy_disease_spec = 0,              
                         dur_R = 1000 * 365,                       
                         dur_bpsv = 1000 * 365,                    
                         dur_spec = 1000 * 365,                    
                         coverage_bpsv = 0.01,                     
                         coverage_spec = 0.01,                     
                         vaccination_rate_bpsv = 0.001,             
                         vaccination_rate_spec = 0.001,             
                         min_age_group_index_priority = 13,        
                         min_age_group_index_non_priority = 4,     
                         runtime = nothing$runtime,
                         seeding_cases = nothing$seeding_cases,
                         output = "full")

nothing_output <- nimue::format(no_vaccine$model_output, 
                                compartments = c("D"),
                            reduce_age = TRUE) %>%
  filter(t > 1, compartment == "deaths") %>%
  group_by(replicate, t) %>%
  filter(t < cutoff_time) %>%
  mutate(value = 1000 * value / nothing$population_size)
nothing_output$scenario <- "anothing"

### Specific vaccine
spec_only <- baseline_scenarios_reduced %>%
  filter(R0 == R0_check, specific_vaccine_start == 250, NPI_int == NPI_int_check, detection_hosp == detection_hosp_check, vaccine_scenario == "specific_only")
spec_vaccine <- run_sars_x(population_size = spec_only$population_size,
                           country = spec_only$country,
                           hosp_bed_capacity = spec_only$hosp_bed_capacity,
                           ICU_bed_capacity = spec_only$ICU_bed_capacity,
                           Rt = spec_only$Rt, 
                           tt_Rt = spec_only$tt_Rt,
                           Tg = spec_only$Tg,
                           IFR = spec_only$IFR,
                           vaccine_scenario = "specific_only", 
                           detection_time = spec_only$detection_time,                      
                           bpsv_start = 100000,                          
                           bpsv_protection_delay = 10000,
                           specific_protection_delay = spec_only$specific_protection_delay, 
                           specific_vaccine_start = spec_only$specific_vaccine_start,                      
                           efficacy_infection_bpsv = 0,           
                           efficacy_disease_bpsv = 0,              
                           efficacy_infection_spec = spec_only$efficacy_infection_spec,           
                           efficacy_disease_spec = spec_only$efficacy_disease_spec,              
                           dur_R = spec_only$dur_R,                       
                           dur_bpsv = spec_only$dur_bpsv,                    
                           dur_spec = spec_only$dur_spec,                    
                           coverage_bpsv = 0.01,                     
                           coverage_spec = spec_only$coverage_spec,                     
                           vaccination_rate_bpsv = 0.001,             
                           vaccination_rate_spec = spec_only$vaccination_rate_spec,             
                           min_age_group_index_priority = spec_only$min_age_group_index_priority,        
                           min_age_group_index_non_priority = spec_only$min_age_group_index_non_priority,     
                           runtime = spec_only$runtime,
                           seeding_cases = spec_only$seeding_cases,
                           output = "full")

spec_only_output <- nimue::format(spec_vaccine$model_output, 
                                  compartments = c("D"),
                                  reduce_age = TRUE) %>%
  filter(t > 1, compartment == "deaths") %>%
  group_by(replicate, t) %>%
  filter(t < cutoff_time) %>%
  mutate(value = 1000 * value / nothing$population_size)
spec_only_output$scenario <- "bspecific_only"

### BPSV
spec_bpsv <- baseline_scenarios_reduced %>%
  filter(R0 == R0_check, specific_vaccine_start == 250, NPI_int == NPI_int_check, detection_hosp == detection_hosp_check, vaccine_scenario == "both_vaccines")
both_vaccines <- run_sars_x(population_size = spec_bpsv$population_size,
                           country = spec_bpsv$country,
                           hosp_bed_capacity = spec_bpsv$hosp_bed_capacity,
                           ICU_bed_capacity = spec_bpsv$ICU_bed_capacity,
                           Rt = spec_bpsv$Rt, 
                           tt_Rt = spec_bpsv$tt_Rt,
                           Tg = spec_bpsv$Tg,
                           IFR = spec_bpsv$IFR,
                           vaccine_scenario = "both_vaccines", 
                           detection_time = spec_bpsv$detection_time,                      
                           bpsv_start = spec_bpsv$bpsv_start,                          
                           bpsv_protection_delay = spec_bpsv$bpsv_protection_delay,
                           specific_protection_delay = spec_bpsv$specific_protection_delay, 
                           specific_vaccine_start = spec_bpsv$specific_vaccine_start,                      
                           efficacy_infection_bpsv = spec_bpsv$efficacy_infection_bpsv,           
                           efficacy_disease_bpsv = spec_bpsv$efficacy_disease_bpsv,              
                           efficacy_infection_spec = spec_bpsv$efficacy_infection_spec,           
                           efficacy_disease_spec = spec_bpsv$efficacy_disease_spec,              
                           dur_R = spec_bpsv$dur_R,                       
                           dur_bpsv = spec_bpsv$dur_bpsv,                    
                           dur_spec = spec_bpsv$dur_spec,                    
                           coverage_bpsv = spec_bpsv$coverage_bpsv,                     
                           coverage_spec = spec_bpsv$coverage_spec,                     
                           vaccination_rate_bpsv = spec_bpsv$vaccination_rate_bpsv,             
                           vaccination_rate_spec = spec_bpsv$vaccination_rate_spec,             
                           min_age_group_index_priority = spec_bpsv$min_age_group_index_priority,        
                           min_age_group_index_non_priority = spec_bpsv$min_age_group_index_non_priority,     
                           runtime = spec_bpsv$runtime,
                           seeding_cases = spec_bpsv$seeding_cases,
                           output = "full")

both_output <- nimue::format(both_vaccines$model_output, 
                             compartments = c("D"),
                             reduce_age = TRUE) %>%
  filter(t > 1, compartment == "deaths") %>%
  group_by(replicate, t) %>%
  filter(t < cutoff_time) %>%
  mutate(value = 1000 * value / nothing$population_size)
both_output$scenario <- "cboth"

overall <- rbind(nothing_output, spec_only_output, both_output)

panel_apt1 <- ggplot(overall, aes(x = t, y = value, col = scenario)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("#948D9B", "#F2C0D7", "#D387AB"),
                      labels = c("No Vaccines", "Disease-Specific Only", "Disease-Specific + BPSV"),
                      name = "250 Days to\nDisease-Specific Vaccine") +
  theme_bw() +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(40, 300)) +
  labs(x = "Time (Days)", y = "Daily Deaths per 1,000 Population")

overall_deaths <- overall %>%
  group_by(scenario) %>%
  summarise(total = sum(value))

panel_apt2 <- ggplot(overall_deaths, aes(x = scenario, y = total, fill = scenario)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#948D9B", "#F2C0D7", "#D387AB"),
                    labels = c("No Vaccines", "Disease-Specific Only", "Disease-Specific + BPSV"),
                    name = "250 Days to\nDisease-Specific Vaccine") +
  scale_x_discrete(labels = c("No\nVaccines", "Disease-\nSpecific\nOnly", "Disease-\nSpecific\n+ BPSV")) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Total Deaths\nper 1,000 Population")

first_half <- cowplot::plot_grid(panel_apt1, panel_apt2, nrow = 1, 
                                 rel_widths = c(2.5, 1), axis = "b", align = "h",
                                 labels = c("A"))

first_half <- panel_apt1 + 
  annotation_custom(
    ggplotGrob(panel_apt2), 
    xmin = 230, xmax = 310, ymin = max(overall$value) / 5, ymax = max(overall$value)) +
  annotation_custom(
    ggplotGrob(inset_detection_plot), 
    xmin = 45, xmax = 140, 
    ymin = max(overall$value) / 2.5, ymax = max(overall$value))

first_half_pt2 <- cowplot::plot_grid(first_half, detection_sensitivity_plot250, nrow = 1, 
                                     align = "h", axis = "bt", 
                                     rel_widths = c(2.5, 1))
# 11.82 x 3.95 dimensions

### Running individual simulations for 100 days

### No vaccines
nothing <- baseline_scenarios_reduced %>%
  filter(R0 == R0_check, specific_vaccine_start == 100, NPI_int == NPI_int_check, detection_hosp == detection_hosp_check, vaccine_scenario == "specific_only")
no_vaccine <- run_sars_x(population_size = nothing$population_size,
                         country = nothing$country,
                         hosp_bed_capacity = nothing$hosp_bed_capacity,
                         ICU_bed_capacity = nothing$ICU_bed_capacity,
                         Rt = nothing$Rt, 
                         tt_Rt = nothing$tt_Rt,
                         Tg = nothing$Tg,
                         IFR = nothing$IFR,
                         vaccine_scenario = "specific_only", 
                         detection_time = 5,                      
                         bpsv_start = 100000,                          
                         bpsv_protection_delay = 10000,
                         specific_protection_delay = 10, 
                         specific_vaccine_start = 600,                      
                         efficacy_infection_bpsv = 0,           
                         efficacy_disease_bpsv = 0,              
                         efficacy_infection_spec = 0,           
                         efficacy_disease_spec = 0,              
                         dur_R = nothing$dur_R,                       
                         dur_bpsv = nothing$dur_bpsv,                    
                         dur_spec = nothing$dur_spec,                    
                         coverage_bpsv = 0.01,                     
                         coverage_spec = 0.01,                     
                         vaccination_rate_bpsv = 0.001,             
                         vaccination_rate_spec = 0.001,             
                         min_age_group_index_priority = nothing$min_age_group_index_priority,        
                         min_age_group_index_non_priority = nothing$min_age_group_index_non_priority,     
                         runtime = nothing$runtime,
                         seeding_cases = nothing$seeding_cases,
                         output = "full")

nothing_output <- nimue::format(no_vaccine$model_output, 
                                compartments = c("D"),
                                reduce_age = TRUE) %>%
  filter(t > 1, compartment == "deaths") %>%
  group_by(replicate, t) %>%
  filter(t < cutoff_time) %>%
  mutate(value = 1000 * value / nothing$population_size)
nothing_output$scenario <- "anothing"

### Specific vaccine
spec_only <- baseline_scenarios_reduced %>%
  filter(R0 == R0_check, specific_vaccine_start == 100, NPI_int == NPI_int_check, detection_hosp == detection_hosp_check, vaccine_scenario == "specific_only")
spec_vaccine <- run_sars_x(population_size = spec_only$population_size,
                           country = spec_only$country,
                           hosp_bed_capacity = spec_only$hosp_bed_capacity,
                           ICU_bed_capacity = spec_only$ICU_bed_capacity,
                           Rt = spec_only$Rt, 
                           tt_Rt = spec_only$tt_Rt,
                           Tg = spec_only$Tg,
                           IFR = spec_only$IFR,
                           vaccine_scenario = "specific_only", 
                           detection_time = spec_only$detection_time,                      
                           bpsv_start = 100000,                          
                           bpsv_protection_delay = 10000,
                           specific_protection_delay = spec_only$specific_protection_delay, 
                           specific_vaccine_start = spec_only$specific_vaccine_start,                      
                           efficacy_infection_bpsv = 0,           
                           efficacy_disease_bpsv = 0,              
                           efficacy_infection_spec = spec_only$efficacy_infection_spec,           
                           efficacy_disease_spec = spec_only$efficacy_disease_spec,              
                           dur_R = spec_only$dur_R,                       
                           dur_bpsv = spec_only$dur_bpsv,                    
                           dur_spec = spec_only$dur_spec,                    
                           coverage_bpsv = 0.01,                     
                           coverage_spec = spec_only$coverage_spec,                     
                           vaccination_rate_bpsv = 0.001,             
                           vaccination_rate_spec = spec_only$vaccination_rate_spec,             
                           min_age_group_index_priority = spec_only$min_age_group_index_priority,        
                           min_age_group_index_non_priority = spec_only$min_age_group_index_non_priority,     
                           runtime = spec_only$runtime,
                           seeding_cases = spec_only$seeding_cases,
                           output = "full")

spec_only_output <- nimue::format(spec_vaccine$model_output, 
                                  compartments = c("D"),
                                  reduce_age = TRUE) %>%
  filter(t > 1, compartment == "deaths") %>%
  group_by(replicate, t) %>%
  filter(t < cutoff_time) %>%
  mutate(value = 1000 * value / nothing$population_size)
spec_only_output$scenario <- "bspecific_only"

### BPSV
spec_bpsv <- baseline_scenarios_reduced %>%
  filter(R0 == R0_check, specific_vaccine_start == 100, NPI_int == NPI_int_check, detection_hosp == detection_hosp_check, vaccine_scenario == "both_vaccines")
both_vaccines <- run_sars_x(population_size = spec_bpsv$population_size,
                            country = spec_bpsv$country,
                            hosp_bed_capacity = spec_bpsv$hosp_bed_capacity,
                            ICU_bed_capacity = spec_bpsv$ICU_bed_capacity,
                            Rt = spec_bpsv$Rt, 
                            tt_Rt = spec_bpsv$tt_Rt,
                            Tg = spec_bpsv$Tg,
                            IFR = spec_bpsv$IFR,
                            vaccine_scenario = "both_vaccines", 
                            detection_time = spec_bpsv$detection_time,                      
                            bpsv_start = spec_bpsv$bpsv_start,                          
                            bpsv_protection_delay = spec_bpsv$bpsv_protection_delay,
                            specific_protection_delay = spec_bpsv$specific_protection_delay, 
                            specific_vaccine_start = spec_bpsv$specific_vaccine_start,                      
                            efficacy_infection_bpsv = spec_bpsv$efficacy_infection_bpsv,           
                            efficacy_disease_bpsv = spec_bpsv$efficacy_disease_bpsv,              
                            efficacy_infection_spec = spec_bpsv$efficacy_infection_spec,           
                            efficacy_disease_spec = spec_bpsv$efficacy_disease_spec,              
                            dur_R = spec_bpsv$dur_R,                       
                            dur_bpsv = spec_bpsv$dur_bpsv,                    
                            dur_spec = spec_bpsv$dur_spec,                    
                            coverage_bpsv = spec_bpsv$coverage_bpsv,                     
                            coverage_spec = spec_bpsv$coverage_spec,                     
                            vaccination_rate_bpsv = spec_bpsv$vaccination_rate_bpsv,             
                            vaccination_rate_spec = spec_bpsv$vaccination_rate_spec,             
                            min_age_group_index_priority = spec_bpsv$min_age_group_index_priority,        
                            min_age_group_index_non_priority = spec_bpsv$min_age_group_index_non_priority,     
                            runtime = spec_bpsv$runtime,
                            seeding_cases = spec_bpsv$seeding_cases,
                            output = "full")

both_output <- nimue::format(both_vaccines$model_output, 
                             compartments = c("D"),
                             reduce_age = TRUE) %>%
  filter(t > 1, compartment == "deaths") %>%
  group_by(replicate, t) %>%
  filter(t < cutoff_time) %>%
  mutate(value = 1000 * value / spec_bpsv$population_size)
both_output$scenario <- "cboth"

overall <- rbind(nothing_output, spec_only_output, both_output)

panel_apt3 <- ggplot(overall, aes(x = t, y = value, col = scenario)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("#948D9B", "#ADBDFF", "#4F96F9"),
                      labels = c("No Vaccines", "Disease-Specific Only", "Disease-Specific + BPSV"),
                      name = "100 Days to\nDisease-Specific Vaccine") +
  theme_bw() +
  coord_cartesian(xlim = c(40, 300)) +
  theme(legend.position = "none") +
  labs(x = "Time (Days)", y = "Daily Deaths per 1,000 Population")

overall_deaths <- overall %>%
  group_by(scenario) %>%
  summarise(total = sum(value))

panel_apt4 <- ggplot(overall_deaths, aes(x = scenario, y = total, fill = scenario)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#948D9B", "#ADBDFF", "#4F96F9"),
                    labels = c("No Vaccines", "Disease-Specific Only", "Disease-Specific + BPSV"),
                    name = "100 Days to\nDisease-Specific Vaccine") +
  scale_x_discrete(labels = c("No\nVaccines", "Disease-\nSpecific\nOnly", "Disease-\nSpecific\n+ BPSV")) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Total Deaths\nper 1,000 Population")

second_half <- cowplot::plot_grid(panel_apt3, panel_apt4, nrow = 1, 
                                  rel_widths = c(2.5, 1), axis = "b", align = "h",
                                  labels = c("B"))

second_half <- panel_apt3 + 
  annotation_custom(
    ggplotGrob(panel_apt4), 
    xmin = 230, xmax = 310, ymin = max(overall$value) / 5, ymax = max(overall$value)) + annotation_custom(
      ggplotGrob(inset_detection_plot), 
      xmin = 45, xmax = 140, 
      ymin = max(overall$value) / 2.5, ymax = max(overall$value))

second_half_pt2 <- cowplot::plot_grid(second_half, detection_sensitivity_plot100, nrow = 1, 
                                      align = "h", axis = "bt", 
                                      rel_widths = c(2.5, 1))

## dimensions = 11.82 x 3.95


# epidemic_curves <- cowplot::plot_grid(second_half, first_half, nrow = 2)
# 
# cowplot::plot_grid(epidemic_curves, detection_sensitivity_plot, ncol = 2, rel_widths = c(2, 1))

# model_outputs3 <- model_outputs2 %>%
#   pivot_longer(cols = c(bpsv_deaths_averted, deaths_spec),
#                names_to = "scenario",
#                values_to = "deaths")
# 
# ggplot(model_outputs3) +
#   geom_bar(aes(x = detection_hosp, y = 1000 * deaths / unique(baseline_scenarios$population_size),
#                fill = interaction(factor(specific_vaccine_start), scenario)), 
#            stat = "identity", size = 1) +
#   geom_bar(aes(x = detection_hosp, y = 1000 * deaths / unique(baseline_scenarios$population_size),
#                fill = interaction(factor(specific_vaccine_start), scenario)),
#            stat = "identity", size = 1) +
#   theme_bw() +
#   facet_grid(specific_vaccine_start ~ .,
#              labeller = as_labeller(c(`100` = "Specific Vaccine in 100 Days", `250` = "Specific Vaccine in 250 Days"))) +
#   guides(col = guide_legend(title = "R0"),
#          fill = "none") +
#   theme(strip.background = element_rect(fill="#F5F5F5"),
#         legend.position = "right")
# R0_check <- 2.5
# detection_hosp_check <- 1
# cutoff_time <- 300
# NPI_int_check <- 7
# 
# check <- scenarios %>%
#   filter(NPI_int == NPI_int_check, R0 == R0_check, specific_vaccine_start == 100, 
#          detection_hosp == detection_hosp_check, vaccine_scenario == "both_vaccines")
# index <- which(colnames(check2) %in% colnames(spec_bpsv))
# check_reduced <- check2[, index]
# 
# x <- gtools::smartbind(check, spec_bpsv)
# NPI_colours <- c("#C64191", "#F0803C", "#0D84A9")
# ggplot(subset(model_outputs2, NPI_int %in% c(7))) +
#   geom_line(aes(x = detection_time, y = 1000 * bpsv_deaths_averted / unique(baseline_scenarios$population_size), 
#                 col = factor(NPI_int)), size = 1) +
#   geom_point(aes(x = detection_time, y = 1000 * bpsv_deaths_averted / unique(baseline_scenarios$population_size), 
#                  fill = factor(NPI_int)), size = 1.5, pch = 21) +
#   theme_bw() +
#   facet_wrap(specific_vaccine_start ~ ., nrow = 2, scales = "free_y", 
#              labeller = as_labeller(c(`100` = "Specific Vaccine in 100 Days", `250` = "Specific Vaccine in 250 Days"))) +
#   labs(x = "Detection Time (Days)", y = "Additional Deaths Averted By BPSV (Per 1,000 Population)") +
#   guides(col = guide_legend(title = "NPI Scenario"),
#          fill = "none") +
#   scale_colour_manual(values = NPI_colours) +
#   scale_fill_manual(values = NPI_colours) +
#   ylim(c(0, NA)) +
#   scale_y_continuous(position = "right") +
#   theme(strip.background = element_rect(fill="#F5F5F5"),
#         legend.position = "bottom")
# 
# 
# ggplot(subset(model_outputs2, NPI_int %in% c(4, 7, 8))) +
#   geom_line(aes(x = detection_time, y = 1000 * deaths_spec / unique(baseline_scenarios$population_size), 
#                 col = factor(NPI_int)), size = 1) +
#   geom_point(aes(x = detection_time, y = 1000 * deaths_spec / unique(baseline_scenarios$population_size), 
#                  fill = factor(NPI_int)), size = 1.5, pch = 21) +
#   theme_bw() +
#   facet_grid(specific_vaccine_start ~ ., #scales = "free_y",
#              labeller = as_labeller(c(`100` = "Specific Vaccine in 100 Days", `250` = "Specific Vaccine in 250 Days"))) +
#   labs(x = "Hospitalisation Detection Threshold", y = "Total Deaths (Per 1,000 Population)") +
#   guides(col = guide_legend(title = "NPI Scenario"),
#          fill = "none") +
#   ylim(c(0, NA)) +
#   theme(strip.background = element_rect(fill="#F5F5F5"),
#         legend.position = "right")
# 
# scale_colour_manual(values = c("#B8336A", "#726DA8", "#42A1B6"))
  