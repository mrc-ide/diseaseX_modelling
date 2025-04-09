# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/ring_vax_simulation.R"))
source(here::here("functions/implement_quarantine.R"))

### Fixed model parameters

### SARS-CoV-1-specific parameters (longer Tg, lower presymptomatic transmission, lower proportion asymptomatic infection)
SC1_generation_time <- function(n) { rgamma(n, shape = 24, rate = 2) } 
SC1_infection_to_onset <- function(n) { rgamma(n, shape = 0.1, rate = 1) } 
SC1_prop_asymptomatic <- 0

### SARS-CoV-2-specific parameters (shorter Tg, higher presymptomatic transmission, higher proportion asymptomatic infection)
SC2_generation_time <- function(n) { rgamma(n, shape = 13.5, rate = 2) } # 6.75 day generation time Gamam distributed (as per Walker et al, Science, 2020)
SC2_infection_to_onset <- function(n) { rgamma(n, shape = 13.5/3, rate = 2) } # ~35% of transmission presymptomatic (per SARS-CoV-2, slightly lower than but roughly aligned with: https://bmjopen.bmj.com/content/11/6/e041240)
SC2_prop_asymptomatic <- 0.15

vaccine_start <- 21 # detection + logistical delay of 3 weeks to vaccination initiation
vaccine_coverage <- 0.8
pop <- 10^10
check_final_size <- 2500
initial_immune <- 0
seeding_cases <- 5

quarantine_time <- function(n) {return(rep(2.6, n))} # from https://pmc.ncbi.nlm.nih.gov/articles/PMC7511527/?utm_source=chatgpt.com#cesec10
prob_quarantine_scan <- c(0, 0.5, 0.7, 0.9) # from https://pmc.ncbi.nlm.nih.gov/articles/PMC7511527/?utm_source=chatgpt.com#cesec10
quarantine_efficacy_scan <- c(0.35, 0.65) # from https://pmc.ncbi.nlm.nih.gov/articles/PMC7511527/?utm_source=chatgpt.com#cesec10

vaccine_efficacy_infection <- 0.35 
vaccine_efficacy_transmission <- 0.35 
vaccine_logistical_delay <- 2
R0_scan <- c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5)
iterations <- 2

## R0 sensitivity analysis (Figure 1B)
fresh_run_R0_sensitivity_analysis <- FALSE
if (fresh_run_R0_sensitivity_analysis) {
  
  ### Setting up R0 scan and the storage for each of the different protection delays
  SC1_storage_nothing <- array(data = NA, dim = c(iterations, length(R0_scan), length(prob_quarantine_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_2weeks <- array(data = NA, dim = c(iterations, length(R0_scan), length(prob_quarantine_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_1week <- array(data = NA, dim = c(iterations, length(R0_scan), length(prob_quarantine_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_2days <- array(data = NA, dim = c(iterations, length(R0_scan), length(prob_quarantine_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_instant <- array(data = NA, dim = c(iterations, length(R0_scan), length(prob_quarantine_scan), length(quarantine_efficacy_scan)))
  
  SC2_storage_nothing <- array(data = NA, dim = c(iterations, length(R0_scan), length(prob_quarantine_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_2weeks <- array(data = NA, dim = c(iterations, length(R0_scan), length(prob_quarantine_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_1week <- array(data = NA, dim = c(iterations, length(R0_scan), length(prob_quarantine_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_2days <- array(data = NA, dim = c(iterations, length(R0_scan), length(prob_quarantine_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_instant <- array(data = NA, dim = c(iterations, length(R0_scan), length(prob_quarantine_scan), length(quarantine_efficacy_scan)))
  
  set.seed(2000)
  seeds <- runif(n = iterations, min = 1, max = 10^9)
  
  for (i in 1:length(R0_scan)) {
    for (j in 1:length(prob_quarantine_scan)) {
      for (k in 1:length(quarantine_efficacy_scan)) {
        for (l in 1:iterations) {
          
          ## 2 weeks delay
          SC1_vacc_2weeks <- ring_vax_bp_sim(offspring = "pois",
                                             mn_offspring = R0_scan[i],
                                             generation_time = SC1_generation_time,
                                             t0 = 0, tf = Inf, population = pop, check_final_size = check_final_size, 
                                             initial_immune = initial_immune,
                                             seeding_cases = seeding_cases, 
                                             seed = seeds[l],
                                             prop_asymptomatic = SC1_prop_asymptomatic,
                                             infection_to_onset = SC1_infection_to_onset,
                                             vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                             vaccine_efficacy_infection = vaccine_efficacy_infection,
                                             vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                             vaccine_logistical_delay = vaccine_logistical_delay,
                                             vaccine_protection_delay = 14,
                                             time_to_quarantine = quarantine_time,
                                             prob_quarantine = prob_quarantine_scan[j],
                                             quarantine_efficacy = quarantine_efficacy_scan[k])
          SC1_storage_vacc_2weeks[l, i, j, k] <- sum(!is.na(SC1_vacc_2weeks$time_infection))
          
          SC2_vacc_2weeks <- ring_vax_bp_sim(offspring = "pois",
                                             mn_offspring = R0_scan[i],
                                             generation_time = SC2_generation_time,
                                             t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
                                             initial_immune = initial_immune,
                                             seeding_cases = seeding_cases, 
                                             seed = seeds[l],
                                             prop_asymptomatic = SC2_prop_asymptomatic,
                                             infection_to_onset = SC2_infection_to_onset,
                                             vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                             vaccine_efficacy_infection = vaccine_efficacy_infection,
                                             vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                             vaccine_logistical_delay = vaccine_logistical_delay,
                                             vaccine_protection_delay = 14,
                                             time_to_quarantine = quarantine_time,
                                             prob_quarantine = prob_quarantine_scan[j],
                                             quarantine_efficacy = quarantine_efficacy_scan[k])
          SC2_storage_vacc_2weeks[l, i, j, k] <- sum(!is.na(SC2_vacc_2weeks$time_infection))
          
          ## 1 week delay
          SC1_vacc_1week <- ring_vax_bp_sim(offspring = "pois",
                                            mn_offspring = R0_scan[i],
                                            generation_time = SC1_generation_time,
                                            t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
                                            initial_immune = initial_immune,
                                            seeding_cases = seeding_cases, 
                                            seed = seeds[l],
                                            prop_asymptomatic = SC1_prop_asymptomatic,
                                            infection_to_onset = SC1_infection_to_onset,
                                            vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                            vaccine_efficacy_infection = vaccine_efficacy_infection,
                                            vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                            vaccine_logistical_delay = vaccine_logistical_delay,
                                            vaccine_protection_delay = 7,
                                            time_to_quarantine = quarantine_time,
                                            prob_quarantine = prob_quarantine_scan[j],
                                            quarantine_efficacy = quarantine_efficacy_scan[k])
          SC1_storage_vacc_1week[l, i, j, k] <- sum(!is.na(SC1_vacc_1week$time_infection))
          
          SC2_vacc_1week <- ring_vax_bp_sim(offspring = "pois",
                                            mn_offspring = R0_scan[i],
                                            generation_time = SC2_generation_time,
                                            t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
                                            initial_immune = initial_immune,
                                            seeding_cases = seeding_cases, 
                                            seed = seeds[l],
                                            prop_asymptomatic = SC2_prop_asymptomatic,
                                            infection_to_onset = SC2_infection_to_onset,
                                            vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                            vaccine_efficacy_infection = vaccine_efficacy_infection,
                                            vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                            vaccine_logistical_delay = vaccine_logistical_delay,
                                            vaccine_protection_delay = 7,
                                            time_to_quarantine = quarantine_time,
                                            prob_quarantine = prob_quarantine_scan[j],
                                            quarantine_efficacy = quarantine_efficacy_scan[k])
          SC2_storage_vacc_1week[l, i, j, k] <- sum(!is.na(SC2_vacc_1week$time_infection))
          
          ## 2 days delay
          SC1_vacc_2days <- ring_vax_bp_sim(offspring = "pois",
                                            mn_offspring = R0_scan[i],
                                            generation_time = SC1_generation_time,
                                            t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
                                            initial_immune = initial_immune,
                                            seeding_cases = seeding_cases, 
                                            seed = seeds[l],
                                            prop_asymptomatic = SC1_prop_asymptomatic,
                                            infection_to_onset = SC1_infection_to_onset,
                                            vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                            vaccine_efficacy_infection = vaccine_efficacy_infection,
                                            vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                            vaccine_logistical_delay = vaccine_logistical_delay,
                                            vaccine_protection_delay = 2,
                                            time_to_quarantine = quarantine_time,
                                            prob_quarantine = prob_quarantine_scan[j],
                                            quarantine_efficacy = quarantine_efficacy_scan[k])
          SC1_storage_vacc_2days[l, i, j, k] <- sum(!is.na(SC1_vacc_2days$time_infection))
          
          SC2_vacc_2days <- ring_vax_bp_sim(offspring = "pois",
                                            mn_offspring = R0_scan[i],
                                            generation_time = SC2_generation_time,
                                            t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
                                            initial_immune = initial_immune,
                                            seeding_cases = seeding_cases, 
                                            seed = seeds[l],
                                            prop_asymptomatic = SC2_prop_asymptomatic,
                                            infection_to_onset = SC2_infection_to_onset,
                                            vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                            vaccine_efficacy_infection = vaccine_efficacy_infection,
                                            vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                            vaccine_logistical_delay = vaccine_logistical_delay,
                                            vaccine_protection_delay = 2,
                                            time_to_quarantine = quarantine_time,
                                            prob_quarantine = prob_quarantine_scan[j],
                                            quarantine_efficacy = quarantine_efficacy_scan[k])
          SC2_storage_vacc_2days[l, i, j, k] <- sum(!is.na(SC2_vacc_2days$time_infection))
          
          ## Instant
          SC1_vacc_instant <- ring_vax_bp_sim(offspring = "pois",
                                              mn_offspring = R0_scan[i],
                                              generation_time = SC1_generation_time,
                                              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
                                              initial_immune = initial_immune,
                                              seeding_cases = seeding_cases, 
                                              seed = seeds[l],
                                              prop_asymptomatic = SC1_prop_asymptomatic,
                                              infection_to_onset = SC1_infection_to_onset,
                                              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                              vaccine_efficacy_infection = vaccine_efficacy_infection,
                                              vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                              vaccine_logistical_delay = vaccine_logistical_delay,
                                              vaccine_protection_delay = 0.01,
                                              time_to_quarantine = quarantine_time,
                                              prob_quarantine = prob_quarantine_scan[j],
                                              quarantine_efficacy = quarantine_efficacy_scan[k])
          SC1_storage_vacc_instant[l, i, j, k] <- sum(!is.na(SC1_vacc_instant$time_infection))
          
          SC2_vacc_instant <- ring_vax_bp_sim(offspring = "pois",
                                              mn_offspring = R0_scan[i],
                                              generation_time = SC2_generation_time,
                                              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
                                              initial_immune = initial_immune,
                                              seeding_cases = seeding_cases, 
                                              seed = seeds[l],
                                              prop_asymptomatic = SC2_prop_asymptomatic,
                                              infection_to_onset = SC2_infection_to_onset,
                                              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                              vaccine_efficacy_infection = vaccine_efficacy_infection,
                                              vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                              vaccine_logistical_delay = vaccine_logistical_delay,
                                              vaccine_protection_delay = 0.01,
                                              time_to_quarantine = quarantine_time,
                                              prob_quarantine = prob_quarantine_scan[j],
                                              quarantine_efficacy = quarantine_efficacy_scan[k])
          SC2_storage_vacc_instant[l, i, j, k] <- sum(!is.na(SC2_vacc_instant$time_infection))
          
          ## Nothing
          SC1_no_vacc <- ring_vax_bp_sim(offspring = "pois",
                                         mn_offspring = R0_scan[i],
                                         generation_time = SC1_generation_time,
                                         t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, initial_immune = initial_immune,
                                         seeding_cases = seeding_cases, 
                                         seed = seeds[l],
                                         prop_asymptomatic = SC1_prop_asymptomatic,
                                         infection_to_onset = SC1_infection_to_onset,
                                         vaccine_start = 1000, vaccine_coverage = vaccine_coverage,
                                         vaccine_efficacy_infection = vaccine_efficacy_infection,
                                         vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                         vaccine_logistical_delay = vaccine_logistical_delay,
                                         vaccine_protection_delay = 100,
                                         time_to_quarantine = quarantine_time,
                                         prob_quarantine = prob_quarantine_scan[j],
                                         quarantine_efficacy = quarantine_efficacy_scan[k])
          SC1_storage_nothing[l, i, j, k] <- sum(!is.na(SC1_no_vacc$time_infection))
          
          SC2_no_vacc <- ring_vax_bp_sim(offspring = "pois",
                                         mn_offspring = R0_scan[i],
                                         generation_time = SC2_generation_time,
                                         t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, initial_immune = initial_immune,
                                         seeding_cases = seeding_cases, 
                                         seed = seeds[l],
                                         prop_asymptomatic = SC2_prop_asymptomatic,
                                         infection_to_onset = SC2_infection_to_onset,
                                         vaccine_start = 1000, vaccine_coverage = vaccine_coverage,
                                         vaccine_efficacy_infection = vaccine_efficacy_infection,
                                         vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                         vaccine_logistical_delay = vaccine_logistical_delay,
                                         vaccine_protection_delay = 100,
                                         time_to_quarantine = quarantine_time,
                                         prob_quarantine = prob_quarantine_scan[j],
                                         quarantine_efficacy = quarantine_efficacy_scan[k])
          SC2_storage_nothing[l, i, j, k] <- sum(!is.na(SC2_no_vacc$time_infection))
          print(paste0("l = ", l))
        }
        print(paste0("k = ", k))
      }
      print(paste0("j = ", j))
    }
    print(paste0("i = ", i))
  }
  
  # Convert each array to a long data frame.
  df_SC1_nothing <- convert_array(arr = SC1_storage_nothing, scenario = "no_vaccination", pathogen = "SARS-CoV-1", R0_scan, prob_quarantine_scan, quarantine_efficacy_scan)
  df_SC2_nothing <- convert_array(SC2_storage_nothing, scenario = "no_vaccination", pathogen = "SARS-CoV-2", R0_scan, prob_quarantine_scan, quarantine_efficacy_scan)
  
  df_SC1_2weeks <- convert_array(SC1_storage_vacc_2weeks, scenario = "2weeks_delay", pathogen = "SARS-CoV-1", R0_scan, prob_quarantine_scan, quarantine_efficacy_scan)
  df_SC2_2weeks <- convert_array(SC2_storage_vacc_2weeks, scenario = "2weeks_delay", pathogen = "SARS-CoV-2", R0_scan, prob_quarantine_scan, quarantine_efficacy_scan)
  
  df_SC1_1week <- convert_array(SC1_storage_vacc_1week, scenario = "1week_delay", pathogen = "SARS-CoV-1", R0_scan, prob_quarantine_scan, quarantine_efficacy_scan)
  df_SC2_1week <- convert_array(SC2_storage_vacc_1week, scenario = "1week_delay", pathogen = "SARS-CoV-2", R0_scan, prob_quarantine_scan, quarantine_efficacy_scan)
  
  df_SC1_2days <- convert_array(SC1_storage_vacc_2days, scenario = "2days_delay", pathogen = "SARS-CoV-1", R0_scan, prob_quarantine_scan, quarantine_efficacy_scan)
  df_SC2_2days <- convert_array(SC2_storage_vacc_2days, scenario = "2days_delay", pathogen = "SARS-CoV-2", R0_scan, prob_quarantine_scan, quarantine_efficacy_scan)
  
  df_SC1_instant <- convert_array(SC1_storage_vacc_instant, scenario = "no_delay", pathogen = "SARS-CoV-1", R0_scan, prob_quarantine_scan, quarantine_efficacy_scan)
  df_SC2_instant <- convert_array(SC2_storage_vacc_instant, scenario = "no_delay", pathogen = "SARS-CoV-2", R0_scan, prob_quarantine_scan, quarantine_efficacy_scan)
  
  # Combine all results into a single data frame.
  overall_bp_df <- bind_rows(
    df_SC1_nothing, df_SC1_2weeks, df_SC1_1week, df_SC1_2days, df_SC1_instant,
    df_SC2_nothing, df_SC2_2weeks, df_SC2_1week, df_SC2_2days, df_SC2_instant
  )
  
  saveRDS(object = overall_bp_df, file = "outputs/Figure1_branchingProcess_Containment/New_Fig1_ringVaccination_paramScan.rds")
} else {
  overall_bp_df <- readRDS("outputs/Figure1_branchingProcess_Containment/New_Fig1_ringVaccination_paramScan.rds")
}

containment_df <- overall_bp_df %>%
  mutate(contained = ifelse(epidemic_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(R0, scenario, pathogen, prob_quarantine, quarantine_efficacy) %>%
  summarise(proportion_contained = sum(contained) / iterations) %>%
  mutate(proportion_contained = ifelse(R0 == 1.00, 1, proportion_contained)) %>%
  mutate(scenario = ifelse(scenario == "no_vaccination", "zno_vaccination", scenario)) %>%
  mutate(scenario = ifelse(scenario == "2days_delay", "bvacc_2days_protectDelay", scenario)) %>%
  mutate(scenario = ifelse(scenario == "1week_delay", "dvacc_1week_protectDelay", scenario)) %>%
  mutate(scenario = ifelse(scenario == "2weeks_delay", "evacc_2weeks_protectDelay", scenario)) %>%
  mutate(scenario = ifelse(scenario == "no_delay", "avacc_no_delay", scenario)) 
containment_df$scenario <- factor(containment_df$scenario, 
                                  levels = c("avacc_no_delay",
                                             "bvacc_2days_protectDelay",
                                             "dvacc_1week_protectDelay", 
                                             "evacc_2weeks_protectDelay",
                                             "zno_vaccination"))
containment_df2 <- containment_df %>%
  arrange(scenario) 
ggplot(subset(containment_df2, pathogen == "SARS-CoV-1"), 
       aes(x = R0, y = 100 * proportion_contained, col = scenario)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_grid(prob_quarantine~quarantine_efficacy) + 
  scale_colour_manual(values = c("#CA2E6B", "#88C5EE", "#236897", "#13496E", "black"), 
                      labels = c("No Delay", "2 Days", "1 Week", "2 Weeks", "No Vaccination"),
                      name = "Vaccine\nProtection\nDelay",
                      guide = guide_legend(reverse = TRUE)) +
  labs(x = "R0", y = "% Outbreaks Contained") +
  theme(strip.background = element_rect(fill = "white"))

ggplot(subset(containment_df2, pathogen == "SARS-CoV-2"), 
       aes(x = R0, y = 100 * proportion_contained, col = scenario)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_grid(prob_quarantine~quarantine_efficacy) + 
  scale_colour_manual(values = c("#CA2E6B", "#88C5EE", "#236897", "#13496E", "black"), 
                      labels = c("No Delay", "2 Days", "1 Week", "2 Weeks", "No Vaccination"),
                      name = "Vaccine\nProtection\nDelay",
                      guide = guide_legend(reverse = TRUE)) +
  labs(x = "R0", y = "% Outbreaks Contained") +
  theme(strip.background = element_rect(fill = "white"))

###### Old Code ######

## Plotting Figure 1B
containment_df <- overall_bp_df %>%
  mutate(contained = ifelse(`Epidemic Size` < (0.9 * check_final_size), 1, 0)) %>%
  group_by(actualR0, R0, scenario, pathogen) %>%
  summarise(proportion_contained = sum(contained) / iterations) %>%
  mutate(proportion_contained = ifelse(actualR0 == 1.00, 1, proportion_contained)) %>%
  mutate(scenario = ifelse(scenario == "no_vaccination", "zno_vaccination", scenario)) %>%
  mutate(scenario = ifelse(scenario == "2days_delay", "bvacc_2days_protectDelay", scenario)) %>%
  mutate(scenario = ifelse(scenario == "1week_delay", "dvacc_1week_protectDelay", scenario)) %>%
  mutate(scenario = ifelse(scenario == "2weeks_delay", "evacc_2weeks_protectDelay", scenario)) %>%
  mutate(scenario = ifelse(scenario == "no_delay", "avacc_no_delay", scenario)) 
containment_df$scenario <- factor(containment_df$scenario, 
                                  levels = c("avacc_no_delay",
                                             "bvacc_2days_protectDelay",
                                             "dvacc_1week_protectDelay", 
                                             "evacc_2weeks_protectDelay",
                                             "zno_vaccination"))
containment_df2 <- containment_df %>%
  arrange(scenario) 
containment_plot <- ggplot(containment_df2, aes(x = actualR0, y = 100 * proportion_contained, col = scenario)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_grid(.~pathogen) + 
  scale_colour_manual(values = c("#CA2E6B", "#88C5EE", "#236897", "#13496E", "black"), 
                      labels = c("No Delay", "2 Days", "1 Week", "2 Weeks", "No Vaccination"),
                      name = "Vaccine\nProtection\nDelay",
                      guide = guide_legend(reverse = TRUE)) +
  labs(x = "R0", y = "% Outbreaks Contained") +
  theme(strip.background = element_rect(fill = "white"))
ggsave(filename = "figures/Figure_1_BranchingProcess/Fig1BC_ContainmentPlot.pdf", plot = containment_plot, 
       width = 8, height = 3.1)

## Vaccination-Related Sensitivity Analyses Heatmaps
fresh_run_vaccination_heatmaps <- FALSE
if (fresh_run_vaccination_heatmaps) {
  
  ## R0 sensitivity analysis for all the runs
  R0_seq <- R0_scan[R0_scan > 1]
  vaccine_protection_delay <- 7
  
  ## Sensitivity Analysis - R0 vs Ratio of Tg to Protection Delay
  Tg_ratio_seq <- seq(1, 4, 0.5)
  storage_R0_TgRatio_sensitivity <- array(data = NA, dim = c(length(R0_seq), length(Tg_ratio_seq), iterations))
  for (i in 1:length(R0_seq)) {
    for (j in 1:length(Tg_ratio_seq)) {
      for (k in 1:iterations) {
        generation_time <- function(n) { rgamma(n, shape = 2 * vaccine_protection_delay * Tg_ratio_seq[j], rate = 2) } 
        infection_to_onset <- function(n) { rgamma(n, shape = (2 * vaccine_protection_delay * Tg_ratio_seq[j])/3, rate = 2) } # keeping proportion of presymptomatic transmission constant as Tg varies
        bp_out <- chain_sim_susc_ring_vacc2(offspring = "pois",
                                            mn_offspring = R0_seq[i],
                                            generation_time = generation_time,
                                            t0 = 0, tf = Inf, pop = pop, 
                                            check_final_size = check_final_size, initial_immune = initial_immune,
                                            seeding_cases = seeding_cases, prop_asymptomatic = prop_asymptomatic,
                                            infection_to_onset = infection_to_onset,
                                            vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                            vaccine_efficacy_infection = vaccine_efficacy_infection,
                                            vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                            vaccine_logistical_delay = vaccine_logistical_delay,
                                            vaccine_protection_delay = vaccine_protection_delay)
        storage_R0_TgRatio_sensitivity[i, j, k] <- sum(!is.na(bp_out$time_infection))
      }
    }
    print(i)
  }
  storage_R0_TgRatio_df <- reshape2::melt(storage_R0_TgRatio_sensitivity)
  storage_R0_TgRatio_df$R0 <- rep(R0_seq, length(Tg_ratio_seq) * iterations)
  storage_R0_TgRatio_df$TgRatio <- rep(Tg_ratio_seq, each = length(R0_seq), times = iterations)
  storage_R0_TgRatio_df$iteration <- rep(1:iterations, each = length(R0_seq) * length(Tg_ratio_seq))
  storage_R0_TgRatio_df <- storage_R0_TgRatio_df %>%
    # select(-(Var1:Var3)) %>%
    rename(final_size = value) %>%
    mutate(contained = ifelse(final_size < (0.9 * check_final_size), 1, 0)) %>%
    group_by(R0, TgRatio) %>%
    summarise(proportion_contained = sum(contained) / n())
  saveRDS(object = storage_R0_TgRatio_df, file = "outputs/Figure1_branchingProcess_Containment/Fig1D_branchingProcess_ringVaccination_R0TgRatioScan.rds")
  
  ## Sensitivity Analysis - R0 vs Vaccine Efficacy
  Tg_ratio_fixed <- 2.5
  generation_time <- function(n) { rgamma(n, shape = 2 * vaccine_protection_delay * Tg_ratio_fixed, rate = 2) } 
  infection_to_onset <- function(n) { rgamma(n, shape = (2 * vaccine_protection_delay * Tg_ratio_fixed)/3, rate = 2) }
  vaccine_efficacy_seq <- seq(0.3, 0.9, 0.1)
  storage_R0_efficacy_sensitivity <- array(data = NA, dim = c(length(R0_seq), length(vaccine_efficacy_seq), iterations))
  for (i in 1:length(R0_seq)) {
    for (j in 1:length(vaccine_efficacy_seq)) {
      for (k in 1:iterations) {
        bp_out <- chain_sim_susc_ring_vacc2(offspring = "pois",
                                            mn_offspring = R0_seq[i],
                                            generation_time = generation_time,
                                            t0 = 0, tf = Inf, pop = pop,
                                            check_final_size = check_final_size, initial_immune = initial_immune,
                                            seeding_cases = seeding_cases, prop_asymptomatic = prop_asymptomatic,
                                            infection_to_onset = infection_to_onset,
                                            vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                            vaccine_efficacy_infection = vaccine_efficacy_seq[j],
                                            vaccine_efficacy_transmission = vaccine_efficacy_seq[j],
                                            vaccine_logistical_delay = vaccine_logistical_delay,
                                            vaccine_protection_delay = vaccine_protection_delay)
        storage_R0_efficacy_sensitivity[i, j, k] <- sum(!is.na(bp_out$time_infection))
      }
    }
    print(i)
  }
  storage_R0_efficacy_df <- reshape2::melt(storage_R0_efficacy_sensitivity)
  storage_R0_efficacy_df$R0 <- rep(R0_seq, length(vaccine_efficacy_seq) * iterations)
  storage_R0_efficacy_df$vaccine_efficacy <- rep(vaccine_efficacy_seq, each = length(R0_seq), times = iterations)
  storage_R0_efficacy_df$iteration <- rep(1:iterations, each = length(R0_seq) * length(vaccine_efficacy_seq))
  storage_R0_efficacy_df <- storage_R0_efficacy_df %>%
    select(-(Var1:Var3)) %>%
    rename(final_size = value) %>%
    mutate(contained = ifelse(final_size < (0.9 * check_final_size), 1, 0)) %>%
    group_by(R0, vaccine_efficacy) %>%
    summarise(proportion_contained = sum(contained) / n())
  saveRDS(object = storage_R0_efficacy_df, file = "outputs/Figure1_branchingProcess_Containment/Fig1E_branchingProcess_ringVaccination_R0EfficacyScan.rds")
  
  ## Sensitivity Analysis - R0 vs Presymptomatic Transmission Proportion
  Tg_ratio_fixed <- 2.5
  generation_time <- function(n) { rgamma(n, shape = 2 * vaccine_protection_delay * Tg_ratio_fixed, rate = 2) } 
  proportion_presymptomatic_seq <- seq(0.1, 0.7, 0.1)
  storage_R0_preSymp_sensitivity <- array(data = NA, dim = c(length(R0_seq), length(proportion_presymptomatic_seq), iterations))
  for (i in 1:length(R0_seq)) {
    for (j in 1:length(proportion_presymptomatic_seq)) {
      for (k in 1:iterations) {
        infection_to_onset <- function(n) { rgamma(n, shape = (2 * vaccine_protection_delay * Tg_ratio_fixed) * proportion_presymptomatic_seq[j], rate = 2) }
        bp_out <- chain_sim_susc_ring_vacc2(offspring = "pois",
                                            mn_offspring = R0_seq[i],
                                            generation_time = generation_time,
                                            t0 = 0, tf = Inf, pop = 10^7,
                                            check_final_size = check_final_size, initial_immune = initial_immune,
                                            seeding_cases = seeding_cases, prop_asymptomatic = prop_asymptomatic,
                                            infection_to_onset = infection_to_onset,
                                            vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                            vaccine_efficacy_infection = vaccine_efficacy_infection,
                                            vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                            vaccine_logistical_delay = vaccine_logistical_delay,
                                            vaccine_protection_delay = vaccine_protection_delay)
        storage_R0_preSymp_sensitivity[i, j, k] <- sum(!is.na(bp_out$time_infection))
      }
    }
    print(i)
  }
  storage_R0_preSymp_df <- reshape2::melt(storage_R0_preSymp_sensitivity)
  storage_R0_preSymp_df$R0 <- rep(R0_seq, length(proportion_presymptomatic_seq) * iterations)
  storage_R0_preSymp_df$prop_preSymp <- rep(proportion_presymptomatic_seq, each = length(R0_seq), times = iterations)
  storage_R0_preSymp_df$iteration <- rep(1:iterations, each = length(R0_seq) * length(proportion_presymptomatic_seq))
  storage_R0_preSymp_df <- storage_R0_preSymp_df %>%
    select(-(Var1:Var3)) %>%
    rename(final_size = value) %>%
    mutate(contained = ifelse(final_size < (0.9 * check_final_size), 1, 0)) %>%
    group_by(R0, prop_preSymp) %>%
    summarise(proportion_contained = sum(contained) / n())
  saveRDS(object = storage_R0_preSymp_df, file = "outputs/Figure1_branchingProcess_Containment/Fig1F_branchingProcess_ringVaccination_R0preSympScan.rds")
  
} else {
  storage_R0_TgRatio_df <- readRDS("outputs/Figure1_branchingProcess_Containment/Fig1D_branchingProcess_ringVaccination_R0TgRatioScan.rds")
  storage_R0_efficacy_df <- readRDS("outputs/Figure1_branchingProcess_Containment/Fig1E_branchingProcess_ringVaccination_R0EfficacyScan.rds")
  storage_R0_preSymp_df <- readRDS("outputs/Figure1_branchingProcess_Containment/Fig1F_branchingProcess_ringVaccination_R0preSympScan.rds")
}

## Creating Vaccination-Related Heatmaps

### Tg Ratio
R0_TgRatio_plot <- ggplot(storage_R0_TgRatio_df, aes(x = R0, y = TgRatio, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0, 1), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  labs(x = "R0",
       y = "Ratio of Tg to Vaccine Protection Delay") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none",
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)
ggsave(file = "figures/Figure_1_BranchingProcess/Fig1C_TgRatio_Sensitivity_Analysis.pdf", plot = R0_TgRatio_plot, width = 2.4, height = 2.4)

### Vaccine Efficacy
R0_vaccine_efficacy_plot <- ggplot(storage_R0_efficacy_df, aes(x = R0, y = 100 * vaccine_efficacy, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0, 1), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  labs(x = "R0",
       y = "Vaccine Efficacy (%)") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none",
        panel.border = element_rect(linetype = "solid", color = "black", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)
ggsave(filename = "figures/Figure_1_BranchingProcess/Fig1D_VaccineEfficacy_Sensitivity_Analysis.pdf", plot = R0_vaccine_efficacy_plot, width = 2.4, height = 2.4)

### Proportion of Presymptomatic Transmission
R0_preSymp_plot <- ggplot(storage_R0_preSymp_df, aes(x = R0, y = 100 * prop_preSymp, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0, 1), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  labs(x = "R0",
       y = "% Presymptomatic Transmission") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none",
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)
ggsave(filename = "figures/Figure_1_BranchingProcess/Fig1E_preSymp_Sensitivity_Analysis.pdf", 
       plot = R0_preSymp_plot, width = 2.4, height = 2.4)

## Heatmaps legend
legend <- R0_TgRatio_plot + theme(legend.position = "bottom")
ggsave(file = "figures/Figure_1_BranchingProcess/Legend_ringVaccinationHeatmap.pdf", plot = legend, width = 2.4, height = 2.4)
