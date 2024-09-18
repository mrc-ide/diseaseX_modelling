# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/branching_process_spatial_vaccination.R"))

### Fixed model parameters

### SARS-CoV-1-specific parameters (longer Tg, lower presymptomatic transmission, lower proportion asymptomatic infection)
SC1_generation_time <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
SC1_infection_to_onset <- function(n) { rgamma(n, shape = 0.1, rate = 1) } ## (negligible assumed currently)
SC1_prop_asymptomatic <- 0
SC1_prob_hosp <- 0.9 # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
SC1_hospitalisation_delay <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/

### SARS-CoV-2-specific parameters (shorter Tg, higher presymptomatic transmission, higher proportion asymptomatic infection)
SC2_generation_time <- function(n) { rgamma(n, shape = 13.5, rate = 2) } # 6.75 day generation time Gamma distributed (as per Walker et al, Science, 2020)
SC2_infection_to_onset <- function(n) { rgamma(n, shape = 13.5/3, rate = 2) } # ~35% of transmission presymptomatic (per SARS-CoV-2, slightly lower than but roughly aligned with: https://bmjopen.bmj.com/content/11/6/e041240)
SC2_prop_asymptomatic <- 0.15
SC2_prob_hosp <- 0.05
SC2_hospitalisation_delay <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/

### Spatial kernel related parameters
mu <- 25
size <- 5
spatial_kernel <- function(n) { rnbinom(n, size = 5, mu = 25) }
spatial_kernel_25 <- 10 * mu # qnbinom(p = 0.25, size = 5, mu = 25)
spatial_kernel_50 <- 50 * mu # qnbinom(p = 0.50, size = 5, mu = 25)
spatial_kernel_75 <- 100 * mu # qnbinom(p = 0.75, size = 5, mu = 25)
spatial_kernel_95 <- 250 * mu # qnbinom(p = 0.95, size = 5, mu = 25)

detection_threshold <- 5
vaccine_coverage <- 0.8
pop <- 10^10
check_final_size <- 2000
initial_immune <- 0
seeding_cases <- 5

vaccine_efficacy_infection <- 0.35 
vaccine_efficacy_transmission <- 0.35 
vaccine_efficacy_disease <- 0.75
vaccine_logistical_delay <- 2
vaccine_protection_delay <- 7
R0_scan <- c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5)
iterations <- 25

## R0 sensitivity analysis (Figure 1B)
fresh_run_R0_sensitivity_analysis <- FALSE
if (fresh_run_R0_sensitivity_analysis) {
  
  ### Setting up R0 scan and the storage for each of the different protection delays
  SC1_storage_nothing <- matrix(nrow = iterations, ncol = length(R0_scan))
  SC1_storage_vacc_25 <- matrix(nrow = iterations, ncol = length(R0_scan))
  SC1_storage_vacc_50 <- matrix(nrow = iterations, ncol = length(R0_scan))
  SC1_storage_vacc_75 <- matrix(nrow = iterations, ncol = length(R0_scan))
  SC1_storage_vacc_95 <- matrix(nrow = iterations, ncol = length(R0_scan))
  
  SC2_storage_nothing <- matrix(nrow = iterations, ncol = length(R0_scan))
  SC2_storage_vacc_25 <- matrix(nrow = iterations, ncol = length(R0_scan))
  SC2_storage_vacc_50 <- matrix(nrow = iterations, ncol = length(R0_scan))
  SC2_storage_vacc_75 <- matrix(nrow = iterations, ncol = length(R0_scan))
  SC2_storage_vacc_95 <- matrix(nrow = iterations, ncol = length(R0_scan))
  
  for (i in 1:length(R0_scan)) {
    for (j in 1:iterations) {
      
      ## 25% of spatial kernel
      SC1_vacc_25 <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                          generation_time = SC1_generation_time,
                                          spatial_kernel = spatial_kernel,
                                          t0 = 0, tf = Inf,
                                          check_final_size = check_final_size,
                                          seeding_cases = seeding_cases,
                                          prop_asymptomatic = SC1_prop_asymptomatic,
                                          prob_hosp = SC1_prob_hosp,
                                          hospitalisation_delay = SC1_hospitalisation_delay,
                                          detection_threshold = detection_threshold,
                                          vaccine_campaign_radius = spatial_kernel_25,
                                          vaccine_coverage = vaccine_coverage,
                                          vaccine_efficacy_infection = vaccine_efficacy_infection,
                                          vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                          vaccine_efficacy_disease = vaccine_efficacy_disease,
                                          vaccine_logistical_delay = vaccine_logistical_delay,
                                          vaccine_protection_delay = vaccine_protection_delay)
      SC1_storage_vacc_25[j, i] <- sum(!is.na(SC1_vacc_25$time_infection))
      
      SC2_vacc_25 <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                          generation_time = SC2_generation_time,
                                          spatial_kernel = spatial_kernel,
                                          t0 = 0, tf = Inf,
                                          check_final_size = check_final_size,
                                          seeding_cases = seeding_cases,
                                          prop_asymptomatic = SC2_prop_asymptomatic,
                                          prob_hosp = SC2_prob_hosp,
                                          hospitalisation_delay = SC2_hospitalisation_delay,
                                          detection_threshold = detection_threshold,
                                          vaccine_campaign_radius = spatial_kernel_25,
                                          vaccine_coverage = vaccine_coverage,
                                          vaccine_efficacy_infection = vaccine_efficacy_infection,
                                          vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                          vaccine_efficacy_disease = vaccine_efficacy_disease,
                                          vaccine_logistical_delay = vaccine_logistical_delay,
                                          vaccine_protection_delay = vaccine_protection_delay)
      SC2_storage_vacc_25[j, i] <- sum(!is.na(SC2_vacc_25$time_infection))
      
      ## 50% of spatial kernel
      SC1_vacc_50 <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                          generation_time = SC1_generation_time,
                                          spatial_kernel = spatial_kernel,
                                          t0 = 0, tf = Inf,
                                          check_final_size = check_final_size,
                                          seeding_cases = seeding_cases,
                                          prop_asymptomatic = SC1_prop_asymptomatic,
                                          prob_hosp = SC1_prob_hosp,
                                          hospitalisation_delay = SC1_hospitalisation_delay,
                                          detection_threshold = detection_threshold,
                                          vaccine_campaign_radius = spatial_kernel_50,
                                          vaccine_coverage = vaccine_coverage,
                                          vaccine_efficacy_infection = vaccine_efficacy_infection,
                                          vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                          vaccine_efficacy_disease = vaccine_efficacy_disease,
                                          vaccine_logistical_delay = vaccine_logistical_delay,
                                          vaccine_protection_delay = vaccine_protection_delay)
      SC1_storage_vacc_50[j, i] <- sum(!is.na(SC1_vacc_50$time_infection))
      
      SC2_vacc_50 <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                          generation_time = SC2_generation_time,
                                          spatial_kernel = spatial_kernel,
                                          t0 = 0, tf = Inf,
                                          check_final_size = check_final_size,
                                          seeding_cases = seeding_cases,
                                          prop_asymptomatic = SC2_prop_asymptomatic,
                                          prob_hosp = SC2_prob_hosp,
                                          hospitalisation_delay = SC2_hospitalisation_delay,
                                          detection_threshold = detection_threshold,
                                          vaccine_campaign_radius = spatial_kernel_50,
                                          vaccine_coverage = vaccine_coverage,
                                          vaccine_efficacy_infection = vaccine_efficacy_infection,
                                          vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                          vaccine_efficacy_disease = vaccine_efficacy_disease,
                                          vaccine_logistical_delay = vaccine_logistical_delay,
                                          vaccine_protection_delay = vaccine_protection_delay)
      SC2_storage_vacc_50[j, i] <- sum(!is.na(SC2_vacc_50$time_infection))
      
      ## 75% of spatial kernel
      SC1_vacc_75 <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                          generation_time = SC1_generation_time,
                                          spatial_kernel = spatial_kernel,
                                          t0 = 0, tf = Inf,
                                          check_final_size = check_final_size,
                                          seeding_cases = seeding_cases,
                                          prop_asymptomatic = SC1_prop_asymptomatic,
                                          prob_hosp = SC1_prob_hosp,
                                          hospitalisation_delay = SC1_hospitalisation_delay,
                                          detection_threshold = detection_threshold,
                                          vaccine_campaign_radius = spatial_kernel_75,
                                          vaccine_coverage = vaccine_coverage,
                                          vaccine_efficacy_infection = vaccine_efficacy_infection,
                                          vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                          vaccine_efficacy_disease = vaccine_efficacy_disease,
                                          vaccine_logistical_delay = vaccine_logistical_delay,
                                          vaccine_protection_delay = vaccine_protection_delay)
      SC1_storage_vacc_75[j, i] <- sum(!is.na(SC1_vacc_75$time_infection))
      
      SC2_vacc_75 <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                          generation_time = SC2_generation_time,
                                          spatial_kernel = spatial_kernel,
                                          t0 = 0, tf = Inf,
                                          check_final_size = check_final_size,
                                          seeding_cases = seeding_cases,
                                          prop_asymptomatic = SC2_prop_asymptomatic,
                                          prob_hosp = SC2_prob_hosp,
                                          hospitalisation_delay = SC2_hospitalisation_delay,
                                          detection_threshold = detection_threshold,
                                          vaccine_campaign_radius = spatial_kernel_75,
                                          vaccine_coverage = vaccine_coverage,
                                          vaccine_efficacy_infection = vaccine_efficacy_infection,
                                          vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                          vaccine_efficacy_disease = vaccine_efficacy_disease,
                                          vaccine_logistical_delay = vaccine_logistical_delay,
                                          vaccine_protection_delay = vaccine_protection_delay)
      SC2_storage_vacc_75[j, i] <- sum(!is.na(SC2_vacc_75$time_infection))
      
      ## 95% of spatial kernel
      SC1_vacc_95 <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                          generation_time = SC1_generation_time,
                                          spatial_kernel = spatial_kernel,
                                          t0 = 0, tf = Inf,
                                          check_final_size = check_final_size,
                                          seeding_cases = seeding_cases,
                                          prop_asymptomatic = SC1_prop_asymptomatic,
                                          prob_hosp = SC1_prob_hosp,
                                          hospitalisation_delay = SC1_hospitalisation_delay,
                                          detection_threshold = detection_threshold,
                                          vaccine_campaign_radius = spatial_kernel_95,
                                          vaccine_coverage = vaccine_coverage,
                                          vaccine_efficacy_infection = vaccine_efficacy_infection,
                                          vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                          vaccine_efficacy_disease = vaccine_efficacy_disease,
                                          vaccine_logistical_delay = vaccine_logistical_delay,
                                          vaccine_protection_delay = vaccine_protection_delay)
      SC1_storage_vacc_95[j, i] <- sum(!is.na(SC1_vacc_95$time_infection))
      
      SC2_vacc_95 <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                          generation_time = SC2_generation_time,
                                          spatial_kernel = spatial_kernel,
                                          t0 = 0, tf = Inf,
                                          check_final_size = check_final_size,
                                          seeding_cases = seeding_cases,
                                          prop_asymptomatic = SC2_prop_asymptomatic,
                                          prob_hosp = SC2_prob_hosp,
                                          hospitalisation_delay = SC2_hospitalisation_delay,
                                          detection_threshold = detection_threshold,
                                          vaccine_campaign_radius = spatial_kernel_95,
                                          vaccine_coverage = vaccine_coverage,
                                          vaccine_efficacy_infection = vaccine_efficacy_infection,
                                          vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                          vaccine_efficacy_disease = vaccine_efficacy_disease,
                                          vaccine_logistical_delay = vaccine_logistical_delay,
                                          vaccine_protection_delay = vaccine_protection_delay)
      SC2_storage_vacc_95[j, i] <- sum(!is.na(SC2_vacc_95$time_infection))
      
      ## Nothing
      SC1_no_vacc <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                          generation_time = SC1_generation_time,
                                          spatial_kernel = spatial_kernel,
                                          t0 = 0, tf = Inf,
                                          check_final_size = check_final_size,
                                          seeding_cases = seeding_cases,
                                          prop_asymptomatic = SC1_prop_asymptomatic,
                                          prob_hosp = SC1_prob_hosp,
                                          hospitalisation_delay = SC1_hospitalisation_delay,
                                          detection_threshold = detection_threshold,
                                          vaccine_campaign_radius = 0.01,
                                          vaccine_coverage = 0.01,
                                          vaccine_efficacy_infection = 0.001,
                                          vaccine_efficacy_transmission = 0.001,
                                          vaccine_efficacy_disease = 0.001,
                                          vaccine_logistical_delay = 1000,
                                          vaccine_protection_delay = 1000)
      SC1_storage_nothing[j, i] <- sum(!is.na(SC1_no_vacc$time_infection))
      
      SC2_no_vacc <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                          generation_time = SC2_generation_time,
                                          spatial_kernel = spatial_kernel,
                                          t0 = 0, tf = Inf,
                                          check_final_size = check_final_size,
                                          seeding_cases = seeding_cases,
                                          prop_asymptomatic = SC2_prop_asymptomatic,
                                          prob_hosp = SC2_prob_hosp,
                                          hospitalisation_delay = SC2_hospitalisation_delay,
                                          detection_threshold = detection_threshold,
                                          vaccine_campaign_radius = 0.01,
                                          vaccine_coverage = 0.01,
                                          vaccine_efficacy_infection = 0.001,
                                          vaccine_efficacy_transmission = 0.001,
                                          vaccine_efficacy_disease = 0.001,
                                          vaccine_logistical_delay = 1000,
                                          vaccine_protection_delay = 1000)
      SC2_storage_nothing[j, i] <- sum(!is.na(SC2_no_vacc$time_infection))
      
      
      if (j %% 25 == 0) {
        print(j)
      }
    }
    print(i)
  }
  
  ## Creating overall dataframe with all the results
  SC1_no_vacc <- data.frame(iteration = 1:iterations, scenario = "no_vaccination", pathogen = "SARS-CoV-1", SC1_storage_nothing)
  colnames(SC1_no_vacc) <- c("iteration", "scenario", "pathogen", paste0("R0=", R0_scan))
  SC2_no_vacc <- data.frame(iteration = 1:iterations, scenario = "no_vaccination", pathogen = "SARS-CoV-2", SC2_storage_nothing)
  colnames(SC2_no_vacc) <- c("iteration", "scenario", "pathogen", paste0("R0=", R0_scan))
  
  SC1_25 <- data.frame(iteration = 1:iterations, scenario = "20_spatial", pathogen = "SARS-CoV-1", SC1_storage_vacc_25)
  colnames(SC1_25) <- c("iteration", "scenario", "pathogen", paste0("R0=", R0_scan))
  SC2_25 <- data.frame(iteration = 1:iterations, scenario = "20_spatial", pathogen = "SARS-CoV-2", SC2_storage_vacc_25)
  colnames(SC2_25) <- c("iteration", "scenario", "pathogen", paste0("R0=", R0_scan))
  
  SC1_50 <- data.frame(iteration = 1:iterations, scenario = "50_spatial", pathogen = "SARS-CoV-1", SC1_storage_vacc_50)
  colnames(SC1_50) <- c("iteration", "scenario", "pathogen", paste0("R0=", R0_scan))
  SC2_50 <- data.frame(iteration = 1:iterations, scenario = "50_spatial", pathogen = "SARS-CoV-2", SC2_storage_vacc_50)
  colnames(SC2_50) <- c("iteration", "scenario", "pathogen", paste0("R0=", R0_scan))
  
  SC1_75 <- data.frame(iteration = 1:iterations, scenario = "75_spatial", pathogen = "SARS-CoV-1", SC1_storage_vacc_75)
  colnames(SC1_75) <- c("iteration", "scenario", "pathogen", paste0("R0=", R0_scan))
  SC2_75 <- data.frame(iteration = 1:iterations, scenario = "75_spatial", pathogen = "SARS-CoV-2", SC2_storage_vacc_75)
  colnames(SC2_75) <- c("iteration", "scenario", "pathogen", paste0("R0=", R0_scan))
  
  SC1_95 <- data.frame(iteration = 1:iterations, scenario = "95_spatial", pathogen = "SARS-CoV-1", SC1_storage_vacc_95)
  colnames(SC1_95) <- c("iteration", "scenario", "pathogen", paste0("R0=", R0_scan))
  SC2_95 <- data.frame(iteration = 1:iterations, scenario = "95_spatial", pathogen = "SARS-CoV-2", SC2_storage_vacc_95)
  colnames(SC2_95) <- c("iteration", "scenario", "pathogen", paste0("R0=", R0_scan))
  
  overall_bp_df <- rbind(SC1_no_vacc, SC2_no_vacc, 
                         SC1_25, SC2_25, 
                         SC1_50, SC2_50,
                         SC1_75, SC2_75,
                         SC1_95, SC2_95) %>%
    pivot_longer(cols = starts_with("R0"), names_to = "R0", values_to = "Epidemic Size") %>%
    mutate(actualR0 = as.numeric(gsub("R0=", "", R0)))
  saveRDS(object = overall_bp_df, file = "outputs/Figure1_branchingProcess_Containment/Fig1HI_branchingProcess_spatialVaccination_R0Scan.rds")
} else {
  overall_bp_df <- readRDS("outputs/Figure1_branchingProcess_Containment/Fig1HI_branchingProcess_spatialVaccination_R0Scan.rds")
}

## Plotting Figure 1B
containment_df <- overall_bp_df %>%
  mutate(contained = ifelse(`Epidemic Size` < (0.9 * check_final_size), 1, 0)) %>%
  group_by(actualR0, R0, scenario, pathogen) %>%
  summarise(proportion_contained = sum(contained) / iterations) %>%
  mutate(proportion_contained = ifelse(actualR0 == 1.00, 1, proportion_contained)) %>%
  mutate(scenario = ifelse(scenario == "no_vaccination", "zno_vaccination", scenario)) %>%
  mutate(scenario = ifelse(scenario == "20_spatial", "d20_spatial_kernelScenario", scenario)) %>%
  mutate(scenario = ifelse(scenario == "50_spatial", "c50_spatial_kernelScenario", scenario)) %>%
  mutate(scenario = ifelse(scenario == "75_spatial", "b75_spatial_kernelScenario", scenario)) %>%
  mutate(scenario = ifelse(scenario == "95_spatial", "a95_spatial_kernelScenario", scenario)) 
containment_df$scenario <- factor(containment_df$scenario, 
                                  levels = c("a95_spatial_kernelScenario",
                                             "b75_spatial_kernelScenario",
                                             "c50_spatial_kernelScenario", 
                                             "d20_spatial_kernelScenario",
                                             "zno_vaccination"))
containment_df2 <- containment_df %>%
  arrange(scenario) 
containment_plot <- ggplot(containment_df2, aes(x = actualR0, y = 100 * proportion_contained, col = scenario)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_grid(.~pathogen) + 
  scale_colour_manual(values = c("#CA2E6B", "#88C5EE", "#236897", "#13496E", "black"), 
                      labels = c("95%", "75%", "50%", "25%", "No Vaccination"),
                      name = "Vaccine\nCampaign\nCoverage",
                      guide = guide_legend(reverse = TRUE)) +
  labs(x = "R0", y = "% Outbreaks Contained") +
  theme(strip.background = element_rect(fill = "white"))
ggsave(filename = "figures/Figure_1_BranchingProcess/Fig1HI_SpatialVaccination_ContainmentPlot.pdf", plot = containment_plot, 
       width = 8, height = 3.1)

## Vaccination-Related Sensitivity Analyses Heatmaps
fresh_run_vaccination_heatmaps <- TRUE
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
  saveRDS(object = storage_R0_TgRatio_df, file = "outputs/Figure1_branchingProcess_Containment/branchingProcess_R0_TgRatio_scan.rds")
  
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
  saveRDS(object = storage_R0_efficacy_df, file = "outputs/Figure1_branchingProcess_Containment/branchingProcess_R0_Efficacy_scan.rds")
  
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
  saveRDS(object = storage_R0_preSymp_df, file = "outputs/Figure1_branchingProcess_Containment/branchingProcess_R0_preSymp_scan.rds")
  
} else {
  storage_R0_TgRatio_df <- readRDS("outputs/Figure1_branchingProcess_Containment/branchingProcess_R0_TgRatio_scan.rds")
  storage_R0_efficacy_df <- readRDS("outputs/Figure1_branchingProcess_Containment/branchingProcess_R0_Efficacy_scan.rds")
  storage_R0_preSymp_df <- readRDS("outputs/Figure1_branchingProcess_Containment/branchingProcess_R0_preSymp_scan.rds")
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
ggsave(file = "figures/Figure_1_BranchingProcess/Figure_1C_TgRatio_Sensitivity_Analysis.pdf", plot = R0_TgRatio_plot, width = 2.4, height = 2.4)

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
ggsave(filename = "figures/Figure_1_BranchingProcess/Figure_1D_VaccineEfficacy_Sensitivity_Analysis.pdf", plot = R0_vaccine_efficacy_plot, width = 2.4, height = 2.4)

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
ggsave(filename = "figures/Figure_1_BranchingProcess/Figure_1E_preSymp_Sensitivity_Analysis.pdf", 
       plot = R0_preSymp_plot, width = 2.4, height = 2.4)

## Heatmaps legend
legend <- R0_TgRatio_plot + theme(legend.position = "bottom")
ggsave(file = "figures/Figure_1_BranchingProcess/Legend_vaccinationHeatmap.pdf", plot = legend, width = 2.4, height = 2.4)


############################ COME BACK TO THIS #######################

## Sensitivity Analysis - R0 vs Vaccine Efficacy MONOCLONAL
R0_seq <- R0_scan[-c(1, 2)]
vaccine_efficacy_seq <- seq(0.1, 1, 0.1)
iterations <- 50
storage_R0_efficacy_sensitivity <- array(data = NA, dim = c(length(R0_seq), length(vaccine_efficacy_seq), iterations))
for (i in 1:length(R0_seq)) {
  for (j in 1:length(vaccine_efficacy_seq)) {
    for (k in 1:iterations) {
      bp_out <- chain_sim_susc_ring_vacc2(offspring = "pois",
                                          mn_offspring = R0_seq[i],
                                          generation_time = generation_time,
                                          t0 = 0, tf = Inf, pop = 10^7, 
                                          check_final_size = check_final_size, initial_immune = 0,
                                          seeding_cases = 5, prop_asymptomatic = 0.2,
                                          infection_to_onset = infection_to_onset,
                                          vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                          vaccine_efficacy_infection = vaccine_efficacy_seq[j],
                                          vaccine_efficacy_transmission = vaccine_efficacy_seq[j],
                                          vaccine_logistical_delay = 2,
                                          vaccine_protection_delay = 0)
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
saveRDS(object = storage_R0_efficacy_df, file = "outputs/branching_process_outputs/branchingProcess_R0_Efficacy_scan_Monoclonal.rds")

storage_R0_efficacy_df <- readRDS("outputs/branching_process_outputs/branchingProcess_R0_Efficacy_scan_Monoclonal.rds")
R0_vaccine_efficacy_plot <- ggplot(storage_R0_efficacy_df, aes(x = R0, y = vaccine_efficacy, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 1), begin = 0.25, end = 1, 
                       name = "Proportion\nContained",
                       direction = -1) +
  scale_fill_distiller(palette = "PuRd", direction = 1,
                       name = "Proportion\nContained") +
  labs(x = "R0",
       y = "Monoclonal Efficacy") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none",
        panel.border = element_rect(linetype = "solid", color = "black", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

