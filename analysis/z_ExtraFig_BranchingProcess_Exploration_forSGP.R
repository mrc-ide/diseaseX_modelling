# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/branching_process.R"))
source(here::here("functions/branching_process2.R"))

## Currently missing:
### Hospitalisations as a trigger for implementing the vaccine (i.e. detection and initiation not based on calendar time but case count)
### Vaccine efficacy waning

### Fixed model parameters
generation_time <- function(n) { rgamma(n, shape = 13.5, rate = 2) } # 6.75 day generation time Gamam distributed (as per Walker et al, Science, 2020)
infection_to_onset <- function(n) { rgamma(n, shape = 13.5/3, rate = 2) } # ~35% of transmission presymptomatic (per SARS-CoV-2, slightly lower than but roughly aligned with: https://bmjopen.bmj.com/content/11/6/e041240)
vaccine_start <- 21 # detection + logistical delay of 3 weeks to vaccination initiation (need to get implementation based on cases not calendar time in)
vaccine_coverage <- 0.8
pop <- 10^10
check_final_size <- 1000
initial_immune <- 0
seeding_cases <- 5
prop_asymptomatic <- 0.15
vaccine_efficacy_infection <- 0.75
vaccine_efficacy_transmission <- 0.55
vaccine_logistical_delay <- 2
R0_scan <- c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5)
iterations <- 250

## R0 sensitivity analysis (Figure 1B)
fresh_run_R0_sensitivity_analysis <- FALSE
if (fresh_run_R0_sensitivity_analysis) {
  
  ## Simulation parameters 
  # hist(generation_time(10^5))
  # hist(infection_to_onset(10^5))
  # median(generation_time(10^5) / infection_to_onset(10^5))
  
  ### Setting up R0 scan and the storage for each of the different protection delays
  storage_nothing <- matrix(nrow = iterations, ncol = length(R0_scan))
  storage_vacc_2weeks <- matrix(nrow = iterations, ncol = length(R0_scan))
  storage_vacc_1week <- matrix(nrow = iterations, ncol = length(R0_scan))
  storage_vacc_Halfweek <- matrix(nrow = iterations, ncol = length(R0_scan))
  storage_vacc_2days <- matrix(nrow = iterations, ncol = length(R0_scan))
  storage_monoclonal <- matrix(nrow = iterations, ncol = length(R0_scan))
  
  for (i in 1:length(R0_scan)) {
    for (j in 1:iterations) {
      
      vacc_2weeks <- chain_sim_susc_ring_vacc2(offspring = "pois",
                                               mn_offspring = R0_scan[i],
                                               generation_time = generation_time,
                                               t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, initial_immune = initial_immune,
                                               seeding_cases = seeding_cases, prop_asymptomatic = prop_asymptomatic,
                                               infection_to_onset = infection_to_onset,
                                               vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                               vaccine_efficacy_infection = vaccine_efficacy_infection,
                                               vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                               vaccine_logistical_delay = vaccine_logistical_delay,
                                               vaccine_protection_delay = 14)
      storage_vacc_2weeks[j, i] <- sum(!is.na(vacc_2weeks$time_infection))
      
      vacc_1week <- chain_sim_susc_ring_vacc2(offspring = "pois",
                                              mn_offspring = R0_scan[i],
                                              generation_time = generation_time,
                                              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, initial_immune = initial_immune,
                                              seeding_cases = seeding_cases, prop_asymptomatic = prop_asymptomatic,
                                              infection_to_onset = infection_to_onset,
                                              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                              vaccine_efficacy_infection = vaccine_efficacy_infection,
                                              vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                              vaccine_logistical_delay = vaccine_logistical_delay,
                                              vaccine_protection_delay = 7)
      storage_vacc_1week[j, i] <- sum(!is.na(vacc_1week$time_infection))
      
      vacc_Halfweek <- chain_sim_susc_ring_vacc2(offspring = "pois",
                                                 mn_offspring = R0_scan[i],
                                                 generation_time = generation_time,
                                                 t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, initial_immune = initial_immune,
                                                 seeding_cases = seeding_cases, prop_asymptomatic = prop_asymptomatic,
                                                 infection_to_onset = infection_to_onset,
                                                 vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                                 vaccine_efficacy_infection = vaccine_efficacy_infection,
                                                 vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                                 vaccine_logistical_delay = vaccine_logistical_delay,
                                                 vaccine_protection_delay = 3.5)
      storage_vacc_Halfweek[j, i] <- sum(!is.na(vacc_Halfweek$time_infection))
      
      vacc_2days <- chain_sim_susc_ring_vacc2(offspring = "pois",
                                              mn_offspring = R0_scan[i],
                                              generation_time = generation_time,
                                              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, initial_immune = initial_immune,
                                              seeding_cases = seeding_cases, prop_asymptomatic = prop_asymptomatic,
                                              infection_to_onset = infection_to_onset,
                                              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                              vaccine_efficacy_infection = vaccine_efficacy_infection,
                                              vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                              vaccine_logistical_delay = vaccine_logistical_delay,
                                              vaccine_protection_delay = 2)
      storage_vacc_2days[j, i] <- sum(!is.na(vacc_2days$time_infection))
      
      monoclonal <- chain_sim_susc_ring_vacc2(offspring = "pois",
                                              mn_offspring = R0_scan[i],
                                              generation_time = generation_time,
                                              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, initial_immune = initial_immune,
                                              seeding_cases = seeding_cases, prop_asymptomatic = prop_asymptomatic,
                                              infection_to_onset = infection_to_onset,
                                              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                              vaccine_efficacy_infection = vaccine_efficacy_infection,
                                              vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                              vaccine_logistical_delay = vaccine_logistical_delay,
                                              vaccine_protection_delay = 0)
      storage_monoclonal[j, i] <- sum(!is.na(monoclonal$time_infection))
      
      no_vacc <- chain_sim_susc_ring_vacc2(offspring = "pois",
                                           mn_offspring = R0_scan[i],
                                           generation_time = generation_time,
                                           t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, initial_immune = initial_immune,
                                           seeding_cases = seeding_cases, prop_asymptomatic = prop_asymptomatic,
                                           infection_to_onset = infection_to_onset,
                                           vaccine_start = 1000, vaccine_coverage = vaccine_coverage,
                                           vaccine_efficacy_infection = vaccine_efficacy_infection,
                                           vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                           vaccine_logistical_delay = vaccine_logistical_delay,
                                           vaccine_protection_delay = 100)
      storage_nothing[j, i] <- sum(!is.na(no_vacc$time_infection))
      if (j %% 25 == 0) {
        print(j)
      }
    }
    print(i)
  }
  
  ## Creating overall dataframe with all the results
  no_vacc <- data.frame(iteration = 1:iterations, scenario = "no_vaccination", storage_nothing)
  colnames(no_vacc) <- c("iteration", "scenario", paste0("R0=", R0_scan))
  
  vacc_2weeks <- data.frame(iteration = 1:iterations, scenario = "vacc_2weeks_protectDelay", storage_vacc_2weeks)
  colnames(vacc_2weeks) <- c("iteration", "scenario", paste0("R0=", R0_scan))
  
  vacc_1week <- data.frame(iteration = 1:iterations, scenario = "vacc_1week_protectDelay", storage_vacc_1week)
  colnames(vacc_1week) <- c("iteration", "scenario", paste0("R0=", R0_scan))
  
  vacc_Halfweek <- data.frame(iteration = 1:iterations, scenario = "vacc_Halfweek_protectDelay", storage_vacc_Halfweek)
  colnames(vacc_Halfweek) <- c("iteration", "scenario", paste0("R0=", R0_scan))
  
  vacc_2days <- data.frame(iteration = 1:iterations, scenario = "vacc_2days_protectDelay", storage_vacc_2days)
  colnames(vacc_2days) <- c("iteration", "scenario", paste0("R0=", R0_scan))
  
  monoclonal <- data.frame(iteration = 1:iterations, scenario = "monoclonal", storage_monoclonal)
  colnames(monoclonal) <- c("iteration", "scenario", paste0("R0=", R0_scan))
  
  overall_bp_df <- rbind(no_vacc, vacc_2weeks, vacc_1week, vacc_Halfweek, vacc_2days, monoclonal) %>%
    pivot_longer(cols = starts_with("R0"), names_to = "R0", values_to = "Epidemic Size") %>%
    mutate(actualR0 = as.numeric(gsub("R0=", "", R0)))
  saveRDS(object = overall_bp_df, file = "outputs/zExtra_SGP_Outputs/branchingProcess_R0_Scan.rds")
} else {
  overall_bp_df <- readRDS("outputs/zExtra_SGP_Outputs/branchingProcess_R0_Scan.rds")
}

## Plotting Figure 1B
containment_df <- overall_bp_df %>%
  mutate(contained = ifelse(`Epidemic Size` < (0.9 * check_final_size), 1, 0)) %>%
  group_by(actualR0, R0, scenario) %>%
  summarise(proportion_contained = sum(contained) / iterations) %>%
  mutate(proportion_contained = ifelse(actualR0 == 1.00, 1, proportion_contained)) %>%
  mutate(scenario = ifelse(scenario == "no_vaccination", "zno_vaccination", scenario)) %>%
  mutate(scenario = ifelse(scenario == "monoclonal", "a_monoclonal", scenario)) %>%
  mutate(scenario = ifelse(scenario == "vacc_2days_protectDelay", "bvacc_2days_protectDelay", scenario)) %>%
  mutate(scenario = ifelse(scenario == "vacc_Halfweek_protectDelay", "cvacc_Halfweek_protectDelay", scenario)) %>%
  mutate(scenario = ifelse(scenario == "vacc_1week_protectDelay", "dvacc_1week_protectDelay", scenario)) %>%
  mutate(scenario = ifelse(scenario == "vacc_2weeks_protectDelay", "evacc_2weeks_protectDelay", scenario)) 
containment_df$scenario <- factor(containment_df$scenario, 
                                  levels = c("a_monoclonal", 
                                             "bvacc_2days_protectDelay",
                                             "cvacc_Halfweek_protectDelay", 
                                             "dvacc_1week_protectDelay", 
                                             "evacc_2weeks_protectDelay",
                                             "zno_vaccination"))
containment_df2 <- containment_df %>%
  arrange(scenario) 
containment_plot <- ggplot(containment_df2, aes(x = actualR0, y = 100 * proportion_contained, col = scenario)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_colour_manual(values = c("#CA2E6B", "#88C5EE", "#3E9BDB", "#236897", "#13496E", "black"), 
                      labels = c("Immediate", "2 Days", "3.5 Days","1 Week", "2 Weeks", "No Vaccination"),
                      name = "Vaccine\nProtection\nDelay",
                      guide = guide_legend(reverse = TRUE)) +
  labs(x = "R0", y = "% Outbreaks Contained")
ggsave(filename = "figures/zExtra_SGP_Figures/B_SGP_ContainmentPlot.pdf", plot = containment_plot, 
       width = 8, height = 3.1)

#########################################################
## Vaccination-Related Sensitivity Analyses Heatmaps
#########################################################
fresh_run_vaccination_heatmaps <- FALSE
generation_time <- function(n) { rgamma(n, shape = 13.5, rate = 2) } # 6.75 day generation time Gamam distributed (as per Walker et al, Science, 2020)
infection_to_onset <- function(n) { rgamma(n, shape = 13.5/3, rate = 2) } # keeping proportion of presymptomatic transmission constant as Tg varies
if (fresh_run_vaccination_heatmaps) {
  
  ## R0 sensitivity analysis for all the runs
  R0_seq <- R0_scan[R0_scan > 1]
  vaccine_protection_delay <- 7
  
  ## Sensitivity Analysis - R0 vs Vaccine Efficacy
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
  saveRDS(object = storage_R0_efficacy_df, file = "outputs/zExtra_SGP_Outputs/C_SGP_vaccination_R0_Efficacy_scan.rds")
  
  ## Sensitivity Analysis - R0 vs Protection Delay
  vaccine_protection_delay_seq <- seq(0, 14, 2)
  storage_R0_delay_sensitivity <- array(data = NA, dim = c(length(R0_seq), length(vaccine_protection_delay_seq), iterations))
  for (i in 1:length(R0_seq)) {
    for (j in 1:length(vaccine_protection_delay_seq)) {
      for (k in 1:iterations) {
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
                                            vaccine_protection_delay = vaccine_protection_delay_seq[j])
        storage_R0_delay_sensitivity[i, j, k] <- sum(!is.na(bp_out$time_infection))
      }
    }
    print(i)
  }
  storage_R0_delay_df <- reshape2::melt(storage_R0_delay_sensitivity)
  storage_R0_delay_df$R0 <- rep(R0_seq, length(vaccine_protection_delay_seq) * iterations)
  storage_R0_delay_df$protection_delay <- rep(vaccine_protection_delay_seq, each = length(R0_seq), times = iterations)
  storage_R0_delay_df$iteration <- rep(1:iterations, each = length(R0_seq) * length(vaccine_protection_delay_seq))
  storage_R0_delay_df <- storage_R0_delay_df %>%
    select(-(Var1:Var3)) %>%
    rename(final_size = value) %>%
    mutate(contained = ifelse(final_size < (0.9 * check_final_size), 1, 0)) %>%
    group_by(R0, protection_delay) %>%
    summarise(proportion_contained = sum(contained) / n())
  saveRDS(object = storage_R0_delay_df, file = "outputs/zExtra_SGP_Outputs/D_SGP_vaccination_R0_protectionDelay_scan.rds")
  
  ## Sensitivity Analysis - R0 vs Presymptomatic Transmission Proportion
  proportion_presymptomatic_seq <- seq(0.1, 0.7, 0.1)
  storage_R0_preSymp_sensitivity <- array(data = NA, dim = c(length(R0_seq), length(proportion_presymptomatic_seq), iterations))
  for (i in 1:length(R0_seq)) {
    for (j in 1:length(proportion_presymptomatic_seq)) {
      for (k in 1:iterations) {
        infection_to_onset <- function(n) { rgamma(n, shape = 13.5 * proportion_presymptomatic_seq[j], rate = 2) }
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
  saveRDS(object = storage_R0_preSymp_df, file = "outputs/zExtra_SGP_Outputs/E_SGP_vaccination_R0_preSymp_scan.rds")
  
} else {
  storage_R0_efficacy_df <- readRDS("outputs/zExtra_SGP_Outputs/C_SGP_vaccination_R0_Efficacy_scan.rds")
  storage_R0_delay_df <- readRDS("outputs/zExtra_SGP_Outputs/D_SGP_vaccination_R0_protectionDelay_scan.rds")
  storage_R0_preSymp_df <- readRDS("outputs/zExtra_SGP_Outputs/E_SGP_vaccination_R0_preSymp_scan.rds")
}

## Creating Vaccination-Related Heatmaps

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
ggsave(filename = "figures/zExtra_SGP_Figures/C_SGP_vaccination_R0_Efficacy_plot.pdf", plot = R0_vaccine_efficacy_plot, width = 2.4, height = 2.4)

### Vaccine Efficacy
R0_vaccine_protectionDelay_plot <- ggplot(storage_R0_delay_df, aes(x = R0, y = protection_delay, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0, 1), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  labs(x = "R0",
       y = "Protection Delay (Days)") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none",
        panel.border = element_rect(linetype = "solid", color = "black", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)
ggsave(filename = "figures/zExtra_SGP_Figures/D_SGP_vaccination_R0_protectionDelay__plot.pdf", plot = R0_vaccine_protectionDelay_plot, width = 2.4, height = 2.4)

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
ggsave(filename = "figures/zExtra_SGP_Figures/E_SGP_vaccination_R0_preSymp_plot.pdf", 
       plot = R0_preSymp_plot, width = 2.4, height = 2.4)

## Heatmaps legend
vaccine_legend <- R0_preSymp_plot + theme(legend.position = "right")
ggsave(file = "figures/zExtra_SGP_Figures/Legend_vaccinationHeatmap.pdf", plot = vaccine_legend, width = 2.4, height = 2.4)


#########################################################
## Monoclonal-Related Sensitivity Analyses Heatmaps
#########################################################
fresh_run_monoclonal_heatmaps <- FALSE
generation_time <- function(n) { rgamma(n, shape = 13.5, rate = 2) } # 6.75 day generation time Gamam distributed (as per Walker et al, Science, 2020)
infection_to_onset <- function(n) { rgamma(n, shape = 13.5/3, rate = 2) } # keeping proportion of presymptomatic transmission constant as Tg varies
if (fresh_run_monoclonal_heatmaps) {
  
  ## R0 sensitivity analysis for all the runs
  R0_seq <- R0_scan[R0_scan > 1]
  vaccine_protection_delay <- 0
  
  ## Sensitivity Analysis - R0 vs Monoclonal Efficacy
  monoclonal_efficacy_seq <- seq(0.3, 0.9, 0.1)
  storage_R0_monoclonal_efficacy_sensitivity <- array(data = NA, dim = c(length(R0_seq), length(monoclonal_efficacy_seq), iterations))
  for (i in 1:length(R0_seq)) {
    for (j in 1:length(monoclonal_efficacy_seq)) {
      for (k in 1:iterations) {
        bp_out <- chain_sim_susc_ring_vacc2(offspring = "pois",
                                            mn_offspring = R0_seq[i],
                                            generation_time = generation_time,
                                            t0 = 0, tf = Inf, pop = pop,
                                            check_final_size = check_final_size, initial_immune = initial_immune,
                                            seeding_cases = seeding_cases, prop_asymptomatic = prop_asymptomatic,
                                            infection_to_onset = infection_to_onset,
                                            vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                            vaccine_efficacy_infection = monoclonal_efficacy_seq[j],
                                            vaccine_efficacy_transmission = monoclonal_efficacy_seq[j],
                                            vaccine_logistical_delay = vaccine_logistical_delay,
                                            vaccine_protection_delay = vaccine_protection_delay)
        storage_R0_monoclonal_efficacy_sensitivity[i, j, k] <- sum(!is.na(bp_out$time_infection))
      }
    }
    print(i)
  }
  storage_R0_monoclonal_efficacy_df <- reshape2::melt(storage_R0_monoclonal_efficacy_sensitivity)
  storage_R0_monoclonal_efficacy_df$R0 <- rep(R0_seq, length(monoclonal_efficacy_seq) * iterations)
  storage_R0_monoclonal_efficacy_df$monoclonal_efficacy <- rep(monoclonal_efficacy_seq, each = length(R0_seq), times = iterations)
  storage_R0_monoclonal_efficacy_df$iteration <- rep(1:iterations, each = length(R0_seq) * length(monoclonal_efficacy_seq))
  storage_R0_monoclonal_efficacy_df <- storage_R0_monoclonal_efficacy_df %>%
    select(-(Var1:Var3)) %>%
    rename(final_size = value) %>%
    mutate(contained = ifelse(final_size < (0.9 * check_final_size), 1, 0)) %>%
    group_by(R0, monoclonal_efficacy) %>%
    summarise(proportion_contained = sum(contained) / n())
  saveRDS(object = storage_R0_monoclonal_efficacy_df, file = "outputs/zExtra_SGP_Outputs/F_SGP_monoclonal_R0_Efficacy_scan.rds")
  
  ## Sensitivity Analysis - R0 vs Presymptomatic Transmission Proportion
  proportion_presymptomatic_seq <- seq(0.1, 0.7, 0.1)
  storage_R0_preSymp_monoclonal_sensitivity <- array(data = NA, dim = c(length(R0_seq), length(proportion_presymptomatic_seq), iterations))
  for (i in 1:length(R0_seq)) {
    for (j in 1:length(proportion_presymptomatic_seq)) {
      for (k in 1:iterations) {
        infection_to_onset <- function(n) { rgamma(n, shape = 13.5 * proportion_presymptomatic_seq[j], rate = 2) }
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
        storage_R0_preSymp_monoclonal_sensitivity[i, j, k] <- sum(!is.na(bp_out$time_infection))
      }
    }
    print(i)
  }
  storage_R0_preSymp_monoclonal_df <- reshape2::melt(storage_R0_preSymp_monoclonal_sensitivity)
  storage_R0_preSymp_monoclonal_df$R0 <- rep(R0_seq, length(proportion_presymptomatic_seq) * iterations)
  storage_R0_preSymp_monoclonal_df$prop_preSymp <- rep(proportion_presymptomatic_seq, each = length(R0_seq), times = iterations)
  storage_R0_preSymp_monoclonal_df$iteration <- rep(1:iterations, each = length(R0_seq) * length(proportion_presymptomatic_seq))
  storage_R0_preSymp_monoclonal_df <- storage_R0_preSymp_monoclonal_df %>%
    select(-(Var1:Var3)) %>%
    rename(final_size = value) %>%
    mutate(contained = ifelse(final_size < (0.9 * check_final_size), 1, 0)) %>%
    group_by(R0, prop_preSymp) %>%
    summarise(proportion_contained = sum(contained) / n())
  saveRDS(object = storage_R0_preSymp_monoclonal_df, file = "outputs/zExtra_SGP_Outputs/G_SGP_monoclonal_R0_preSymp_scan.rds")
  
} else {
  storage_R0_monoclonal_efficacy_df <- readRDS("outputs/zExtra_SGP_Outputs/F_SGP_monoclonal_R0_Efficacy_scan.rds")
  storage_R0_preSymp_monoclonal_df <- readRDS("outputs/zExtra_SGP_Outputs/G_SGP_monoclonal_R0_preSymp_scan.rds")
}

## Creating Monoclonal-Related Heatmaps

### Monoclonal Efficacy
R0_monoclonal_efficacy_plot <- ggplot(storage_R0_monoclonal_efficacy_df, aes(x = R0, y = 100 * monoclonal_efficacy, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 1), begin = 0.25, end = 1, 
                       name = "Proportion\nContained",
                       direction = -1) +
  labs(x = "R0",
       y = "Monoclonal Efficacy (%)") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none",
        panel.border = element_rect(linetype = "solid", color = "black", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)
ggsave(filename = "figures/zExtra_SGP_Figures/F_SGP_monoclonal_R0_Efficacy_plot.pdf", plot = R0_monoclonal_efficacy_plot, width = 2.4, height = 2.4)

### Proportion of Presymptomatic Transmission
R0_monoclonal_preSymp_plot <- ggplot(storage_R0_preSymp_monoclonal_df, aes(x = R0, y = 100 * prop_preSymp, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 1), begin = 0.25, end = 1, 
                       name = "Proportion\nContained",
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
ggsave(filename = "figures/zExtra_SGP_Figures/G_SGP_monoclonal_R0_preSymp_plot.pdf", plot = R0_monoclonal_preSymp_plot, width = 2.4, height = 2.4)

## Heatmaps legend
monoclonal_legend <- R0_monoclonal_preSymp_plot + theme(legend.position = "right")
ggsave(file = "figures/zExtra_SGP_Figures/Legend_monoclonalHeatmap.pdf", plot = monoclonal_legend, width = 2.4, height = 2.4)

