# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/branching_process_spatial_vaccination.R"))

### Epi parameters
# SC1_generation_time <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
# SC1_infection_to_onset <- function(n) { rgamma(n, shape = 0.1, rate = 1) } ## Ask Azra for values (negligible assumed currently)
# SC1_prop_asymptomatic <- 0
# SC1_prob_hosp <- 0.95 # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
# SC1_hospitalisation_delay <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/

SC1_generation_time <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
SC1_infection_to_onset <- function(n) { rgamma(n, shape = 0.1, rate = 1) } ## Ask Azra for values (negligible assumed currently)
SC1_prop_asymptomatic <- 0
SC1_prob_hosp <- 1 # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
SC1_hospitalisation_delay <- function(n) { 2 } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/

### Spatial kernel parameters
mu <- 1
size <- 5
spatial_kernel <- function(n) { rnbinom(n, size = 5, mu = mu) }
spatial_kernel <- function(n) { 0 }
vaccination_radius <- 1000 * mu

### Vaccine related parameters
vaccine_coverage <- 1
vaccine_efficacy_infection <- 0.5 
vaccine_efficacy_transmission <- 0.5 
vaccine_efficacy_disease <- 0
vaccine_logistical_delay <- 2
vaccine_protection_delay <- 2

### Other parameters
pop <- 10^10
check_final_size <- 250
initial_immune <- 0
seeding_cases <- 3

set.seed(10)
test <- spatial_bp_geog_vacc(mn_offspring = 6,
                             generation_time = SC1_generation_time,
                             spatial_kernel = spatial_kernel,
                             t0 = 0, tf = Inf,
                             check_final_size = check_final_size,
                             seeding_cases = seeding_cases,
                             prop_asymptomatic = SC1_prop_asymptomatic,
                             prob_hosp = SC1_prob_hosp,
                             hospitalisation_delay = SC1_hospitalisation_delay,
                             detection_threshold = 5,
                             vaccine_campaign_radius = vaccination_radius,
                             vaccine_coverage = vaccine_coverage,
                             vaccine_efficacy_infection = vaccine_efficacy_infection,
                             vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                             vaccine_efficacy_disease = vaccine_efficacy_disease,
                             vaccine_logistical_delay = vaccine_logistical_delay,
                             vaccine_protection_delay = vaccine_protection_delay)
sum(!is.na(test$time_infection))

x <- test %>%
  select(ancestor, n_offspring, n_offspring_new, n_offspring_new_new, hospitalised, time_infection, 
         vaccinated, time_vaccinated, vaccinated_before_infection, time_protected, protected_before_infection)
table(x$vaccinated)

hist(x$n_offspring)
hist(x$n_offspring_new)
hist(x$n_offspring_new_new)

## note that vax campaign can never be triggered until all infections are generated from index case.

mn_offspring = 6
generation_time = SC1_generation_time
spatial_kernel = spatial_kernel
t0 = 0
tf = Inf
check_final_size = check_final_size
seeding_cases = seeding_cases
prop_asymptomatic = SC1_prop_asymptomatic
prob_hosp = SC1_prob_hosp
hospitalisation_delay = SC1_hospitalisation_delay
detection_threshold = 5
vaccine_campaign_radius = vaccination_radius
vaccine_coverage = vaccine_coverage
vaccine_efficacy_infection = vaccine_efficacy_infection
vaccine_efficacy_transmission = vaccine_efficacy_transmission
vaccine_efficacy_disease = vaccine_efficacy_disease
vaccine_logistical_delay = vaccine_logistical_delay
vaccine_protection_delay = vaccine_protection_delay



### Sensitivity analysis parameters
R0_scan <- c(0.75, 1, 1.5)
surveillance_scan <- c(1, 10)
iterations <- 15

storage <- array(data = NA, dim = c(iterations, length(R0_scan), length(surveillance_scan)))
for (i in 1:length(R0_scan)) {
  for (j in 1:length(surveillance_scan)) {
    for (k in 1:iterations) {
      
      temp <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                   generation_time = SC1_generation_time,
                                   spatial_kernel = spatial_kernel,
                                   t0 = 0, tf = Inf,
                                   check_final_size = check_final_size,
                                   seeding_cases = seeding_cases,
                                   prop_asymptomatic = SC1_prop_asymptomatic,
                                   prob_hosp = SC1_prob_hosp,
                                   hospitalisation_delay = SC1_hospitalisation_delay,
                                   detection_threshold = surveillance_scan[j],
                                   vaccine_campaign_radius = vaccination_radius,
                                   vaccine_coverage = vaccine_coverage,
                                   vaccine_efficacy_infection = vaccine_efficacy_infection,
                                   vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                   vaccine_efficacy_disease = vaccine_efficacy_disease,
                                   vaccine_logistical_delay = vaccine_logistical_delay,
                                   vaccine_protection_delay = vaccine_protection_delay)
      storage[k, i, j] <- sum(!is.na(temp$time_infection))
  
      print(paste0("i = ", i, ", j = ", j, ", k = ", k))
    }
  }
}


x <- reshape2::melt(storage)
colnames(x) <- c("iteration", "R0", "surveillance", "outbreak_size")
test <- x %>%
  mutate(contained = ifelse(outbreak_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(R0, surveillance) %>%
  summarise(proportion_contained = sum(contained) / iterations)

ggplot(test, aes(x = R0, y = proportion_contained, col = factor(surveillance))) +
  geom_line()
  
## Creating overall dataframe with all the results
SC1_no_vacc <- data.frame(iteration = 1:iterations, scenario = "no_vaccination", pathogen = "SARS-CoV-1", SC1_storage_nothing)
colnames(SC1_no_vacc) <- c("iteration", "scenario", "pathogen", paste0("R0=", R0_scan))
SC2_no_vacc <- data.frame(iteration = 1:iterations, scenario = "no_vaccination", pathogen = "SARS-CoV-2", SC2_storage_nothing)
colnames(SC2_no_vacc) <- c("iteration", "scenario", "pathogen", paste0("R0=", R0_scan))


overall_bp_df <- rbind(SC1_no_vacc, SC2_no_vacc, 
                       SC1_25, SC2_25, 
                       SC1_50, SC2_50,
                       SC1_75, SC2_75,
                       SC1_95, SC2_95) %>%
  pivot_longer(cols = starts_with("R0"), names_to = "R0", values_to = "Epidemic Size") %>%
  mutate(actualR0 = as.numeric(gsub("R0=", "", R0)))


## Plotting Figure 1B
containment_df <-  
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
ggsave(filename = "figures/Figure_1_BranchingProcess/NEW_Figure_1B_ContainmentPlot_SpatialVaccination.pdf", plot = containment_plot, 
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

