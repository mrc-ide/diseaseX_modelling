# Load required libraries
source(here::here("main.R"))
library(parallel)
library(tictoc)

# Load required functions
source(here::here("functions/spatial_vax_simulation.R"))
source(here::here("functions/implement_quarantine.R"))
source(here::here("functions/helper_functions.R"))
source(here::here("functions/time_to_nth_infection.R"))

### Fixed Model Parameters

### Quarantine time (based on Kucharksi et al: https://pmc.ncbi.nlm.nih.gov/articles/PMC7511527/?utm_source=chatgpt.com#cesec10)
days <- 1:6
p <- c(0.20, 0.20, 0.20, 0.20, 0.20, 0.20)
sampled_quarantine_times <- sample(days, size = 1000, replace = TRUE, prob = p)
fit_gamma_fdplus <- fitdist(sampled_quarantine_times, "gamma")
quarantine_time_closure <- function(quarantine_time_shape, quarantine_time_rate) {
  function(n) {
    rgamma(n, shape = quarantine_time_shape, rate = quarantine_time_rate)
  }
}
quarantine_time <- quarantine_time_closure(quarantine_time_shape = 3.342, #eval(fit_gamma_fdplus$estimate["shape"]
                                           quarantine_time_rate = 0.951) # eval(fit_gamma_fdplus$estimate["rate"])
prob_quarantine_contact_traced <- 0.47 # from the same article as above - ~53% traced (via app-like) * 90% adhering. 
                                       # assuming tracing is via app-like hence this is lower than vaccine coverage (also contingent on tracing)
prob_quarantine_symptoms <- 0.9        # from the same article as above 

### SC1 parameters
SC1_generation_time <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
SC1_infection_to_onset <- function(n) { rgamma(n, shape = 0.1, rate = 1) } # (negligible assumed currently)
SC1_prop_asymptomatic <- 0
SC1_prob_hosp <- 0.9 # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
SC1_hospitalisation_delay <- function(n) { 12 } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/

### SC2 parameters
SC2_generation_time <- function(n) { rgamma(n, shape = 13.5, rate = 2) } # 6.75 day generation time Gamam distributed (as per Walker et al, Science, 2020)
SC2_infection_to_onset <- function(n) { rgamma(n, shape = 13.5/3, rate = 2) } # ~35% of transmission presymptomatic (per SARS-CoV-2, slightly lower than but roughly aligned with: https://bmjopen.bmj.com/content/11/6/e041240)
SC2_prop_asymptomatic <- 0.15
SC2_prob_hosp <- 0.05
SC2_hospitalisation_delay <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/

### Spatial kernel parameters
mu <- 10
size <- 4
spatial_kernel <- function(n) { rnbinom(n, size = 4, mu = 10) }

### Vaccine related parameters
vaccine_coverage <- 0.8
vaccine_efficacy_disease <- 0.95
vaccine_logistical_delay <- 2
vaccine_protection_delay <- 7

### Other parameters
pop <- 10^10
check_final_size <- 2500
initial_immune <- 0
seeding_cases <- 3
iterations <- 100

#########################################################################
## R0 sensitivity analysis (Figure 2B)
#########################################################################
R0_scan <- c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5)
surveillance_scan <- c(1, 10, 25, 50, 75)
vaccine_efficacy_infection_scan <- c(0.35, 0.75) 
vaccine_efficacy_transmission_scan <- c(0.35, 0.5) 
spatial_ratio_scan <- c(50, 100)
quarantine_efficacy_scan <- c(0, 0.35, 0.65) # from the same article as above
length(R0_scan) * length(surveillance_scan) * length(vaccine_efficacy_infection_scan) * length(spatial_ratio_scan) * length(quarantine_efficacy_scan) * (50/ (60 * 60)) 

fresh_run_R0_sensitivity_analysis <- FALSE
tic()
n <- 2000
if (fresh_run_R0_sensitivity_analysis) {
  
  ## Setting up the cluster to support the parallel runs
  no_cores <- min(iterations, 10)
  cl <- makeCluster(no_cores)
  set.seed(2000)
  seeds <- runif(n = iterations, min = 1, max = 10^9)
  clusterExport(cl, list("mu", "R0_scan", "SC1_generation_time", "spatial_kernel", "calculate_Reff", "calculate_R0",
                         "check_final_size", "seeding_cases", "SC1_prop_asymptomatic", "time_to_nth_infection",
                         "SC1_prob_hosp", "SC1_hospitalisation_delay", "surveillance_scan", "n",
                         "SC1_infection_to_onset", "SC2_infection_to_onset", "implement_quarantine",
                         "vaccine_coverage", "vaccine_efficacy_infection_scan", "vaccine_efficacy_transmission_scan",
                         "vaccine_efficacy_disease", "vaccine_logistical_delay", "vaccine_protection_delay",
                         "spatial_ratio_scan", "spatial_vax_bp_sim", "spatial_calc", "seeds", "pop",
                         "SC2_generation_time", "SC2_prop_asymptomatic", "SC2_prob_hosp", "SC2_hospitalisation_delay",
                         "quarantine_time", "quarantine_efficacy_scan", "prob_quarantine_symptoms", "prob_quarantine_contact_traced", "quarantine_efficacy_scan"))
  clusterEvalQ(cl, {
    library(dplyr) 
    library(tidyr)
  })
  
  ## Running to generate all the scenarios that have the vaccine
  SC1_storage <- array(data = NA, dim = c(iterations, length(R0_scan), length(surveillance_scan), length(spatial_ratio_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan), 4))
  SC2_storage <- array(data = NA, dim = c(iterations, length(R0_scan), length(surveillance_scan), length(spatial_ratio_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan), 4))
  for (i in 1:length(R0_scan)) {
    for (j in 1:length(surveillance_scan)) {
      for (k in 1:length(spatial_ratio_scan)) {
        for (l in 1:length(vaccine_efficacy_infection_scan)) {
          for (m in 1:length(quarantine_efficacy_scan)) {
            
            # Setup parallel processing for the iterations
            # tic()
            clusterExport(cl, list("i", "j", "k", "l", "m"))
            results <- parLapply(cl, 1:iterations, function(n) {
              
              this_seed <- seeds[n]
              vaccination_radius <- spatial_ratio_scan[k] * mu
              
              # SARS-CoV-1 Pathogen Archetype
              SC1_temp <- spatial_vax_bp_sim(offspring = "pois",
                                             mn_offspring = R0_scan[i],
                                             generation_time = SC1_generation_time,
                                             spatial_kernel = spatial_kernel,
                                             t0 = 0, tf = Inf,
                                             initial_immune = 0,
                                             check_final_size = check_final_size,
                                             seeding_cases = seeding_cases,
                                             prop_asymptomatic = SC1_prop_asymptomatic,
                                             infection_to_onset = SC1_infection_to_onset,
                                             prob_hosp = SC1_prob_hosp,
                                             hospitalisation_delay = SC1_hospitalisation_delay,
                                             detection_threshold = surveillance_scan[j],
                                             vaccine_campaign_radius = vaccination_radius,
                                             vaccine_coverage = vaccine_coverage,
                                             vaccine_efficacy_infection = vaccine_efficacy_infection_scan[l],
                                             vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[l],
                                             vaccine_efficacy_disease = vaccine_efficacy_disease,
                                             vaccine_logistical_delay = vaccine_logistical_delay,
                                             vaccine_protection_delay = vaccine_protection_delay,
                                             seed = this_seed,
                                             population = pop,
                                             time_to_quarantine = quarantine_time,
                                             prob_quarantine_contact_traced = prob_quarantine_contact_traced,
                                             prob_quarantine_symptoms = prob_quarantine_symptoms,
                                             quarantine_efficacy = quarantine_efficacy_scan[m])
              SC1_count <- sum(!is.na(SC1_temp$time_infection))
              SC1_time_to_n <- time_to_nth_infection(tdf = SC1_temp, n = n)[[1]]
              SC1_Reff <- calculate_Reff(SC1_temp, "spatial_vax")
              SC1_R0 <- calculate_R0(SC1_temp)
              
              # SARS-CoV-2 Pathogen Archetype
              SC2_temp <- spatial_vax_bp_sim(offspring = "pois",
                                             mn_offspring = R0_scan[i],
                                             generation_time = SC2_generation_time,
                                             spatial_kernel = spatial_kernel,
                                             t0 = 0, tf = Inf,
                                             initial_immune = 0,
                                             check_final_size = check_final_size,
                                             seeding_cases = seeding_cases,
                                             prop_asymptomatic = SC2_prop_asymptomatic,
                                             infection_to_onset = SC2_infection_to_onset,
                                             prob_hosp = SC2_prob_hosp,
                                             hospitalisation_delay = SC2_hospitalisation_delay,
                                             detection_threshold = surveillance_scan[j],
                                             vaccine_campaign_radius = vaccination_radius,
                                             vaccine_coverage = vaccine_coverage,
                                             vaccine_efficacy_infection = vaccine_efficacy_infection_scan[l],
                                             vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[l],
                                             vaccine_efficacy_disease = vaccine_efficacy_disease,
                                             vaccine_logistical_delay = vaccine_logistical_delay,
                                             vaccine_protection_delay = vaccine_protection_delay,
                                             seed = this_seed,
                                             population = pop,
                                             time_to_quarantine = quarantine_time,
                                             prob_quarantine_contact_traced = prob_quarantine_contact_traced,
                                             prob_quarantine_symptoms = prob_quarantine_symptoms,
                                             quarantine_efficacy = quarantine_efficacy_scan[m])
              SC2_count <- sum(!is.na(SC2_temp$time_infection))
              SC2_time_to_n <- time_to_nth_infection(tdf = SC2_temp, n = n)[[1]]
              SC2_Reff <- calculate_Reff(SC2_temp, "spatial_vax")
              SC2_R0 <- calculate_R0(SC2_temp)
              
              list(SC1_count = SC1_count, SC2_count = SC2_count,
                   SC1_time_to_n = SC1_time_to_n, SC2_time_to_n = SC2_time_to_n,
                   SC1_Reff = SC1_Reff, SC2_Reff = SC2_Reff,
                   SC1_R0 = SC1_R0, SC2_R0 = SC2_R0)
            })
            
            # Extract results and store them in the respective storage arrays
            for (n in 1:iterations) {
              SC1_storage[n, i, j, k, l, m, 1] <- results[[n]]$SC1_count
              SC2_storage[n, i, j, k, l, m, 1] <- results[[n]]$SC2_count
              
              SC1_storage[n, i, j, k, l, m, 2] <- results[[n]]$SC1_time_to_n
              SC2_storage[n, i, j, k, l, m, 2] <- results[[n]]$SC2_time_to_n
              
              SC1_storage[n, i, j, k, l, m, 3] <- results[[n]]$SC1_Reff
              SC2_storage[n, i, j, k, l, m, 3] <- results[[n]]$SC2_Reff
              
              SC1_storage[n, i, j, k, l, m, 4] <- results[[n]]$SC1_R0
              SC2_storage[n, i, j, k, l, m, 4] <- results[[n]]$SC2_R0
            }
            print(paste0("i = ", i, ", j = ", j, ", k = ", k, ", l = " , l, ", m = ", m))
            # toc()
          }
        }
      }
    }
  }
  
  ## Running to generate the scenarios that don't have the vaccine
  SC1_no_vaccine_storage <- array(data = NA, dim = c(iterations, length(R0_scan), 1, 1, 1, length(quarantine_efficacy_scan), 4))
  SC2_no_vaccine_storage <- array(data = NA, dim = c(iterations, length(R0_scan), 1, 1, 1, length(quarantine_efficacy_scan), 4))
  for (i in 1:length(R0_scan)) {
    for (j in 1:length(quarantine_efficacy_scan)) {
      
      # Setup parallel processing for the iterations
      clusterExport(cl, list("i", "j"))
      results <- parLapply(cl, 1:iterations, function(k) {
        
        # Setting the seed
        this_seed <- seeds[k]
        
        # SARS-CoV-1 Pathogen Archetype
        SC1_temp <- spatial_vax_bp_sim(offspring = "pois",
                                       mn_offspring = R0_scan[i],
                                       generation_time = SC1_generation_time,
                                       spatial_kernel = spatial_kernel,
                                       t0 = 0, tf = Inf,
                                       check_final_size = check_final_size,
                                       seeding_cases = seeding_cases,
                                       prop_asymptomatic = SC1_prop_asymptomatic,
                                       infection_to_onset = SC1_infection_to_onset,
                                       prob_hosp = SC1_prob_hosp,
                                       hospitalisation_delay = SC1_hospitalisation_delay,
                                       detection_threshold = 1000,
                                       vaccine_campaign_radius = 1,
                                       vaccine_coverage = 0,
                                       vaccine_efficacy_infection = 0,
                                       vaccine_efficacy_transmission = 0,
                                       vaccine_efficacy_disease = 0,
                                       vaccine_logistical_delay = 100,
                                       vaccine_protection_delay = 100,
                                       seed = this_seed,
                                       initial_immune = 0,
                                       population = pop,
                                       time_to_quarantine = quarantine_time,
                                       prob_quarantine_contact_traced = prob_quarantine_contact_traced,
                                       prob_quarantine_symptoms = prob_quarantine_symptoms,
                                       quarantine_efficacy = quarantine_efficacy_scan[j])
        SC1_count <- sum(!is.na(SC1_temp$time_infection))
        SC1_time_to_n <- time_to_nth_infection(tdf = SC1_temp, n = n)[[1]]
        SC1_Reff <- calculate_Reff(SC1_temp, "spatial_vax")
        SC1_R0 <- calculate_R0(SC1_temp)
        
        # SARS-CoV-2 Pathogen Archetype
        SC2_temp <- spatial_vax_bp_sim(mn_offspring = R0_scan[i],
                                       generation_time = SC2_generation_time,
                                       spatial_kernel = spatial_kernel,
                                       t0 = 0, tf = Inf,
                                       check_final_size = check_final_size,
                                       seeding_cases = seeding_cases,
                                       prop_asymptomatic = SC2_prop_asymptomatic,
                                       infection_to_onset = SC2_infection_to_onset,
                                       prob_hosp = SC2_prob_hosp,
                                       hospitalisation_delay = SC2_hospitalisation_delay,
                                       detection_threshold = 1000,
                                       vaccine_campaign_radius = 1,
                                       vaccine_coverage = 0,
                                       vaccine_efficacy_infection = 0,
                                       vaccine_efficacy_transmission = 0,
                                       vaccine_efficacy_disease = 0,
                                       vaccine_logistical_delay = 100,
                                       vaccine_protection_delay = 100,
                                       seed = this_seed,
                                       initial_immune = 0,
                                       population = pop,
                                       time_to_quarantine = quarantine_time,
                                       prob_quarantine_contact_traced = prob_quarantine_contact_traced,
                                       prob_quarantine_symptoms = prob_quarantine_symptoms,
                                       quarantine_efficacy = quarantine_efficacy_scan[j])
        SC2_count <- sum(!is.na(SC2_temp$time_infection))
        SC2_time_to_n <- time_to_nth_infection(tdf = SC2_temp, n = n)[[1]]
        SC2_Reff <- calculate_Reff(SC2_temp, "spatial_vax")
        SC2_R0 <- calculate_R0(SC2_temp)
        
        list(SC1_count = SC1_count, SC2_count = SC2_count,
             SC1_time_to_n = SC1_time_to_n, SC2_time_to_n = SC2_time_to_n,
             SC1_Reff = SC1_Reff, SC2_Reff = SC2_Reff,
             SC1_R0 = SC1_R0, SC2_R0 = SC2_R0)
        
      })

      # Extract results and store them in the respective storage arrays
      for (k in 1:iterations) {
        SC1_no_vaccine_storage[k, i, 1, 1, 1, j, 1] <- results[[n]]$SC1_count
        SC2_no_vaccine_storage[k, i, 1, 1, 1, j, 1] <- results[[n]]$SC2_count

        SC1_no_vaccine_storage[k, i, 1, 1, 1, j, 2] <- results[[n]]$SC1_time_to_n
        SC2_no_vaccine_storage[k, i, 1, 1, 1, j, 2] <- results[[n]]$SC2_time_to_n
        
        SC1_no_vaccine_storage[k, i, 1, 1, 1, j, 3] <- results[[n]]$SC1_Reff
        SC2_no_vaccine_storage[k, i, 1, 1, 1, j, 3] <- results[[n]]$SC2_Reff
        
        SC1_no_vaccine_storage[k, i, 1, 1, 1, j, 4] <- results[[n]]$SC1_R0
        SC2_no_vaccine_storage[k, i, 1, 1, 1, j, 4] <- results[[n]]$SC2_R0
      }
      
      print(paste0("i = ", i, ", j = ", j))
    }
  }
  stopCluster(cl)
  
  outcome_names <- c("epidemic_size", "time_to_n", "Reff", "R0")
  
  # SARS-CoV-1 Results Processing
  
  ## With vaccine
  SC1_vaccine_reshaped <- reshape2::melt(SC1_storage)
  colnames(SC1_vaccine_reshaped) <- c("iteration", "R0", "surveillance", "spatial_ratio", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  SC1_vaccine_reshaped <- SC1_vaccine_reshaped %>%
    mutate(iteration = as.integer(iteration),
           input_R0 = R0_scan[R0],
           surveillance = surveillance_scan[surveillance],
           spatial_ratio = spatial_ratio_scan[spatial_ratio], 
           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy], 
           vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy], 
           quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
           outcome = outcome_names[outcome]) %>% 
    dplyr::select(iteration, input_R0, -R0, surveillance, spatial_ratio, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
                  quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
    pivot_wider(names_from = "outcome", values_from = value)
  SC1_vaccine_reshaped$vaccine <- "yes_vaccine"
  
  ## Without vaccine
  SC1_no_vaccine_reshaped <- reshape2::melt(SC1_no_vaccine_storage)
  colnames(SC1_no_vaccine_reshaped) <- c("iteration", "R0", "surveillance", "spatial_ratio", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  SC1_no_vaccine_reshaped <- SC1_no_vaccine_reshaped %>%
    mutate(iteration = as.integer(iteration),
           input_R0 = R0_scan[R0],
           surveillance = surveillance_scan[surveillance],
           spatial_ratio = spatial_ratio_scan[spatial_ratio], 
           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy], 
           vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy], 
           quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
           outcome = outcome_names[outcome]) %>% 
    dplyr::select(iteration, input_R0, -R0, surveillance, spatial_ratio, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
                  quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
    pivot_wider(names_from = "outcome", values_from = value)
  SC1_no_vaccine_reshaped$vaccine <- "no_vaccine"
  
  SC1_reshaped2 <- rbind(SC1_vaccine_reshaped, SC1_no_vaccine_reshaped)
  SC1_reshaped2$pathogen <- "SARS-CoV-1"
  saveRDS(SC1_reshaped2, "outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_SC1_paramScan.rds")
  
  # SARS-CoV-2 Results Processing
  
  ## With vaccine
  SC2_vaccine_reshaped <- reshape2::melt(SC2_storage)
  colnames(SC2_vaccine_reshaped) <- c("iteration", "R0", "surveillance", "spatial_ratio", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  SC2_vaccine_reshaped <- SC2_vaccine_reshaped %>%
    mutate(iteration = as.integer(iteration),
           input_R0 = R0_scan[R0],
           surveillance = surveillance_scan[surveillance],
           spatial_ratio = spatial_ratio_scan[spatial_ratio], 
           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy], 
           vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy], 
           quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
           outcome = outcome_names[outcome]) %>% 
    dplyr::select(iteration, input_R0, -R0, surveillance, spatial_ratio, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
                  quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
    pivot_wider(names_from = "outcome", values_from = value)
  SC2_vaccine_reshaped$vaccine <- "yes_vaccine"
  
  ## Without vaccine
  SC2_no_vaccine_reshaped <- reshape2::melt(SC2_no_vaccine_storage)
  colnames(SC2_no_vaccine_reshaped) <- c("iteration", "R0", "surveillance", "spatial_ratio", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  SC2_no_vaccine_reshaped <- SC2_no_vaccine_reshaped %>%
    mutate(iteration = as.integer(iteration),
           input_R0 = R0_scan[R0],
           surveillance = surveillance_scan[surveillance],
           spatial_ratio = spatial_ratio_scan[spatial_ratio], 
           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy], 
           vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy], 
           quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
           outcome = outcome_names[outcome]) %>% 
    dplyr::select(iteration, input_R0, -R0, surveillance, spatial_ratio, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
                  quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
    pivot_wider(names_from = "outcome", values_from = value)
  SC2_no_vaccine_reshaped$vaccine <- "no_vaccine"
  
  SC2_reshaped2 <- rbind(SC2_vaccine_reshaped, SC2_no_vaccine_reshaped)
  SC2_reshaped2$pathogen <- "SARS-CoV-2"
  saveRDS(SC2_reshaped2, "outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_SC2_paramScan.rds")
  
} else {
  
  SC1_reshaped2 <- readRDS("outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_SC1_paramScan.rds")
  SC2_reshaped2 <- readRDS("outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_SC2_paramScan.rds")
  
}
toc()

## Plotting the proportion of outbreaks controlled
SC1_no_vaccine_higher_efficacy_dummy_df <- SC1_reshaped2 %>%
  filter(vaccine == "no_vaccine") %>%
  mutate(vaccine_efficacy_infection = 0.75)
SC2_no_vaccine_higher_efficacy_dummy_df <- SC2_reshaped2 %>%
  filter(vaccine == "no_vaccine") %>%
  mutate(vaccine_efficacy_infection = 0.75)

overall_spatial_vax_df <- rbind(SC1_reshaped2, SC2_reshaped2, SC1_no_vaccine_higher_efficacy_dummy_df, SC2_no_vaccine_higher_efficacy_dummy_df) %>%
  filter(surveillance != 100) %>%
  filter(spatial_ratio == 50) %>%
  mutate(contained = ifelse(epidemic_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(input_R0, pathogen, quarantine_efficacy) %>%
  mutate(time_to_n_relative = time_to_n / time_to_n[vaccine == "no_vaccine"]) %>%
  group_by(input_R0, pathogen, vaccine, vaccine_efficacy_infection, quarantine_efficacy, surveillance, spatial_ratio) %>%
  summarise(proportion_contained = sum(contained) / iterations,
            R0_actual = mean(input_R0, na.rm = TRUE),
            Reff_mean = median(Reff, na.rm = TRUE),
            Reff_lower = quantile(Reff, 0.1, na.rm = TRUE),
            Reff_upper = quantile(Reff, 0.9, na.rm = TRUE),
            time_to_n = mean(time_to_n, na.rm = TRUE),
            time_to_n_relative = mean(time_to_n_relative, na.rm = TRUE)) %>%
  mutate(proportion_contained = ifelse(input_R0 == 1.00, 1, proportion_contained)) 

combo_tbl  <- expand_grid(pathogen =  c("SARS-CoV-2", "SARS-CoV-1"), quarantine_efficacy = c(0, 0.65))
combo_plot <- combo_tbl %>% 
  mutate(p = map2(pathogen, quarantine_efficacy, ~ make_stacked_plot_Reff_spatialvax(overall_spatial_vax_df, .x, .y, 0.35, c(1, 2))))

Fig1BCDE <- plot_grid(plotlist = list(combo_plot$p[[3]], combo_plot$p[[1]], combo_plot$p[[4]], combo_plot$p[[2]]),
                      nrow = length(c(0, 0.65)), ncol  = 2,
                      labels = c("B", "C", "D", "E"), label_size = 10)

## Plotting Supplementary Figure looking at time to epidemic threshold
overall_spatial_vax_df$vaccine_quarantine_elision <- paste0("Vaccine Efficacy = ", overall_spatial_vax_df$vaccine_efficacy_infection, "\nQuarantine Efficacy = ", overall_spatial_vax_df$quarantine_efficacy)
time_to_n_plot <- ggplot(subset(overall_spatial_vax_df, quarantine_efficacy != 0.35),
                         aes(x = input_R0, y = time_to_n_relative, col = interaction(vaccine, factor(surveillance)))) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_grid(vaccine_quarantine_elision~pathogen) + 
  scale_colour_manual(labels = c("No\nVaccine", paste0(surveillance_scan, " Hosp.")),
                      values = palette,
                      name = "Surveillance\nThreshold\nTrigger") +
  labs(x = "R0", y = "Fold Increase in Time to Epidemic Threshold") +
  theme(strip.background = element_rect(fill = "white"))
ggsave(plot = time_to_n_plot, filename = "figures/Figure_1_BranchingProcess/FigS2_ParamScan_timetoN.pdf", height = 8.5, width = 8)

#########################################################################
## Parameter scan sensitivity analyses (Figure 2C-E)
#########################################################################
fresh_run_vaccination_heatmaps <- TRUE
tic()
if (fresh_run_vaccination_heatmaps) {
  
  # Generating the seeds
  set.seed(2000)
  seeds <- runif(n = iterations, min = 1, max = 10^9)
  
  ## Fixed values for the runs
  R0_seq <- R0_scan[R0_scan > 1]
  surveillance_threshold_fixed <- 10
  spatial_ratio_fixed <- 50
  vaccine_efficacy_infection_scan <- c(0, 0.35, 0.75)
  vaccine_efficacy_transmission_scan <- c(0, 0.35, 0.5)
  outcome_names <- c("epidemic_size", "time_to_n", "Reff", "R0")
  time_to_n_indicator <- 2000
  
  #######################################################################
  ## Sensitivity Analysis - R0 vs Spatial Vax Radius
  #######################################################################
  
  ## Parameter scan arguments
  spatial_ratio_scan_full <- c(1, 10, 25, 50, 75, 100)
  storage_R0_SpatialRadius_sensitivity <- array(data = NA, dim = c(iterations, length(R0_seq), length(spatial_ratio_scan_full), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan), 4))
  
  ## Setting up the cluster to support the parallel runs
  no_cores <- min(iterations, 10)
  cl <- makeCluster(no_cores)
  set.seed(2000)
  seeds <- runif(n = iterations, min = 1, max = 10^9)
  clusterExport(cl, list("mu", "R0_seq", "SC1_generation_time", "spatial_kernel", "calculate_Reff", "calculate_R0", "surveillance_threshold_fixed",
                         "check_final_size", "seeding_cases", "SC1_prop_asymptomatic", "time_to_nth_infection", "spatial_ratio_scan_full",
                         "SC1_prob_hosp", "SC1_hospitalisation_delay", "surveillance_scan", "n",
                         "SC1_infection_to_onset", "SC2_infection_to_onset", "implement_quarantine",
                         "vaccine_coverage", "vaccine_efficacy_infection_scan", "vaccine_efficacy_transmission_scan",
                         "vaccine_efficacy_disease", "vaccine_logistical_delay", "vaccine_protection_delay",
                         "spatial_ratio_scan", "spatial_vax_bp_sim", "spatial_calc", "seeds", "pop",
                         "SC2_generation_time", "SC2_prop_asymptomatic", "SC2_prob_hosp", "SC2_hospitalisation_delay",
                         "quarantine_time", "quarantine_efficacy_scan", "prob_quarantine_symptoms", "prob_quarantine_contact_traced", "quarantine_efficacy_scan"))
  clusterEvalQ(cl, {
    library(dplyr) 
    library(tidyr)
  })
  
  ## Running the simulations
  for (i in 1:length(R0_seq)) {
    for (j in 1:length(spatial_ratio_scan_full)) {
      for (k in 1:length(vaccine_efficacy_infection_scan)) {
        for (l in 1:length(quarantine_efficacy_scan)) {
          
          # Setup parallel processing for the iterations
          clusterExport(cl, list("i", "j", "k", "l"))
          results <- parLapply(cl, 1:iterations, function(m) {
            
            # SARS-CoV-2 Pathogen Archetype
            SC2_temp <- spatial_vax_bp_sim(offspring = "pois",
                                           mn_offspring = R0_seq[i],
                                           generation_time = SC2_generation_time,
                                           spatial_kernel = spatial_kernel,
                                           t0 = 0, tf = Inf,
                                           initial_immune = 0,
                                           check_final_size = check_final_size,
                                           seeding_cases = seeding_cases,
                                           prop_asymptomatic = SC2_prop_asymptomatic,
                                           infection_to_onset = SC2_infection_to_onset,
                                           prob_hosp = SC2_prob_hosp,
                                           hospitalisation_delay = SC2_hospitalisation_delay,
                                           detection_threshold = surveillance_threshold_fixed,
                                           vaccine_campaign_radius = spatial_ratio_scan_full[j] * mu,
                                           vaccine_coverage = vaccine_coverage,
                                           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[k],
                                           vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[k],
                                           vaccine_efficacy_disease = vaccine_efficacy_disease,
                                           vaccine_logistical_delay = vaccine_logistical_delay,
                                           vaccine_protection_delay = vaccine_protection_delay,
                                           seed = seeds[m],
                                           population = pop,
                                           time_to_quarantine = quarantine_time,
                                           prob_quarantine_contact_traced = prob_quarantine_contact_traced,
                                           prob_quarantine_symptoms = prob_quarantine_symptoms,
                                           quarantine_efficacy = quarantine_efficacy_scan[l])
            count <- sum(!is.na(SC2_temp$time_infection))
            time_to_n <- time_to_nth_infection(tdf = SC2_temp, n = time_to_n_indicator)[[1]]
            Reff <- calculate_Reff(SC2_temp, "spatial_vax")
            R0 <- calculate_R0(SC2_temp)
            
            list(count = count, time_to_n = time_to_n, Reff = Reff, R0 = R0)})
          
          # Extract results and store them in the respective storage arrays
          for (n in 1:iterations) {
            storage_R0_SpatialRadius_sensitivity[n, i, j, k, l, 1] <- results[[n]]$count
            storage_R0_SpatialRadius_sensitivity[n, i, j, k, l, 2] <- results[[n]]$time_to_n
            storage_R0_SpatialRadius_sensitivity[n, i, j, k, l, 3] <- results[[n]]$Reff
            storage_R0_SpatialRadius_sensitivity[n, i, j, k, l, 4] <- results[[n]]$R0
          }
          print(paste0("i = ", i, ", j = ", j, ", k = ", k, ", l = " , l))
        }
      }
    }
  }
  stopCluster(cl)
  
  ## Processing the simulations
  reshaped_R0_SpatialRadius_sensitivity <- reshape2::melt(storage_R0_SpatialRadius_sensitivity)
  colnames(reshaped_R0_SpatialRadius_sensitivity) <- c("iteration", "R0", "spatial_ratio", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  reshaped_R0_SpatialRadius_sensitivity <- reshaped_R0_SpatialRadius_sensitivity %>%
    mutate(iteration = as.integer(iteration),
           input_R0 = R0_seq[R0],
           spatial_ratio = spatial_ratio_scan_full[spatial_ratio], 
           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy], 
           vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy], 
           quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
           outcome = outcome_names[outcome]) %>% 
    dplyr::select(iteration, input_R0, -R0, spatial_ratio, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
                  quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
    pivot_wider(names_from = "outcome", values_from = value)
  saveRDS(object = reshaped_R0_SpatialRadius_sensitivity, file = "outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_R0SpatialScan.rds")
  
  #######################################################################
  ## Sensitivity Analysis - R0 vs Vaccine Efficacy
  #######################################################################
  
  ## Parameter scan arguments
  vaccine_efficacy_scan_full <- seq(0.3, 0.9, 0.1)
  storage_R0_VaccineEff_sensitivity <- array(data = NA, dim = c(iterations, length(R0_seq), length(vaccine_efficacy_scan_full), length(quarantine_efficacy_scan), 4))
  
  ## Setting up the cluster to support the parallel runs
  no_cores <- min(iterations, 10)
  cl <- makeCluster(no_cores)
  set.seed(2000)
  seeds <- runif(n = iterations, min = 1, max = 10^9)
  clusterExport(cl, list("mu", "R0_seq", "SC1_generation_time", "spatial_kernel", "calculate_Reff", "calculate_R0", "surveillance_threshold_fixed", 
                         "check_final_size", "seeding_cases", "SC1_prop_asymptomatic", "time_to_nth_infection", "spatial_ratio_fixed",
                         "SC1_prob_hosp", "SC1_hospitalisation_delay", "surveillance_scan", "n",
                         "SC1_infection_to_onset", "SC2_infection_to_onset", "implement_quarantine",
                         "vaccine_coverage", "vaccine_efficacy_scan_full", "vaccine_efficacy_scan_full",
                         "vaccine_efficacy_disease", "vaccine_logistical_delay", "vaccine_protection_delay",
                         "spatial_ratio_scan", "spatial_vax_bp_sim", "spatial_calc", "seeds", "pop",
                         "SC2_generation_time", "SC2_prop_asymptomatic", "SC2_prob_hosp", "SC2_hospitalisation_delay",
                         "quarantine_time", "quarantine_efficacy_scan", "prob_quarantine_symptoms", "prob_quarantine_contact_traced", "quarantine_efficacy_scan"))
  clusterEvalQ(cl, {
    library(dplyr) 
    library(tidyr)
  })
  
  ## Running the simulations
  for (i in 1:length(R0_seq)) {
    for (k in 1:length(vaccine_efficacy_scan_full)) {
      for (l in 1:length(quarantine_efficacy_scan)) {
        
        # Setup parallel processing for the iterations
        clusterExport(cl, list("i", "k", "l"))
        results <- parLapply(cl, 1:iterations, function(m) {
          
          # SARS-CoV-2 Pathogen Archetype
          SC2_temp <- spatial_vax_bp_sim(offspring = "pois",
                                         mn_offspring = R0_seq[i],
                                         generation_time = SC2_generation_time,
                                         spatial_kernel = spatial_kernel,
                                         t0 = 0, tf = Inf,
                                         initial_immune = 0,
                                         check_final_size = check_final_size,
                                         seeding_cases = seeding_cases,
                                         prop_asymptomatic = SC2_prop_asymptomatic,
                                         infection_to_onset = SC2_infection_to_onset,
                                         prob_hosp = SC2_prob_hosp,
                                         hospitalisation_delay = SC2_hospitalisation_delay,
                                         detection_threshold = surveillance_threshold_fixed,
                                         vaccine_campaign_radius = spatial_ratio_fixed * mu,
                                         vaccine_coverage = vaccine_coverage,
                                         vaccine_efficacy_infection = vaccine_efficacy_scan_full[k],
                                         vaccine_efficacy_transmission = vaccine_efficacy_scan_full[k],
                                         vaccine_efficacy_disease = vaccine_efficacy_disease,
                                         vaccine_logistical_delay = vaccine_logistical_delay,
                                         vaccine_protection_delay = vaccine_protection_delay,
                                         seed = seeds[m],
                                         population = pop,
                                         time_to_quarantine = quarantine_time,
                                         prob_quarantine_contact_traced = prob_quarantine_contact_traced,
                                         prob_quarantine_symptoms = prob_quarantine_symptoms,
                                         quarantine_efficacy = quarantine_efficacy_scan[l])
          count <- sum(!is.na(SC2_temp$time_infection))
          time_to_n <- time_to_nth_infection(tdf = SC2_temp, n = time_to_n_indicator)[[1]]
          Reff <- calculate_Reff(SC2_temp, "spatial_vax")
          R0 <- calculate_R0(SC2_temp)
          
          list(count = count, time_to_n = time_to_n, Reff = Reff, R0 = R0)})
        
        # Extract results and store them in the respective storage arrays
        for (n in 1:iterations) {
          storage_R0_VaccineEff_sensitivity[n, i, k, l, 1] <- results[[n]]$count
          storage_R0_VaccineEff_sensitivity[n, i, k, l, 2] <- results[[n]]$time_to_n
          storage_R0_VaccineEff_sensitivity[n, i, k, l, 3] <- results[[n]]$Reff
          storage_R0_VaccineEff_sensitivity[n, i, k, l, 4] <- results[[n]]$R0
        }
        print(paste0("i = ", i, ", k = ", k, ", l = " , l))
      }
    }
  }
  stopCluster(cl)
  
  ## Processing the simulations
  reshaped_R0_VaccineEff_sensitivity <- reshape2::melt(storage_R0_VaccineEff_sensitivity)
  colnames(reshaped_R0_VaccineEff_sensitivity) <- c("iteration", "R0", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  reshaped_R0_VaccineEff_sensitivity <- reshaped_R0_VaccineEff_sensitivity %>%
    mutate(iteration = as.integer(iteration),
           input_R0 = R0_seq[R0],
           vaccine_efficacy_infection = vaccine_efficacy_scan_full[vaccine_efficacy], 
           vaccine_efficacy_transmission = vaccine_efficacy_scan_full[vaccine_efficacy], 
           quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
           outcome = outcome_names[outcome]) %>% 
    dplyr::select(iteration, input_R0, -R0, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
                  quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
    pivot_wider(names_from = "outcome", values_from = value)
  saveRDS(object = reshaped_R0_VaccineEff_sensitivity, file = "outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_R0Efficacy.rds")
  
  #######################################################################
  ## Sensitivity Analysis - R0 vs Surveillance Threshold
  #######################################################################
  
  ## Parameter scan arguments
  storage_R0_SurvThreshold_sensitivity <- array(data = NA, dim = c(iterations, length(R0_seq), length(surveillance_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan), 4))
  
  ## Setting up the cluster to support the parallel runs
  no_cores <- min(iterations, 10)
  cl <- makeCluster(no_cores)
  set.seed(2000)
  seeds <- runif(n = iterations, min = 1, max = 10^9)
  clusterExport(cl, list("mu", "R0_seq", "SC1_generation_time", "spatial_kernel", "calculate_Reff", "calculate_R0", "surveillance_threshold_fixed",
                         "check_final_size", "seeding_cases", "SC1_prop_asymptomatic", "time_to_nth_infection", "spatial_ratio_scan_full",
                         "SC1_prob_hosp", "SC1_hospitalisation_delay", "surveillance_scan", "n", "spatial_ratio_fixed",
                         "SC1_infection_to_onset", "SC2_infection_to_onset", "implement_quarantine",
                         "vaccine_coverage", "vaccine_efficacy_infection_scan", "vaccine_efficacy_transmission_scan",
                         "vaccine_efficacy_disease", "vaccine_logistical_delay", "vaccine_protection_delay",
                         "spatial_ratio_scan", "spatial_vax_bp_sim", "spatial_calc", "seeds", "pop",
                         "SC2_generation_time", "SC2_prop_asymptomatic", "SC2_prob_hosp", "SC2_hospitalisation_delay",
                         "quarantine_time", "quarantine_efficacy_scan", "prob_quarantine_symptoms", "prob_quarantine_contact_traced", "quarantine_efficacy_scan"))
  clusterEvalQ(cl, {
    library(dplyr) 
    library(tidyr)
  })
  
  ## Running the simulations
  for (i in 1:length(R0_seq)) {
    for (j in 1:length(surveillance_scan)) {
      for (k in 1:length(vaccine_efficacy_infection_scan)) {
        for (l in 1:length(quarantine_efficacy_scan)) {
          
          # Setup parallel processing for the iterations
          clusterExport(cl, list("i", "j", "k", "l"))
          results <- parLapply(cl, 1:iterations, function(m) {
            
            # SARS-CoV-2 Pathogen Archetype
            SC2_temp <- spatial_vax_bp_sim(offspring = "pois",
                                           mn_offspring = R0_seq[i],
                                           generation_time = SC2_generation_time,
                                           spatial_kernel = spatial_kernel,
                                           t0 = 0, tf = Inf,
                                           initial_immune = 0,
                                           check_final_size = check_final_size,
                                           seeding_cases = seeding_cases,
                                           prop_asymptomatic = SC2_prop_asymptomatic,
                                           infection_to_onset = SC2_infection_to_onset,
                                           prob_hosp = SC2_prob_hosp,
                                           hospitalisation_delay = SC2_hospitalisation_delay,
                                           detection_threshold = surveillance_scan[j],
                                           vaccine_campaign_radius = spatial_ratio_fixed * mu,
                                           vaccine_coverage = vaccine_coverage,
                                           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[k],
                                           vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[k],
                                           vaccine_efficacy_disease = vaccine_efficacy_disease,
                                           vaccine_logistical_delay = vaccine_logistical_delay,
                                           vaccine_protection_delay = vaccine_protection_delay,
                                           seed = seeds[m],
                                           population = pop,
                                           time_to_quarantine = quarantine_time,
                                           prob_quarantine_contact_traced = prob_quarantine_contact_traced,
                                           prob_quarantine_symptoms = prob_quarantine_symptoms,
                                           quarantine_efficacy = quarantine_efficacy_scan[l])
            count <- sum(!is.na(SC2_temp$time_infection))
            time_to_n <- time_to_nth_infection(tdf = SC2_temp, n = time_to_n_indicator)[[1]]
            Reff <- calculate_Reff(SC2_temp, "spatial_vax")
            R0 <- calculate_R0(SC2_temp)
            
            list(count = count, time_to_n = time_to_n, Reff = Reff, R0 = R0)})
          
          # Extract results and store them in the respective storage arrays
          for (n in 1:iterations) {
            storage_R0_SurvThreshold_sensitivity[n, i, j, k, l, 1] <- results[[n]]$count
            storage_R0_SurvThreshold_sensitivity[n, i, j, k, l, 2] <- results[[n]]$time_to_n
            storage_R0_SurvThreshold_sensitivity[n, i, j, k, l, 3] <- results[[n]]$Reff
            storage_R0_SurvThreshold_sensitivity[n, i, j, k, l, 4] <- results[[n]]$R0
          }
          print(paste0("i = ", i, ", j = ", j, ", k = ", k, ", l = " , l))
        }
      }
    }
  }
  stopCluster(cl)
  
  ## Processing the simulations
  reshaped_R0_SurvThreshold_sensitivity <- reshape2::melt(storage_R0_SurvThreshold_sensitivity)
  colnames(reshaped_R0_SurvThreshold_sensitivity) <- c("iteration", "R0", "surveillance_threshold", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  reshaped_R0_SurvThreshold_sensitivity <- reshaped_R0_SurvThreshold_sensitivity %>%
    mutate(iteration = as.integer(iteration),
           input_R0 = R0_seq[R0],
           surveillance_threshold = surveillance_scan[surveillance_threshold], 
           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy], 
           vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy], 
           quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
           outcome = outcome_names[outcome]) %>% 
    dplyr::select(iteration, input_R0, -R0, surveillance_threshold, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
                  quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
    pivot_wider(names_from = "outcome", values_from = value)
  saveRDS(object = reshaped_R0_SurvThreshold_sensitivity, file = "outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_R0SurvThreshold.rds")
  
} else {
  
  reshaped_R0_SpatialRadius_sensitivity <- readRDS(file = "outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_R0SpatialScan.rds")
  reshaped_R0_VaccineEff_sensitivity <- readRDS(file = "outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_R0Efficacy.rds")
  reshaped_R0_SurvThreshold_sensitivity <- readRDS(file = "outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_R0SurvThreshold.rds")
  
}
toc()

## need to re run these with the new vaccine efficacy scans (that include 0) in them
## Creating Vaccination-Related Heatmaps

### R0_SpatialRadius
# R0_SpatialRadius_df <- reshaped_R0_SpatialRadius_sensitivity %>% 
#   ungroup() %>% 
#   mutate(contained  = (epidemic_size < 0.9 * check_final_size),
#          time_to_n_2 = ifelse(is.na(time_to_n),
#                               max(time_to_n, na.rm = TRUE),
#                               time_to_n)) %>% 
#   filter(input_R0 == 1.25, spatial_ratio == 100, quarantine_efficacy == 0) %>%
#   group_by(iteration, input_R0, spatial_ratio, quarantine_efficacy) %>%  # ← add iteration
#   mutate(
#     time_to_n_relative   = time_to_n  / time_to_n[vaccine_efficacy_infection < 0.01],
#     time_to_n_2_relative = time_to_n_2 / time_to_n_2[vaccine_efficacy_infection < 0.01])


R0_SpatialRadius_df <- reshaped_R0_SpatialRadius_sensitivity %>%
  ungroup() %>%
  mutate(contained = ifelse(epidemic_size < (0.9 * check_final_size), 1, 0)) %>%
  mutate(time_to_n_2 = ifelse(is.na(time_to_n), max(time_to_n, na.rm = TRUE), time_to_n)) %>%
  group_by(iteration, input_R0, spatial_ratio, quarantine_efficacy) %>%
  mutate(time_to_n_check = time_to_n[vaccine_efficacy_infection == 0],
         time_to_n_relative = time_to_n / time_to_n[vaccine_efficacy_infection == 0],
         time_to_n_2_relative = time_to_n_2 / time_to_n_2[vaccine_efficacy_infection == 0]) %>%
  ungroup() %>%
  group_by(input_R0, spatial_ratio, vaccine_efficacy_infection, quarantine_efficacy) %>%
  summarise(proportion_contained = sum(contained) / iterations,
            avg_time_to_n = mean(time_to_n, na.rm = TRUE),
            avg_time_to_n_2 = mean(time_to_n_2, na.rm = TRUE),
            avg_time_to_n_relative = mean(time_to_n_relative, na.rm = TRUE),
            avg_time_to_n_2_relative = mean(time_to_n_2_relative, na.rm = TRUE),
            avg_R0 = mean(R0),
            avg_Reff = mean(Reff)) %>%
  mutate(avg_time_to_n_2_relative_plot = ifelse(proportion_contained < 0.9, avg_time_to_n_2_relative, NA))






         
### Old plotting figures

# Results Plotting
overall <- rbind(SC1_reshaped2, SC2_reshaped2) %>%
  mutate(contained = ifelse(outbreak_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(R0, surveillance, spatial_ratio, vaccine_efficacy_infection, quarantine_efficacy, vaccine, pathogen) %>%
  summarise(proportion_contained = sum(contained) / iterations,
            avg_time_to_n = mean(time_to_n, na.rm = TRUE),
            avg_R0 = mean(R0),
            avg_Reff = mean(Reff))

### Plotting the output
vaccine_efficacy_index <- which(vaccine_efficacy_infection_scan == 0.35)
spatial_ratio_index <- which(spatial_ratio_scan == 50)
surveillance_scan_index <- which(surveillance_scan == 10)
surveillance_scan_exclude <- which(surveillance_scan == 25)

overall2 <- overall %>%
  filter(vaccine_efficacy == vaccine_efficacy_index | vaccine_efficacy == 0 ,
         spatial_ratio == spatial_ratio_index | spatial_ratio == 0,
         surveillance != surveillance_scan_exclude) %>%
  mutate(R0 = ifelse(vaccine == "no_vaccine", R0, R0_scan[R0]))

fig1HI <- ggplot(overall2, aes(x = R0, y = 100 * proportion_contained, col = factor(surveillance))) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_grid(. ~ pathogen,
             labeller = as_labeller(c(`SARS-CoV-1` = "SARS-CoV-1",
                                      `SARS-CoV-2` = "SARS-CoV-2"))) +
  scale_colour_manual(labels = c("No\nVaccine", paste0(surveillance_scan, " Hosp.")),
                      values = c("#474747", c("#E3AFCB", "#D474A4", "#B52F7B", "#9C105A", "#6B0045")),
                      name = "Surveillance\nThreshold\nTrigger") +
  labs(y = "% Outbreaks Contained") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))
ggsave(filename = "figures/Figure_1_BranchingProcess/Fig1HI_SpatialVaccination_ContainmentPlot.pdf", 
       plot = fig1c, 
       width = 8, height = 3.1)

R0_surveillance <- overall %>%
  mutate(R0 = ifelse(vaccine == "no_vaccine", R0, R0_scan[R0])) %>%
  filter(pathogen == "SARS-CoV-2" & 
           spatial_ratio == spatial_ratio_index & 
           vaccine_efficacy == vaccine_efficacy_index & 
           R0 > 1 & 
           vaccine == "yes_vaccine")
R0_surveillance_plot <- ggplot(R0_surveillance, aes(x = R0, y = surveillance , fill = 100 * proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 100), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  scale_y_continuous(breaks = 1:length(unique(R0_surveillance$surveillance)), labels = surveillance_scan) +
  scale_x_continuous(breaks = c(1.5, 2.0, 2.5), labels = c(1.5, 2, 2.5)) +
  labs(x = "R0",
       y = "Surveillance Threshold") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE) +
  theme(legend.position = "none")
ggsave(filename = "figures/Figure_1_BranchingProcess/Fig1L_SurvThresh_SpatialVax.pdf", 
       plot = R0_surveillance_plot, width = 2.5, height = 2.34)


R0_spatial <- overall %>%
  mutate(R0 = ifelse(vaccine == "no_vaccine", R0, R0_scan[R0])) %>%
  filter(pathogen == "SARS-CoV-2" & 
           surveillance == surveillance_scan_index & 
           vaccine_efficacy == vaccine_efficacy_index & 
           R0 > 1 & 
           vaccine == "yes_vaccine")
R0_spatial_plot <- ggplot(R0_spatial, aes(x = R0, y = spatial_ratio, fill = 100 * proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 100), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  scale_y_continuous(breaks = 1:length(unique(R0_spatial$spatial_ratio)), labels = spatial_ratio_scan) +
  scale_x_continuous(breaks = c(1.5, 2.0, 2.5), labels = c(1.5, 2, 2.5)) +
  labs(x = "R0",
       y = "Ratio Spatial Vax Radius") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE) +
  theme(legend.position = "none")
ggsave(filename = "figures/Figure_1_BranchingProcess/Fig1J_SpatRatio_SpatialVax.pdf", 
       plot = R0_spatial_plot, width = 2.5, height = 2.34)



vaccine_efficacy_central_index <- which(vaccine_efficacy_infection_scan == 0.35)
R0_efficacy <- overall %>%
  mutate(R0 = ifelse(vaccine == "no_vaccine", R0, R0_scan[R0])) %>%
  filter(pathogen == "SARS-CoV-2" & 
           surveillance == surveillance_scan_index & 
           spatial_ratio == spatial_ratio_index & 
           R0 > 1 & 
           vaccine == "yes_vaccine" &
           vaccine_efficacy != vaccine_efficacy_central_index) %>%
  mutate(vaccine_efficacy = ifelse(vaccine_efficacy > 2, vaccine_efficacy - 1, vaccine_efficacy))
R0_efficacy_plot <- ggplot(R0_efficacy, aes(x = R0, y = vaccine_efficacy, fill = 100 * proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 100), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  scale_y_continuous(breaks = c(2, 4, 6), 
                     labels = c(40, 60, 80)) +
  scale_x_continuous(breaks = c(1.5, 2.0, 2.5), labels = c(1.5, 2, 2.5)) +
  labs(x = "R0",
       y = "Vaccine Efficacy (%)") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE) +
  theme(legend.position = "none")
ggsave(filename = "figures/Figure_1_BranchingProcess/Fig1K_Efficacy_SpatialVax.pdf", 
       plot = R0_efficacy_plot, width = 2.5, height = 2.34)


heatmap_legend <- cowplot::get_legend(R0_spatial_plot)
heatmaps <- cowplot::plot_grid(R0_spatial_plot + theme(legend.position = "none"), 
                               R0_efficacy_plot + theme(legend.position = "none"),
                               R0_surveillance_plot + theme(legend.position = "none"), 
                               heatmap_legend,
                               ncol = 4, rel_widths = c(1, 1, 1, 0.3))

legend <- R0_spatial_plot + theme(legend.position = "bottom")
ggsave(file = "figures/Figure_1_BranchingProcess/LegendSpatialVaxHeatmap.pdf", plot = legend, width = 2.4, height = 2.4)
