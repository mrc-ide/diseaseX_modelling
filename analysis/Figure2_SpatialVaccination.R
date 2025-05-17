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
SC2_generation_time <- function(n) { rgamma(n, shape = 13.5, rate = 2) } # 6.75 day generation time Gamma distributed (as per Walker et al, Science, 2020)
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
time_to_n_indicator <- check_final_size * 0.9
initial_immune <- 0
seeding_cases <- 3
iterations <- 200

#########################################################################
## R0 sensitivity analysis (Figure 2B)
#########################################################################
R0_scan <- c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5)
surveillance_scan <- c(1, 10, 25, 50, 75, 10000)
vaccine_efficacy_infection_scan <- c(0.35, 0.75) 
vaccine_efficacy_transmission_scan <- c(0.35, 0.5) 
spatial_ratio_scan <- c(50) # , 100)
quarantine_efficacy_scan <- c(0, 0.65) # 0.35, 0.65) # from the same article as above
length(R0_scan) * length(surveillance_scan) * length(vaccine_efficacy_infection_scan) * length(spatial_ratio_scan) * length(quarantine_efficacy_scan)

fresh_run_R0_sensitivity_analysis <- TRUE
if (fresh_run_R0_sensitivity_analysis) {
  
  ## Setting up the cluster to support the parallel runs
  no_cores <- min(iterations, 10)
  cl <- makeCluster(no_cores)
  set.seed(2000)
  seeds <- runif(n = iterations, min = 1, max = 10^9)
  clusterExport(cl, list("mu", "R0_scan", "SC1_generation_time", "spatial_kernel", "calculate_Reff", "calculate_R0",
                         "check_final_size", "seeding_cases", "SC1_prop_asymptomatic", "time_to_nth_infection",
                         "SC1_prob_hosp", "SC1_hospitalisation_delay", "surveillance_scan", "time_to_n_indicator",
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
              SC1_time_to_n <- time_to_nth_infection(tdf = SC1_temp, n = time_to_n_indicator)[[1]]
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
              SC2_time_to_n <- time_to_nth_infection(tdf = SC2_temp, n = time_to_n_indicator)[[1]]
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
  # SC1_no_vaccine_storage <- array(data = NA, dim = c(iterations, length(R0_scan), 1, 1, 1, length(quarantine_efficacy_scan), 4))
  # SC2_no_vaccine_storage <- array(data = NA, dim = c(iterations, length(R0_scan), 1, 1, 1, length(quarantine_efficacy_scan), 4))
  # for (i in 1:length(R0_scan)) {
  #   for (j in 1:length(quarantine_efficacy_scan)) {
  #     
  #     # Setup parallel processing for the iterations
  #     clusterExport(cl, list("i", "j"))
  #     results <- parLapply(cl, 1:iterations, function(k) {
  #       
  #       # Setting the seed
  #       this_seed <- seeds[k]
  #       
  #       # SARS-CoV-1 Pathogen Archetype
  #       SC1_temp <- spatial_vax_bp_sim(offspring = "pois",
  #                                      mn_offspring = R0_scan[i],
  #                                      generation_time = SC1_generation_time,
  #                                      spatial_kernel = spatial_kernel,
  #                                      t0 = 0, tf = Inf,
  #                                      check_final_size = check_final_size,
  #                                      seeding_cases = seeding_cases,
  #                                      prop_asymptomatic = SC1_prop_asymptomatic,
  #                                      infection_to_onset = SC1_infection_to_onset,
  #                                      prob_hosp = SC1_prob_hosp,
  #                                      hospitalisation_delay = SC1_hospitalisation_delay,
  #                                      detection_threshold = 1000,
  #                                      vaccine_campaign_radius = 1,
  #                                      vaccine_coverage = 0,
  #                                      vaccine_efficacy_infection = 0,
  #                                      vaccine_efficacy_transmission = 0,
  #                                      vaccine_efficacy_disease = 0,
  #                                      vaccine_logistical_delay = 100,
  #                                      vaccine_protection_delay = 100,
  #                                      seed = this_seed,
  #                                      initial_immune = 0,
  #                                      population = pop,
  #                                      time_to_quarantine = quarantine_time,
  #                                      prob_quarantine_contact_traced = prob_quarantine_contact_traced,
  #                                      prob_quarantine_symptoms = prob_quarantine_symptoms,
  #                                      quarantine_efficacy = quarantine_efficacy_scan[j])
  #       SC1_count <- sum(!is.na(SC1_temp$time_infection))
  #       SC1_time_to_n <- time_to_nth_infection(tdf = SC1_temp, n = time_to_n_indicator)[[1]]
  #       SC1_Reff <- calculate_Reff(SC1_temp, "spatial_vax")
  #       SC1_R0 <- calculate_R0(SC1_temp)
  #       
  #       # SARS-CoV-2 Pathogen Archetype
  #       SC2_temp <- spatial_vax_bp_sim(mn_offspring = R0_scan[i],
  #                                      generation_time = SC2_generation_time,
  #                                      spatial_kernel = spatial_kernel,
  #                                      t0 = 0, tf = Inf,
  #                                      check_final_size = check_final_size,
  #                                      seeding_cases = seeding_cases,
  #                                      prop_asymptomatic = SC2_prop_asymptomatic,
  #                                      infection_to_onset = SC2_infection_to_onset,
  #                                      prob_hosp = SC2_prob_hosp,
  #                                      hospitalisation_delay = SC2_hospitalisation_delay,
  #                                      detection_threshold = 1000,
  #                                      vaccine_campaign_radius = 1,
  #                                      vaccine_coverage = 0,
  #                                      vaccine_efficacy_infection = 0,
  #                                      vaccine_efficacy_transmission = 0,
  #                                      vaccine_efficacy_disease = 0,
  #                                      vaccine_logistical_delay = 100,
  #                                      vaccine_protection_delay = 100,
  #                                      seed = this_seed,
  #                                      initial_immune = 0,
  #                                      population = pop,
  #                                      time_to_quarantine = quarantine_time,
  #                                      prob_quarantine_contact_traced = prob_quarantine_contact_traced,
  #                                      prob_quarantine_symptoms = prob_quarantine_symptoms,
  #                                      quarantine_efficacy = quarantine_efficacy_scan[j])
  #       SC2_count <- sum(!is.na(SC2_temp$time_infection))
  #       SC2_time_to_n <- time_to_nth_infection(tdf = SC2_temp, n = time_to_n_indicator)[[1]]
  #       SC2_Reff <- calculate_Reff(SC2_temp, "spatial_vax")
  #       SC2_R0 <- calculate_R0(SC2_temp)
  #       
  #       list(SC1_count = SC1_count, SC2_count = SC2_count,
  #            SC1_time_to_n = SC1_time_to_n, SC2_time_to_n = SC2_time_to_n,
  #            SC1_Reff = SC1_Reff, SC2_Reff = SC2_Reff,
  #            SC1_R0 = SC1_R0, SC2_R0 = SC2_R0)
  #       
  #     })
  # 
  #     # Extract results and store them in the respective storage arrays
  #     for (k in 1:iterations) {
  #       SC1_no_vaccine_storage[k, i, 1, 1, 1, j, 1] <- results[[n]]$SC1_count
  #       SC2_no_vaccine_storage[k, i, 1, 1, 1, j, 1] <- results[[n]]$SC2_count
  # 
  #       SC1_no_vaccine_storage[k, i, 1, 1, 1, j, 2] <- results[[n]]$SC1_time_to_n
  #       SC2_no_vaccine_storage[k, i, 1, 1, 1, j, 2] <- results[[n]]$SC2_time_to_n
  #       
  #       SC1_no_vaccine_storage[k, i, 1, 1, 1, j, 3] <- results[[n]]$SC1_Reff
  #       SC2_no_vaccine_storage[k, i, 1, 1, 1, j, 3] <- results[[n]]$SC2_Reff
  #       
  #       SC1_no_vaccine_storage[k, i, 1, 1, 1, j, 4] <- results[[n]]$SC1_R0
  #       SC2_no_vaccine_storage[k, i, 1, 1, 1, j, 4] <- results[[n]]$SC2_R0
  #     }
  #     
  #     print(paste0("i = ", i, ", j = ", j))
  #   }
  # }
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
  # SC1_vaccine_reshaped$vaccine <- "yes_vaccine"
  
  ## Without vaccine
  # SC1_no_vaccine_reshaped <- reshape2::melt(SC1_no_vaccine_storage)
  # colnames(SC1_no_vaccine_reshaped) <- c("iteration", "R0", "surveillance", "spatial_ratio", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  # SC1_no_vaccine_reshaped <- SC1_no_vaccine_reshaped %>%
  #   mutate(iteration = as.integer(iteration),
  #          input_R0 = R0_scan[R0],
  #          surveillance = surveillance_scan[surveillance],
  #          spatial_ratio = spatial_ratio_scan[spatial_ratio], 
  #          vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy], 
  #          vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy], 
  #          quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
  #          outcome = outcome_names[outcome]) %>% 
  #   dplyr::select(iteration, input_R0, -R0, surveillance, spatial_ratio, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
  #                 quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
  #   pivot_wider(names_from = "outcome", values_from = value)
  # SC1_no_vaccine_reshaped$vaccine <- "no_vaccine"
  
  SC1_reshaped2 <- rbind(SC1_vaccine_reshaped) # , SC1_no_vaccine_reshaped)
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
  # SC2_vaccine_reshaped$vaccine <- "yes_vaccine"
  
  ## Without vaccine
  # SC2_no_vaccine_reshaped <- reshape2::melt(SC2_no_vaccine_storage)
  # colnames(SC2_no_vaccine_reshaped) <- c("iteration", "R0", "surveillance", "spatial_ratio", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  # SC2_no_vaccine_reshaped <- SC2_no_vaccine_reshaped %>%
  #   mutate(iteration = as.integer(iteration),
  #          input_R0 = R0_scan[R0],
  #          surveillance = surveillance_scan[surveillance],
  #          spatial_ratio = spatial_ratio_scan[spatial_ratio], 
  #          vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy], 
  #          vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy], 
  #          quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
  #          outcome = outcome_names[outcome]) %>% 
  #   dplyr::select(iteration, input_R0, -R0, surveillance, spatial_ratio, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
  #                 quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
  #   pivot_wider(names_from = "outcome", values_from = value)
  # SC2_no_vaccine_reshaped$vaccine <- "no_vaccine"
  
  SC2_reshaped2 <- rbind(SC2_vaccine_reshaped) # , SC2_no_vaccine_reshaped)
  SC2_reshaped2$pathogen <- "SARS-CoV-2"
  saveRDS(SC2_reshaped2, "outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_SC2_paramScan.rds")
  
} else {
  
  SC1_reshaped2 <- readRDS("outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_SC1_paramScan.rds")
  SC2_reshaped2 <- readRDS("outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_SC2_paramScan.rds")
  
}

## Plotting the proportion of outbreaks controlled
overall_spatial_vax_df <- rbind(SC1_reshaped2, SC2_reshaped2) %>%
  mutate(surveillance = ifelse(surveillance == 1, "1 Hosp.", surveillance)) %>%
  mutate(surveillance = ifelse(surveillance == 10, "10 Hosp.", surveillance)) %>%
  mutate(surveillance = ifelse(surveillance == 25, "25 Hosp", surveillance)) %>%
  mutate(surveillance = ifelse(surveillance == 50, "50 Hosp", surveillance)) %>%
  mutate(surveillance = ifelse(surveillance == 75, "75 Hosp", surveillance)) %>%
  mutate(surveillance = ifelse(surveillance == 10000, "zno_vaccination", surveillance))

overall_spatial_df <- overall_spatial_vax_df %>%
  filter(spatial_ratio == 50) %>%
  mutate(contained = ifelse(epidemic_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(input_R0, pathogen, quarantine_efficacy) %>%
  mutate(time_to_n_relative = time_to_n / time_to_n[surveillance == "zno_vaccination"]) %>%
  group_by(input_R0, pathogen, vaccine_efficacy_infection, quarantine_efficacy, surveillance, spatial_ratio) %>%
  summarise(proportion_contained = sum(contained) / n(),
            R0_actual = mean(input_R0, na.rm = TRUE),
            Reff_mean = median(Reff, na.rm = TRUE),
            Reff_lower = quantile(Reff, 0.1, na.rm = TRUE),
            Reff_upper = quantile(Reff, 0.9, na.rm = TRUE), 
            time_to_n = mean(time_to_n, na.rm = TRUE),
            time_to_n_relative = mean(time_to_n_relative, na.rm = TRUE)) %>%
  filter(vaccine_efficacy_infection == 0.35)
combo_tbl  <- expand_grid(pathogen =  c("SARS-CoV-2", "SARS-CoV-1"), quarantine_efficacy = c(0, 0.65))
combo_plot <- combo_tbl %>% 
  mutate(p = map2(pathogen, quarantine_efficacy, ~ make_stacked_plot_Reff_spatialvax(overall_spatial_df, .x, .y, 0.35, c(1, 2))))

Fig1BCDE <- plot_grid(plotlist = list(combo_plot$p[[3]], combo_plot$p[[1]], combo_plot$p[[4]], combo_plot$p[[2]]),
                      nrow = length(c(0, 0.65)), ncol  = 2,
                      labels = c("B", "C", "D", "E"), label_size = 10)

## Plotting Supplementary Figure looking at time to epidemic threshold
overall_spatial_vax_df2 <- rbind(SC1_reshaped2, SC2_reshaped2) %>%
  mutate(surveillance = ifelse(surveillance == 1, "1 Hosp.", surveillance)) %>%
  mutate(surveillance = ifelse(surveillance == 10, "10 Hosp.", surveillance)) %>%
  mutate(surveillance = ifelse(surveillance == 25, "25 Hosp", surveillance)) %>%
  mutate(surveillance = ifelse(surveillance == 50, "50 Hosp", surveillance)) %>%
  mutate(surveillance = ifelse(surveillance == 75, "75 Hosp", surveillance)) %>%
  mutate(surveillance = ifelse(surveillance == 10000, "zno_vaccination", surveillance))

overall_spatial_df2 <- overall_spatial_vax_df2 %>%
  filter(spatial_ratio == 50) %>%
  mutate(contained = ifelse(epidemic_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(input_R0, pathogen, quarantine_efficacy) %>%
  mutate(time_to_n_relative = time_to_n / time_to_n[surveillance == "zno_vaccination"]) %>%
  group_by(input_R0, pathogen, vaccine_efficacy_infection, quarantine_efficacy, surveillance, spatial_ratio) %>%
  summarise(proportion_contained = sum(contained) / n(),
            R0_actual = mean(input_R0, na.rm = TRUE),
            Reff_mean = median(Reff, na.rm = TRUE),
            Reff_lower = quantile(Reff, 0.1, na.rm = TRUE),
            Reff_upper = quantile(Reff, 0.9, na.rm = TRUE), 
            time_to_n = mean(time_to_n, na.rm = TRUE),
            time_to_n_relative = mean(time_to_n_relative, na.rm = TRUE)) 

## Plotting Supplementary Figure looking at time to epidemic threshold
palette <- c("#E3AFCB", "#D474A4", "#B52F7B", "#9C105A", "#6B0045", "#474747")
overall_spatial_df2$vaccine_quarantine_elision <- paste0("Vaccine Efficacy = ", overall_spatial_df2$vaccine_efficacy_infection, "\nQuarantine Efficacy = ", overall_spatial_df2$quarantine_efficacy)
overall_spatial_df3 <- overall_spatial_df2 %>%
  filter(input_R0 > 1) %>%
  group_by(vaccine_quarantine_elision, pathogen, surveillance) %>%
  mutate(time_to_n_relative_mean2 = ifelse(is.na(time_to_n_relative), max(time_to_n_relative, na.rm = TRUE), time_to_n_relative))

SC2_rectangle_df <- overall_spatial_df3 %>% 
  filter(surveillance != 10000,
         quarantine_efficacy != 0.35, 
         pathogen == "SARS-CoV-2") %>%
  arrange(input_R0) %>%                                   
  group_by(vaccine_quarantine_elision, surveillance) %>%
  filter(proportion_contained != 1) %>%
  slice(1) %>%
  mutate(xmin = 0.75, xmax = input_R0-0.25, ymin = -Inf, ymax =  Inf)

SC2_time_to_n_plot <- ggplot(subset(overall_spatial_df2, quarantine_efficacy != 0.35 &
                                      pathogen == "SARS-CoV-2" & 
                                      surveillance != 10000),
                             aes(x = input_R0, y = time_to_n_relative, col = factor(surveillance))) +
  geom_line() +
  geom_rect_pattern(data = SC2_rectangle_df, 
                    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = ymax), 
                    inherit.aes = FALSE, fill = "grey80", pattern = "stripe",
                    pattern_angle = 45, pattern_size = 0.4, pattern_density = 0.01,
                    pattern_spacing = 0.1, pattern_colour = "black", 
                    alpha = 1) +
  geom_point(shape = 21, fill  = "white", size  = 8, stroke = 1.1) +
  geom_text(aes(label = scales::percent(proportion_contained, accuracy = 3)),
            size = 2.8, vjust = -2, hjust = 0.5) +
  geom_text(aes(label = paste0(round(time_to_n_relative, 1), "x")),
            size = 2.8, vjust = 0.5, hjust = 0.5, col = "black") +
  theme_bw() +
  facet_grid(vaccine_quarantine_elision ~ surveillance,
             labeller = labeller(scenario = c(`1` = "1 Hosp.",
                                              `10` = "10 Hosp.",
                                              `25` = "25 Hosp.",
                                              `50` = "50 Hosp.",
                                              `75` = "75 Hosp."))) + 
  scale_colour_manual(labels = c(paste0(surveillance_scan, " Hosp."), "No\nVaccine"),
                      values = palette,
                      name = "Surveillance\nThreshold\nTrigger") +
  labs(x = "R0", y = "Fold Increase in Time to Epidemic Threshold") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(size = 7)) +
  lims(y = c(0, 7), x = c(0.75, 2.65))

SC1_rectangle_df <- overall_spatial_df3 %>% 
  filter(surveillance != 10000,
         quarantine_efficacy != 0.35, 
         pathogen == "SARS-CoV-1") %>%
  arrange(R0) %>%                                   
  group_by(vaccine_quarantine_elision, surveillance) %>%
  filter(proportion_contained != 1) %>%
  slice(1) %>%
  mutate(xmin = 0.75, xmax = R0-0.25, ymin = -Inf, ymax =  Inf)

SC1_time_to_n_plot <- ggplot(subset(overall_spatial_df2, quarantine_efficacy != 0.35 &
                                      pathogen == "SARS-CoV-1" & 
                                      surveillance != 10000),
                             aes(x = input_R0, y = time_to_n_relative, col = factor(surveillance))) +
  geom_line() +
  geom_rect_pattern(data = SC1_rectangle_df, 
                    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = ymax), 
                    inherit.aes = FALSE, fill = "grey80", pattern = "stripe",
                    pattern_angle = 45, pattern_size = 0.4, pattern_density = 0.01,
                    pattern_spacing = 0.1, pattern_colour = "black", 
                    alpha = 1) +
  geom_point(shape = 21, fill  = "white", size  = 8, stroke = 1.1) +
  geom_text(aes(label = scales::percent(proportion_contained, accuracy = 3)),
            size = 2.8, vjust = -2, hjust = 0.5) +
  geom_text(aes(label = paste0(round(time_to_n_relative, 1), "x")),
            size = 2.8, vjust = 0.5, hjust = 0.5, col = "black") +
  theme_bw() +
  facet_grid(vaccine_quarantine_elision ~ surveillance,
             labeller = labeller(scenario = c(`1` = "1 Hosp.",
                                              `10` = "10 Hosp.",
                                              `25` = "25 Hosp.",
                                              `50` = "50 Hosp.",
                                              `75` = "75 Hosp."))) + 
  scale_colour_manual(labels = c(paste0(surveillance_scan, " Hosp."), "No\nVaccine"),
                      values = palette,
                      name = "Surveillance\nThreshold\nTrigger") +
  labs(x = "R0", y = "Fold Increase in Time to Epidemic Threshold") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(size = 7)) +
  lims(y = c(0, 7), x = c(0.75, 2.65))

time_to_n_plot <- cowplot::plot_grid(SC1_time_to_n_plot, SC2_time_to_n_plot, labels = c("A", "B"), nrow = 2)
ggsave(plot = time_to_n_plot, filename = "figures/Figure_2_SpatialVaccination/FigS2_ParamScan_timetoN.pdf", height = 8.5, width = 8)

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

  #######################################################################
  ## Sensitivity Analysis - R0 vs Spatial Vax Radius
  #######################################################################
  
  ## Parameter scan arguments
  spatial_ratio_scan_full <- c(1, 10, 25, 50, 75, 100)
  storage_R0_SpatialRadius_sensitivity <- array(data = NA, dim = c(iterations, length(R0_seq), length(spatial_ratio_scan_full), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan), 5))
  
  ## Setting up the cluster to support the parallel runs
  no_cores <- min(iterations, 10)
  cl <- makeCluster(no_cores)
  set.seed(2000)
  seeds <- runif(n = iterations, min = 1, max = 10^9)
  clusterExport(cl, list("mu", "R0_seq", "SC1_generation_time", "spatial_kernel", "calculate_Reff", "calculate_R0", "surveillance_threshold_fixed",
                         "check_final_size", "seeding_cases", "SC1_prop_asymptomatic", "time_to_nth_infection", "spatial_ratio_scan_full",
                         "SC1_prob_hosp", "SC1_hospitalisation_delay", "surveillance_scan", "time_to_n_indicator",
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
            
            list(count = count, time_to_n = time_to_n, Reff = Reff, R0 = R0, radius = spatial_ratio_scan_full[j])})
          
          # Extract results and store them in the respective storage arrays
          for (n in 1:iterations) {
            storage_R0_SpatialRadius_sensitivity[n, i, j, k, l, 1] <- results[[n]]$count
            storage_R0_SpatialRadius_sensitivity[n, i, j, k, l, 2] <- results[[n]]$time_to_n
            storage_R0_SpatialRadius_sensitivity[n, i, j, k, l, 3] <- results[[n]]$Reff
            storage_R0_SpatialRadius_sensitivity[n, i, j, k, l, 4] <- results[[n]]$R0
            storage_R0_SpatialRadius_sensitivity[n, i, j, k, l, 5] <- results[[n]]$radius
          }
          print(paste0("i = ", i, ", j = ", j, ", k = ", k, ", l = " , l))
        }
      }
    }
  }
  stopCluster(cl)
  
  ## Processing the simulations
  outcome_names_test <- c("epidemic_size", "time_to_n", "Reff", "R0", "radius")
  reshaped_R0_SpatialRadius_sensitivity <- reshape2::melt(storage_R0_SpatialRadius_sensitivity)
  colnames(reshaped_R0_SpatialRadius_sensitivity) <- c("iteration", "R0", "spatial_ratio", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  reshaped_R0_SpatialRadius_sensitivity <- reshaped_R0_SpatialRadius_sensitivity %>%
    mutate(iteration = as.integer(iteration),
           input_R0 = R0_seq[R0],
           spatial_ratio = spatial_ratio_scan_full[spatial_ratio], 
           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy], 
           vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[vaccine_efficacy], 
           quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
           outcome = outcome_names_test[outcome]) %>% 
    dplyr::select(iteration, input_R0, -R0, spatial_ratio, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
                  quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
    pivot_wider(names_from = "outcome", values_from = value)
  saveRDS(object = reshaped_R0_SpatialRadius_sensitivity, file = "outputs/Figure1_branchingProcess_Containment/Fig2_spatialVaccination_R0SpatialScan.rds")

  #######################################################################
  ## Sensitivity Analysis - R0 vs Vaccine Efficacy
  #######################################################################
  
  ## Parameter scan arguments
  vaccine_efficacy_scan_full <- c(0, seq(0.3, 0.9, 0.1))
  storage_R0_VaccineEff_sensitivity <- array(data = NA, dim = c(iterations, length(R0_seq), length(vaccine_efficacy_scan_full), length(quarantine_efficacy_scan), 4))
  
  ## Setting up the cluster to support the parallel runs
  no_cores <- min(iterations, 10)
  cl <- makeCluster(no_cores)
  set.seed(2000)
  seeds <- runif(n = iterations, min = 1, max = 10^9)
  clusterExport(cl, list("mu", "R0_seq", "SC1_generation_time", "spatial_kernel", "calculate_Reff", "calculate_R0", "surveillance_threshold_fixed", 
                         "check_final_size", "seeding_cases", "SC1_prop_asymptomatic", "time_to_nth_infection", "spatial_ratio_fixed",
                         "SC1_prob_hosp", "SC1_hospitalisation_delay", "surveillance_scan", "time_to_n_indicator",
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
  surveillance_scan_full <- c(1, 10, 25, 50, 75, 100)
  storage_R0_SurvThreshold_sensitivity <- array(data = NA, dim = c(iterations, length(R0_seq), length(surveillance_scan_full), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan), 4))
    
  ## Setting up the cluster to support the parallel runs
  no_cores <- min(iterations, 10)
  cl <- makeCluster(no_cores)
  set.seed(2000)
  seeds <- runif(n = iterations, min = 1, max = 10^9)
  clusterExport(cl, list("mu", "R0_seq", "SC1_generation_time", "spatial_kernel", "calculate_Reff", "calculate_R0", "surveillance_threshold_fixed",
                         "check_final_size", "seeding_cases", "SC1_prop_asymptomatic", "time_to_nth_infection", "spatial_ratio_scan_full",
                         "SC1_prob_hosp", "SC1_hospitalisation_delay", "surveillance_scan_full", "spatial_ratio_fixed",
                         "SC1_infection_to_onset", "SC2_infection_to_onset", "implement_quarantine", "time_to_n_indicator",
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
    for (j in 1:length(surveillance_scan_full)) {
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
                                           detection_threshold = surveillance_scan_full[j],
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
           surveillance_threshold = surveillance_scan_full[surveillance_threshold], 
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

# Creating Vaccination-Related Heatmaps

### R0_SpatialRadius
R0_SpatialRadius_df <- reshaped_R0_SpatialRadius_sensitivity %>%
  ungroup() %>%
  mutate(contained = ifelse(epidemic_size < (0.9 * check_final_size), 1, 0)) %>%
  mutate(time_to_n_2 = ifelse(is.na(time_to_n), max(time_to_n, na.rm = TRUE), time_to_n)) %>%
  group_by(iteration, input_R0, spatial_ratio) %>%
  mutate(time_to_n_check = time_to_n[vaccine_efficacy_infection == 0 & quarantine_efficacy == 0],
         time_to_n_relative = time_to_n / time_to_n[vaccine_efficacy_infection == 0 & quarantine_efficacy == 0],
         time_to_n_2_relative = time_to_n_2 / time_to_n_2[vaccine_efficacy_infection == 0 & quarantine_efficacy == 0]) %>%
  ungroup() %>%
  group_by(input_R0, spatial_ratio, vaccine_efficacy_infection, quarantine_efficacy) %>%
  summarise(proportion_contained = sum(contained) / n(),
            avg_time_to_n = mean(time_to_n, na.rm = TRUE),
            avg_time_to_n_2 = mean(time_to_n_2, na.rm = TRUE),
            avg_time_to_n_relative = mean(time_to_n_relative, na.rm = TRUE),
            avg_time_to_n_2_relative = mean(time_to_n_2_relative, na.rm = TRUE),
            avg_time_to_n_check = mean(time_to_n_check, na.rm = TRUE),
            avg_R0 = mean(R0),
            avg_Reff = mean(Reff)) %>%
  mutate(avg_time_to_n_2_relative_plot = ifelse(proportion_contained < 0.9, avg_time_to_n_2_relative, NA))

### Vaccine Efficacy
R0_efficacy_df <- reshaped_R0_VaccineEff_sensitivity %>%
  ungroup() %>%
  mutate(contained = ifelse(epidemic_size < (0.9 * check_final_size), 1, 0)) %>%
  mutate(time_to_n_2 = ifelse(is.na(time_to_n), max(time_to_n, na.rm = TRUE), time_to_n)) %>%
  group_by(iteration, input_R0) %>%
  mutate(time_to_n_check = time_to_n[vaccine_efficacy_infection == 0 & quarantine_efficacy == 0],
         time_to_n_relative = time_to_n / time_to_n[vaccine_efficacy_infection == 0 & quarantine_efficacy == 0],
         time_to_n_2_relative = time_to_n_2 / time_to_n_2[vaccine_efficacy_infection == 0 & quarantine_efficacy == 0]) %>%
  ungroup() %>%
  group_by(input_R0, vaccine_efficacy_infection, quarantine_efficacy) %>%
  summarise(proportion_contained = sum(contained) / n(),
            avg_time_to_n = mean(time_to_n, na.rm = TRUE),
            avg_time_to_n_2 = mean(time_to_n_2, na.rm = TRUE),
            avg_time_to_n_relative = mean(time_to_n_relative, na.rm = TRUE),
            avg_time_to_n_2_relative = mean(time_to_n_2_relative, na.rm = TRUE),
            avg_R0 = mean(R0),
            avg_Reff = mean(Reff)) %>%
  mutate(avg_time_to_n_2_relative_plot = ifelse(proportion_contained < 0.9, avg_time_to_n_2_relative, NA))

### Surveillance Threshold
R0_SurvThresh_df <- reshaped_R0_SurvThreshold_sensitivity %>%
  ungroup() %>%
  mutate(contained = ifelse(epidemic_size < (0.9 * check_final_size), 1, 0)) %>%
  mutate(time_to_n_2 = ifelse(is.na(time_to_n), max(time_to_n, na.rm = TRUE), time_to_n)) %>%
  group_by(iteration, input_R0, surveillance_threshold) %>%
  mutate(time_to_n_check = time_to_n[vaccine_efficacy_infection == 0 & quarantine_efficacy == 0],
         time_to_n_relative = time_to_n / time_to_n[vaccine_efficacy_infection == 0 & quarantine_efficacy == 0],
         time_to_n_2_relative = time_to_n_2 / time_to_n_2[vaccine_efficacy_infection == 0 & quarantine_efficacy == 0]) %>%
  ungroup() %>%
  group_by(input_R0, surveillance_threshold, vaccine_efficacy_infection, quarantine_efficacy) %>%
  summarise(proportion_contained = sum(contained) / n(),
            avg_time_to_n = mean(time_to_n, na.rm = TRUE),
            avg_time_to_n_2 = mean(time_to_n_2, na.rm = TRUE),
            avg_time_to_n_relative = mean(time_to_n_relative, na.rm = TRUE),
            avg_time_to_n_2_relative = mean(time_to_n_2_relative, na.rm = TRUE),
            avg_R0 = mean(R0),
            avg_Reff = mean(Reff)) %>%
  mutate(avg_time_to_n_2_relative_plot = ifelse(proportion_contained < 0.9, avg_time_to_n_2_relative, NA))

############################################################
### Parameter Scan Figure Plots
############################################################

############################################################
### R0 / Spatial Vax Radius Parameter Scans
############################################################
main_contained_R0_SpatialRadius_plot <- ggplot(subset(R0_SpatialRadius_df, quarantine_efficacy != 0.35 & vaccine_efficacy_infection == 0.35), aes(x = input_R0, y = factor(spatial_ratio), fill = 100 * proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 100), begin = 0.175, end = 1, name = "Proportion\nContained",
                       oob = scales::squish,
                       direction = -1) +
  labs(x = "R0", y = "Ratio Spatial Vax Radius") +
  facet_grid(quarantine_efficacy ~ .,
             labeller = labeller(vaccine_efficacy_infection = c(`0` = "No Vaccine", 
                                                                `0.35` = "Vaccine Efficacy = 35%", 
                                                                `0.75` = "Vaccine Efficacy = 75%"),
                                 quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

SI_contained_R0_SpatialRadius_plot <- ggplot(subset(R0_SpatialRadius_df, quarantine_efficacy != 0.35), aes(x = input_R0, y = factor(spatial_ratio), fill = 100 * proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 100), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  labs(x = "R0", y = "Ratio Spatial Vax Radius") +
  facet_grid(quarantine_efficacy ~ vaccine_efficacy_infection,
             labeller = labeller(vaccine_efficacy_infection = c(`0` = "No Vaccine", 
                                                                `0.35` = "Vaccine Efficacy = 35%", 
                                                                `0.75` = "Vaccine Efficacy = 75%"),
                                 quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text.y = element_text(size = 8),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

SI_Reff_R0_SpatialRadius_plot <- ggplot(subset(R0_SpatialRadius_df, quarantine_efficacy != 0.35), aes(x = input_R0, y = factor(spatial_ratio), fill = 100 * (1 - (avg_Reff / avg_R0)))) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 70), begin = 0.175, end = 1,
                       name = "% Red.\nin R",
                       oob = scales::squish,
                       direction = 1) +
  labs(x = "R0", y = "Ratio Spatial Vax Radius") +
  facet_grid(quarantine_efficacy ~ vaccine_efficacy_infection,
             labeller = labeller(vaccine_efficacy_infection = c(`0` = "No Vaccine", 
                                                                `0.35` = "Vaccine Efficacy = 35%", 
                                                                `0.75` = "Vaccine Efficacy = 75%"),
                                 quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text.y = element_text(size = 8),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

SI_timetoN_R0_SpatialRadius_plot_part1 <- ggplot(R0_SpatialRadius_df, aes(x = input_R0, y = factor(spatial_ratio), fill = avg_time_to_n_relative,
                 alpha = 1 - proportion_contained)) +
  geom_tile(colour = "grey") +
  scale_fill_viridis_c(option = "rocket", limits = c(0.8, max(R0_SpatialRadius_df$avg_time_to_n_relative, na.rm = TRUE)), begin = 0.175, end = 1,
                       oob = scales::squish,
                       name = "Fold Increase\nin Time to Epidemic\nThreshold",
                       direction = 1) +
  geom_tile(data = filter(R0_SpatialRadius_df, proportion_contained < 0.5), aes(x = input_R0, y = factor(spatial_ratio)),
            fill = NA, colour = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = filter(R0_SpatialRadius_df, proportion_contained < 0.5), aes(x = input_R0, y = factor(spatial_ratio),
                                                                                 label = paste0(sprintf("%.1f", avg_time_to_n_relative), "x")),
            colour = "white", size = 3, inherit.aes = FALSE) +
  scale_alpha(name = "% Outbreaks\nNot Contained",
              limits = c(0, 1)) +
  labs(x = "R0", y = "Ratio Spatial Vax Radius") +
  facet_grid(quarantine_efficacy ~ vaccine_efficacy_infection,
             labeller = labeller(vaccine_efficacy_infection = c(`0` = "No Vaccine",
                                                                `0.35` = "Vaccine Efficacy = 35%", 
                                                                `0.75` = "Vaccine Efficacy = 75%"),
                                 quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(size = 8),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE) 

R0_SpatialRadius_legend_df <- expand.grid(colour_val = seq(0.9, max(R0_SpatialRadius_df$avg_time_to_n_relative, na.rm = TRUE), length.out = 10),
                                          alpha_val  = seq(0, 1,   length.out = 10))
SI_timetoN_R0_SpatialRadius_plot_part2 <- ggplot(R0_SpatialRadius_legend_df, aes(x = 100 * alpha_val, y = colour_val, fill = colour_val, alpha = alpha_val)) +
  geom_raster() +
  scale_fill_viridis_c(option = "rocket", limits = c(.9, max(R0_SpatialRadius_df$avg_time_to_n_relative, na.rm = TRUE)),
                       begin = 0.175, end = 1, oob = scales::squish, name = "Fold Increase\n(Time to threshold)") +
  scale_alpha(range = c(0, 1), name  = "% Outbreaks\nNot Contained") +
  geom_tile(data = filter(R0_SpatialRadius_legend_df, alpha_val > 0.5), aes(x = 100 * alpha_val, y = colour_val),
            fill = NA, colour = "black", inherit.aes = FALSE) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("100%", "75%", "50%", "25%", "0%")) +
  guides(fill  = "none", alpha = "none") +
  labs(y = "Fold Increas in Time to\nEpidemic Threshold", x = "% Outbreaks Controlled") +
  coord_cartesian(expand = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.background = element_blank(), panel.grid = element_blank())

R0_SpatialRadius_legend <- cowplot::plot_grid(NULL, SI_timetoN_R0_SpatialRadius_plot_part2, NULL, nrow = 3, rel_heights = c(1, 2, 1))

SI_timetoN_R0_SpatialRadius_plot <- cowplot::plot_grid(SI_timetoN_R0_SpatialRadius_plot_part1, R0_SpatialRadius_legend, nrow = 1, rel_widths = c(3, 1))

SI_R0_SpatialRadius_top_two_thirds <- cowplot::plot_grid(SI_contained_R0_SpatialRadius_plot, SI_Reff_R0_SpatialRadius_plot, nrow = 2,
                                                   align = "v", axis = "r", labels = c("A", "B"))
SI_R0_SpatialRadius_overall <- cowplot::plot_grid(SI_R0_SpatialRadius_top_two_thirds, SI_timetoN_R0_SpatialRadius_plot, nrow = 2, rel_heights = c(2, 1))
ggsave(file = "figures/Figure_2_SpatialVaccination/FigS2_R0SpatialRadius_overall.pdf", plot = SI_R0_SpatialRadius_overall, width = 8, height = 11)

############################################################
### R0 / Vaccine Efficacy Parameter Scans
############################################################
main_contained_R0_efficacy_plot <- ggplot(subset(R0_efficacy_df, quarantine_efficacy != 0.35 & vaccine_efficacy_infection != 0.0), aes(x = input_R0, y = 100 * vaccine_efficacy_infection, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 1), begin = 0.175, end = 1, name = "Proportion\nContained",
                       oob = scales::squish,
                       direction = -1) +
  labs(x = "R0", y = "Vaccine Efficacy (%)") +
  facet_grid(quarantine_efficacy ~ .,
             labeller = labeller(quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text.y = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

SI_Reff_R0_efficacy_plot <- ggplot(subset(R0_efficacy_df, quarantine_efficacy != 0.35 & vaccine_efficacy_infection != 0.0), aes(x = input_R0, y = 100 * vaccine_efficacy_infection, fill = 100 * (1 - (avg_Reff / avg_R0)))) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 70), begin = 0.175, end = 1,
                       name = "% Red.\nin R",
                       oob = scales::squish,
                       direction = 1) +
  labs(x = "R0", y = "Vaccine Efficacy (%)") +
  facet_grid(. ~ quarantine_efficacy,
             labeller = labeller(quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

SI_timetoN_R0_efficacy_plot_part1 <- ggplot(subset(R0_efficacy_df, vaccine_efficacy_infection != 0.0),
                                            aes(x = input_R0, y = 100 * vaccine_efficacy_infection, fill = avg_time_to_n_relative, alpha = 100 * (1 - proportion_contained))) +
  geom_tile(colour = "grey") +
  scale_fill_viridis_c(option = "rocket", limits = c(0.8, max(R0_efficacy_df$avg_time_to_n_relative, na.rm = TRUE)), begin = 0.175, end = 1,
                       # oob = scales::squish,
                       name = "Fold Increase\nin Time to Epidemic\nThreshold",
                       direction = 1) +
  geom_tile(data = filter(subset(R0_efficacy_df, vaccine_efficacy_infection != 0.0), proportion_contained < 0.5), aes(x = input_R0, y = 100 * vaccine_efficacy_infection),
            fill = NA, colour = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = filter(subset(R0_efficacy_df, vaccine_efficacy_infection != 0.0), proportion_contained < 0.5), aes(x = input_R0, y = 100 * vaccine_efficacy_infection,
                                                                                                                      label = paste0(sprintf("%.1f", avg_time_to_n_relative), "x")),
            colour = "white", size = 3, inherit.aes = FALSE) +
  scale_alpha(name = "% Outbreaks\nNot Contained",
              limits = c(0, 100)) +
  labs(x = "R0", y = "Vaccine Efficacy (%)") +
  facet_grid(. ~ quarantine_efficacy,
             labeller = labeller(vaccine_efficacy_infection = c(`0` = "No Vaccine",
                                                                `0.35` = "Vaccine Efficacy = 35%", 
                                                                `0.75` = "Vaccine Efficacy = 75%"),
                                 quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(size = 8),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE) 

R0_Efficacy_legend_df <- expand.grid(colour_val = seq(0.9, max(subset(R0_efficacy_df, vaccine_efficacy_infection != 0.0)$avg_time_to_n_relative, na.rm = TRUE), length.out = 10),
                                     alpha_val  = seq(0, 1,   length.out = 10))
SI_timetoN_R0_Efficacy_plot_part2 <- ggplot(R0_Efficacy_legend_df, aes(x = 100 * alpha_val, y = colour_val, fill = colour_val, alpha = alpha_val)) +
  geom_raster() +
  scale_fill_viridis_c(option = "rocket", limits = c(.9, max(subset(R0_efficacy_df, vaccine_efficacy_infection != 0.0)$avg_time_to_n_relative, na.rm = TRUE)),
                       begin = 0.175, end = 1, oob = scales::squish, name = "Fold Increase\n(Time to threshold)") +
  scale_alpha(range = c(0, 1), name  = "% Outbreaks\nNot Contained") +
  geom_tile(data = filter(R0_Efficacy_legend_df, alpha_val > 0.5), aes(x = 100 * alpha_val, y = colour_val),
            fill = NA, colour = "black", inherit.aes = FALSE) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("100%", "75%", "50%", "25%", "0%")) +
  guides(fill  = "none", alpha = "none") +
  labs(y = "Fold Increas in Time to\nEpidemic Threshold", x = "% Outbreaks Controlled") +
  coord_cartesian(expand = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.background = element_blank(), panel.grid = element_blank())

R0_efficacy_legend <- cowplot::plot_grid(NULL, SI_timetoN_R0_Efficacy_plot_part2, NULL, nrow = 3, rel_heights = c(1, 2, 1))

SI_timetoN_R0_efficacy_plot <- cowplot::plot_grid(SI_timetoN_R0_efficacy_plot_part1, R0_efficacy_legend, nrow = 1, rel_widths = c(3, 1))

SI_R0_efficacy_overall <- cowplot::plot_grid(SI_Reff_R0_efficacy_plot, SI_timetoN_R0_efficacy_plot, nrow = 2, rel_heights = c(1, 1), labels = c("A", "B"))
ggsave(file = "figures/Figure_2_SpatialVaccination/FigS2_R0Efficacy_overall.pdf", plot = SI_R0_efficacy_overall, width = 8, height = 6)

############################################################
### R0 / Surveillance Threshold Parameter Scans
############################################################
main_contained_R0_SurvThresh_plot <- ggplot(subset(R0_SurvThresh_df, quarantine_efficacy != 0.35 & vaccine_efficacy_infection == 0.35), aes(x = input_R0, y = factor(surveillance_threshold), fill = 100 * proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 100), begin = 0.175, end = 1, name = "Proportion\nContained",
                       oob = scales::squish,
                       direction = -1) +
  labs(x = "R0", y = "Surveillance\nThreshold") +
  facet_grid(quarantine_efficacy ~ .,
             labeller = labeller(vaccine_efficacy_infection = c(`0` = "No Vaccine", 
                                                                `0.35` = "Vaccine Efficacy = 35%", 
                                                                `0.75` = "Vaccine Efficacy = 75%"),
                                 quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

SI_contained_R0_SurvThresh_plot <- ggplot(subset(R0_SurvThresh_df, quarantine_efficacy != 0.35), aes(x = input_R0, y = factor(surveillance_threshold), fill = 100 * proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 100), begin = 0.175, end = 1, name = "Proportion\nContained",
                       oob = scales::squish,
                       direction = -1) +
  labs(x = "R0", y = "Surveillance Threshold") +
  facet_grid(quarantine_efficacy ~ vaccine_efficacy_infection,
             labeller = labeller(vaccine_efficacy_infection = c(`0` = "No Vaccine", 
                                                                `0.35` = "Vaccine Efficacy = 35%", 
                                                                `0.75` = "Vaccine Efficacy = 75%"),
                                 quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

SI_Reff_R0_SurvThresh_plot <- ggplot(subset(R0_SurvThresh_df, quarantine_efficacy != 0.35), aes(x = input_R0, y = factor(surveillance_threshold), fill = 100 * (1 - (avg_Reff / avg_R0)))) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 70), begin = 0.175, end = 1,
                       oob = scales::squish,
                       name = "% Red.\nin R",
                       direction = 1) +
  labs(x = "R0", y = "Surveillance Threshold") +
  facet_grid(quarantine_efficacy ~ vaccine_efficacy_infection,
             labeller = labeller(vaccine_efficacy_infection = c(`0` = "No Vaccine", 
                                                                `0.35` = "Vaccine Efficacy = 35%", 
                                                                `0.75` = "Vaccine Efficacy = 75%"),
                                 quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

SI_timetoN_R0_SurvThresh_plot_part1 <- ggplot(subset(R0_SurvThresh_df, quarantine_efficacy != 0.35), 
                                           aes(x = input_R0, y = factor(surveillance_threshold), fill = avg_time_to_n_relative, alpha = 100 * (1 - proportion_contained))) +
  geom_tile(colour = "grey") +
  scale_fill_viridis_c(option = "rocket", limits = c(0.8, max(R0_SurvThresh_df$avg_time_to_n_relative, na.rm = TRUE)), begin = 0.175, end = 1,
                       oob = scales::squish,
                       name = "Fold Increase\nin Time to Epidemic\nThreshold",
                       direction = 1) +
  geom_tile(data = filter(R0_SurvThresh_df, proportion_contained < 0.5), aes(x = input_R0, y = factor(surveillance_threshold)),
            fill = NA, colour = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = filter(R0_SurvThresh_df, proportion_contained < 0.5), aes(x = input_R0, y = factor(surveillance_threshold),
                                                                             label = paste0(sprintf("%.1f", avg_time_to_n_relative), "x")),
            colour = "white", size = 3, inherit.aes = FALSE) +
  scale_alpha(name = "% Outbreaks\nNot Contained",
              limits = c(0, 100)) +
  labs(x = "R0", y = "Surveillance Threshold") +
  facet_grid(quarantine_efficacy ~ vaccine_efficacy_infection,
             labeller = labeller(vaccine_efficacy_infection = c(`0` = "No Vaccine",
                                                                `0.35` = "Vaccine Efficacy = 35%", 
                                                                `0.75` = "Vaccine Efficacy = 75%"),
                                 quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        legend.position = "none",
        strip.text.y = element_text(size = 8),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

R0_SurvThresh_legend_df <- expand.grid(colour_val = seq(0.9, max(R0_SurvThresh_df$avg_time_to_n_relative, na.rm = TRUE), length.out = 10),
                                       alpha_val  = seq(0, 1,   length.out = 10))
SI_timetoN_R0_SurvThresh_plot_part2 <- ggplot(R0_SurvThresh_legend_df, aes(x = 100 * alpha_val, y = colour_val, fill = colour_val, alpha = alpha_val)) +
  geom_raster() +
  scale_fill_viridis_c(option = "rocket", limits = c(.9, max(R0_SurvThresh_df$avg_time_to_n_relative, na.rm = TRUE)),
                       begin = 0.175, end = 1, oob = scales::squish, name = "Fold Increase\n(Time to threshold)") +
  scale_alpha(range = c(0, 1), name  = "% Outbreaks\nNot Contained") +
  geom_tile(data = filter(R0_SurvThresh_legend_df, alpha_val > 0.5), aes(x = 100 * alpha_val, y = colour_val),
            fill = NA, colour = "black", inherit.aes = FALSE) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("100%", "75%", "50%", "25%", "0%")) +
  guides(fill  = "none", alpha = "none") +
  labs(y = "Fold Increas in Time to\nEpidemic Threshold", x = "% Outbreaks Controlled") +
  coord_cartesian(expand = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        panel.background = element_blank(), panel.grid = element_blank())

R0_SurvThresh_legend <- cowplot::plot_grid(NULL, SI_timetoN_R0_SurvThresh_plot_part2, NULL, nrow = 3, rel_heights = c(1, 2, 1))

SI_timetoN_R0_SurvThresh_plot <- cowplot::plot_grid(SI_timetoN_R0_SurvThresh_plot_part1, R0_SurvThresh_legend, nrow = 1, rel_widths = c(3, 1))

SI_R0_SurvThresh_top_two_thirds <- cowplot::plot_grid(SI_contained_R0_SurvThresh_plot, SI_Reff_R0_SurvThresh_plot, nrow = 2,
                                                         align = "v", axis = "r", labels = c("A", "B"))
SI_R0_SurvThresh_overall <- cowplot::plot_grid(SI_R0_SurvThresh_top_two_thirds, SI_timetoN_R0_SurvThresh_plot, nrow = 2, rel_heights = c(2, 1))
ggsave(file = "figures/Figure_2_SpatialVaccination/FigS2_R0SurvThresh_overall.pdf", plot = SI_R0_SurvThresh_overall, width = 8, height = 11)

## Overall main figure
Fig1FGH <- cowplot::plot_grid(main_contained_R0_SpatialRadius_plot + theme(legend.position = "none"),
                              main_contained_R0_efficacy_plot + theme(legend.position = "none"), 
                              main_contained_R0_SurvThresh_plot + theme(legend.position = "none"),
                              nrow = 1,
                              labels = c("F", "G", "H"), rel_widths = c(1, 1, 1.08))
overall_figure2 <- cowplot::plot_grid(Fig1BCDE, Fig1FGH, nrow = 2, rel_heights = c(1.25, 1))
ggsave(file = "figures/Figure_2_SpatialVaccination/Fig2_Overall.pdf", plot = overall_figure2, width = 8, height = 9.5)



### Old Code
# SI_timetoN_R0_SpatialRadius_plot <- ggplot(subset(R0_SpatialRadius_df, quarantine_efficacy != 0.35 & vaccine_efficacy_infection != 0),
#                                            aes(x = input_R0, y = factor(spatial_ratio), fill = avg_time_to_n_relative)) +
#   geom_tile(colour = "black") +
#   scale_fill_viridis_c(option = "rocket", limits = c(0.9, 1.5), begin = 0.175, end = 1,
#                        breaks = c(1, 1.1, 1.2, 1.3, 1.4, 1.5),
#                        oob = scales::squish,
#                        name = "Fold Increase\nin Time to Epidemic\nThreshold",
#                        direction = 1) +
#   scale_alpha(name = "% Outbreaks\nNot Contained") +
#   labs(x = "R0", y = "Ratio Spatial Vax Radius") +
#   facet_grid(vaccine_efficacy_infection ~ quarantine_efficacy,
#              labeller = labeller(vaccine_efficacy_infection = c(`0.35` = "Vaccine Efficacy = 35%", 
#                                                                 `0.75` = "Vaccine Efficacy = 75%"),
#                                  quarantine_efficacy = c(`0`   = "No Quarantine", 
#                                                          `0.35` = "Quarantine Efficacy = 35%", 
#                                                          `0.65`   = "Quarantine Efficacy = 65%"))) +  
#   theme(axis.text = element_text(angle = 0),
#         plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 12),
#         panel.background = element_blank(),
#         strip.background = element_rect(fill = "white", colour = "black"),
#         panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
#   coord_cartesian(expand = FALSE) 