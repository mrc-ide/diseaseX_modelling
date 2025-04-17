# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/ring_vax_simulation.R"))
source(here::here("functions/implement_quarantine.R"))
source(here::here("functions/time_to_nth_infection.R"))
source(here::here("functions/helper_functions.R"))

### Fixed model parameters

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
quarantine_efficacy_scan <- c(0, 0.35, 0.65) # from the same article as above

### SARS-CoV-1-specific parameters (longer Tg, lower presymptomatic transmission, lower proportion asymptomatic infection)
SC1_generation_time <- function(n) { rgamma(n, shape = 24, rate = 2) } 
SC1_infection_to_onset <- function(n) { rgamma(n, shape = 0.1, rate = 1) } 
SC1_prop_asymptomatic <- 0

### SARS-CoV-2-specific parameters (shorter Tg, higher presymptomatic transmission, higher proportion asymptomatic infection)
SC2_generation_time <- function(n) { rgamma(n, shape = 13.5, rate = 2) } # 6.75 day generation time Gamam distributed (as per Walker et al, Science, 2020)
SC2_infection_to_onset <- function(n) { rgamma(n, shape = 13.5/3, rate = 2) } # ~35% of transmission presymptomatic (per SARS-CoV-2, slightly lower than but roughly aligned with: https://bmjopen.bmj.com/content/11/6/e041240)
SC2_prop_asymptomatic <- 0.15
SC2_isolation_Tg_fraction <- unname((fit_gamma_fdplus$estimate["shape"] / fit_gamma_fdplus$estimate["rate"]) / (13.5 / 2)) # isolation time is what fraction of generation time on average

### Vaccine-related parameters
vaccine_start <- 21 # detection + logistical delay of 3 weeks to vaccination initiation
vaccine_coverage <- 0.8
vaccine_efficacy_infection_scan <- c(0.35, 0.75) 
vaccine_efficacy_transmission_scan <- c(0.35, 0.5) 
vaccine_logistical_delay <- 2

### Other parameters
pop <- 10^10
check_final_size <- 2500
initial_immune <- 0
seeding_cases <- 5
iterations <- 20
R0_scan <- c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5)

### Setting up the cluster for parallel running
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

## R0 sensitivity analysis (Figure 1B)
fresh_run_R0_sensitivity_analysis <- TRUE
n <- 2000
if (fresh_run_R0_sensitivity_analysis) {
  
  ### Final Epidemic Size - Setting up R0 scan and the storage for each of the different protection delays
  SC1_storage_nothing <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_2weeks <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_1week <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_2days <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_instant <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  
  SC2_storage_nothing <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_2weeks <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_1week <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_2days <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_instant <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  
  ### Time to nth Case - Setting up R0 scan and the storage for each of the different protection delays
  SC1_storage_nothing_nth_day <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_2weeks_nth_day <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_1week_nth_day <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_2days_nth_day <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_instant_nth_day <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  
  SC2_storage_nothing_nth_day <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_2weeks_nth_day <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_1week_nth_day <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_2days_nth_day <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_instant_nth_day <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  
  ### R0 - Setting up R0 scan and the storage for each of the different protection delays
  SC1_storage_nothing_R0 <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_2weeks_R0 <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_1week_R0 <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_2days_R0 <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_instant_R0 <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  
  SC2_storage_nothing_R0 <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_2weeks_R0 <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_1week_R0 <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_2days_R0 <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_instant_R0 <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  
  ### Reff - Setting up R0 scan and the storage for each of the different protection delays
  SC1_storage_nothing_Reff <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_2weeks_Reff <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_1week_Reff <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_2days_Reff <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC1_storage_vacc_instant_Reff <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  
  SC2_storage_nothing_Reff <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_2weeks_Reff <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_1week_Reff <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_2days_Reff <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  SC2_storage_vacc_instant_Reff <- array(data = NA, dim = c(iterations, length(R0_scan), length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan)))
  
  set.seed(2000)
  seeds <- runif(n = iterations, min = 1, max = 10^9)
  
  for (i in 1:length(R0_scan)) {
    
    for (j in 1:length(vaccine_efficacy_infection_scan)) {
    
      for (k in 1:length(quarantine_efficacy_scan)) {
        
        out_list <- foreach(l = seq_len(iterations), # .combine = 'list',
                            # .multicombine = TRUE, 
                            .export = c("ring_vax_bp_sim", "time_to_nth_infection"),
                            .packages = c("dplyr", "tidyr")) %dopar% {
            
            # Grab the seed for this iteration
            this_seed <- seeds[l]
            
            # ----------------------------------------------------------------
            # SC1_vacc_2weeks
            # ----------------------------------------------------------------
            SC1_vacc_2weeks <- ring_vax_bp_sim(
              offspring = "pois",
              mn_offspring = R0_scan[i],
              generation_time = SC1_generation_time,
              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
              initial_immune = initial_immune,
              seeding_cases = seeding_cases, 
              seed = this_seed,
              prop_asymptomatic = SC1_prop_asymptomatic,
              infection_to_onset = SC1_infection_to_onset,
              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
              vaccine_efficacy_infection = vaccine_efficacy_infection_scan[j],
              vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[j],
              vaccine_logistical_delay = vaccine_logistical_delay,
              vaccine_protection_delay = 14,
              time_to_quarantine = quarantine_time,
              prob_quarantine_contact_traced = prob_quarantine_contact_traced,
              prob_quarantine_symptoms = prob_quarantine_symptoms,
              quarantine_efficacy = quarantine_efficacy_scan[k]
            )
            sc1_2w_size <- sum(!is.na(SC1_vacc_2weeks$time_infection))
            sc1_2w_to_n <- time_to_nth_infection(tdf = SC1_vacc_2weeks, n = n)[[1]]
            sc1_2w_Reff <- calculate_Reff(SC1_vacc_2weeks)
            sc1_2w_R0 <- calculate_R0(SC1_vacc_2weeks)
            
            # ----------------------------------------------------------------
            # SC2_vacc_2weeks
            # ----------------------------------------------------------------
            SC2_vacc_2weeks <- ring_vax_bp_sim(
              offspring = "pois",
              mn_offspring = R0_scan[i],
              generation_time = SC2_generation_time,
              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
              initial_immune = initial_immune,
              seeding_cases = seeding_cases, 
              seed = this_seed,
              prop_asymptomatic = SC2_prop_asymptomatic,
              infection_to_onset = SC2_infection_to_onset,
              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
              vaccine_efficacy_infection = vaccine_efficacy_infection_scan[j],
              vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[j],
              vaccine_logistical_delay = vaccine_logistical_delay,
              vaccine_protection_delay = 14,
              time_to_quarantine = quarantine_time,
              prob_quarantine_contact_traced = prob_quarantine_contact_traced,
              prob_quarantine_symptoms = prob_quarantine_symptoms,
              quarantine_efficacy = quarantine_efficacy_scan[k]
            )
            sc2_2w_size   <- sum(!is.na(SC2_vacc_2weeks$time_infection))
            sc2_2w_to_n    <- time_to_nth_infection(tdf = SC2_vacc_2weeks, n = n)[[1]]
            sc2_2w_Reff <- calculate_Reff(SC2_vacc_2weeks)
            sc2_2w_R0 <- calculate_R0(SC2_vacc_2weeks)
            
            # ----------------------------------------------------------------
            # SC1_vacc_1week
            # ----------------------------------------------------------------
            SC1_vacc_1week <- ring_vax_bp_sim(
              offspring = "pois",
              mn_offspring = R0_scan[i],
              generation_time = SC1_generation_time,
              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
              initial_immune = initial_immune,
              seeding_cases = seeding_cases, 
              seed = this_seed,
              prop_asymptomatic = SC1_prop_asymptomatic,
              infection_to_onset = SC1_infection_to_onset,
              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
              vaccine_efficacy_infection = vaccine_efficacy_infection_scan[j],
              vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[j],
              vaccine_logistical_delay = vaccine_logistical_delay,
              vaccine_protection_delay = 7,
              time_to_quarantine = quarantine_time,
              prob_quarantine_contact_traced = prob_quarantine_contact_traced,
              prob_quarantine_symptoms = prob_quarantine_symptoms,
              quarantine_efficacy = quarantine_efficacy_scan[k]
            )
            sc1_1w_size <- sum(!is.na(SC1_vacc_1week$time_infection))
            sc1_1w_to_n  <- time_to_nth_infection(tdf = SC1_vacc_1week, n = n)[[1]]
            sc1_1w_Reff <- calculate_Reff(SC1_vacc_1week)
            sc1_1w_R0 <- calculate_R0(SC1_vacc_1week)
            
            # ----------------------------------------------------------------
            # SC2_vacc_1week
            # ----------------------------------------------------------------
            SC2_vacc_1week <- ring_vax_bp_sim(
              offspring = "pois",
              mn_offspring = R0_scan[i],
              generation_time = SC2_generation_time,
              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
              initial_immune = initial_immune,
              seeding_cases = seeding_cases, 
              seed = this_seed,
              prop_asymptomatic = SC2_prop_asymptomatic,
              infection_to_onset = SC2_infection_to_onset,
              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
              vaccine_efficacy_infection = vaccine_efficacy_infection_scan[j],
              vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[j],
              vaccine_logistical_delay = vaccine_logistical_delay,
              vaccine_protection_delay = 7,
              time_to_quarantine = quarantine_time,
              prob_quarantine_contact_traced = prob_quarantine_contact_traced,
              prob_quarantine_symptoms = prob_quarantine_symptoms,
              quarantine_efficacy = quarantine_efficacy_scan[k]
            )
            sc2_1w_size <- sum(!is.na(SC2_vacc_1week$time_infection))
            sc2_1w_to_n  <- time_to_nth_infection(tdf = SC2_vacc_1week, n = n)[[1]]
            sc2_1w_Reff <- calculate_Reff(SC2_vacc_1week)
            sc2_1w_R0 <- calculate_R0(SC2_vacc_1week)
            
            # ----------------------------------------------------------------
            # SC1_vacc_2days
            # ----------------------------------------------------------------
            SC1_vacc_2days <- ring_vax_bp_sim(
              offspring = "pois",
              mn_offspring = R0_scan[i],
              generation_time = SC1_generation_time,
              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
              initial_immune = initial_immune,
              seeding_cases = seeding_cases, 
              seed = this_seed,
              prop_asymptomatic = SC1_prop_asymptomatic,
              infection_to_onset = SC1_infection_to_onset,
              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
              vaccine_efficacy_infection = vaccine_efficacy_infection_scan[j],
              vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[j],
              vaccine_logistical_delay = vaccine_logistical_delay,
              vaccine_protection_delay = 2,
              time_to_quarantine = quarantine_time,
              prob_quarantine_contact_traced = prob_quarantine_contact_traced,
              prob_quarantine_symptoms = prob_quarantine_symptoms,
              quarantine_efficacy = quarantine_efficacy_scan[k]
            )
            sc1_2d_size <- sum(!is.na(SC1_vacc_2days$time_infection))
            sc1_2d_to_n  <- time_to_nth_infection(tdf = SC1_vacc_2days, n = n)[[1]]
            sc1_2d_Reff <- calculate_Reff(SC1_vacc_2days)
            sc1_2d_R0 <- calculate_R0(SC1_vacc_2days)
          
            # ----------------------------------------------------------------
            # SC2_vacc_2days
            # ----------------------------------------------------------------
            SC2_vacc_2days <- ring_vax_bp_sim(
              offspring = "pois",
              mn_offspring = R0_scan[i],
              generation_time = SC2_generation_time,
              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
              initial_immune = initial_immune,
              seeding_cases = seeding_cases, 
              seed = this_seed,
              prop_asymptomatic = SC2_prop_asymptomatic,
              infection_to_onset = SC2_infection_to_onset,
              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
              vaccine_efficacy_infection = vaccine_efficacy_infection_scan[j],
              vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[j],
              vaccine_logistical_delay = vaccine_logistical_delay,
              vaccine_protection_delay = 2,
              time_to_quarantine = quarantine_time,
              prob_quarantine_contact_traced = prob_quarantine_contact_traced,
              prob_quarantine_symptoms = prob_quarantine_symptoms,
              quarantine_efficacy = quarantine_efficacy_scan[k]
            )
            sc2_2d_size <- sum(!is.na(SC2_vacc_2days$time_infection))
            sc2_2d_to_n  <- time_to_nth_infection(tdf = SC2_vacc_2days, n = n)[[1]]
            sc2_2d_Reff <- calculate_Reff(SC2_vacc_2days)
            sc2_2d_R0 <- calculate_R0(SC2_vacc_2days)
            
            # ----------------------------------------------------------------
            # SC1_vacc_instant
            # ----------------------------------------------------------------
            SC1_vacc_instant <- ring_vax_bp_sim(
              offspring = "pois",
              mn_offspring = R0_scan[i],
              generation_time = SC1_generation_time,
              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
              initial_immune = initial_immune,
              seeding_cases = seeding_cases, 
              seed = this_seed,
              prop_asymptomatic = SC1_prop_asymptomatic,
              infection_to_onset = SC1_infection_to_onset,
              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
              vaccine_efficacy_infection = vaccine_efficacy_infection_scan[j],
              vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[j],
              vaccine_logistical_delay = vaccine_logistical_delay,
              vaccine_protection_delay = 0.01,
              time_to_quarantine = quarantine_time,
              prob_quarantine_contact_traced = prob_quarantine_contact_traced,
              prob_quarantine_symptoms = prob_quarantine_symptoms,
              quarantine_efficacy = quarantine_efficacy_scan[k]
            )
            sc1_inst_size <- sum(!is.na(SC1_vacc_instant$time_infection))
            sc1_inst_to_n  <- time_to_nth_infection(tdf = SC1_vacc_instant, n = n)[[1]]
            sc1_inst_Reff <- calculate_Reff(SC1_vacc_instant)
            sc1_inst_R0 <- calculate_R0(SC1_vacc_instant)
            
            # ----------------------------------------------------------------
            # SC2_vacc_instant
            # ----------------------------------------------------------------
            SC2_vacc_instant <- ring_vax_bp_sim(
              offspring = "pois",
              mn_offspring = R0_scan[i],
              generation_time = SC2_generation_time,
              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
              initial_immune = initial_immune,
              seeding_cases = seeding_cases, 
              seed = this_seed,
              prop_asymptomatic = SC2_prop_asymptomatic,
              infection_to_onset = SC2_infection_to_onset,
              vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
              vaccine_efficacy_infection = vaccine_efficacy_infection_scan[j],
              vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[j],
              vaccine_logistical_delay = vaccine_logistical_delay,
              vaccine_protection_delay = 0.01,
              time_to_quarantine = quarantine_time,
              prob_quarantine_contact_traced = prob_quarantine_contact_traced,
              prob_quarantine_symptoms = prob_quarantine_symptoms,
              quarantine_efficacy = quarantine_efficacy_scan[k]
            )
            sc2_inst_size <- sum(!is.na(SC2_vacc_instant$time_infection))
            sc2_inst_to_n  <- time_to_nth_infection(tdf = SC2_vacc_instant, n = n)[[1]]
            sc2_inst_Reff <- calculate_Reff(SC2_vacc_instant)
            sc2_inst_R0 <- calculate_R0(SC2_vacc_instant)
            
            # ----------------------------------------------------------------
            # SC1_no_vacc
            # ----------------------------------------------------------------
            SC1_no_vacc <- ring_vax_bp_sim(
              offspring = "pois",
              mn_offspring = R0_scan[i],
              generation_time = SC1_generation_time,
              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
              initial_immune = initial_immune,
              seeding_cases = seeding_cases, 
              seed = this_seed,
              prop_asymptomatic = SC1_prop_asymptomatic,
              infection_to_onset = SC1_infection_to_onset,
              vaccine_start = 1000, vaccine_coverage = vaccine_coverage,
              vaccine_efficacy_infection = vaccine_efficacy_infection_scan[j],
              vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[j],
              vaccine_logistical_delay = vaccine_logistical_delay,
              vaccine_protection_delay = 100,
              time_to_quarantine = quarantine_time,
              prob_quarantine_contact_traced = prob_quarantine_contact_traced,
              prob_quarantine_symptoms = prob_quarantine_symptoms,
              quarantine_efficacy = quarantine_efficacy_scan[k]
            )
            sc1_nothing_size <- sum(!is.na(SC1_no_vacc$time_infection))
            sc1_nothing_to_n  <- time_to_nth_infection(tdf = SC1_no_vacc, n = n)[[1]]
            sc1_nothing_Reff <- calculate_Reff(SC1_no_vacc)
            sc1_nothing_R0 <- calculate_R0(SC1_no_vacc)
            
            # ----------------------------------------------------------------
            # SC2_no_vacc
            # ----------------------------------------------------------------
            SC2_no_vacc <- ring_vax_bp_sim(
              offspring = "pois",
              mn_offspring = R0_scan[i],
              generation_time = SC2_generation_time,
              t0 = 0, tf = Inf, pop = pop, check_final_size = check_final_size, 
              initial_immune = initial_immune,
              seeding_cases = seeding_cases, 
              seed = this_seed,
              prop_asymptomatic = SC2_prop_asymptomatic,
              infection_to_onset = SC2_infection_to_onset,
              vaccine_start = 1000, vaccine_coverage = vaccine_coverage,
              vaccine_efficacy_infection = vaccine_efficacy_infection_scan[j],
              vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[j],
              vaccine_logistical_delay = vaccine_logistical_delay,
              vaccine_protection_delay = 100,
              time_to_quarantine = quarantine_time,
              prob_quarantine_contact_traced = prob_quarantine_contact_traced,
              prob_quarantine_symptoms = prob_quarantine_symptoms,
              quarantine_efficacy = quarantine_efficacy_scan[k]
            )
            sc2_nothing_size <- sum(!is.na(SC2_no_vacc$time_infection))
            sc2_nothing_to_n  <- time_to_nth_infection(tdf = SC2_no_vacc, n = n)[[1]]
            sc2_nothing_Reff <- calculate_Reff(SC2_no_vacc)
            sc2_nothing_R0 <- calculate_R0(SC2_no_vacc)
            
            # Return a small named list with everything needed
            x <- list(sc1_2w_size = sc1_2w_size, sc1_2w_to_n = sc1_2w_to_n, sc1_2w_Reff = sc1_2w_Reff, sc1_2w_R0 = sc1_2w_R0, 
                 sc2_2w_size = sc2_2w_size, sc2_2w_to_n = sc2_2w_to_n, sc2_2w_Reff = sc2_2w_Reff, sc2_2w_R0 = sc2_2w_R0,
                 sc1_1w_size = sc1_1w_size, sc1_1w_to_n = sc1_1w_to_n, sc1_1w_Reff = sc1_1w_Reff, sc1_1w_R0 = sc1_1w_R0,
                 sc2_1w_size = sc2_1w_size, sc2_1w_to_n = sc2_1w_to_n, sc2_1w_Reff = sc2_1w_Reff, sc2_1w_R0 = sc2_1w_R0,
                 sc1_2d_size = sc1_2d_size, sc1_2d_to_n = sc1_2d_to_n, sc1_2d_Reff = sc1_2d_Reff, sc1_2d_R0 = sc1_2d_R0,
                 sc2_2d_size = sc2_2d_size, sc2_2d_to_n = sc2_2d_to_n, sc2_2d_Reff = sc2_2d_Reff, sc2_2d_R0 = sc2_2d_R0,
                 sc1_inst_size = sc1_inst_size, sc1_inst_to_n = sc1_inst_to_n, sc1_inst_Reff = sc1_inst_Reff, sc1_inst_R0 = sc1_inst_R0,
                 sc2_inst_size = sc2_inst_size, sc2_inst_to_n = sc2_inst_to_n, sc2_inst_Reff = sc2_inst_Reff, sc2_inst_R0 = sc2_inst_R0,
                 sc1_nothing_size = sc1_nothing_size, sc1_nothing_to_n = sc1_nothing_to_n, sc1_nothing_Reff = sc1_nothing_Reff, sc1_nothing_R0 = sc1_nothing_R0,
                 sc2_nothing_size = sc2_nothing_size, sc2_nothing_to_n = sc2_nothing_to_n, sc2_nothing_Reff = sc2_nothing_Reff, sc2_nothing_R0 = sc2_nothing_R0)
            return(x)
        }
        
        for (l in seq_len(iterations)) {
          tmp <- out_list[[l]]
          SC1_storage_vacc_2weeks[l, i, j, k] <- tmp$sc1_2w_size
          SC1_storage_vacc_2weeks_nth_day[l, i, j, k] <- tmp$sc1_2w_to_n
          SC1_storage_vacc_2weeks_R0[l, i, j, k] <- tmp$sc1_2w_R0
          SC1_storage_vacc_2weeks_Reff[l, i, j, k] <- tmp$sc1_2w_Reff
          
          SC2_storage_vacc_2weeks[l, i, j, k] <- tmp$sc2_2w_size
          SC2_storage_vacc_2weeks_nth_day[l, i, j, k] <- tmp$sc2_2w_to_n
          SC2_storage_vacc_2weeks_R0[l, i, j, k] <- tmp$sc2_2w_R0
          SC2_storage_vacc_2weeks_Reff[l, i, j, k] <- tmp$sc2_2w_Reff
          
          SC1_storage_vacc_1week[l, i, j, k] <- tmp$sc1_1w_size
          SC1_storage_vacc_1week_nth_day[l, i, j, k] <- tmp$sc1_1w_to_n
          SC1_storage_vacc_1week_R0[l, i, j, k] <- tmp$sc1_1w_R0
          SC1_storage_vacc_1week_Reff[l, i, j, k] <- tmp$sc1_1w_Reff
          
          SC2_storage_vacc_1week[l, i, j, k] <- tmp$sc2_1w_size
          SC2_storage_vacc_1week_nth_day[l, i, j, k]  <- tmp$sc2_1w_to_n
          SC2_storage_vacc_1week_R0[l, i, j, k] <- tmp$sc2_1w_R0
          SC2_storage_vacc_1week_Reff[l, i, j, k] <- tmp$sc2_1w_Reff
          
          SC1_storage_vacc_2days[l, i, j, k] <- tmp$sc1_2d_size
          SC1_storage_vacc_2days_nth_day[l, i, j, k]  <- tmp$sc1_2d_to_n
          SC1_storage_vacc_2days_R0[l, i, j, k] <- tmp$sc1_2d_R0
          SC1_storage_vacc_2days_Reff[l, i, j, k] <- tmp$sc1_2d_Reff
          
          SC2_storage_vacc_2days[l, i, j, k] <- tmp$sc2_2d_size
          SC2_storage_vacc_2days_nth_day[l, i, j, k]  <- tmp$sc2_2d_to_n
          SC2_storage_vacc_2days_R0[l, i, j, k] <- tmp$sc2_2d_R0
          SC2_storage_vacc_2days_Reff[l, i, j, k] <- tmp$sc2_2d_Reff
          
          SC1_storage_vacc_instant[l, i, j, k] <- tmp$sc1_inst_size
          SC1_storage_vacc_instant_nth_day[l, i, j, k]<- tmp$sc1_inst_to_n
          SC1_storage_vacc_instant_R0[l, i, j, k] <- tmp$sc1_inst_R0
          SC1_storage_vacc_instant_Reff[l, i, j, k] <- tmp$sc1_inst_Reff
          
          SC2_storage_vacc_instant[l, i, j, k] <- tmp$sc2_inst_size
          SC2_storage_vacc_instant_nth_day[l, i, j, k]<- tmp$sc2_inst_to_n
          SC2_storage_vacc_instant_R0[l, i, j, k] <- tmp$sc2_inst_R0
          SC2_storage_vacc_instant_Reff[l, i, j, k] <- tmp$sc2_inst_Reff
          
          SC1_storage_nothing[l, i, j, k] <- tmp$sc1_nothing_size
          SC1_storage_nothing_nth_day[l, i, j, k] <- tmp$sc1_nothing_to_n
          SC1_storage_nothing_R0[l, i, j, k] <- tmp$sc1_nothing_R0
          SC1_storage_nothing_Reff[l, i, j, k] <- tmp$sc1_nothing_Reff
          
          SC2_storage_nothing[l, i, j, k] <- tmp$sc2_nothing_size
          SC2_storage_nothing_nth_day[l, i, j, k] <- tmp$sc2_nothing_to_n
          SC2_storage_nothing_R0[l, i, j, k] <- tmp$sc2_nothing_R0
          SC2_storage_nothing_Reff[l, i, j, k] <- tmp$sc2_nothing_Reff
          
        }
        print(paste0("k = ", k))
      }
      print(paste0("j = ", j))
    }
    print(paste0("i = ", i))
  }
  stopCluster(cl)
  
  # Convert each array to a long data frame - epidemic size
  df_SC1_nothing <- convert_array(arr = SC1_storage_nothing, scenario = "no_vaccination", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_nothing <- convert_array(SC2_storage_nothing, scenario = "no_vaccination", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_2weeks <- convert_array(SC1_storage_vacc_2weeks, scenario = "2weeks_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_2weeks <- convert_array(SC2_storage_vacc_2weeks, scenario = "2weeks_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_1week <- convert_array(SC1_storage_vacc_1week, scenario = "1week_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_1week <- convert_array(SC2_storage_vacc_1week, scenario = "1week_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_2days <- convert_array(SC1_storage_vacc_2days, scenario = "2days_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_2days <- convert_array(SC2_storage_vacc_2days, scenario = "2days_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_instant <- convert_array(SC1_storage_vacc_instant, scenario = "no_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_instant <- convert_array(SC2_storage_vacc_instant, scenario = "no_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)

  # Convert each array to a long data frame - epidemic size
  df_SC1_nothing_time_to_n <- convert_array_time_to_n(arr = SC1_storage_nothing_nth_day, scenario = "no_vaccination", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_nothing_time_to_n <- convert_array_time_to_n(SC2_storage_nothing_nth_day, scenario = "no_vaccination", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_2weeks_time_to_n <- convert_array_time_to_n(SC1_storage_vacc_2weeks_nth_day, scenario = "2weeks_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_2weeks_time_to_n <- convert_array_time_to_n(SC2_storage_vacc_2weeks_nth_day, scenario = "2weeks_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_1week_time_to_n <- convert_array_time_to_n(SC1_storage_vacc_1week_nth_day, scenario = "1week_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_1week_time_to_n <- convert_array_time_to_n(SC2_storage_vacc_1week_nth_day, scenario = "1week_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_2days_time_to_n <- convert_array_time_to_n(SC1_storage_vacc_2days_nth_day, scenario = "2days_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_2days_time_to_n <- convert_array_time_to_n(SC2_storage_vacc_2days_nth_day, scenario = "2days_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_instant_time_to_n <- convert_array_time_to_n(SC1_storage_vacc_instant_nth_day, scenario = "no_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_instant_time_to_n <- convert_array_time_to_n(SC2_storage_vacc_instant_nth_day, scenario = "no_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  
  # Convert each array to a long data frame - Reff
  df_SC1_nothing_Reff <- convert_array_Reff(arr = SC1_storage_nothing_Reff, scenario = "no_vaccination", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_nothing_Reff <- convert_array_Reff(SC2_storage_nothing_Reff, scenario = "no_vaccination", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_2weeks_Reff <- convert_array_Reff(SC1_storage_vacc_2weeks_Reff, scenario = "2weeks_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_2weeks_Reff <- convert_array_Reff(SC2_storage_vacc_2weeks_Reff, scenario = "2weeks_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_1week_Reff <- convert_array_Reff(SC1_storage_vacc_1week_Reff, scenario = "1week_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_1week_Reff <- convert_array_Reff(SC2_storage_vacc_1week_Reff, scenario = "1week_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_2days_Reff <- convert_array_Reff(SC1_storage_vacc_2days_Reff, scenario = "2days_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_2days_Reff <- convert_array_Reff(SC2_storage_vacc_2days_Reff, scenario = "2days_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_instant_Reff <- convert_array_Reff(SC1_storage_vacc_instant_Reff, scenario = "no_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_instant_Reff <- convert_array_Reff(SC2_storage_vacc_instant_Reff, scenario = "no_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  
  # Convert each array to a long data frame - R0
  df_SC1_nothing_R0 <- convert_array_R0(arr = SC1_storage_nothing_R0, scenario = "no_vaccination", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_nothing_R0 <- convert_array_R0(SC2_storage_nothing_R0, scenario = "no_vaccination", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_2weeks_R0 <- convert_array_R0(SC1_storage_vacc_2weeks_R0, scenario = "2weeks_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_2weeks_R0 <- convert_array_R0(SC2_storage_vacc_2weeks_R0, scenario = "2weeks_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_1week_R0 <- convert_array_R0(SC1_storage_vacc_1week_R0, scenario = "1week_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_1week_R0 <- convert_array_R0(SC2_storage_vacc_1week_R0, scenario = "1week_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_2days_R0 <- convert_array_R0(SC1_storage_vacc_2days_R0, scenario = "2days_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_2days_R0 <- convert_array_R0(SC2_storage_vacc_2days_R0, scenario = "2days_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC1_instant_R0 <- convert_array_R0(SC1_storage_vacc_instant_R0, scenario = "no_delay", pathogen = "SARS-CoV-1", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  df_SC2_instant_R0 <- convert_array_R0(SC2_storage_vacc_instant_R0, scenario = "no_delay", pathogen = "SARS-CoV-2", R0_scan, vaccine_efficacy_transmission_scan, vaccine_efficacy_infection_scan, quarantine_efficacy_scan)
  
  # Combine all results into a single data frame.
  overall_bp_df <- bind_rows(
    df_SC1_nothing, df_SC1_2weeks, df_SC1_1week, df_SC1_2days, df_SC1_instant,
    df_SC2_nothing, df_SC2_2weeks, df_SC2_1week, df_SC2_2days, df_SC2_instant) %>%
    left_join(bind_rows(
      df_SC1_nothing_time_to_n, df_SC1_2weeks_time_to_n, df_SC1_1week_time_to_n, df_SC1_2days_time_to_n, df_SC1_instant_time_to_n,
      df_SC2_nothing_time_to_n, df_SC2_2weeks_time_to_n, df_SC2_1week_time_to_n, df_SC2_2days_time_to_n, df_SC2_instant_time_to_n), 
      by = c("iteration", "scenario", "pathogen", "R0", "vaccine_efficacy_infection", "vaccine_efficacy_transmission", "quarantine_efficacy")) %>%
    left_join(bind_rows(
      df_SC1_nothing_Reff, df_SC1_2weeks_Reff, df_SC1_1week_Reff, df_SC1_2days_Reff, df_SC1_instant_Reff,
      df_SC2_nothing_Reff, df_SC2_2weeks_Reff, df_SC2_1week_Reff, df_SC2_2days_Reff, df_SC2_instant_Reff), 
      by = c("iteration", "scenario", "pathogen", "R0", "vaccine_efficacy_infection", "vaccine_efficacy_transmission", "quarantine_efficacy")) %>%
    left_join(bind_rows(
      df_SC1_nothing_R0, df_SC1_2weeks_R0, df_SC1_1week_R0, df_SC1_2days_R0, df_SC1_instant_R0,
      df_SC2_nothing_R0, df_SC2_2weeks_R0, df_SC2_1week_R0, df_SC2_2days_R0, df_SC2_instant_R0), 
      by = c("iteration", "scenario", "pathogen", "R0", "vaccine_efficacy_infection", "vaccine_efficacy_transmission", "quarantine_efficacy"))
  
  saveRDS(object = overall_bp_df,
          file = "outputs/Figure1_branchingProcess_Containment/Fig1_ringVaccination_paramScan.rds")
} else {
  overall_bp_df <- readRDS("outputs/Figure1_branchingProcess_Containment/Fig1_ringVaccination_paramScan.rds")
}

## Plotting the proportion of outbreaks controlled
containment_df <- overall_bp_df %>%
  mutate(contained = ifelse(epidemic_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(R0, pathogen, vaccine_efficacy_infection, quarantine_efficacy) %>%
  mutate(time_to_n_relative = time_to_n / time_to_n[scenario == "no_vaccination"]) %>%
  group_by(R0, scenario, pathogen, vaccine_efficacy_infection, quarantine_efficacy) %>%
  summarise(proportion_contained = sum(contained) / iterations,
            R0_actual = mean(R0_actual, na.rm = TRUE),
            Reff_mean = median(Reff, na.rm = TRUE),
            Reff_lower = quantile(Reff, 0.1, na.rm = TRUE),
            Reff_upper = quantile(Reff, 0.9, na.rm = TRUE),
            time_to_n = mean(time_to_n, na.rm = TRUE),
            time_to_n_relative = mean(time_to_n_relative, na.rm = TRUE)) %>%
  mutate(proportion_contained = ifelse(R0 == 1.00, 1, proportion_contained)) %>%
  mutate(scenario = ifelse(scenario == "no_vaccination", "zno_vaccination", scenario)) %>%
  mutate(scenario = ifelse(scenario == "2days_delay", "bvacc_2days_protectDelay", scenario)) %>%
  mutate(scenario = ifelse(scenario == "1week_delay", "dvacc_1week_protectDelay", scenario)) %>%
  mutate(scenario = ifelse(scenario == "2weeks_delay", "evacc_2weeks_protectDelay", scenario)) %>%
  mutate(scenario = ifelse(scenario == "no_delay", "avacc_no_delay", scenario))
containment_df$scenario <- factor(containment_df$scenario, 
                                  levels = c("avacc_no_delay", "bvacc_2days_protectDelay", "dvacc_1week_protectDelay", 
                                             "evacc_2weeks_protectDelay", "zno_vaccination"))
containment_df2 <- containment_df %>%
  arrange(scenario) 

## Plotting Fig1B - Reff and Proportion of Outbreaks Contained
combo_tbl  <- expand_grid(pathogen =  c("SARS-CoV-2", "SARS-CoV-1"), quarantine_efficacy = c(0, 0.65))
combo_plot <- combo_tbl %>% 
  mutate(p = map2(pathogen, quarantine_efficacy, ~ make_stacked_plot_Reff(containment_df2, .x, .y, c(1, 2))))

Fig1BCDE <- plot_grid(plotlist = list(combo_plot$p[[2]], combo_plot$p[[1]], combo_plot$p[[4]], combo_plot$p[[3]]),
                      nrow = length(c(0, 0.65)), ncol  = length(pathogens),
                      labels = c("B", "C", "D", "E"), label_size = 10)

## Plotting Supplementary Figure looking at time to epidemic threshold
containment_df2$vaccine_quarantine_elision <- paste0("Vaccine Effiacy = ", containment_df2$vaccine_efficacy_infection, "\nQuarantine Effiacy = ", containment_df2$quarantine_efficacy)
time_to_n_plot <- ggplot(subset(containment_df2, quarantine_efficacy != 0.35),
       aes(x = R0, y = time_to_n_relative, col = scenario)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_grid(vaccine_quarantine_elision~pathogen) + 
  scale_colour_manual(values = c("#CA2E6B", "#88C5EE", "#236897", "#13496E", "black"),
                      labels = c("No Delay", "2 Days", "1 Week", "2 Weeks", "No Vaccination"),
                      name = "Vaccine\nProtection\nDelay",
                      guide = guide_legend(reverse = TRUE)) +
  labs(x = "R0", y = "Fold Increase in Time to Epidemic Threshold") +
  theme(strip.background = element_rect(fill = "white"))
ggsave(plot = time_to_n_plot, filename = "figures/Figure_1_BranchingProcess/FigS1_ParamScan_timetoN.pdf", height = 8.5, width = 8)

####################################################################################################################################
## Vaccination-Related Sensitivity Analyses Heatmaps
### Note that old paper results were with 0.75 vaccine efficacy against infection and 0.5 against onwards transmission. 
### Have updated that here.
####################################################################################################################################
fresh_run_vaccination_heatmaps <- TRUE
if (fresh_run_vaccination_heatmaps) {
  
  set.seed(2000)
  seeds <- runif(n = iterations, min = 1, max = 10^9)
  outcome_names <- c("epidemic_size", "time_to_n", "Reff", "R0")
  
  ## R0 sensitivity analysis for all the runs
  R0_seq <- R0_scan[R0_scan > 1]
  vaccine_protection_delay <- 7
  vaccine_efficacy_infection_scan <- c(0, 0.35, 0.75)
  vaccine_efficacy_transmission_scan <- c(0, 0.35, 0.75)
  
  #######################################################################
  ## Sensitivity Analysis - R0 vs Ratio of Tg to Protection Delay
  #######################################################################
  num_cores <- min(iterations, parallel::detectCores() - 1)
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  Tg_ratio_seq <- seq(1, 4, 0.5)
  storage_R0_TgRatio_sensitivity <- array(data = NA, dim = c(iterations, length(R0_seq), length(Tg_ratio_seq), 
                                                             length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan), 4))
  for (i in 1:length(R0_seq)) {
    for (j in 1:length(Tg_ratio_seq)) {
      for (m in seq_along(vaccine_efficacy_infection_scan)) {
        for (k in 1:length(quarantine_efficacy_scan)) {
          
          # Parallelize over l (iterations) with foreach
          out_list <- foreach(
            l = seq_len(iterations),
            .export       = c("ring_vax_bp_sim", "quarantine_time_closure"),
            .packages     = c("stats", "dplyr", "tidyr")  # 'stats' for rgamma, if needed
          ) %dopar% {
            
            generation_time <- function(n) { rgamma(n, shape = 2 * vaccine_protection_delay * Tg_ratio_seq[j], rate = 2) } 
            infection_to_onset <- function(n) { rgamma(n, shape = (2 * vaccine_protection_delay * Tg_ratio_seq[j])/3, rate = 2) } # keeping proportion of presymptomatic transmission constant as Tg varies
            quarantine_time <- quarantine_time_closure(quarantine_time_shape = SC2_isolation_Tg_fraction * Tg_ratio_seq[j] * vaccine_protection_delay, quarantine_time_rate = 1) # keeping proportion isolating over time the same
            bp_out <- ring_vax_bp_sim(offspring = "pois",
                                      mn_offspring = R0_seq[i],
                                      generation_time = generation_time,
                                      seed = seeds[l],
                                      t0 = 0, tf = Inf, population = pop, 
                                      check_final_size = check_final_size, initial_immune = initial_immune,
                                      seeding_cases = seeding_cases, prop_asymptomatic = SC2_prop_asymptomatic,
                                      infection_to_onset = infection_to_onset,
                                      vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                      vaccine_efficacy_infection = vaccine_efficacy_infection_scan[m],
                                      vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[m],
                                      vaccine_logistical_delay = vaccine_logistical_delay,
                                      vaccine_protection_delay = vaccine_protection_delay,
                                      time_to_quarantine = quarantine_time,
                                      prob_quarantine_contact_traced = prob_quarantine_contact_traced,
                                      prob_quarantine_symptoms = prob_quarantine_symptoms,
                                      quarantine_efficacy = quarantine_efficacy_scan[k])
            size <- sum(!is.na(bp_out$time_infection))
            to_n <- time_to_nth_infection(tdf = bp_out, n = n)[[1]]
            Reff <- calculate_Reff(bp_out)
            R0 <- calculate_R0(bp_out)
            
            x <- list(size = size, to_n = to_n, Reff = Reff, R0 = R0)
            return(x)
          }
          
          for (l in seq_len(iterations)) {
            tmp <- out_list[[l]]
            storage_R0_TgRatio_sensitivity[l, i, j, m, k, 1] <- tmp$size
            storage_R0_TgRatio_sensitivity[l, i, j, m, k, 2] <- tmp$to_n
            storage_R0_TgRatio_sensitivity[l, i, j, m, k, 3] <- tmp$Reff
            storage_R0_TgRatio_sensitivity[l, i, j, m, k, 4] <- tmp$R0
          }
        }
      }
    }
    print(paste("Finished i =", i, "of", length(R0_seq)))
  }
  stopCluster(cl)
  
  ## Processing the simulations
  reshaped_R0_TgRatio_sensitivity <- reshape2::melt(storage_R0_TgRatio_sensitivity)
  colnames(reshaped_R0_TgRatio_sensitivity) <- c("iteration", "R0", "TgRatio", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  reshaped_R0_TgRatio_sensitivity <- reshaped_R0_TgRatio_sensitivity %>%
    mutate(iteration = as.integer(iteration),
           input_R0 = R0_seq[R0],
           Tg_Ratio = Tg_ratio_seq[TgRatio],
           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy], 
           vaccine_efficacy_transmission = vaccine_efficacy_infection_scan[vaccine_efficacy], 
           quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
           outcome = outcome_names[outcome]) %>% 
    dplyr::select(iteration, input_R0, -R0, Tg_Ratio, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
                  quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
    pivot_wider(names_from = "outcome", values_from = value)
  saveRDS(object = reshaped_R0_TgRatio_sensitivity, file = "outputs/Figure1_branchingProcess_Containment/Fig1_ringVaccination_R0TgRatio.rds")
  storage_R0_TgRatio_df <- reshaped_R0_TgRatio_sensitivity
  
  #######################################################################
  ## Sensitivity Analysis - R0 vs Vaccine Efficacy
  #######################################################################
  num_cores <- parallel::detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  Tg_ratio_fixed <- 2.5
  generation_time <- function(n) { rgamma(n, shape = 2 * vaccine_protection_delay * Tg_ratio_fixed, rate = 2) }
  infection_to_onset <- function(n) { rgamma(n, shape = (2 * vaccine_protection_delay * Tg_ratio_fixed)/3, rate = 2) }
  vaccine_efficacy_seq <- c(0, seq(0.3, 0.9, 0.1))
  storage_R0_efficacy_sensitivity <- array(data = NA, dim = c(iterations, length(R0_seq), length(vaccine_efficacy_seq), 
                                                              length(quarantine_efficacy_scan), 4))
  for (i in 1:length(R0_seq)) {
    for (j in 1:length(vaccine_efficacy_seq)) {
      for (k in 1:length(quarantine_efficacy_scan)) {
        
        # Parallelize over l (iterations) with foreach
        out_list <- foreach(
          l = seq_len(iterations),
          .export       = c("SC2_isolation_Tg_fraction", "ring_vax_bp_sim", "quarantine_time_closure"),
          .packages     = c("stats", "dplyr", "tidyr")  # 'stats' for rgamma, if needed
        ) %dopar% {
          
          quarantine_time <- quarantine_time_closure(quarantine_time_shape = SC2_isolation_Tg_fraction * Tg_ratio_fixed * vaccine_protection_delay, quarantine_time_rate = 1) # keeping proportion isolating over time the same
          bp_out <- ring_vax_bp_sim(offspring = "pois",
                                    mn_offspring = R0_seq[i],
                                    generation_time = generation_time,
                                    seed = seeds[l],
                                    t0 = 0, tf = Inf, population = pop,
                                    check_final_size = check_final_size, initial_immune = initial_immune,
                                    seeding_cases = seeding_cases, prop_asymptomatic = SC2_prop_asymptomatic,
                                    infection_to_onset = infection_to_onset,
                                    vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                    vaccine_efficacy_infection = vaccine_efficacy_seq[j],
                                    vaccine_efficacy_transmission = vaccine_efficacy_seq[j],
                                    vaccine_logistical_delay = vaccine_logistical_delay,
                                    vaccine_protection_delay = vaccine_protection_delay,
                                    time_to_quarantine = quarantine_time,
                                    prob_quarantine_contact_traced = prob_quarantine_contact_traced,
                                    prob_quarantine_symptoms = prob_quarantine_symptoms,
                                    quarantine_efficacy = quarantine_efficacy_scan[k])
          size <- sum(!is.na(bp_out$time_infection))
          to_n <- time_to_nth_infection(tdf = bp_out, n = n)[[1]]
          Reff <- calculate_Reff(bp_out)
          R0 <- calculate_R0(bp_out)
          
          x <- list(size = size, to_n = to_n, Reff = Reff, R0 = R0)
          return(x)
        }
        
        for (l in seq_len(iterations)) {
          tmp <- out_list[[l]]
          storage_R0_efficacy_sensitivity[l, i, j, k, 1] <- tmp$size
          storage_R0_efficacy_sensitivity[l, i, j, k, 2] <- tmp$to_n
          storage_R0_efficacy_sensitivity[l, i, j, k, 3] <- tmp$Reff
          storage_R0_efficacy_sensitivity[l, i, j, k, 4] <- tmp$R0
        }
      }
    }
    print(paste("Finished i =", i, "of", length(R0_seq)))
  }
  stopCluster(cl)
  
  ## Processing the simulations
  reshaped_R0_efficacy_sensitivity <- reshape2::melt(storage_R0_efficacy_sensitivity)
  colnames(reshaped_R0_efficacy_sensitivity) <- c("iteration", "R0", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  reshaped_R0_efficacy_sensitivity <- reshaped_R0_efficacy_sensitivity %>%
    mutate(iteration = as.integer(iteration),
           input_R0 = R0_seq[R0],
           vaccine_efficacy_infection = vaccine_efficacy_seq[vaccine_efficacy], 
           vaccine_efficacy_transmission = vaccine_efficacy_seq[vaccine_efficacy], 
           quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
           outcome = outcome_names[outcome]) %>% 
    dplyr::select(iteration, input_R0, -R0, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
                  quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
    pivot_wider(names_from = "outcome", values_from = value)
  saveRDS(object = reshaped_R0_efficacy_sensitivity, file = "outputs/Figure1_branchingProcess_Containment/Fig1_ringVaccination_R0Efficacy.rds")
  storage_R0_efficacy_df <- reshaped_R0_efficacy_sensitivity
    
  #######################################################################
  ## Sensitivity Analysis - R0 vs Presymptomatic Transmission Proportion
  #######################################################################
  num_cores <- parallel::detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  Tg_ratio_fixed <- 2.5
  generation_time <- function(n) { rgamma(n, shape = 2 * vaccine_protection_delay * Tg_ratio_fixed, rate = 2) } 
  proportion_presymptomatic_seq <- seq(0.1, 0.7, 0.1)
  storage_R0_preSymp_sensitivity <- array(data = NA, dim = c(iterations, length(R0_seq), length(proportion_presymptomatic_seq), 
                                                             length(vaccine_efficacy_infection_scan), length(quarantine_efficacy_scan), 4))
  for (i in 1:length(R0_seq)) {
    for (j in 1:length(proportion_presymptomatic_seq)) {
      for (m in seq_along(vaccine_efficacy_infection_scan)) {
        for (k in 1:length(quarantine_efficacy_scan)) {
          
          # Parallelize over l (iterations) with foreach
          out_list <- foreach(
            l = seq_len(iterations),
            .export       = c("SC2_isolation_Tg_fraction", "ring_vax_bp_sim", "quarantine_time_closure"),
            .packages     = c("stats", "dplyr", "tidyr")  # 'stats' for rgamma, if needed
          ) %dopar% {
            
            quarantine_time <- quarantine_time_closure(quarantine_time_shape = SC2_isolation_Tg_fraction * Tg_ratio_fixed * vaccine_protection_delay, quarantine_time_rate = 1) # keeping proportion isolating over time the same
            infection_to_onset <- function(n) { rgamma(n, shape = (2 * vaccine_protection_delay * Tg_ratio_fixed) * proportion_presymptomatic_seq[j], rate = 2) }
            bp_out <- ring_vax_bp_sim(offspring = "pois",
                                      mn_offspring = R0_seq[i],
                                      generation_time = generation_time,
                                      seed = seeds[l],
                                      t0 = 0, tf = Inf, pop = 10^10,
                                      check_final_size = check_final_size, initial_immune = initial_immune,
                                      seeding_cases = seeding_cases, prop_asymptomatic = SC2_prop_asymptomatic,
                                      infection_to_onset = infection_to_onset,
                                      vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                      vaccine_efficacy_infection = vaccine_efficacy_infection_scan[m],
                                      vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[m],
                                      vaccine_logistical_delay = vaccine_logistical_delay,
                                      vaccine_protection_delay = vaccine_protection_delay,
                                      time_to_quarantine = quarantine_time,
                                      prob_quarantine_contact_traced = prob_quarantine_contact_traced,
                                      prob_quarantine_symptoms = prob_quarantine_symptoms,
                                      quarantine_efficacy = quarantine_efficacy_scan[k])
            size <- sum(!is.na(bp_out$time_infection))
            to_n <- time_to_nth_infection(tdf = bp_out, n = n)[[1]]
            Reff <- calculate_Reff(bp_out)
            R0 <- calculate_R0(bp_out)
            
            x <- list(size = size, to_n = to_n, Reff = Reff, R0 = R0)
            return(x)
          }
          
          for (l in seq_len(iterations)) {
            tmp <- out_list[[l]]
            storage_R0_preSymp_sensitivity[l, i, j, m, k, 1] <- tmp$size
            storage_R0_preSymp_sensitivity[l, i, j, m, k, 2] <- tmp$to_n
            storage_R0_preSymp_sensitivity[l, i, j, m, k, 3] <- tmp$Reff
            storage_R0_preSymp_sensitivity[l, i, j, m, k, 4] <- tmp$R0
          }
        }
      }
    }
    print(paste("Finished i =", i, "of", length(R0_seq)))
  }
  stopCluster(cl)
  
  ## Processing the simulations
  reshaped_R0_preSymp_sensitivity <- reshape2::melt(storage_R0_preSymp_sensitivity)
  colnames(reshaped_R0_preSymp_sensitivity) <- c("iteration", "R0", "preSymp", "vaccine_efficacy", "quarantine_efficacy", "outcome", "value")
  reshaped_R0_preSymp_sensitivity <- reshaped_R0_preSymp_sensitivity %>%
    mutate(iteration = as.integer(iteration),
           input_R0 = R0_seq[R0],
           prop_preSymp = proportion_presymptomatic_seq[preSymp],
           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[vaccine_efficacy], 
           vaccine_efficacy_transmission = vaccine_efficacy_infection_scan[vaccine_efficacy], 
           quarantine_efficacy = quarantine_efficacy_scan[quarantine_efficacy],
           outcome = outcome_names[outcome]) %>% 
    dplyr::select(iteration, input_R0, -R0, prop_preSymp, vaccine_efficacy_infection, vaccine_efficacy_transmission, 
                  quarantine_efficacy, outcome, value, -vaccine_efficacy) %>%
    pivot_wider(names_from = "outcome", values_from = value)
  saveRDS(object = reshaped_R0_preSymp_sensitivity, file = "outputs/Figure1_branchingProcess_Containment/Fig1_ringVaccination_R0preSymp.rds")
  storage_R0_preSymp_df <- reshaped_R0_preSymp_sensitivity
  
} else {
  storage_R0_TgRatio_df <- readRDS("outputs/Figure1_branchingProcess_Containment/Fig1_ringVaccination_R0TgRatio.rds")
  storage_R0_efficacy_df <- readRDS("outputs/Figure1_branchingProcess_Containment/Fig1_ringVaccination_R0Efficacy.rds")
  storage_R0_preSymp_df <- readRDS("outputs/Figure1_branchingProcess_Containment/Fig1_ringVaccination_R0preSymp.rds")
}

## Creating Vaccination-Related Heatmaps

### Tg Ratio
R0_TgRatio_df <- storage_R0_TgRatio_df %>%
  ungroup() %>%
  mutate(contained = ifelse(epidemic_size < (0.9 * check_final_size), 1, 0)) %>%
  mutate(time_to_n_2 = ifelse(is.na(time_to_n), max(time_to_n, na.rm = TRUE), time_to_n)) %>%
  group_by(input_R0, Tg_Ratio, quarantine_efficacy) %>%
  mutate(time_to_n_relative = time_to_n / time_to_n[vaccine_efficacy_infection == 0],
         time_to_n_2_relative = time_to_n_2 / time_to_n_2[vaccine_efficacy_infection == 0]) %>%
  ungroup() %>%
  group_by(input_R0, Tg_Ratio, vaccine_efficacy_infection, quarantine_efficacy) %>%
  summarise(proportion_contained = sum(contained) / iterations,
            avg_time_to_n = mean(time_to_n, na.rm = TRUE),
            avg_time_to_n_2 = mean(time_to_n_2, na.rm = TRUE),
            avg_time_to_n_relative = mean(time_to_n_relative, na.rm = TRUE),
            avg_time_to_n_2_relative = mean(time_to_n_2_relative, na.rm = TRUE),
            avg_R0 = mean(R0),
            avg_Reff = mean(Reff)) %>%
  mutate(avg_time_to_n_2_relative_plot = ifelse(proportion_contained < 0.9, avg_time_to_n_2_relative, NA))

### Vaccine Efficacy
R0_efficacy_df <- storage_R0_efficacy_df %>%
  mutate(contained = ifelse(epidemic_size < (0.9 * check_final_size), 1, 0)) %>%
  mutate(time_to_n_2 = ifelse(is.na(time_to_n), max(time_to_n, na.rm = TRUE), time_to_n)) %>%
  group_by(input_R0, quarantine_efficacy) %>%
  mutate(time_to_n_relative = time_to_n / time_to_n[vaccine_efficacy_infection == 0],
         time_to_n_2_relative = time_to_n_2 / time_to_n_2[vaccine_efficacy_infection == 0]) %>%
  ungroup() %>%
  group_by(input_R0, vaccine_efficacy_infection, quarantine_efficacy) %>%
  summarise(proportion_contained = sum(contained) / iterations,
            avg_time_to_n = mean(time_to_n, na.rm = TRUE),
            avg_time_to_n_2 = mean(time_to_n_2, na.rm = TRUE),
            avg_time_to_n_relative = mean(time_to_n_relative, na.rm = TRUE),
            avg_time_to_n_2_relative = mean(time_to_n_2_relative, na.rm = TRUE),
            avg_R0 = mean(R0),
            avg_Reff = mean(Reff)) %>%
  mutate(avg_time_to_n_2_relative_plot = ifelse(proportion_contained < 0.9, avg_time_to_n_2_relative, NA))

### Proportion of Presymptomatic Transmission
R0_preSymp_df <- storage_R0_preSymp_df %>%
  mutate(contained = ifelse(epidemic_size < (0.9 * check_final_size), 1, 0)) %>%
  mutate(time_to_n_2 = ifelse(is.na(time_to_n), max(time_to_n, na.rm = TRUE), time_to_n)) %>%
  group_by(input_R0, prop_preSymp, quarantine_efficacy) %>%
  mutate(time_to_n_relative = time_to_n / time_to_n[vaccine_efficacy_infection == 0],
         time_to_n_2_relative = time_to_n_2 / time_to_n_2[vaccine_efficacy_infection == 0]) %>%
  ungroup() %>%
  group_by(input_R0, prop_preSymp, vaccine_efficacy_infection, quarantine_efficacy) %>%
  summarise(proportion_contained = sum(contained) / iterations,
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
### R0 / Tg Ratio Parameter Scans
############################################################
main_contained_R0_TgRatio_plot <- ggplot(subset(R0_TgRatio_df, quarantine_efficacy != 0.35 & vaccine_efficacy_infection == 0.35), aes(x = input_R0, y = Tg_Ratio, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0, 1), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  labs(x = "R0", y = "Ratio of Tg to Vaccine Protection Delay") +
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

SI_contained_R0_TgRatio_plot <- ggplot(subset(R0_TgRatio_df, quarantine_efficacy != 0.35), aes(x = input_R0, y = Tg_Ratio, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0, 1), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  labs(x = "R0", y = "Ratio of Tg to Vaccine Protection Delay") +
  facet_grid(vaccine_efficacy_infection ~ quarantine_efficacy,
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
ggsave(file = "figures/Figure_1_BranchingProcess/FigS1_R0TgRatio_ParameterScan_ContainedPlot.pdf", plot = SI_contained_R0_TgRatio_plot, width = 9, height = 7.2)

SI_Reff_R0_TgRatio_plot <- ggplot(subset(R0_TgRatio_df, quarantine_efficacy != 0.35), aes(x = input_R0, y = Tg_Ratio, fill = 100 * (1 - (avg_Reff / avg_R0)))) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0, 70), begin = 0.175, end = 1,
                       name = "% Red.\nin R0",
                       direction = 1) +
  labs(x = "R0",
       y = "Ratio of Tg to Vaccine Protection Delay") +
  facet_grid(vaccine_efficacy_infection ~ quarantine_efficacy,
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
ggsave(file = "figures/Figure_1_BranchingProcess/FigS1_R0TgRatio_ParameterScan_ReffPlot.pdf", plot = SI_Reff_R0_TgRatio_plot, width = 9, height = 7.2)

SI_timetoN_R0_TgRatio_plot <- ggplot(subset(R0_TgRatio_df, quarantine_efficacy != 0.35 & vaccine_efficacy_infection != 0), 
       aes(x = input_R0, y = Tg_Ratio, fill = avg_time_to_n_relative, alpha = 100 * (1 - proportion_contained))) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0.9, 1.5), begin = 0.175, end = 1,
                       breaks = c(1, 1.1, 1.2, 1.3, 1.4, 1.5),
                       name = "Fold Increase\nin Time to Epidemic\nThreshold",
                       direction = 1) +
  scale_alpha(name = "% Outbreaks\nNot Contained") +
  labs(x = "R0", y = "Ratio of Tg to Vaccine Protection Delay") +
  facet_grid(vaccine_efficacy_infection ~ quarantine_efficacy,
             labeller = labeller(vaccine_efficacy_infection = c(`0.35` = "Vaccine Efficacy = 35%", 
                                                                `0.75` = "Vaccine Efficacy = 75%"),
                                 quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)
ggsave(file = "figures/Figure_1_BranchingProcess/FigS1_R0TgRatio_ParameterScan_timetoNPlot.pdf", plot = SI_timetoN_R0_TgRatio_plot, width = 9, height = 4.8)


############################################################
### R0 / Vaccine Efficacy Parameter Scans
############################################################
main_contained_R0_efficacy_plot <- ggplot(subset(R0_efficacy_df, quarantine_efficacy != 0.35 & vaccine_efficacy_infection != 0.0), aes(x = input_R0, y = 100 * vaccine_efficacy_infection, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0, 1), begin = 0.175, end = 1, name = "Proportion\nContained",
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
  scale_fill_viridis_c(option = "mako", limits = c(20, 70), begin = 0.175, end = 1,
                       name = "% Red.\nin R0",
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
ggsave(file = "figures/Figure_1_BranchingProcess/FigS1_R0VaxEff_ParameterScan_ReffPlot.pdf", plot = SI_Reff_R0_efficacy_plot, width = 6, height = 2.4)

SI_timetoN_R0_efficacy_plot <- ggplot(subset(R0_efficacy_df, quarantine_efficacy != 0.35 & vaccine_efficacy_infection != 0), 
                                     aes(x = input_R0, y = 100 * vaccine_efficacy_infection, fill = avg_time_to_n_relative, alpha = 100 * (1 - proportion_contained))) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0.9, 1.5), begin = 0.175, end = 1,
                       breaks = c(1, 1.1, 1.2, 1.3, 1.4, 1.5),
                       name = "Fold Increase\nin Time to Epidemic\nThreshold",
                       direction = 1) +
  scale_alpha(name = "% Outbreaks\nNot Contained") +
  labs(x = "R0", y = "Vaccine Efficacy (%)") +
  facet_grid(. ~ quarantine_efficacy,
             labeller = labeller(quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)
ggsave(file = "figures/Figure_1_BranchingProcess/FigS1_R0TgRatio_ParameterScan_timetoNPlot.pdf", plot = SI_timetoN_R0_efficacy_plot, width = 6, height = 2.4)

#################################################################
### R0 / Proportion Presymptomatic Transmission Parameter Scans
#################################################################
main_contained_R0_preSymp_plot <- ggplot(subset(R0_preSymp_df, quarantine_efficacy != 0.35 & vaccine_efficacy_infection == 0.35), aes(x = input_R0, y = 100 * prop_preSymp, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0, 1), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  labs(x = "R0", y = "% Presymptomatic\nTransmission") +
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
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

SI_contained_R0_preSymp_plot <- ggplot(subset(R0_preSymp_df, quarantine_efficacy != 0.35), aes(x = input_R0, y = 100 * prop_preSymp, fill = proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0, 1), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  labs(x = "R0", y = "% Presymptomatic\nTransmission") +
  facet_grid(vaccine_efficacy_infection ~ quarantine_efficacy,
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
ggsave(file = "figures/Figure_1_BranchingProcess/FigS1_R0preSymp_ParameterScan_ContainedPlot.pdf", plot = SI_contained_R0_preSymp_plot, width = 9, height = 7.2)

SI_Reff_R0_preSymp_plot <- ggplot(subset(R0_preSymp_df, quarantine_efficacy != 0.35), aes(x = input_R0, y = 100 * prop_preSymp, fill = 100 * (1 - (avg_Reff / avg_R0)))) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0, 70), begin = 0.175, end = 1,
                       name = "% Red.\nin R0",
                       direction = 1) +
  labs(x = "R0", y = "% Presymptomatic\nTransmission") +
  facet_grid(vaccine_efficacy_infection ~ quarantine_efficacy,
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
ggsave(file = "figures/Figure_1_BranchingProcess/FigS1_R0preSymp_ParameterScan_ReffPlot.pdf", plot = SI_Reff_R0_preSymp_plot, width = 9, height = 7.2)

SI_timetoN_R0_preSymp_plot <- ggplot(subset(R0_preSymp_df, quarantine_efficacy != 0.35 & vaccine_efficacy_infection != 0), 
                                     aes(x = input_R0, y = 100 * prop_preSymp, fill = avg_time_to_n_relative, alpha = 100 * (1 - proportion_contained))) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "mako", limits = c(0.9, 1.5), begin = 0.175, end = 1,
                       breaks = c(1, 1.1, 1.2, 1.3, 1.4, 1.5),
                       name = "Fold Increase\nin Time to Epidemic\nThreshold",
                       direction = 1) +
  scale_alpha(name = "% Outbreaks\nNot Contained") +
  labs(x = "R0", y = "% Presymptomatic\nTransmission") +
  facet_grid(vaccine_efficacy_infection ~ quarantine_efficacy,
             labeller = labeller(vaccine_efficacy_infection = c(`0.35` = "Vaccine Efficacy = 35%", 
                                                                `0.75` = "Vaccine Efficacy = 75%"),
                                 quarantine_efficacy = c(`0`   = "No Quarantine", 
                                                         `0.35` = "Quarantine Efficacy = 35%", 
                                                         `0.65`   = "Quarantine Efficacy = 65%"))) +  
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)
ggsave(file = "figures/Figure_1_BranchingProcess/FigS1_R0preSymp_ParameterScan_timetoNPlot.pdf", plot = SI_timetoN_R0_preSymp_plot, width = 9, height = 4.8)

## Overall figure
Fig1FGH <- cowplot::plot_grid(main_contained_R0_TgRatio_plot + theme(legend.position = "none"),
                              main_contained_R0_efficacy_plot + theme(legend.position = "none"), 
                              main_contained_R0_preSymp_plot + theme(legend.position = "none"),
                              nrow = 1,
                              labels = c("F", "G", "H"), rel_widths = c(1, 1, 1.12))
overall_figure1 <- cowplot::plot_grid(Fig1BCDE, Fig1FGH, nrow = 2, rel_heights = c(1.25, 1))
ggsave(file = "figures/Figure_1_BranchingProcess/Fig1_Overall.pdf", plot = overall_figure1, width = 8, height = 9.5)


