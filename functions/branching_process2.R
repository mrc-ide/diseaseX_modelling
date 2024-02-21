chain_sim_susc_ring_vacc2 <- function(offspring = c("pois", "nbinom"), mn_offspring, disp_offspring, 
                                     generation_time, t0 = 0, tf = Inf, pop, check_final_size, initial_immune = 0,
                                     seeding_cases, prop_asymptomatic, infection_to_onset,                                 
                                     vaccine_start, vaccine_coverage, vaccine_efficacy_infection,                         
                                     vaccine_efficacy_transmission, vaccine_logistical_delay, vaccine_protection_delay) {
  
  offspring <- match.arg(offspring)
  if (offspring == "pois") {
    if (!missing(disp_offspring)) {
      warning("argument disp_offspring not used for\n poisson offspring distribution.")
    }
    ### NOTE THAT THIS IS MISSING THE SUSCEPTIBLE DEPLETION CURRENTLY. DID THIS BECAUSE SOME WEIRDNESS AROUND R0 = 1 GIVING ONGOING EPIDEMICS.
    # offspring_fun <- function(n, susc) {
    #   truncdist::rtrunc(n, spec = "pois", lambda = mn_offspring * susc/pop, b = susc)
    # }
    offspring_fun <- function(n, susc) {
      rpois(n, lambda = mn_offspring)
    }
  }
  else if (offspring == "nbinom") {
    if (disp_offspring <= 1) {
      stop("Offspring distribution 'nbinom' requires argument\n disp_offspring > 1. Use 'pois' if there is no overdispersion.")
    }
    offspring_fun <- function(n, susc) {
      new_mn <- mn_offspring * susc/pop
      size <- new_mn/(disp_offspring - 1)
      truncdist::rtrunc(n, spec = "nbinom", b = susc, mu = new_mn, size = size)
    }
  }
  
  # Pre-allocate a dataframe with the maximum size
  max_cases <- check_final_size
  tdf <- data.frame(
    id = integer(max_cases),
    ancestor = integer(max_cases),
    generation = integer(max_cases),
    time_infection = NA_real_,
    time_onset = numeric(max_cases),
    vaccinated = integer(max_cases),
    time_vaccinated = numeric(max_cases),
    vaccinated_before_infection = integer(max_cases),
    vaccinated_after_infection = integer(max_cases),
    time_protected = numeric(max_cases),
    protected_before_infection = integer(max_cases),
    protected_after_infection = integer(max_cases),
    asymptomatic = integer(max_cases),
    n_offspring = integer(max_cases),
    n_offspring_new = integer(max_cases),
    offspring_generated = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Initialize the dataframe with the seeding cases
  tdf[1:seeding_cases, ] <- data.frame(
    id = seq_len(seeding_cases),
    ancestor = NA_integer_,
    generation = 1L,
    time_infection = t0 + seq(from = 0, to = 0.01, length.out = seeding_cases),
    time_onset = NA,
    vaccinated = 0,
    time_vaccinated = NA,
    vaccinated_before_infection = NA,
    vaccinated_after_infection = NA,
    time_protected = NA,
    protected_before_infection = NA,
    protected_after_infection = NA,
    asymptomatic = integer(seeding_cases),
    n_offspring = NA_integer_,
    n_offspring_new = NA_integer_,
    offspring_generated = FALSE
  )
  
  susc <- pop - initial_immune - 1L
  time_infection_index <- t0

  # while (any(tdf$time_infection[!tdf$offspring_generated & !is.na(tdf$time_infection)] <= tf) & susc > 0) {
  ### NOTE: check_final_size doesn't really work when you're trying to get incidence over time
  ###       because at the check_final_size mark you might have lots of infections you haven't generated offspring 
  ###       for, which in turn will go on to generate offspring that will contribute to the incidence within the
  ###       timeframe you're simulating. It's incomplete in that regard
  while (any(tdf$time_infection[!tdf$offspring_generated & !is.na(tdf$time_infection)] <= tf) & susc > 0 & nrow(tdf) <= check_final_size) {
    
    time_infection_index <- min(tdf$time_infection[tdf$offspring_generated == 0 & !is.na(tdf$time_infection)])              # Note: Is not an issue in practice, but I don't think this is currently set up to handle >= 2 infections with same infection time currently
    idx <- which(tdf$time_infection == time_infection_index & !tdf$offspring_generated)[1] # get the id of the earliest unsimulated infection
    id_parent <- tdf$id[idx]                                                               # parent of the earliest unsimulated infection
    t_parent <- tdf$time_infection[idx]                                                    # infection time of the earliest unsimulated infection
    gen_parent <- tdf$generation[idx]                                                      # generation of the earliest unsimulated infection
    current_max_id <- max(tdf$id)                                                          # total number of infections in the dataframe currently (so we can figure out how to label the new infections)
    index_vaccinated <- tdf$vaccinated[idx]                                                # whether or not the index case (the "parent") is vaccinated
    time_vaccinated <- tdf$time_vaccinated[idx]                                            # when the index case (the "parent") was vaccinated
    time_protected <- tdf$time_protected[idx]                                              # when the index case (the "parent") was protected
    onset_time_index_case <- infection_to_onset(n = 1)                                     # generate the time from infection to symptom onset for the index case
    tdf$time_onset[idx] <- onset_time_index_case                                           # --
    index_asymptomatic <- tdf$asymptomatic[idx]                                            # whether or not the index case (the "parent") is asymptomatic (influences whether contacts get ring vaccinated or not)
    
    n_offspring <- offspring_fun(1, susc) 
    tdf$n_offspring[idx] <- n_offspring
    
    if (index_vaccinated == 1) {
      if (time_protected < time_infection_index) {
        n_offspring <- sum(rbinom(n = n_offspring, size = 1, prob = 1 - vaccine_efficacy_transmission))
      }
    }
    if (n_offspring %% 1 > 0) { # Checking offspring function is correctly returning integers
      stop("Offspring distribution must return integers")
    }
    tdf$n_offspring_new[idx] <- n_offspring
    tdf$offspring_generated[idx] <- TRUE
    
    if (n_offspring > 0) {
      
      new_times <- generation_time(n_offspring)
      if (any(new_times < 0)) {
        stop("Generation times must be >= 0.")
      }
      if (time_infection_index < vaccine_start) {
        asymptomatic <- rbinom(n = n_offspring, size = 1, prob = prop_asymptomatic)
        tdf[(current_max_id+1):(current_max_id+n_offspring), "id"] <- c(current_max_id + seq_len(n_offspring))
        tdf[(current_max_id+1):(current_max_id+n_offspring), "ancestor"] <- id_parent
        tdf[(current_max_id+1):(current_max_id+n_offspring), "generation"] <- gen_parent + 1L
        tdf[(current_max_id+1):(current_max_id+n_offspring), "time_infection"] <- new_times + t_parent
        tdf[(current_max_id+1):(current_max_id+n_offspring), "time_onset"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "vaccinated"] <- 0
        tdf[(current_max_id+1):(current_max_id+n_offspring), "time_vaccinated"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "vaccinated_before_infection"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "vaccinated_after_infection"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "time_protected"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "protected_before_infection"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "protected_after_infection"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "asymptomatic"] <- asymptomatic
        tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring_new"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "offspring_generated"] <- FALSE
      } else {
        if (index_asymptomatic == 1) {
          asymptomatic <- rbinom(n = n_offspring, size = 1, prob = prop_asymptomatic)
          tdf[(current_max_id+1):(current_max_id+n_offspring), "id"] <- c(current_max_id + seq_len(n_offspring))
          tdf[(current_max_id+1):(current_max_id+n_offspring), "ancestor"] <- id_parent
          tdf[(current_max_id+1):(current_max_id+n_offspring), "generation"] <- gen_parent + 1L
          tdf[(current_max_id+1):(current_max_id+n_offspring), "time_infection"] <- new_times + t_parent
          tdf[(current_max_id+1):(current_max_id+n_offspring), "time_onset"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "vaccinated"] <- 0
          tdf[(current_max_id+1):(current_max_id+n_offspring), "time_vaccinated"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "vaccinated_before_infection"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "vaccinated_after_infection"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "time_protected"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "protected_before_infection"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "protected_after_infection"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "asymptomatic"] <- asymptomatic
          tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring_new"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "offspring_generated"] <- FALSE
        } else {
          
          are_they_vaccinated <- vector(mode = "integer", length = n_offspring) ## vector of whether secondary infections get vaccinated
          did_they_have_potential_protection <- vector(mode = "integer", length = n_offspring)
          infection_retained <- vector(mode = "integer", length = n_offspring)  ## vector of whether secondary infections get retained (i.e. not prevented by ring vaccination)
          infection_retained[1:length(infection_retained)] <- 1                 ## default to infections being retained; and then flow through below to see if they get removed
          time_to_secondary_vaccination <- onset_time_index_case + vaccine_logistical_delay                    ## Time between index case infected and secondary cases ring vaccinated
          time_to_secondary_vaccination_protection <- time_to_secondary_vaccination + vaccine_protection_delay ## Time between index case infected and secondary cases protected with vaccination
          for (i in 1:n_offspring) {
            if (time_to_secondary_vaccination <= new_times[i]) { # if infection occurs AFTER (potential) vaccination
              are_they_vaccinated[i] <- rbinom(n = 1, size = 1, prob = vaccine_coverage)
              did_they_have_potential_protection[i] <- ifelse(time_to_secondary_vaccination_protection <= new_times[i], 1, 0)
              were_they_protected <- rbinom(n = 1, size = 1, prob = vaccine_efficacy_infection * did_they_have_potential_protection[i])
              infection_retained[i] <- ifelse(are_they_vaccinated[i] == 1 & were_they_protected == 1, 0, 1)
            } else { # if infection occurs BEFORE vaccination
              are_they_vaccinated[i] <- 0 ## eliding together "unvaccinated" and "vaccinated after infection occurs"
              infection_retained[i] <- 1 ## note that implicitly here we're implicitly "saying" these folks aren't vaccinated.
            }
          }
          new_n_offspring <- sum(infection_retained)
          new_new_times <- new_times[which(infection_retained == 1)]
          vaccinated <- ifelse(are_they_vaccinated[which(infection_retained == 1)] == 0, 0, 1)                                            # of the retained infections, which are vaccinated
          time_vaccinated <- ifelse(are_they_vaccinated == 1, time_to_secondary_vaccination, NA)[which(infection_retained == 1)]           # of the retained infections, when are they vaccinated (relative to infection time of index)
          vaccinated_before_infection <- ifelse(are_they_vaccinated[which(infection_retained == 1)] == 0, NA, 1)                           # (currently we combine all individuals not vaccinated and vaccinated after infection into "unvaccinated", so all vaccinated individuals necessarily got vaccinated before infection)
          time_protected <- ifelse(are_they_vaccinated == 1, time_to_secondary_vaccination_protection, NA)[which(infection_retained == 1)] # of the retained infections, when are they protected (relative to infection time of index)
          protected_before_infection <- ifelse(is.na(time_protected), NA, ifelse(time_protected <= new_new_times, 1, 0))
          protected_after_infection <- ifelse(is.na(time_protected), NA, ifelse(time_protected > new_new_times, 1, 0))
          asymptomatic <- rbinom(n = new_n_offspring, size = 1, prob = prop_asymptomatic)
          if (new_n_offspring != 0) {
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "id"] <- c(current_max_id + seq_len(new_n_offspring))
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "ancestor"] <- id_parent
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "generation"] <- gen_parent + 1L
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "time_infection"] <- new_new_times + t_parent
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "time_onset"] <- NA
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "vaccinated"] <- vaccinated
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "time_vaccinated"] <- time_vaccinated + t_parent
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "vaccinated_before_infection"] <- vaccinated_before_infection
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "vaccinated_after_infection"] <- ifelse(is.na(vaccinated_before_infection), NA, 0)
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "time_protected"] <- time_protected + t_parent
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "protected_before_infection"] <- protected_before_infection
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "protected_after_infection"] <- protected_after_infection
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "asymptomatic"] <- asymptomatic
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring"] <- NA
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring_new"] <- NA
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "offspring_generated"] <- FALSE
          }
        }
      }
    }
    susc <- susc - n_offspring
  }
  tdf <- tdf[tdf$time_infection <= tf, ]
  tdf <- tdf[order(tdf$time_infection, tdf$id), ]
  tdf$offspring_generated <- NULL
  tdf$abs_time_onset <- tdf$time_onset + tdf$time_infection
  return(tdf)
}

## Testing out the branching process
## Loading required libraries
# library(tictoc); library(profvis)
# 
# ## Checking both versions produce the same results
# generation_time <- function(n) { rgamma(n, shape = 12, rate = 2) }
# infection_to_onset <- function(n) { rgamma(n, shape = 3, rate = 2) }
# iterations <- 25
# R0_scan <- c(0.75, 1, 1.25, 1.5, 1.75, 2, 3)
# storage_vacc <- matrix(nrow = iterations, ncol = length(R0_scan))
# storage_vacc2 <- matrix(nrow = iterations, ncol = length(R0_scan))
# for (i in 1:length(R0_scan)) {
#   for (j in 1:iterations) {
# 
#     test_vacc <- chain_sim_susc_ring_vacc(offspring = "pois",
#                                           mn_offspring = R0_scan[i],
#                                           generation_time = generation_time,
#                                           t0 = 0, tf = Inf, pop = 10^8, check_final_size = 1000, initial_immune = 0,
#                                           seeding_cases = 5, prop_asymptomatic = 0,
#                                           infection_to_onset = infection_to_onset,
#                                           vaccine_start = 5, vaccine_coverage = 0.8,
#                                           vaccine_efficacy_infection = 0.75,
#                                           vaccine_efficacy_transmission = 0.75,
#                                           vaccine_logistical_delay = 2,
#                                           vaccine_protection_delay = 2)
#     storage_vacc[j, i] <- nrow(test_vacc)
#     
#     test_vacc2 <- chain_sim_susc_ring_vacc2(offspring = "pois",
#                                             mn_offspring = R0_scan[i],
#                                             generation_time = generation_time,
#                                             t0 = 0, tf = Inf, pop = 10^8, check_final_size = 1000, initial_immune = 0,
#                                             seeding_cases = 5, prop_asymptomatic = 0,
#                                             infection_to_onset = infection_to_onset,
#                                             vaccine_start = 5, vaccine_coverage = 0.8,
#                                             vaccine_efficacy_infection = 0.75,
#                                             vaccine_efficacy_transmission = 0.75,
#                                             vaccine_logistical_delay = 2,
#                                             vaccine_protection_delay = 2)
#     storage_vacc2[j, i] <- sum(!is.na(test_vacc2$time_infection))
# 
#   }
#   print(i)
# }
# plot(R0_scan, apply(storage_vacc, 2, median), type = "l", xlab = "R0", ylab = "Final Epidemic Size (Capped at 1,000)")
# lines(R0_scan, apply(storage_vacc2, 2, median), col = "blue")
# 
# 
# ## Checking out relative speeds
# tic()
# storage_vacc <- matrix(nrow = iterations, ncol = length(R0_scan))
# for (i in 1:length(R0_scan)) {
#   for (j in 1:iterations) {
#     test_vacc <- chain_sim_susc_ring_vacc(offspring = "pois",
#                                           mn_offspring = R0_scan[i],
#                                           generation_time = generation_time,
#                                           t0 = 0, tf = Inf, pop = 10^8, check_final_size = 1000, initial_immune = 0,
#                                           seeding_cases = 5, prop_asymptomatic = 0,
#                                           infection_to_onset = infection_to_onset,
#                                           vaccine_start = 5, vaccine_coverage = 0.8,
#                                           vaccine_efficacy_infection = 0.75,
#                                           vaccine_efficacy_transmission = 0.75,
#                                           vaccine_logistical_delay = 2,
#                                           vaccine_protection_delay = 2)
#     storage_vacc[j, i] <- nrow(test_vacc)
#   }
#   print(i)
# }
# toc()
# 
# tic()
# storage_vacc2 <- matrix(nrow = iterations, ncol = length(R0_scan))
# for (i in 1:length(R0_scan)) {
#   for (j in 1:iterations) {
#     test_vacc2 <- chain_sim_susc_ring_vacc2(offspring = "pois",
#                                             mn_offspring = R0_scan[i],
#                                             generation_time = generation_time,
#                                             t0 = 0, tf = Inf, pop = 10^8, check_final_size = 1000, initial_immune = 0,
#                                             seeding_cases = 5, prop_asymptomatic = 0,
#                                             infection_to_onset = infection_to_onset,
#                                             vaccine_start = 5, vaccine_coverage = 0.8,
#                                             vaccine_efficacy_infection = 0.75,
#                                             vaccine_efficacy_transmission = 0.75,
#                                             vaccine_logistical_delay = 2,
#                                             vaccine_protection_delay = 2)
#     storage_vacc2[j, i] <- sum(!is.na(test_vacc2$time_infection))
#   }
#   print(i)
# }
# toc()

# tic()
# test <- chain_sim_susc_ring_vacc(offspring = "pois",
#                                  mn_offspring = 3.5,
#                                  generation_time = generation_time,
#                                  t0 = 0, tf = Inf, pop = 10^8, check_final_size = 1500, initial_immune = 0,
#                                  seeding_cases = 5, prop_asymptomatic = 0,
#                                  infection_to_onset = infection_to_onset,
#                                  vaccine_start = 5, vaccine_coverage = 1,
#                                  vaccine_efficacy_infection = 0.85,
#                                  vaccine_efficacy_transmission = 0.85,
#                                  vaccine_logistical_delay = 3,
#                                  vaccine_protection_delay = 1)
# toc()
# # 
# tic()
# profvis({test <- chain_sim_susc_ring_vacc2(offspring = "pois",
#                                            mn_offspring = 3.5,
#                                            generation_time = generation_time,
#                                            t0 = 0, tf = Inf, pop = 10^8, check_final_size = 1500, initial_immune = 0,
#                                            seeding_cases = 5, prop_asymptomatic = 0,
#                                            infection_to_onset = infection_to_onset,
#                                            vaccine_start = 5, vaccine_coverage = 1,
#                                            vaccine_efficacy_infection = 0.85,
#                                            vaccine_efficacy_transmission = 0.85,
#                                            vaccine_logistical_delay = 3,
#                                            vaccine_protection_delay = 1)})
# toc()






