ring_vax_bp_sim <- function(## Transmission Parameters
                            offspring = c("pois", "nbinom"),   # offspring distribution 
                            mn_offspring,                      # mean of the offspring distribution
                            disp_offspring,                    # overdispersion of the offspring distribution (if negative binomial)
                            
                            ## Natural History Parameters
                            generation_time,                   # generation time distribution
                            prop_asymptomatic,                 # probability of being asymptomatic
                            infection_to_onset,                # time from infection to symptom onset distribution
                             
                            ## Vaccine-Related Parameters
                            vaccine_start,                     # time at which the vaccine becomes available
                            vaccine_coverage,                  # probability that each eligible individual gets vaccinated
                            vaccine_efficacy_infection,        # vaccine efficacy against infection
                            vaccine_efficacy_transmission,     # reduction in transmissibility of breakthrough infections in vaccinated individuals
                            vaccine_logistical_delay,          # delay between symptom onset and vaccination of contacts (and contacts of contacts) occurring
                            vaccine_protection_delay,          # delay between vaccination and protection developing
                            
                            ## Quarantine Related Parameters
                            time_to_quarantine,                # time delay between trigger and quarantining (either symptoms of the primary infection, or symptoms in the secondary infection)
                            prob_quarantine,                   # probability that the individual successfully isolates
                            quarantine_efficacy,               # effectiveness of the quarantine at reducing onwards transmission

                            ## Miscellaneous Parameters
                            t0 = 0,                            # simulation starting time
                            tf = Inf,                          # final simulation time
                            population,                        # total population size
                            check_final_size,                  # maximum number of infections to simulate
                            initial_immune = 0,                # initial number of individuals who are immune
                            seeding_cases,                     # number of seeding cases to start the epidemic with
                            seed                               # stochastic seed
                            ) {
  
  ## Setting the seed
  set.seed(seed)
  
  ## Setting up the number of susceptibles
  susc <- population - initial_immune
  
  ## Setting up the offspring distribution
  offspring <- match.arg(offspring)
  if (offspring == "pois") {
    offspring_fun <- function(n, susc) {
      rpois(n, lambda = mn_offspring * susc/population)
    }
  } else if (offspring == "nbinom") {
    if (disp_offspring <= 1) {
      stop("Offspring distribution 'nbinom' requires argument\n disp_offspring > 1. Use 'pois' if there is no overdispersion.")
    }
    offspring_fun <- function(n, susc) {
      new_mn <- mn_offspring * susc/population
      size <- new_mn/(disp_offspring - 1)
      truncdist::rtrunc(n, spec = "nbinom", b = susc, mu = new_mn, size = size)
    }
  } else {
    stop("offspring specification is wrong")
  }
  
  ## Checking the quarantine trigger
  # Checking likelihood distributions are valid
  quarantine_trigger <- c("contact_tracing", "symptom_onset")
  if (!all(misc$likelihood_distribution %in% valid_distributions)) {
    stop("misc$likelihood_distribution can only contain 'poisson' or 'negative_binomial'.")
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
    quarantined = integer(check_final_size),                   
    time_quarantined_absolute = NA_real_,
    n_offspring = integer(max_cases),
    n_offspring_new = integer(max_cases),
    n_offspring_quarantine = integer(max_cases),
    n_offspring_post_pruning = integer(max_cases),
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
    quarantined = integer(seeding_cases),
    time_quarantined_absolute = NA,
    n_offspring = NA_integer_,
    n_offspring_new = NA_integer_,
    n_offspring_quarantine = NA_integer_,
    n_offspring_post_pruning = NA_integer_,
    offspring_generated = FALSE
  )
  
  time_infection_index <- t0

  ## While we haven't hit the simulation cap size (check_final_size) and any infections exist where we have not yet generated the requisite offspring, 
  ## continue to generate infections
  while ((any(is.na(tdf$n_offspring)) & nrow(tdf) <= check_final_size & susc > 0)) {
    
    ## Getting the timings of the earliest/oldest infection we haven't yet generated infections for - this is the "INDEX INFECTION"
    time_infection_index <- min(tdf$time_infection[tdf$offspring_generated == 0 & !is.na(tdf$time_infection)]) 
    idx <- which(tdf$time_infection == time_infection_index & !tdf$offspring_generated)[1] # get the id of the earliest unsimulated infection
    previous_parent_id <- tdf$ancestor[idx]
    previous_parent_idx <- which(tdf$id == previous_parent_id)
    id_parent <- tdf$id[idx]                                                               # id of the earliest unsimulated infection
    t_parent <- tdf$time_infection[idx]                                                    # infection time of the earliest unsimulated infection
    gen_parent <- tdf$generation[idx]                                                      # generation of the earliest unsimulated infection
    current_max_id <- max(tdf$id)                                                          # total number of infections in the dataframe currently (so we can figure out how to label the new infections)
    index_vaccinated <- tdf$vaccinated[idx]                                                # whether or not the index case (the "parent") is vaccinated
    time_vaccinated <- tdf$time_vaccinated[idx]                                            # when the index case (the "parent") was vaccinated
    time_protected <- tdf$time_protected[idx]                                              # when the index case (the "parent") was protected
    onset_time_index_case <- infection_to_onset(n = 1)                                     # generate the time from infection to symptom onset for the index case
    tdf$time_onset[idx] <- onset_time_index_case                                           # --
    index_asymptomatic <- tdf$asymptomatic[idx]                                            # whether or not the index case (the "parent") is asymptomatic (influences whether contacts get ring vaccinated or not)
    ### CHECK THAT T_PARENT AND time_infection_index ARE THE SAME - THEY SHOULD BE I THINK???
    
    ## Calculating whether or not the individual isolates/quarantines
    ### note that here I swap from using "parent" to refer to "index" and instead as the prior infector. Need to sort this at some point.
    ### (it's just syntax, the code is actually right/doing the right thing)
    
    ## Extracting parent information
    parent_infection_time <- tdf$time_infection[previous_parent_idx];
    parent_asymptomatic <- tdf$asymptomatic[previous_parent_idx];
    parent_onset_time <- ifelse(parent_asymptomatic == 0, tdf$time_onset[previous_parent_idx], NA)
    
    ## If the parent is symptomatic, this triggers quarantine in secondary infections relative to timing of symptoms in parent
    if (parent_asymptomatic == 0) {
      
      index_quarantine <- rbinom(n = 1, size = 1, prob = prob_quarantine)                                     # whether or not the index infection isolates
      index_quarantine_time <- ifelse(index_quarantine == 1, onset_to_quarantine(n = 1), NA)                  # if the infection isolates, how soon after symptom onset they do so
      tdf$quarantined[idx] <- index_quarantine                                                                # adding quarantine indicator to storage dataframe
      absolute_quarantine_time <- parent_infection_time + parent_onset_time + index_quarantine_time           # adding quarantine time in absolute calendar time to the storage dataframe
      tdf$time_quarantined_absolute[idx] <- absolute_quarantine_time

    ## If the parent is asymptomatic but index is symptomatic, this triggers quarantine relative to timing in index
    } else if (parent_asymptomatic == 1 & index_asymptomatic == 0) {
      
      index_quarantine <- rbinom(n = 1, size = 1, prob = prob_quarantine)                                        # whether or not the index infection isolates
      index_quarantine_time <- ifelse(index_quarantine == 1, onset_to_quarantine(n = 1), NA)                     # if the infection isolates, how soon after symptom onset they do so
      tdf$quarantined[idx] <- index_quarantine                                                                   # adding quarantine indicator to storage dataframe
      absolute_quarantine_time <- time_infection_index + onset_time_index_case + index_quarantine_time           # adding quarantine time in absolute calendar time to the storage dataframe
      tdf$time_quarantined_absolute[idx] <- absolute_quarantine_time
    
    ## If both the parent and index are asymptomatic, no quarantining can possibly happen
    } else {
      index_quarantine <- 0
      index_quarantine_time <- NA
      absolute_quarantine_time <- NA
    }
    
    ## Calculating time to vaccinate the contacts of this index infection
    time_to_secondary_vaccination <- onset_time_index_case + vaccine_logistical_delay                    ## Time between index case infected and secondary cases ring vaccinated
    time_to_secondary_vaccination_protection <- time_to_secondary_vaccination + vaccine_protection_delay ## Time between index case infected and secondary cases protected with vaccination

    # Generating offspring for this infection
    n_offspring <- offspring_fun(1, susc) 
    tdf$n_offspring[idx] <- n_offspring
    tdf$offspring_generated[idx] <- TRUE
    
    ## If infection was previously vaccinated and is a breakthrough infection, account for reduced transmissibility
    if (index_vaccinated == 1 & (time_protected < time_infection_index) & n_offspring != 0) {
      n_offspring <- sum(rbinom(n = n_offspring, size = 1, prob = 1 - vaccine_efficacy_transmission))
    }
    tdf$n_offspring_new[idx] <- n_offspring
    new_times <- generation_time(n_offspring)
    absolute_new_times <- time_infection_index + new_times
    
    ## If index infection quarantines, reduce secondary infections - note that quarantine only occurs if infection has symptoms. Only do this if there are offspring to avert.
    if (index_quarantine == 1 & n_offspring != 0) {
      
      # Implement quarantining for index infection
      index_n_offspring <- implement_quarantine(symptom_onset_time = 0, # have wrapped symptom onset time into absolute_quarantine time,
                                                quarantine_time = absolute_quarantine_time,
                                                n_offspring = n_offspring,
                                                offspring_infection_times = absolute_new_times,
                                                quarantine_efficacy = quarantine_efficacy)
      
      # Updating number offspring, their infection times and characteristics to reflect removals due to quarantining
      n_offspring <- index_n_offspring$updated_n_offspring
      new_times <- index_n_offspring$updated_infection_times

    }
    tdf$n_offspring_quarantine[idx] <- n_offspring
    
    ## If there are any offspring remaining, implement ring vaccination
    if (n_offspring > 0) {
      
      ## If vaccine is available but the infection is asymptomatic, no ring vaccination
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
        tdf[(current_max_id+1):(current_max_id+n_offspring), "quarantined"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "time_quarantined_absolute"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring_new"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring_quarantine"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring_post_pruning"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "offspring_generated"] <- FALSE
        
      } else {
        
        ## If vaccine is not available yet, no ring vaccination
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
          tdf[(current_max_id+1):(current_max_id+n_offspring), "quarantined"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "time_quarantined_absolute"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring_new"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring_quarantine"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring_post_pruning"] <- NA
          tdf[(current_max_id+1):(current_max_id+n_offspring), "offspring_generated"] <- FALSE
          
        ## If vaccine is available and the infection is symptomatic, implement ring vaccination
        } else {
          
          ## Creating vectors to store whether or not secondary infections are successfully ring-vaccinated, protected and/or prevented
          are_they_vaccinated <- vector(mode = "integer", length = n_offspring) ## vector of whether secondary infections get vaccinated
          did_they_have_potential_protection <- vector(mode = "integer", length = n_offspring)
          infection_retained <- vector(mode = "integer", length = n_offspring)  ## vector of whether secondary infections get retained (i.e. not prevented by ring vaccination)
          infection_retained[1:length(infection_retained)] <- 1                 ## default to infections being retained; and then flow through below to see if they get removed
          
          ## Looping over secondary infections and evaluating whether they're vaccinated, protected and/or prevented
          for (i in 1:n_offspring) {
            
            # if infection occurs AFTER (potential) ring-vaccination
            if (time_to_secondary_vaccination <= new_times[i]) { 
              are_they_vaccinated[i] <- rbinom(n = 1, size = 1, prob = vaccine_coverage)
              did_they_have_potential_protection[i] <- ifelse(time_to_secondary_vaccination_protection <= new_times[i], 1, 0)
              were_they_protected <- rbinom(n = 1, size = 1, prob = vaccine_efficacy_infection * did_they_have_potential_protection[i])
              infection_retained[i] <- ifelse(are_they_vaccinated[i] == 1 & were_they_protected == 1, 0, 1)
            
            # If infection would otherwise occur BEFORE ring-vaccination
            } else { 
              are_they_vaccinated[i] <- 0 ## eliding together "unvaccinated" and "vaccinated after infection occurs"
              infection_retained[i] <- 1 ## note that implicitly here we're implicitly "saying" these folks aren't vaccinated.
            }
          }
          new_n_offspring <- sum(infection_retained)
          tdf$n_offspring_post_pruning[idx] <- new_n_offspring
          
          ## Pruning the secondary infections after applying ring-vaccination
          retained_index <- which(infection_retained == 1)                                                                             # which secondary infections were NOT averted by ring-vaccination and thus are retained for inclusion in the dataframe
          new_new_times <- new_times[retained_index]
          vaccinated <- ifelse(are_they_vaccinated[retained_index] == 0, 0, 1)                                             # of the retained infections, which are vaccinated
          time_vaccinated <- ifelse(are_they_vaccinated[retained_index] == 1, time_to_secondary_vaccination, NA)           # of the retained infections, when are they vaccinated (relative to infection time of index)
          vaccinated_before_infection <- ifelse(are_they_vaccinated[retained_index] == 0, NA, 1)                           # (currently we combine all individuals not vaccinated and vaccinated after infection into "unvaccinated", so all vaccinated individuals necessarily got vaccinated before infection)
          time_protected <- ifelse(are_they_vaccinated[retained_index] == 1, time_to_secondary_vaccination_protection, NA) # of the retained infections, when are they protected (relative to infection time of index)
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
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "quarantined"] <- NA
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "time_quarantined_absolute"] <- NA
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring"] <- NA
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring_new"] <- NA
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring_quarantine"] <- NA
            tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring_post_pruning"] <- NA
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

