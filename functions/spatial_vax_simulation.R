## Function to calculate the geographical coordinates of infection offspring given parent coordinates
spatial_calc <- function(parent_x_coord, 
                         parent_y_coord,
                         n_offspring, 
                         spatial_kernel) {
  
  
  ## Drawing spatial - note that currently time and distance/direction are completely uncorrelated atm
  distance <- spatial_kernel(n_offspring)
  direction <- runif(n_offspring, 0, 2*pi)
  direction_degrees <- 360 * direction / (2 * pi) 
  subtract <- direction_degrees %/% 90
  direction_degrees_sub90 <- direction_degrees - 90 * subtract
  direction_radians_sub90 <- 2 * pi * direction_degrees_sub90 / 360
  opposite <- distance * sin(direction_radians_sub90)
  adjacent <- sqrt(distance^2 - opposite^2)
  
  ## Creating temporary dataframe
  tdf <- data.frame(distance = NA_real_,
                    x_coordinate = NA_real_,
                    y_coordinate = NA_real_,
                    overall_distance = NA_real_,
                    stringsAsFactors = FALSE)
  
  for (i in 1:n_offspring) {
    if(direction_degrees[i] < 90) {
      tdf[i, "x_coordinate"] <- parent_x_coord + opposite[i]
      tdf[i, "y_coordinate"] <- parent_y_coord + adjacent[i]
      tdf[i, "distance"] <- distance[i]
      tdf[i, "overall_distance"] <- sqrt(tdf[i, "x_coordinate"]^2 + tdf[i, "y_coordinate"]^2)
    } else if (direction_degrees[i] >= 90 & direction_degrees[i] < 180) {
      tdf[i, "x_coordinate"] <- parent_x_coord + opposite[i]
      tdf[i, "y_coordinate"] <- parent_y_coord - adjacent[i]
      tdf[i, "distance"] <- distance[i]
      tdf[i, "overall_distance"] <- sqrt(tdf[i, "x_coordinate"]^2 + tdf[i, "y_coordinate"]^2)
    } else if (direction_degrees[i] >= 180 & direction_degrees[i] < 270) {
      tdf[i, "x_coordinate"] <- parent_x_coord - opposite[i]
      tdf[i, "y_coordinate"] <- parent_y_coord - adjacent[i]
      tdf[i, "distance"] <- distance[i]
      tdf[i, "overall_distance"] <- sqrt(tdf[i, "x_coordinate"]^2 + tdf[i, "y_coordinate"]^2)
    } else {
      tdf[i, "x_coordinate"] <- parent_x_coord - opposite[i]
      tdf[i, "y_coordinate"] <- parent_y_coord + adjacent[i]
      tdf[i, "distance"] <- distance[i]
      tdf[i, "overall_distance"] <- sqrt(tdf[i, "x_coordinate"]^2 + tdf[i, "y_coordinate"]^2)
    }
  }
  return(tdf)
}

## spatial vaccination branching process
spatial_vax_bp_sim <- function(## Transmission Parameters
                               offspring = c("pois", "nbinom"),   # offspring distribution 
                               mn_offspring,                      # mean of the offspring distribution
                               disp_offspring,                    # overdispersion of the offspring distribution (if negative binomial)
                               spatial_kernel,                    # spatial kernel for onwards transmission (distance between infections)
                               
                               ## Natural History Parameters
                               generation_time,                   # generation time distribution
                               prop_asymptomatic,                 # probability of being asymptomatic
                               infection_to_onset,                # time from infection to symptom onset distribution
                               prob_hosp,                         # probability of an infected individual being hospitalised
                               hospitalisation_delay,             # delay between becoming infected and becoming hospitalised
                               
                               ## Vaccine-Related Parameters
                               vaccine_campaign_radius,           # radius of the spatial vaccination campaign
                               vaccine_coverage,                  # probability that each eligible individual gets vaccinated
                               vaccine_efficacy_infection,        # vaccine efficacy against infection
                               vaccine_efficacy_transmission,     # reduction in transmissibility of breakthrough infections in vaccinated individuals
                               vaccine_logistical_delay,          # delay between symptom onset and vaccination of contacts (and contacts of contacts) occurring
                               vaccine_protection_delay,          # delay between vaccination and protection developing
                               vaccine_efficacy_disease,          # vaccine efficacy against disease
                               detection_threshold,               # hospitalisation threshold at which detection occurs
                               
                               ## Quarantine Related Parameters
                               time_to_quarantine,                # time delay between trigger and quarantining (either symptoms of the primary infection, or symptoms in the secondary infection)
                               prob_quarantine_contact_traced,    # probability that the individual successfully isolates given they're a contact of an infection
                               prob_quarantine_symptoms,          # probability that the individual successfully isolates given they have symptoms
                               quarantine_efficacy,               # effectiveness of the quarantine at reducing onwards transmission
                               
                               ## Miscellaneous Parameters
                               t0 = 0, 
                               tf = Inf, 
                               population,
                               check_final_size,
                               initial_immune,
                               seeding_cases,
                               seed
                               ) {
  
  ## Setting the seed
  set.seed(seed)
  
  ## Setting up the number of susceptibles
  susc <- population - initial_immune
  
  ## Setting up the offspring distribution
  ### Note that the modification for susceptible depletion can sometimes give some slightly off results around R=1
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
  
  ## Pre-allocate a dataframe with the maximum number of individuals to be simulated
  max_cases <- check_final_size
  tdf <- data.frame(
    id = integer(max_cases),
    ancestor = integer(max_cases),
    generation = integer(max_cases),
    time_infection = NA_real_,
    time_onset = numeric(max_cases),
    distance = NA_real_,
    x_coordinate = NA_real_,
    y_coordinate = NA_real_,
    overall_distance = NA_real_,
    hospitalised = NA_integer_, 
    vaccinated = integer(max_cases),
    time_vaccinated = numeric(max_cases),
    vaccinated_before_infection = integer(max_cases),
    time_protected = numeric(max_cases),
    protected_before_infection = integer(max_cases),
    asymptomatic = integer(max_cases),
    quarantined = integer(check_final_size),
    time_quarantined_absolute = NA_real_,
    n_offspring = integer(max_cases),
    n_offspring_new = integer(max_cases),
    n_offspring_quarantine = integer(max_cases),
    n_offspring_new_new = integer(max_cases),
    offspring_generated = FALSE,
    stringsAsFactors = FALSE)
  
  ## Initialize the dataframe with the seeding cases
  tdf[1:seeding_cases, ] <- data.frame(
    id = seq_len(seeding_cases),
    ancestor = NA_integer_,
    generation = 1L,
    time_infection = t0 + seq(from = 0, to = 0.01, length.out = seeding_cases),
    time_onset = NA,
    distance = 0,
    x_coordinate = 0,
    y_coordinate = 0,
    overall_distance = 0,
    hospitalised = rbinom(n = seeding_cases, size = 1, prob = prob_hosp),
    vaccinated = 0,
    time_vaccinated = NA,
    vaccinated_before_infection = NA,
    time_protected = NA,
    protected_before_infection = NA,
    asymptomatic = integer(seeding_cases),
    quarantined = integer(seeding_cases),                   
    time_quarantined_absolute = NA_real_,
    n_offspring = NA_integer_,
    n_offspring_new = NA_integer_,
    n_offspring_quarantine = NA_integer_,
    n_offspring_new_new = NA_integer_,
    offspring_generated = FALSE)
  
  time_infection_index <- t0
  
  ## While we haven't hit the simulation cap size (check_final_size) and any infections exist where we have not yet generated the requisite offspring, 
  ## continue to generate infections
  while ((any(is.na(tdf$n_offspring)) & nrow(tdf) <= check_final_size & susc > 0)) {
    
    ## Getting the timings of the earliest/oldest infection we haven't yet generated infections for - this is the "INDEX INFECTION"
    time_infection_index <- min(tdf$time_infection[tdf$offspring_generated == 0 & !is.na(tdf$time_infection)])
    idx <- which(tdf$time_infection == time_infection_index & !tdf$offspring_generated)[1] # get the id of the earliest unsimulated infection
    previous_parent_id <- tdf$ancestor[idx]
    previous_parent_idx <- which(tdf$id == previous_parent_id)
    id_parent <- tdf$id[idx]                                                               # parent of the earliest unsimulated infection
    t_parent <- tdf$time_infection[idx]                                                    # infection time of the earliest unsimulated infection
    gen_parent <- tdf$generation[idx]                                                      # generation of the earliest unsimulated infection
    current_max_id <- max(tdf$id, na.rm = TRUE)                                            # total number of infections in the dataframe currently (so we can figure out how to label the new infections
    index_vaccinated <- tdf$vaccinated[idx]                                                # whether or not the index case (the "parent") is vaccinated
    time_vaccinated <- tdf$time_vaccinated[idx]                                            # when the index case (the "parent") was vaccinated
    time_protected <- tdf$time_protected[idx]                                              # when the index case (the "parent") was protected
    onset_time_index_case <- infection_to_onset(n = 1)                                     # generate the time from infection to symptom onset for the index case
    index_asymptomatic <- tdf$asymptomatic[idx]                                            # whether or not the index case (the "parent") is asymptomatic (influences whether contacts get ring vaccinated or not)
    tdf$time_onset[idx] <- ifelse(index_asymptomatic == 0, onset_time_index_case, NA)      # --
    total_hospitalised <- sum(tdf$hospitalised, na.rm = TRUE)                              # total number of individuals hospitalised (used as spatial vaccination trigger)
    
    ## Calculating whether or not the individual isolates/quarantines
    ### note that here I swap from using "parent" to refer to "index" and instead as the prior infector. Need to sort this at some point.
    ### (it's just syntax, the code is actually right/doing the right thing)
    ## Extracting parent information
    if (identical(previous_parent_idx, integer(0))) { ## for the seeding cases
      parent_infection_time <- 0
      parent_asymptomatic <- 1
      parent_onset_time <- 0
    } else { ## for everyone else
      parent_infection_time <- tdf$time_infection[previous_parent_idx];
      parent_asymptomatic <- tdf$asymptomatic[previous_parent_idx];
      parent_onset_time <- ifelse(parent_asymptomatic == 0, tdf$time_onset[previous_parent_idx], NA)
    }
    
    ## If the parent is symptomatic, this triggers quarantine in secondary infections relative to timing of symptoms in parent
    if (parent_asymptomatic == 0) {
      
      index_quarantine <- rbinom(n = 1, size = 1, prob = prob_quarantine_contact_traced)                      # whether or not the index infection isolates
      index_quarantine_time <- ifelse(index_quarantine == 1, time_to_quarantine(n = 1), NA)                   # if the infection isolates, how soon after symptom onset they do so
      tdf$quarantined[idx] <- index_quarantine                                                                # adding quarantine indicator to storage dataframe
      absolute_quarantine_time <- parent_infection_time + parent_onset_time + index_quarantine_time           # adding quarantine time in absolute calendar time to the storage dataframe
      tdf$time_quarantined_absolute[idx] <- absolute_quarantine_time
      
      ## If the parent is asymptomatic but index is symptomatic, this triggers quarantine relative to timing in index
    } else if (parent_asymptomatic == 1 & index_asymptomatic == 0) {
      
      index_quarantine <- rbinom(n = 1, size = 1, prob = prob_quarantine_symptoms)                               # whether or not the index infection isolates
      index_quarantine_time <- ifelse(index_quarantine == 1, time_to_quarantine(n = 1), NA)                      # if the infection isolates, how soon after symptom onset they do so
      tdf$quarantined[idx] <- index_quarantine                                                                   # adding quarantine indicator to storage dataframe
      absolute_quarantine_time <- time_infection_index + onset_time_index_case + index_quarantine_time           # adding quarantine time in absolute calendar time to the storage dataframe
      tdf$time_quarantined_absolute[idx] <- absolute_quarantine_time
      
      ## If both the parent and index are asymptomatic, no quarantining can possibly happen
    } else {
      index_quarantine <- 0
      index_quarantine_time <- NA
      absolute_quarantine_time <- NA
    }
    
    # Generating offspring for this infection
    n_offspring <- offspring_fun(1, susc) 
    tdf$n_offspring[idx] <- n_offspring
    tdf$offspring_generated[idx] <- TRUE
    
    ## If infection was previously vaccinated and is a breakthrough infection, account for reduced transmissibility
    if (index_vaccinated == 1) {
      if (!is.na(time_protected) & time_protected < time_infection_index) {
        n_offspring <- sum(rbinom(n = n_offspring, size = 1, prob = 1 - vaccine_efficacy_transmission))
      }
    }
    tdf$n_offspring_new[idx] <- n_offspring
    new_times <- generation_time(n_offspring)
    absolute_new_times <- time_infection_index + new_times
    
    ## If index infection quarantines, reduce secondary infections - note that quarantine only occurs if infection has symptoms. Only do this if there are offspring to avert.
    if (index_quarantine == 1 & index_asymptomatic == 0 & n_offspring != 0) {
      
      # Implement quarantining for index infection
      index_n_offspring <- implement_quarantine(symptom_onset_time = 0,  # have wrapped symptom onset time into absolute_quarantine time below
                                                quarantine_time = absolute_quarantine_time,
                                                n_offspring = n_offspring,
                                                offspring_infection_times = absolute_new_times,
                                                quarantine_efficacy = quarantine_efficacy)
      
      # Updating number offspring, their infection times and characteristics to reflect removals due to quarantining
      n_offspring <- index_n_offspring$updated_n_offspring
      new_times <- index_n_offspring$updated_infection_times - time_infection_index
      
    }
    tdf$n_offspring_quarantine[idx] <- n_offspring
    
    ## If there are any offspring remaining, implement spatial vaccination
    if (n_offspring > 0) {
      
      ## Generating location of each newly infected individual
      new_locations <- spatial_calc(parent_x_coord = tdf$x_coordinate[idx],
                                    parent_y_coord = tdf$y_coordinate[idx],
                                    n_offspring = n_offspring,
                                    spatial_kernel = spatial_kernel)
      
      ## Infections occur before vaccine starts
      if (total_hospitalised < detection_threshold) {
        tdf[(current_max_id+1):(current_max_id+n_offspring), "id"] <- c(current_max_id + seq_len(n_offspring))
        tdf[(current_max_id+1):(current_max_id+n_offspring), "ancestor"] <- id_parent
        tdf[(current_max_id+1):(current_max_id+n_offspring), "generation"] <- gen_parent + 1L
        tdf[(current_max_id+1):(current_max_id+n_offspring), "time_infection"] <- new_times + t_parent
        tdf[(current_max_id+1):(current_max_id+n_offspring), "vaccinated"] <- 0
        tdf[(current_max_id+1):(current_max_id+n_offspring), "hospitalised"] <- rbinom(n = n_offspring, size = 1, prob = prob_hosp)
        tdf[(current_max_id+1):(current_max_id+n_offspring), "time_vaccinated"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "vaccinated_before_infection"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "time_protected"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "protected_before_infection"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "asymptomatic"] <- rbinom(n = n_offspring, size = 1, prob = prop_asymptomatic)
        tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring_new"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring_quarantine"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring_new_new"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "offspring_generated"] <- FALSE
        tdf[(current_max_id+1):(current_max_id+n_offspring), "x_coordinate"] <- new_locations$x_coordinate
        tdf[(current_max_id+1):(current_max_id+n_offspring), "y_coordinate"] <- new_locations$y_coordinate
        tdf[(current_max_id+1):(current_max_id+n_offspring), "distance"] <- new_locations$distance
        tdf[(current_max_id+1):(current_max_id+n_offspring), "overall_distance"] <- new_locations$overall_distance
        tdf[(current_max_id+1):(current_max_id+n_offspring), "quarantined"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "time_quarantined_absolute"] <- NA
        tdf$n_offspring_new_new[idx] <- n_offspring
        
      } else {

        ### Calculating the time when the nth hospitalisation gets admitted and they trigger the campaign
        temp <- tdf[tdf$hospitalised == 1 & !is.na(tdf$hospitalised), ]     # getting all the hospitalised individuals
        time_trigger_infection <- temp$time_infection[detection_threshold]  # getting the time of infection of the trigger individual
        x_trigger_infection <- temp$x_coordinate[detection_threshold]       # x-coord of trigger individual
        y_trigger_infection <- temp$y_coordinate[detection_threshold]       # y-coord of trigger individual
        
        ## Creating storage for various quantities relating to vaccination
        time_to_secondary_vaccination <- time_trigger_infection + hospitalisation_delay(1) + vaccine_logistical_delay  ## Time when geographical vaccination campaign starts and people vaccinated
        time_to_secondary_vaccination_protection <- time_to_secondary_vaccination + vaccine_protection_delay           ## Time between  geographical vaccination campaign starting and people protected by the vaccination
        
        are_they_vaccinated <- vector(mode = "integer", length = n_offspring)                ## vector of whether secondary infections get vaccinated
        did_they_have_potential_protection <- vector(mode = "integer", length = n_offspring) ## vector of whether vaccinnated AND protection has developed
        were_they_protected <- vector(mode = "integer", length = n_offspring)                ## did vaccination successfully protect them
        infection_retained <- vector(mode = "integer", length = n_offspring)                 ## vector of whether secondary infections get retained (i.e. not prevented by vaccination)
        infection_retained[1:length(infection_retained)] <- 1                                ## default to infections being retained; and then flow through below to see if they get removed
        
        ## Looping over secondary infections and evaluating whether they're vaccinated, protected and/or prevented
        for (i in 1:n_offspring) {
          
          ## Calculating the distance the new infection is from the trigger infection
          distance_from_trigger_x <- abs(new_locations$x_coordinate[i] - x_trigger_infection)
          distance_from_trigger_y <- abs(new_locations$y_coordinate[i] - y_trigger_infection)
          distance_from_trigger <- sqrt(distance_from_trigger_x^2 + distance_from_trigger_y^2)
          
          ## If they're within the radius and can be potentially vaccinated before they would otherwise be infected
          if (time_to_secondary_vaccination <= (t_parent + new_times[i]) & distance_from_trigger <= vaccine_campaign_radius) { # if infection occurs AFTER (potential) vaccination
            are_they_vaccinated[i] <- rbinom(n = 1, size = 1, prob = vaccine_coverage)
            did_they_have_potential_protection[i] <- ifelse(are_they_vaccinated[i] == 1 & time_to_secondary_vaccination_protection <= (t_parent + new_times[i]), 1, 0)
            were_they_protected[i] <- rbinom(n = 1, size = 1, prob = vaccine_efficacy_infection * did_they_have_potential_protection[i])
            infection_retained[i] <- ifelse(are_they_vaccinated[i] == 1 & were_they_protected[i] == 1, 0, 1)
          
          # if infection occurs BEFORE vaccination can occur OR the infection is too far away from the case that triggers vaccination
          } else { 
            are_they_vaccinated[i] <- 0 ## eliding together "unvaccinated" and "vaccinated after infection occurs"
            infection_retained[i] <- 1 ## note that implicitly here we're "saying" these folks aren't vaccinated.
          }
        }
        new_n_offspring <- sum(infection_retained)
        tdf$n_offspring_new_new[idx] <- new_n_offspring
        
        ## Adding the infections that are retained into our overall dataframe 
        if (new_n_offspring != 0) {
          
          ## Defining index for which infections we keep (i.e. those vaccination doesn't successfully protect)
          retained_infections_index <- which(infection_retained == 1)
          
          ## Subsetting the times and locations of retained infections
          new_new_times <- new_times[retained_infections_index]                                                                      # The time of infection for retained infections
          new_new_locations <- new_locations[retained_infections_index, ]                                                            # The locations of retained infections
          
          ## Getting vaccination status and timing of vaccination relative to infection
          vaccinated <- ifelse(are_they_vaccinated[retained_infections_index] == 1, 1, 0)                                            # of the retained infections, which are vaccinated
          time_vaccinated <- ifelse(are_they_vaccinated[retained_infections_index] == 1, time_to_secondary_vaccination, NA)          # of the retained infections, when are they vaccinated (relative to infection time of index)
          vaccinated_before_infection <- ifelse(are_they_vaccinated[retained_infections_index] == 0, NA, 
                                          ifelse(are_they_vaccinated[retained_infections_index] == 1 & 
                                                   time_to_secondary_vaccination <= (t_parent + new_new_times),
                                                1, 0))
          time_protected <- ifelse(did_they_have_potential_protection[retained_infections_index] == 1, time_to_secondary_vaccination_protection, NA)
          protected_before_infection <- ifelse(are_they_vaccinated[retained_infections_index] == 0, NA, 
                                          ifelse(are_they_vaccinated[which(infection_retained == 1)] == 1 & 
                                                 time_to_secondary_vaccination_protection <= (t_parent + new_new_times), 
                                               1, 0))
          
          ## Calculating whether individuals are hospitalised
          hospitalised_vec <- vector(mode = "integer", length = new_n_offspring)
          for (i in 1:length(protected_before_infection)) {
            if (is.na(protected_before_infection[i]) | (!is.na(protected_before_infection[i]) & protected_before_infection[i] == 0)) {
              hospitalised_vec[i] <- rbinom(n = 1, size = 1, prob = prob_hosp)
            } else {
              hospitalised_vec[i] <- rbinom(n = 1, size = 1, prob = prob_hosp * (1 - vaccine_efficacy_disease))
            }
          }
          
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "id"] <- c(current_max_id + seq_len(new_n_offspring))
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "ancestor"] <- id_parent
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "generation"] <- gen_parent + 1L
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "time_infection"] <- new_new_times + t_parent
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "vaccinated"] <- vaccinated
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "hospitalised"] <- hospitalised_vec
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "time_vaccinated"] <- time_vaccinated
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "vaccinated_before_infection"] <- vaccinated_before_infection
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "time_protected"] <- time_protected
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "protected_before_infection"] <- protected_before_infection
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "asymptomatic"] <- rbinom(n = new_n_offspring, size = 1, prob = prop_asymptomatic)
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring"] <- NA
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring_new"] <- NA
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring_quarantine"] <- NA
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring_new_new"] <- NA
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "offspring_generated"] <- FALSE
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "x_coordinate"] <- new_new_locations$x_coordinate
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "y_coordinate"] <- new_new_locations$y_coordinate
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "distance"] <- new_new_locations$distance
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "overall_distance"] <- new_new_locations$overall_distance
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring"] <- NA
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring_new"] <- NA
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring_quarantine"] <- NA
          tdf[(current_max_id+1):(current_max_id+new_n_offspring), "n_offspring_post_pruning"] <- NA
        }
      }
    } else {
      tdf$n_offspring_new_new[idx] <- 0
    }
    susc <- susc - n_offspring
  }
  tdf <- tdf[tdf$time_infection <= tf, ]
  tdf <- tdf[order(tdf$time_infection, tdf$id), ] # moved to outside the bracket, check that doesn't mess anything up (I don't think so)
  return(tdf)
}
