## Branching process simulation with ring vaccination
### Modified version of chain_sim_susc from bpmodels (https://github.com/epiverse-trace/bpmodels) to include ring vaccination.
### NOTE: Still potentially need to add:
###     1) Probability of a contact not being identified (similar to Kucharski et al' Ebola paper in EID)
###     2) vaccine_start being derived based on cumulative number of cases detected
###
### Function Arguments:
###   offspring = offspring distribution 
###   mn_offspring = mean of the offspring distribution
##    disp_offspring = overdispersion of the offspring distribution (negative binomial only)
###   generation_time = function returning draw from the generation time distribution
###   t0 = starting time
###   tf = final timepoint (won't go to this if check_final_size is reached first)
###   pop = population
###   check_final_size = number of infected individuals the branching process will simulate before stopping
###   initial_immune = number of individuals initially immune to the pathogen
###   seeding cases = starting number of infections
###   prop_asymptomatic = proportion of infections that are asymptomatic
###   infection_to_onset = function returning draw from infection to symptom onset distribution
###   vaccine_start = time at which vaccine becomes available
###   vaccine_coverage = proportion of eligible individuals who receive the vaccination
###   vaccine_efficacy_infection = vaccine efficacy against infection 
###   vaccine_efficacy_transmission = vaccine efficacy against onwards transmission in vaccinated individuals who are infected
###   vaccine_logistical_delay = delay between vaccine becoming available and people being vaccinated
###   vaccine_protection_delay = delay between receiving vaccination and becoming protected

chain_sim_susc_ring_vacc <- function(offspring = c("pois", "nbinom"), mn_offspring, disp_offspring, 
                                     generation_time, t0 = 0, tf = Inf, pop, check_final_size, initial_immune = 0,
                                     seeding_cases, prop_asymptomatic, infection_to_onset,                                 
                                     vaccine_start, vaccine_coverage, vaccine_efficacy_infection,                         
                                     vaccine_efficacy_transmission, vaccine_logistical_delay, vaccine_protection_delay) {
  
  ## Checking correct input of offspring distributions and setting up the offspring function
  offspring <- match.arg(offspring)
  if (offspring == "pois") {
    if (!missing(disp_offspring)) {
      warning("argument disp_offspring not used for\n poisson offspring distribution.")
    }
    offspring_fun <- function(n, susc) {
      truncdist::rtrunc(n, spec = "pois", lambda = mn_offspring * susc/pop, b = susc)
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
  
  ## Creating table to store branching process tree
  tdf <- data.frame(id = seq_along(1:seeding_cases),                                              # Unique Infection ID
                    ancestor = NA_integer_,                                                       # Who infected this person
                    generation = 1L,                                                              # Which generation of onwards infections are they in
                    time_infection = t0 + seq(from = 0, to = 0.01, length.out = seeding_cases),   # Time of infection (note the small seq() added is to avoid non-unique infection times for starting infections)
                    offspring_generated = FALSE,                                                  # Dummy variable indicating whether secondary infections have been generated for this infection
                    time_onset = NA,                                                              # Symptom onset time relative to time of infection
                    vaccinated = 0,                                                               # Whether or not the individual receives a vaccination
                    time_vaccinated = NA,                                                         # The time at which the individual is vaccinated (NA if doesn't get vaccination)
                    time_protected = NA,                                                          # The time at which the individual receives protection from the vaccination
                    asymptomatic = 0)                                                             # Whether or not the individual's infection is asymptomatic or not
                    ## note that currently time_onset is relative to time of infection (time) currently; all other times are relative to start of the outbreak 
  
  ## Setting up the simulation initial conditions (initial time, initial size of susceptible population etc)
  susc <- pop - initial_immune - 1L
  t <- t0
  
  ## Running the branching process - iterating over all infected people and drawing i) the number they infect (from offspring dist.);
  ##                                                                           and ii) the timings of these infections (from generation time dist.)
  
  ### While loop keeps generating new infections whilst:
  ####   1) There are still susceptible people left;
  ####   2) The latest time in the timeframe is < tdf AND total number of infected people is less than check_final_size
  while (any(tdf$time_infection[!tdf$offspring_generated] <= tf) & nrow(tdf) <= check_final_size & susc > 0) {
    
    ## Find the earliest infection in the transmission tree whose infections have not yet been calculated (the "parent" of the secondary infections)
    time_infection_index <- min(tdf$time_infection[!tdf$offspring_generated])              # Note: Is not an issue in practice, but I don't think this is currently set up to handle >= 2 infections with same infection time currently
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
    
    ## Generating number of people infected by the index case (the "parent")
    n_offspring <- offspring_fun(1, susc) 
    
    ### Note: We modify and reduce transmissibility of the index case, if they were vaccinated and protection has arisen
    ###       by the time they got infected (i.e. vaccine protection was insufficient but it had developed).
    ###       We make assumption that you only get this reduction in transmission (mediated by vaccine_efficacy_transmission) if
    ###       you were protected by the time you were infected (and protection failed). If you had been vaccinated and then were
    ###       infected before protection developing, we assume no impact on transmission (i.e. you derived no benefit from the vaccine).
    if (index_vaccinated == 1 & (time_protected < time_infection_index)) {  
      n_offspring <- sum(rbinom(n = n_offspring, size = 1, prob = 1 - vaccine_efficacy_transmission))
    }
    if (n_offspring %% 1 > 0) { # Checking offspring function is correctly returning integers
      stop("Offspring distribution must return integers")
    }
    tdf$offspring_generated[idx] <- TRUE
    
    ## Assigning times of these infections
    if (n_offspring > 0) {
      new_times <- generation_time(n_offspring)
      if (any(new_times < 0)) {
        stop("Generation times must be >= 0.")
      }
      
      ## Modelling Ring Vaccination - All eligible/willing (i.e. it's a function of coverage) secondary cases of an index case receive the vaccination following onset 
      ##                              of symptoms (after a small logistical delay) in the index case (it follows that secondary cases of asymptomatic index cases do not get the vaccination).
      ##                              Vaccinated individuals develop protection some time after receiving the vaccination. If the individual develops
      ##                              protection BEFORE the theoretical time they would have been infected index case, then the infection is potentially
      ##                              averted (with this decided probabilistically based on the vaccine efficacy against infection). If the individual doesn't receive the
      ##                              vaccine, they are not protected. If the individual receives the vaccination, but are infected before protection arises, they 
      ##                              are not protected. If the individual receives the vaccination, protection develops and they are infected anyway, we assume that 
      ##                              the vaccine contributes to reduced transmissibility (with this decided probabilistically based on the vaccine efficacy against transmission).
      ## Note: Things still to do:
      ### 1) Think we're possibly missing the folks currently who get infected before the vaccination takes place but who are asymptomatic and so infection is missed and are vaccinated anyway
      ###    I don't think this makes any difference in practice (as they would never get any benefit of vaccination) but might be worth making this explicit
      ### 2) Check there's nothing weird going on with index case asymptomatic and who's getting vaccinated in a given transmission tree 
      ###
      ### Note: By definition, those included in the generated dataframe (transmission tree) are ONLY those infected, and specifically will be either:
      ###   i) those unvaccinated
      ###  ii) those who got vaccinated but too late (therefore no benefit of vaccination);
      ### iii) those who got vaccinated but it failed to protect (in which case, reduction in transmission due to vaccination is relevant and taken care of above)
      
      ## If the vaccine isn't yet available, no-one gets vaccinated and all new infections are added into the tree
      if (t < vaccine_start) {
        asymptomatic <- rbinom(n = n_offspring, size = 1, prob = prop_asymptomatic)
        new_df <- data.frame(id = current_max_id + seq_len(n_offspring), 
                             ancestor = id_parent,
                             generation = gen_parent + 1L, 
                             time = new_times + t_parent, 
                             offspring_generated = FALSE,
                             time_onset = NA,
                             vaccinated = 0,
                             time_vaccinated = NA,
                             time_protected = NA,
                             asymptomatic = asymptomatic)
        tdf <- rbind(tdf, new_df)
        
      ## If the vaccine is available
      } else {
        
        ## If the index case is asymptomatic, all secondary cases arising from them won't be vaccinated because the index infection won't be detected (asymptomatic)
        if (index_asymptomatic == 1) {
          asymptomatic <- rbinom(n = n_offspring, size = 1, prob = prop_asymptomatic)
          new_df <- data.frame(id = current_max_id + seq_len(n_offspring), 
                               ancestor = id_parent,
                               generation = gen_parent + 1L, 
                               time = new_times + t_parent, 
                               offspring_generated = FALSE,
                               time_onset = NA,
                               vaccinated = 0,
                               time_vaccinated = NA,
                               time_protected = NA,
                               asymptomatic = asymptomatic)
          tdf <- rbind(tdf, new_df)
          
        ## If index case is symptomatic, then the contacts are vaccinated and the potential secondary infections get ring vaccinated
        ## Below, we calculate when they are vaccinated and whether they develop protection in time to prevent infection (i.e. does protection
        ## arise before the time of their infection by the index case - we calculate all of this below.
        } else {
          
          are_they_vaccinated <- vector(mode = "integer", length = n_offspring) ## vector of whether secondary infections get vaccinated
          infection_retained <- vector(mode = "integer", length = n_offspring)  ## vector of whether secondary infections get retained (i.e. not prevented by ring vaccination)
          infection_retained[1:length(infection_retained)] <- 1                 ## default to infections being retained; and then flow through below to see if they get removed
          
          ## Calculating time to vaccination protection
          time_to_secondary_vaccination <- onset_time_index_case + vaccine_logistical_delay                    ## Time between index case infected and secondary cases ring vaccinated
          time_to_secondary_vaccination_protection <- time_to_secondary_vaccination + vaccine_protection_delay ## Time between index case infected and secondary cases protected with vaccination
          
          ## Checking whether vaccination occurs before or after the secondary cases have been generated; and if so, whether or not the infection is prevented by the vaccine
          for (i in 1:n_offspring) {
            if (time_to_secondary_vaccination <= new_times[i]) {
              are_they_vaccinated[i] <- rbinom(n = 1, size = 1, prob = vaccine_coverage)
              did_they_have_potential_protection <- ifelse(time_to_secondary_vaccination_protection <= new_times[i], 1, 0)
              were_they_protected <- rbinom(n = 1, size = 1, prob = vaccine_efficacy_infection * did_they_have_potential_protection)
              infection_retained[i] <- ifelse(are_they_vaccinated[i] == 1 & were_they_protected == 1, 0, 1)
            } else {
              infection_retained[i] <- 1
            }
          }
          
          ## Calculating number of infections which actually occur despite vaccination, and creating relevant inputs to dataframe for these successful secondary infections
          new_n_offspring <- sum(infection_retained)
          new_new_times <- new_times[which(infection_retained == 1)]
          time_vaccinated <- ifelse(are_they_vaccinated == 1, time_to_secondary_vaccination, NA)
          time_protected <- ifelse(are_they_vaccinated == 1, time_to_secondary_vaccination_protection, NA)
          asymptomatic <- rbinom(n = new_n_offspring, size = 1, prob = prop_asymptomatic)
          
          ## Checking whether there are any infections to add to the table
          if (new_n_offspring != 0) {
            new_df <- data.frame(id = current_max_id + seq_len(new_n_offspring), 
                                 ancestor = id_parent,
                                 generation = gen_parent + 1L, 
                                 time = new_new_times + t_parent, 
                                 offspring_generated = FALSE,
                                 time_onset = NA,
                                 vaccinated = are_they_vaccinated[which(infection_retained == 1)],
                                 time_vaccinated = time_vaccinated[which(infection_retained == 1)] + t_parent,
                                 time_protected = time_protected[which(infection_retained == 1)] + t_parent,
                                 asymptomatic = asymptomatic)
            tdf <- rbind(tdf, new_df)
          }
        }
      }
    }
    susc <- susc - n_offspring
  }
  
  ## Bookkeeping for end of function - ordering the dataframe by time of infection, id etc
  tdf <- tdf[tdf$time <= tf, ]
  tdf <- tdf[order(tdf$time, tdf$id), ]
  tdf$offspring_generated <- NULL
  tdf$abs_time_onset <- tdf$time_onset + tdf$time
  
  return(tdf)
}

## Testing and running the model

chain_sim_susc_ring_vacc(offspring = "pois"
                         mn_offspring = 20, 
                         disp_offspring = , 
                                     generation_time, t0 = 0, tf = Inf, pop, check_final_size, initial_immune = 0,
                                     seeding_cases, prop_asymptomatic, infection_to_onset,                                 
                                     vaccine_start, vaccine_coverage, vaccine_efficacy_infection,                         
                                     vaccine_efficacy_transmission, vaccine_logistical_delay, vaccine_protection_delay) 


offspring <- "pois"
mn_offspring <- 20
generation_time <- function(n) {  rgamma(n, shape = 14, rate = 2) } # generation_time <- function(n) {  return(rep(7, n)) } 
t0 <- 0
tf <- Inf
pop <- 10^8
check_final_size <- 10000
initial_immune <- 0
infection_to_onset <- function(n) { rgamma(n, shape = 1, rate = 2) } # infection_to_onset <- function(n) { return(rep(2, n)) } 
vaccine_start <- 1
vaccine_coverage <- 1
vaccine_efficacy_infection <- 0.5
vaccine_efficacy_transmission <- 0
vaccine_logistical_delay <- 1
vaccine_protection_delay <- 1
seeding_cases <- 1
prop_asymptomatic <- 0



x <- tdf %>%
  group_by(ancestor) %>%
  summarise(n =n ())

par(mfrow = c(2, 1))
hist(x$n)
hist(offspring_fun(n = 5000, susc = pop))

mean(x$n)
mean(offspring_fun(n = 5000, susc = pop))



## Tracking IDs of infected individuals in sub-trees
# traverse_tree <- function(check_id, df, counted) {
#   
#   # Get the ids of individuals directly infected by the current individual
#   children <- df %>% 
#     filter(ancestor == check_id) %>% 
#     pull(id)
#   
#   # Remove children that have already been counted
#   children <- setdiff(children, counted)
#   
#   # If there are no children left, return empty list
#   if (length(children) == 0) {
#     return(list())
#   }
#   
#   # Update the list of counted individuals
#   counted <- c(counted, children)
#   
#   # Recursively call this function on all the children and collect the ids
#   infected_individuals <- c(children, unlist(lapply(children, function(x) traverse_tree(x, df, counted))))
#   
#   # Return the list of infected individuals
#   return(infected_individuals)
# }
# ## (Lightly) Modified version of chain_sim_susc from bpmodels (https://github.com/epiverse-trace/bpmodels)
# chain_sim_susc <- function (offspring = c("pois", "nbinom"), mn_offspring, disp_offspring, 
#                             generation_time, t0 = 0, tf = Inf, pop, check_final_size, initial_immune = 0) {
#   
#   ## Checking correct input of offspring distributions and setting up the offspring function
#   offspring <- match.arg(offspring)
#   if (offspring == "pois") {
#     if (!missing(disp_offspring)) {
#       warning("argument disp_offspring not used for\n poisson offspring distribution.")
#     }
#     offspring_fun <- function(n, susc) {
#       truncdist::rtrunc(n, spec = "pois", lambda = mn_offspring * susc/pop, b = susc)
#     }
#   }
#   else if (offspring == "nbinom") {
#     if (disp_offspring <= 1) {
#       stop("Offspring distribution 'nbinom' requires argument\n disp_offspring > 1. Use 'pois' if there is no overdispersion.")
#     }
#     offspring_fun <- function(n, susc) {
#       new_mn <- mn_offspring * susc/pop
#       size <- new_mn/(disp_offspring - 1)
#       truncdist::rtrunc(n, spec = "nbinom", b = susc, mu = new_mn, size = size)
#     }
#   }
#   
#   ## Creating table to store branching process tree
#   tdf <- data.frame(id = 1L, ancestor = NA_integer_, generation = 1L, time = t0, offspring_generated = FALSE)
#   
#   ## Setting up the simulation initial conditions (initial time, initial size of susceptible population etc)
#   susc <- pop - initial_immune - 1L
#   t <- t0
#   
#   ## Running the branching process - iterating over all infected people and drawing i) the number they infect (from offspring dist.);
#   ## and ii) the timings of these infections (from generation time dist.)
#   
#   ### whilst there are still susceptible people and the latest time in the timeframe is < tdf OR number of infected people is less than check_final_size
#   while ((any(tdf$time[!tdf$offspring_generated] <= tf) | nrow(tdf <= check_final_size)) & susc > 0) {
#     
#     ## Find the earliest infection in the transmission tree whose infections have not yet been calculated
#     t <- min(tdf$time[!tdf$offspring_generated])
#     idx <- which(tdf$time == t & !tdf$offspring_generated)[1] # does this ever get stuck with multiple people this is true for (e.g. infected @ exact same time?) - don't think so in practice
#     id_parent <- tdf$id[idx]
#     t_parent <- tdf$time[idx]
#     gen_parent <- tdf$generation[idx]
#     current_max_id <- max(tdf$id)
#     
#     ## Generating number of people infected by this infection
#     n_offspring <- offspring_fun(1, susc)
#     if (n_offspring %% 1 > 0) { # Checking offspring function is acting correctly (returning integers)
#       stop("Offspring distribution must return integers")
#     }
#     tdf$offspring_generated[idx] <- TRUE
#     
#     ## Assigning times of these infections (if the index infection does indeed infect anyone)
#     if (n_offspring > 0) {
#       new_times <- generation_time(n_offspring)
#       if (any(new_times < 0)) {
#         stop("Serial interval must be >= 0.")
#       }
#       new_df <- data.frame(id = current_max_id + seq_len(n_offspring), 
#                            time = new_times + t_parent, ancestor = id_parent, 
#                            generation = gen_parent + 1L, offspring_generated = FALSE)
#       tdf <- rbind(tdf, new_df)
#     }
#     susc <- susc - n_offspring
#   }
#   
#   ## Bookkeeping for end of function - ordering the dataframe by time of infection, id etc
#   tdf <- tdf[tdf$time <= tf, ]
#   tdf <- tdf[order(tdf$time, tdf$id), ]
#   tdf$offspring_generated <- NULL
#   
#   return(tdf)
# }

