## (Lightly) Modified version of chain_sim_susc from bpmodels (https://github.com/epiverse-trace/bpmodels)
chain_sim_susc <- function (offspring = c("pois", "nbinom"), mn_offspring, disp_offspring, 
                            generation_time, t0 = 0, tf = Inf, pop, check_final_size, initial_immune = 0) {
  
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
  tdf <- data.frame(id = 1L, ancestor = NA_integer_, generation = 1L, time = t0, offspring_generated = FALSE)
  
  ## Setting up the simulation initial conditions (initial time, initial size of susceptible population etc)
  susc <- pop - initial_immune - 1L
  t <- t0
  
  ## Running the branching process - iterating over all infected people and drawing i) the number they infect (from offspring dist.);
  ## and ii) the timings of these infections (from generation time dist.)
  
  ### whilst there are still susceptible people and the latest time in the timeframe is < tdf OR number of infected people is less than check_final_size
  while ((any(tdf$time[!tdf$offspring_generated] <= tf) | nrow(tdf <= check_final_size)) & susc > 0) {
    
    ## Find the earliest infection in the transmission tree whose infections have not yet been calculated
    t <- min(tdf$time[!tdf$offspring_generated])
    idx <- which(tdf$time == t & !tdf$offspring_generated)[1] # does this ever get stuck with multiple people this is true for (e.g. infected @ exact same time?) - don't think so in practice
    id_parent <- tdf$id[idx]
    t_parent <- tdf$time[idx]
    gen_parent <- tdf$generation[idx]
    current_max_id <- max(tdf$id)
    
    ## Generating number of people infected by this infection
    n_offspring <- offspring_fun(1, susc)
    if (n_offspring %% 1 > 0) { # Checking offspring function is acting correctly (returning integers)
      stop("Offspring distribution must return integers")
    }
    tdf$offspring_generated[idx] <- TRUE
    
    ## Assigning times of these infections (if the index infection does indeed infect anyone)
    if (n_offspring > 0) {
      new_times <- generation_time(n_offspring)
      if (any(new_times < 0)) {
        stop("Serial interval must be >= 0.")
      }
      new_df <- data.frame(id = current_max_id + seq_len(n_offspring), 
                           time = new_times + t_parent, ancestor = id_parent, 
                           generation = gen_parent + 1L, offspring_generated = FALSE)
      tdf <- rbind(tdf, new_df)
    }
    susc <- susc - n_offspring
  }
  
  ## Bookkeeping for end of function - ordering the dataframe by time of infection, id etc
  tdf <- tdf[tdf$time <= tf, ]
  tdf <- tdf[order(tdf$time, tdf$id), ]
  tdf$offspring_generated <- NULL
  
  return(tdf)
}

## (More Heavily) Modified version of chain_sim_susc from bpmodels (https://github.com/epiverse-trace/bpmodels)
## to include ring vaccination.
## Still to add in:
### cases missed (per Kucharski)
### >1 seeding cases
### reduced onwards transmission

offspring <- "pois"
mn_offspring <- 6
generation_time <- function(n) { rgamma(n, shape = 20, rate = 2) }
t0 <- 0
tf <- Inf
pop <- 10^7
check_final_size <- 5000
initial_immune <- 0
infection_to_onset <- function(n) { rgamma(n, shape = 5, rate = 2) }
vaccine_start <- 5
vaccine_coverage <- 0.85
vaccine_efficacy_infection <- 0.75
vaccine_logistical_delay <- 1
vaccine_protection_delay <- 1
seeding_cases <- 5

0.75 * 0.85

chain_sim_susc_ring_vacc <- function (offspring = c("pois", "nbinom"), mn_offspring, disp_offspring, 
                                      generation_time, t0 = 0, tf = Inf, pop, check_final_size, initial_immune = 0,
                                      seeding_cases, infection_to_onset,                                             ## epi properties relevant to vaccination
                                      vaccine_start, vaccine_coverage, vaccine_efficacy_infection,                   ## vaccine properties i
                                      vaccine_logistical_delay, vaccine_protection_delay) {                          ## vaccine properties ii
  
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
  tdf <- data.frame(id = seq_along(1:seeding_cases), ancestor = NA_integer_, generation = 1L, time = t0 + seq(from = 0, to = 0.01, length.out = seeding_cases), 
                    offspring_generated = FALSE, vaccinated = 0, time_vaccinated = NA)
  
  ## Setting up the simulation initial conditions (initial time, initial size of susceptible population etc)
  susc <- pop - initial_immune - 1L
  t <- t0
  
  ## Running the branching process - iterating over all infected people and drawing i) the number they infect (from offspring dist.);
  ## and ii) the timings of these infections (from generation time dist.)
  
  ### whilst there are still susceptible people and the latest time in the timeframe is < tdf OR number of infected people is less than check_final_size
  while (any(tdf$time[!tdf$offspring_generated] <= tf) & nrow(tdf) <= check_final_size & susc > 0) {
    
    ## Find the earliest infection in the transmission tree whose infections have not yet been calculated
    t <- min(tdf$time[!tdf$offspring_generated])
    idx <- which(tdf$time == t & !tdf$offspring_generated)[1] # does this ever get stuck with multiple people this is true for (e.g. infected @ exact same time?) - don't think so in practice
    id_parent <- tdf$id[idx]
    t_parent <- tdf$time[idx]
    gen_parent <- tdf$generation[idx]
    current_max_id <- max(tdf$id)
    index_vaccinated <- tdf$vaccinated[idx]
    
    ## Generating number of people infected by this infection
    n_offspring <- offspring_fun(1, susc)
    if (n_offspring %% 1 > 0) { # Checking offspring function is acting correctly (returning integers)
      stop("Offspring distribution must return integers")
    }
    tdf$offspring_generated[idx] <- TRUE
    
    ## Assigning times of these infections (if the index infection does indeed infect anyone)
    if (n_offspring > 0) {
      new_times <- generation_time(n_offspring)
      if (any(new_times < 0)) {
        stop("Generation times must be >= 0.")
      }
      
      ## Modelling ring vaccination
      
      ## Is the vaccine available?
      if (t >= vaccine_start) {

        are_they_vaccinated <- vector(mode = "integer", length = n_offspring) ## vector of whether secondary infections get vaccinated
        infection_retained <- vector(mode = "integer", length = n_offspring)  ## vector of whether secondary infections get retained (i.e. not prevented by ring vaccination)
        infection_retained[1:length(infection_retained)] <- 1                 ## default to infections being retained; and then flow through below to see if they get removed
        
        ## Is the index case being considered vaccinated?

        ## Calculating time to vaccination protection for index case
        onset_time <- infection_to_onset(n = 1)  ## time to onset for index case
        time_to_secondary_vaccination <- onset_time + vaccine_logistical_delay                               ## Time between index case infected and secondary cases ring vaccinated
        time_to_secondary_vaccination_protection <- time_to_secondary_vaccination + vaccine_protection_delay ## Time between index case infected and secondary cases protected with vaccination
        
        ## Checking whether vaccination occurs before or after the secondary cases have been generated;
        ## and if so, whether or not the infection is prevented by the vaccine
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
        
        new_n_offspring <- sum(infection_retained)
        new_new_times <- new_times[which(infection_retained == 1)]
        time_vaccinated <- ifelse(are_they_vaccinated == 1, time_to_secondary_vaccination, NA)
        
        ## Checking whether there are any infections to add to the table
        if (new_n_offspring != 0) {
          new_df <- data.frame(id = current_max_id + seq_len(new_n_offspring), 
                               ancestor = id_parent,
                               generation = gen_parent + 1L, 
                               time = new_new_times + t_parent, 
                               offspring_generated = FALSE,
                               vaccinated = are_they_vaccinated[which(infection_retained == 1)],
                               time_vaccinated = time_vaccinated[which(infection_retained == 1)] + t_parent)
          tdf <- rbind(tdf, new_df)
        }
        
      } else {
        
        new_df <- data.frame(id = current_max_id + seq_len(n_offspring), 
                             ancestor = id_parent,
                             generation = gen_parent + 1L, 
                             time = new_times + t_parent, 
                             offspring_generated = FALSE,
                             vaccinated = 0,
                             time_vaccinated = NA)
        tdf <- rbind(tdf, new_df)
      }
    }
    susc <- susc - n_offspring
  }
  
  ## Bookkeeping for end of function - ordering the dataframe by time of infection, id etc
  tdf <- tdf[tdf$time <= tf, ]
  tdf <- tdf[order(tdf$time, tdf$id), ]
  tdf$offspring_generated <- NULL
  
  return(tdf)
}

x <- tdf %>%
  group_by(ancestor) %>%
  summarise(n =n ())

par(mfrow = c(2, 1))
hist(x$n)
hist(offspring_fun(n = 5000, susc = 10^7))

mean(x$n)
mean(offspring_fun(n = 5000, susc = 10^7))



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
