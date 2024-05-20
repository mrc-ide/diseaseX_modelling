# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/branching_process_spatial_vaccination.R"))

### SC1 parameters
SC1_generation_time <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
SC1_infection_to_onset <- function(n) { rgamma(n, shape = 0.1, rate = 1) } ## Ask Azra for values (negligible assumed currently)
SC1_prop_asymptomatic <- 0
SC1_prob_hosp <- 0.9 # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
SC1_hospitalisation_delay <- function(n) { 12 } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/

### Other parameters
pop <- 10^10
check_final_size <- 20000
initial_immune <- 0
seeding_cases <- 3
R0_scan <- c(1)
iterations <- 200

SC1_no_vaccination <- data.frame(expand_grid(R0_scan, iterations = 1:iterations), outbreak_size = NA_real_, vaccine = "no_vaccine")
counter <- 1
for (i in 1:length(R0_scan)) {
  for (j in 1:iterations) {
    
    # SARS-CoV-1 Pathogen Archetype
    SC1_temp <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                     generation_time = SC1_generation_time,
                                     spatial_kernel = spatial_kernel,
                                     t0 = 0, tf = Inf,
                                     check_final_size = check_final_size,
                                     seeding_cases = seeding_cases,
                                     prop_asymptomatic = SC1_prop_asymptomatic,
                                     prob_hosp = SC1_prob_hosp,
                                     hospitalisation_delay = SC1_hospitalisation_delay,
                                     detection_threshold = 1000,
                                     vaccine_campaign_radius = 1,
                                     vaccine_coverage = 0,
                                     vaccine_efficacy_infection = 0,
                                     vaccine_efficacy_transmission = 0,
                                     vaccine_efficacy_disease = 0,
                                     vaccine_logistical_delay = 100,
                                     vaccine_protection_delay = 100)
    SC1_count <- sum(!is.na(SC1_temp$time_infection))
    SC1_no_vaccination[counter, "outbreak_size"] <- SC1_count
    counter <- counter + 1
    print(j)
  }
  print(i)
}

SC1_no_vaccination <- SC1_no_vaccination %>%
  select(iterations, R0_scan, everything()) %>%
  rename(iteration = iterations, R0 = R0_scan)
SC1_no_vaccination$vaccine <- "no_vaccine"

overall <- SC1_no_vaccination %>%
  mutate(contained = ifelse(outbreak_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(R0, vaccine) %>%
  summarise(proportion_contained = sum(contained) / iterations)
