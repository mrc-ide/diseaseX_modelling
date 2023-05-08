### Load functions #############################################################
source("R/functions_multivaccine.R")

### Specify parameters ###############################################################

target_pop <- 1e6
income_group <- c("HIC","UMIC","LMIC","LIC")
R0 <- c(1.5,2.5,4) 
Rt1 <- R0    # Rt at first NPI introduction
Rt2 <- R0   # Rt at lifting
timing1 <- 30      # introduction of NPIs
timing2 <- 365     # lifting of NPIs
ifr_scaling <-  c(0.1,0.5,1,2) # IFR in UK population with uniform attack rate - can take values 0.1, 0.5, 2
coverage <- 0.8      # proportion of the population vaccinated
efficacy_infection_v1 <- 0.5 
efficacy_disease_v1 <-  0.8 
efficacy_infection_v2 <- 0.8
efficacy_disease_v2 <-  0.9 

## SAFIR: 5% per week in HIC/UMIC, 2% per week in LMIC/LIC
## translates to 16 weeks (112 days) for HIC/UMIC, 40 weeks (280 days) for LMIC/LIC for 80% coverage
vaccination_rate <- 0.02  # vaccination rate per week as percentage of population
#vaccine_period <- 112 # days over which to vaccine target coverage
duration_R <- 365 # duration of infection-induced immunity
duration_V <- 365 # duration of vaccine-induced immunity
dur_vacc_delay <- 1 # mean duration from vaccination to protection
seeding_cases <- 100 # define as the number of cases at first sequencing - will need to explore
#vaccine_start <- c(10,20,50,100,365*2) #c(10,20,50,100,365) # days after start of the epidemic - will depend on seeding_cases
runtime <- 365*2
vaccine_1_start <- c(10, 20, 30, runtime)
vaccine_2_start <- c(100, runtime) # note that this needs to be after v1 has been completed - at 2% per week this takes 50-60 days in HIC
lower_priority <- 14  #65+
lower_vaccine <- 4 #15+ 
two_vaccines <- c(0,1)



#################### Scenario table
scenarios <- expand_grid(target_pop = target_pop,
                         income_group = income_group,
                         R0 = R0,
                         Rt1 = Rt1,
                         Rt2 = Rt2,
                         timing1 = timing1,
                         timing2 = timing2,
                         ifr_scaling = ifr_scaling,
                         coverage = coverage,
                         efficacy_infection_v1 = efficacy_infection_v1,
                         efficacy_disease_v1 = efficacy_disease_v1,
                         efficacy_infection_v2 = efficacy_infection_v2,
                         efficacy_disease_v2 = efficacy_disease_v2,
                         vaccination_rate = vaccination_rate,
                         duration_R = duration_R,
                         duration_V = duration_V,
                         dur_vacc_delay = dur_vacc_delay,
                         seeding_cases = seeding_cases,
                         vaccine_1_start = vaccine_1_start,
                         vaccine_2_start = vaccine_2_start,
                         lower_priority = lower_priority,
                         lower_vaccine = lower_vaccine,
                         two_vaccines = two_vaccines,
                        runtime = runtime) %>%
#  filter((efficacy_infection ==0.6 & efficacy_disease==0.8 & (vaccine_start <100 | vaccine_start == 365*2)) | (efficacy_infection==0.8 & efficacy_disease==0.9 & vaccine_start >=100)) %>%
  filter(Rt1 == R0 & Rt2 == R0) %>%
  filter((vaccine_1_start <= 50 & vaccine_2_start == 100 & two_vaccines == 1) | (vaccine_1_start ==runtime & vaccine_2_start == 100 & two_vaccines == 0) 
                | (vaccine_1_start == runtime & vaccine_2_start == runtime ))

nrow(scenarios)

#### Run the model #############################################################
plan(multicore, workers = 4)
system.time({out <- future_pmap(scenarios, run_scenario, .progress = TRUE)})

#### Format output #############################################################
out_combine <- format_out(out,scenarios)

### Save output ################################################################
saveRDS(out_combine, "output/1_unmitigated_mv_scenarios.rds")
