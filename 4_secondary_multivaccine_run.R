## Runs scenarios for secondary country 
## Vary introduction times - 20 days, 50 days or 80 days after spillover in emergence country
## Start sarbecovirus vaccine 30 days after spillover event (i.e. assuming early detection in emergence country)
## Detect circulating virus 30 days after importation - put in place NPIs
## Introduce SARS-X vaccine 130 days after spillover in emergence country
## Lift restrictions when 80% coverage of 60+ population achieved or after 30 days for counterfactual

## Load VSL estimates ####
vsl <- read.csv("data/vsly.csv")
median_hospital_days <-  7.5


### Load functions #############################################################
source("R/functions_multivaccine_secondary.R")

### Specify parameters ###############################################################

target_pop <- 1e6
income_group <-  c("HIC","UMIC","LMIC","LIC")
R0 <-  c(1.5,2.5,4) 
Rt1 <- 1.0 # before actual introduction
Rt1a <- R0 # from introduction
Rt1b <- 1.1 # NPIs introduced
Rt2 <- R0   # Rt at lifting
timing1 <- 1  # keep low level of infection until actual importation
timing1a <- c(20, 50, 80, 110, 140, 170) # this is the time of importation
timing1b <- c(50, 80, 110, 140, 170, 200) # this is the time that lockdown - 30 days after importation
timing2 <-  c(80, 110, 140, 170, 200, 230)   # lifting of NPIs - default of 30 day lockdown for no vaccine, code calculates timing 2 when elderly population vaccinatedifr_scaling <- 1.0 #c(0.1,0.5,1,2) # IFR in UK population with uniform attack rate - can take values 0.1, 0.5, 2
ifr_scaling <- c(0.1,0.5,1,2) # IFR in UK population with uniform attack rate - can take values 0.1, 0.5, 2
coverage <- 0.8      # proportion of the population vaccinated
efficacy_infection_v1 <- 0.35 
efficacy_disease_v1 <-  0.8 
efficacy_infection_v2 <- 0.55
efficacy_disease_v2 <-  0.9 


## SAFIR: 5% per week in HIC/UMIC, 2% per week in LMIC/LIC
## translates to 16 weeks (112 days) for HIC/UMIC, 40 weeks (280 days) for LMIC/LIC for 80% coverage
vaccination_rate <- 0.02  # vaccination rate per week as percentage of population

duration_R <- 5*365 # duration of infection-induced immunity
duration_V <- 5*365 # duration of vaccine-induced immunity
dur_vacc_delay <- 14 # mean duration from vaccination to protection
seeding_cases <- 1 # define as the number of cases at first identification 
runtime <- 500
vaccine_1_start <- c(30, runtime) #c(30, 60, 90, runtime)
## assume that sequencing occurs 7 days after detection, and SARS-X vaccine available 100 days after sequencing
vaccine_2_start <- c(137, runtime) #c(130, 160, 190, runtime) # note that this needs to be after v1 has been completed - at 2% per week this takes 50-60 days in HIC
lower_priority <- 14  #60+
lower_vaccine <- 4 #15+ 
two_vaccines <- c(0,1)


#################### Scenario table
scenarios <- expand_grid(target_pop = target_pop,
                         income_group = income_group,
                         R0 = R0,
                         Rt1 = Rt1,
                         Rt1a = Rt1a,
                         Rt1b = Rt1b,
                         Rt2 = Rt2,
                         timing1 = timing1,
                         timing1a = timing1a,
                         timing1b = timing1b,
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
  filter(Rt2 == R0) %>%
  filter(
    ((two_vaccines==1) &  (vaccine_2_start == 107 + vaccine_1_start)) |
      ((two_vaccines==0) & (vaccine_1_start == runtime)) |
      ((vaccine_1_start == runtime) & (vaccine_2_start == runtime))
  ) %>%
  filter((timing1b == timing1a + 30) & (timing2 == timing1b + 30)) %>%
  filter(Rt1a == R0)
 

nrow(scenarios)

#### Run the model #############################################################
plan(multicore, workers = 2)
system.time({out <- future_pmap(scenarios, run_scenario_2, .progress = TRUE)})

#### Format output #############################################################
scenarios <- scenarios %>%
  select (- timing2)
out_combine <- format_out(out,scenarios)

### Save output ################################################################
saveRDS(out_combine, "output/4_secondary_mv_scenarios.rds")
