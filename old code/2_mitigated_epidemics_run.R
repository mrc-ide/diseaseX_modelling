### Load functions #############################################################
source("R/functions.R")

### Specify parameters ###############################################################

target_pop <- 1e6
income_group <- c("HIC","UMIC","LMIC","LIC")
R0 <- c(1.5,2.5,4) 
Rt1 <- 1.1    # Rt at first NPI introduction
Rt2 <- R0   # Rt at lifting
timing1 <- 10  # introduction of NPIs
timing2 <- c(112+10,112+20,112+50,112+100,365)    # lifting of NPIs - 112 days to 80% coverage 
ifr_scaling <- c(0.1,0.5,1,2) # IFR in UK population with uniform attack rate - can take values 0.1, 0.5, 2
coverage <- 0.8      # proportion of the population vaccinated
vaccine_coverage_mat <- "Elderly"   # this is the prioritisation matrix 
efficacy_infection <- c(0.6,0.8)
efficacy_disease <- c(0.8,0.9)

## SAFIR: 5% per week in HIC/UMIC, 2% per week in LMIC/LIC
## translates to 16 weeks (112 days) for HIC/UMIC, 40 weeks (280 days) for LMIC/LIC for 80% coverage
vaccine_period <- 112 # days over which to vaccine target coverage
duration_R <- 365 # duration of infection-induced immunity
duration_V <- 365 # duration of vaccine-induced immunity
dur_vacc_delay <- 21 # mean duration from vaccination to protection
seeding_cases <- 100 # define as the number of cases at first sequencing - will need to explore
vaccine_start <- c(10,20,50,100,365*2) #c(10,20,50,100) # days after start of the epidemic - will depend on seeding_cases
runtime <- 365*2


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
                         vaccine_coverage_mat = vaccine_coverage_mat,
                         efficacy_infection = efficacy_infection,
                         efficacy_disease = efficacy_disease,
                         vaccine_period = vaccine_period,
                         duration_R = duration_R,
                         duration_V = duration_V,
                         dur_vacc_delay = dur_vacc_delay,
                         seeding_cases = seeding_cases,
                         vaccine_start = vaccine_start,
                         runtime = runtime) %>%
  filter((efficacy_infection ==0.6 & efficacy_disease==0.8 & (vaccine_start <100 | vaccine_start == 2*365)) | (efficacy_infection==0.8 & efficacy_disease==0.9 & vaccine_start >=100)) %>%
#  filter(efficacy_infection ==0.6 & efficacy_disease==0.8) %>%
  filter(Rt2 == R0) %>%
  filter((timing2 == 112 + vaccine_start) | (vaccine_start==365*2 & timing2==365))
  

#scenarios <- scenarios %>%
 #  select(t_start, R0, Rt1, Rt2, timing1, timing2, coverage, mode, efficacy, hs_constraints, age_target, income_group, duration_R, duration_V, dur_vacc_delay, vaccine_period, vaccine_start, timing1, timing2, immunosenescence, seeding_cases)

nrow(scenarios)

#### Run the model #############################################################
plan(multicore, workers = 2)
system.time({out <- future_pmap(scenarios, run_scenario, .progress = TRUE)})

#### Format output #############################################################
out_combine <- format_out(out,scenarios)

### Save output ################################################################
saveRDS(out_combine, "output/2_mitigated_scenarios.rds")
