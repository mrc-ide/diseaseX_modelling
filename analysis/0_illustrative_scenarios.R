# Notes:
## 1) Double check IFR calc (fine to just alter prob_hosp?)
## 2) Are we fine assuming unlimited healthcare capacity?

# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("CEPI_DiseaseX_Modelling/R/functions_multivaccine.R"))

### Defining Central Scenarios

#### Demographic Parameters
target_pop <- 1e6                                                       
raw_pop <- squire::get_population(country = "Argentina")$n                   
standard_pop <- round(raw_pop * target_pop / sum(raw_pop))               
mm <- squire::get_mixing_matrix(country = "Argentina")            

#### Healthcare Parameters
hosp_bed_capacity <- 100000000                                         
ICU_bed_capacity <- 100000000                                         

### Vaccination Parameters
efficacy_infection_v1 <- 0.35              # vaccine efficacy against infection - BNPC
efficacy_disease_v1 <-  0.8                # vaccine efficacy against disease - BNPC
efficacy_infection_v2 <- 0.55              # vaccine efficacy against infection - specific vaccine
efficacy_disease_v2 <-  0.9                # vaccine efficacy against disease - specific vaccine
duration_R <- 1000 * 365                   # duration of infection-induced immunity - assumed large as not considering waning (currently)
duration_V <- 1000 * 365                   # duration of vaccine-induced immunity for both vaccines - assumed large as not considering waning (currently) - also need to make vaccine specific
dur_vacc_delay <- 0.5                      # mean duration from vaccination to protection
vaccine_1_start <- 50                      # time when BNPC starts being distributed
vaccine_2_start <- 100                     # note that this needs to be after BPSV i.e. v1 for 60+ has been completed to avoid having double vax rate (i.e. v1 BPSV and v2 spec for 60+ happening at once)
coverage <- 0.8                            # proportion of the population vaccinated
vaccination_rate <- 0.04                   # vaccination rate per week as percentage of population
v1_daily_doses <- vaccination_rate * sum(standard_pop) / 7
v2_daily_doses <- vaccination_rate * sum(standard_pop) / 7  
lower_priority <- 13                       # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
lower_vaccine <- 4                         # index of the youngest age group that *receives* vaccines (4 = 15+)
priority_age_groups <- lower_priority:17   # creating the index list for priority age groups  
vaccination_age_groups <- lower_vaccine:17 # creating the index list for vaccinable age-groups
vaccine_coverage_mat <-                    # matrix (col = age group, row = prioritisation step) of how to assign vaccines to different age-groups according to priority
  matrix(c(rep(0, 17 - length(priority_age_groups)), 
           rep(coverage, length(priority_age_groups)), 
           rep(0, 17 - length(vaccination_age_groups)), 
           rep(coverage, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)
second_dose_delay <- 0.5                     # controls how many days after "1st dose" people receive second dose; see here: https://github.com/mrc-ide/squire.page/blob/main/inst/odin/nimue_booster.R#L427-L430

elderly_pop_to_vaccinate <- sum(standard_pop[priority_age_groups]) * coverage # 60+s receive primary (BNPCV) and booster (diseaseX-specific); under 60s receive just primary (diseaseX-specific)
time_to_coverage_v1 <- ceiling(elderly_pop_to_vaccinate / v1_daily_doses)     # calculated so that we only vaccinate 60+s with the stockpiled BNPSC
minimum_spec_development_time_allowed <- (time_to_coverage_v1 + second_dose_delay + dur_vacc_delay)
if (minimum_spec_development_time_allowed > (vaccine_2_start - vaccine_1_start)) {
  stop("Virus-specific vaccine developed too soon given speed of vaccination campaign")
}
time_to_coverage_v2_elderly <- ceiling(elderly_pop_to_vaccinate/v2_daily_doses)
primary_doses <- c(0, v1_daily_doses, 0, v2_daily_doses) # final quantity here is doses after final quantity in tt_primary_doses
tt_primary_doses <- c(0, vaccine_1_start, round(vaccine_1_start + time_to_coverage_v1), vaccine_2_start + time_to_coverage_v2_elderly) # consider what it means when this stops too soon 

booster_doses <- c(0, v2_daily_doses, 0) 
tt_booster_doses <- c(0, vaccine_2_start, vaccine_2_start + time_to_coverage_v2_elderly) # note that if this stops too soon, you'll end up with elderly age-groups not being fully boosted - WHY IS THIS???


#setup as a matrix giving VE per age group
ve_v1 <- list(infection = efficacy_infection_v1, disease = efficacy_disease_v1)
ve_v2 <- list(infection = efficacy_infection_v2, disease = efficacy_disease_v2)
ve_i_booster_elderly <- c(ve_v1$infection, ve_v1$infection, ve_v1$infection, 0, ve_v2$infection, ve_v2$infection, 0) 
ve_i_booster_non_elderly <- c(ve_v2$infection, ve_v2$infection, ve_v2$infection, 0, ve_v2$infection, ve_v2$infection, 0)
vaccine_efficacy_infection <- matrix(c( rep(ve_i_booster_non_elderly, 17 - length(priority_age_groups)), 
                                        rep(ve_i_booster_elderly, length(priority_age_groups)) ), ncol = 17)

#repeat for disease 
ve_i_booster_elderly <- c(ve_v1$disease, ve_v1$disease, ve_v1$disease, 0, ve_v2$disease, ve_v2$disease, 0) 
ve_i_booster_non_elderly <- c(ve_v2$disease, ve_v2$disease, ve_v2$disease, 0, ve_v2$disease, ve_v2$disease, 0)
vaccine_efficacy_disease <- matrix(c( rep(ve_i_booster_non_elderly, 17 - length(priority_age_groups)), 
                                      rep(ve_i_booster_elderly, length(priority_age_groups)) ), ncol = 17)

vaccine_booster_follow_up_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(1, length(priority_age_groups)))
vaccine_booster_initial_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(1, length(priority_age_groups)))


# Central things to vary each time we vary another quantity
specific_vaccine_availability <- c(100, 300)
R0 <-  c(1.5, 3) 
Tg <- c(7.5, 15)
current_Tg <- squire.page:::durs_booster$dur_IMild + squire.page:::durs_booster$dur_ICase
Tg_ratio <- Tg / current_Tg
TgVary_dur_IMild <- Tg_ratio * squire.page:::durs_booster$dur_IMild
TgVary_dur_ICase <- Tg_ratio * squire.page:::durs_booster$dur_ICase

raw_pop <- squire::get_population(country = "Argentina")$n  # getting representative country (needs to change from Argentina to median UMIC)                 
standard_pop <- round(raw_pop * 10^7 / sum(raw_pop))        # standardised population size
mm <- squire::get_mixing_matrix(country = "Argentina")      # mixing matrix for that country
contact_rates <- apply(mm, 1, sum)
contact_rates <- c(contact_rates, "80+" = unname(contact_rates[16]))
pop_contact_weighting <- standard_pop * contact_rates
non_severe_deaths <- nimue:::probs$prob_hosp * (1 - nimue:::probs$prob_severe) * nimue:::probs$prob_non_severe_death_treatment
severe_deaths <- nimue:::probs$prob_hosp * nimue:::probs$prob_severe * nimue:::probs$prob_severe_death_treatment
raw_IFR <- 100 * sum(((non_severe_deaths + severe_deaths) * pop_contact_weighting/sum(pop_contact_weighting))) 
target_IFR <- c(0.75, 2)
IFR_scaling_factor <- target_IFR / raw_IFR
Rt_change_unmitigated <- 1                                  # multiplicative factor for R0 in unmitigated scenario
Rt_change_non_mandated <- 0.7                               # multiplicative factor for R0 in non-mandated restrictions scenario                            
Rt_change_stringent <- 0.4                                  # multiplicative factor for R0 in stringent restrictions scenario    
Rt_change <- c(Rt_change_unmitigated, Rt_change_non_mandated, Rt_change_stringent) 
##### Note that the above might need to be changed depending on which R0s we use

### consider potentially just implementing a bunch of restrictions to bring R << 1 in all cases except for
### a specific NPIs figure where we look at this in more detail

# don't forget to run each of these with and without broad vaccine
baseline_scenarios <- expand.grid(R0 = R0,
                                  Rt_change = Rt_change, 
                                  IFR_Scaling_Factor = IFR_scaling_factor,
                                  Tg_Scaling_Factor = Tg_ratio,
                                  specific_vaccine_availability = specific_vaccine_availability,
                                  broad_vaccine_availability = c("Yes", "No"))
baseline_scenarios$coverage <- 0.8                # proportion of the population vaccinated
baseline_scenarios$efficacy_infection_v1 <- 0.8   # proportion of the population vaccinated
baseline_scenarios$efficacy_disease_v1 <-  0.8    # broad vaccine disease efficacy
baseline_scenarios$efficacy_infection_v2 <- 0.55  # specific vaccine infection efficacy
baseline_scenarios$efficacy_disease_v2 <-  0.9    # specific vaccine disease efficacy
baseline_scenarios$hosp_bed_capacity <- 100000000 # hospital bed capacity
baseline_scenarios$ICU_bed_capacity <- 100000000  # ICU bed capacity
baseline_scenarios$vaccination_rate <- 0.04       # vaccination rate per week as percentage of population, 4% per week (20 weeks, 140 days) 

runtime <- 600
tictoc::tic()
r1 <- run_booster( 
  time_period = runtime,                                                     # time to run the model for
  population = standard_pop,                                                 # population to be simulated
  contact_matrix_set = mm,                                                   # mixing matrix
  R0 = R0,                                                                   # basic reproduction number
  tt_R0 = 0,                                                                 # timing of changes in R0 to mimic NPIs
  hosp_bed_capacity = hosp_bed_capacity,                                     # hospital bed capacity
  ICU_bed_capacity = ICU_bed_capacity,                                       # ICU bed capacity
  prob_hosp = prob_hosp,                                                     # probability of requiring hospitalisation
  prob_severe = prob_severe,                                                 # probability of requiring ICU stay
  prob_non_severe_death_no_treatment = prob_non_severe_death_treatment,      # prob death given non-severe disease
  prob_severe_death_treatment = prob_severe_death_treatment,                 # prob death given severe disease
  dur_R = duration_R,                                                        # duration of immunity following infection
  seeding_cases = seeding_cases,                                             # number of cases outbreak starts with
  vaccine_coverage_mat = vaccine_coverage_mat,                               # matrix (col = age group, row = prioritisation step) of how to assign vaccines to different age-groups according to priority
  vaccine_efficacy_infection = vaccine_efficacy_infection,                   # vaccine efficacy against infection (matrix where row = age group, column = vaccine/dose)
  vaccine_efficacy_disease = vaccine_efficacy_disease,                       # vaccine efficacy against disease (matrix where row = age group, column = vaccine/dose)
                                                                             # --> columns = (first dose, second dose, second dose, waned, booster, booster, waned) 
                                                                             # --> note in this hacky version, for age > priority group, first vaccine = BNPSC, second = diseaseX; for age < priority group, only have first vaccine which is diseaseX 
  second_dose_delay = second_dose_delay,                                     # delay between receipt of first and second dose (ignored here as we're modelling single dose vaccines) - second doses are assigned
                                                                             # automatically in squire.page, at second_dose_delay days after the primary dose was delivered
  dur_V = rep(duration_V/2, 4),                                              # duration of vaccination-derived immunity (currently not vaccine specific) 
                                                                             # --> first two elements are for primary vaccine, elements 3 and 4 are for booster, half duration in each compartment giving erlang-2 waning
                                                                             # --> gets converted into gamma_vaccine by https://github.com/mrc-ide/squire.page/blob/073ec7d107419d95aabfaba3735ffdc2609195a4/R/booster_model.R#L575-L584
  protection_delay_rate = 1/dur_vacc_delay,                                  # Delay between full initial vaccine course and protection (see run_booster -> parameters_booster -> apply_dose_delay_booster in squire.page)
                                                                             # Combines with "second_dose_delay" above to generate the amount of time between first dose receipt and second dose receipt (under assumption)
                                                                             # that individuals are protected immediately after receiving second dose. I.e. way of aggregating both delay between doses and delay between dose and protection. 
  primary_doses = primary_doses,                                             # Primary dose (i.e. first dose, with second dose then implicitly generated under the hood) vaccination rate 
  tt_primary_doses = tt_primary_doses,                                       # Timings for variable primary dose vaccination rate 
  booster_doses = booster_doses,                                             # Booster dose vaccination rate 
  tt_booster_doses = tt_booster_doses,                                       # Timings for variable booster dose vaccination rate 
  vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,   # Vector describing which age-groups are eligible for follow-up boosters (re-boostering, basically)
  vaccine_booster_initial_coverage = vaccine_booster_initial_coverage)       # Vector describing which age-groups are eligible for initial boosters
tictoc::toc()

