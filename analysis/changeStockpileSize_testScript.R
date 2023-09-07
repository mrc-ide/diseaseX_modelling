# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))

# Parameters for model running
Rt <- 2
tt_Rt <- 0
IFR <- 1
country <- "Argentina"
population_size <- 10^10
hosp_bed_capacity <- 10^10                                         
ICU_bed_capacity <- 10^10   
Tg <- 5.5
detection_time <- 10
bpsv_start <- 1
bpsv_protection_delay <- 1
specific_vaccine_start <- 200
specific_protection_delay <- 7                 
efficacy_infection_bpsv <- 0.35                
efficacy_disease_bpsv <- 0.5
efficacy_infection_spec <- 0.55             
efficacy_disease_spec <- 0.9                   
dur_R <- 365000000                             
dur_bpsv <- 365000000                          
dur_spec <- 365000000
vaccination_rate <- 0.035
min_age_group_index_priority <- 13
min_age_group_index_non_priority <- 4
runtime <- 1000
seeding_cases <- 50

## Adapting Generation Time
TgVary <- scale_generation_time(target_Tg = Tg)
TgVary_dur_IMild <- TgVary$dur_IMild
TgVary_dur_ICase <- TgVary$dur_ICase

# Adjust prob_hosp
prob_hosp <- scale_IFR(country = country, population_size = population_size, target_IFR = IFR)

# Setting up vaccination stuff

## Setting vaccination rates and coverages
bpsv_stockpile_coverage <- 0.4
coverage <- 0.9
priority_age_groups <- min_age_group_index_priority:17  
vaccination_age_groups <- min_age_group_index_non_priority:17 
vaccine_coverage_mat <- matrix(c(rep(0, 17 - length(priority_age_groups)), rep(coverage, length(priority_age_groups)), 
                                 rep(0, 17 - length(vaccination_age_groups)),  rep(coverage, length(vaccination_age_groups))), ncol = 17, byrow = TRUE)

## Generating vaccine dosing schedules
standard_pop <- generate_standard_pop(country = country, population_size = population_size)
daily_doses <- vaccination_rate * population_size / 7
priority_age_groups <- min_age_group_index_priority:17

### New structure is that we compress the duration of the bpsv campaign so that we only achieve bpsv_stockpile_coverage - overall coverage (applied to specific vaccine) is still applied to vaccine_coverage_matrix
###   we then change the vaccine_efficacy matrix just before specific campaign starts, primary dose efficacy for elderly goes from bpsv -> specific
###   we also turn off second doses so that everyone who gets vaccinated with specific primary doesn't get a second dose
###   we then start the specific vaccination campaign and
###      1) boost all the elderly that got vaccinated with bpsv
###      2) give primary doses (now specific vaccine efficacy profile) to all elderly who *didn't* get the bpsv
###      3) give primary doses to everyone else
###   I think 2) and 3) should be handled automatically within the optimisation steps, but let's see

elderly_pop_to_vaccinate_bpsv <- ceiling(sum(standard_pop[priority_age_groups]) * bpsv_stockpile_coverage)
time_to_coverage_bpsv <- ceiling(elderly_pop_to_vaccinate_bpsv/daily_doses)

elderly_pop_to_vaccinate_spec <- ceiling(sum(standard_pop[priority_age_groups]) * coverage)
time_to_coverage_spec <- ceiling(elderly_pop_to_vaccinate_spec/daily_doses)

primary_doses <- 
  c(rep(0, detection_time),                     ## time between epidemic start and detection
    rep(0, bpsv_start),                         ## time between detection and initiation of BPSV campaign
    rep(0, bpsv_protection_delay),              ## time between initiation of BPSV campaign and people first being protected by that first dose
    rep(daily_doses, time_to_coverage_bpsv),    ## protection (if any) emerges in BPSV-vaccinated primary vaccinated folks
    rep(0, specific_vaccine_start - time_to_coverage_bpsv - bpsv_protection_delay - bpsv_start), # specific vaccine campaign starts specific_vaccine_start days after detection
    rep(0, time_to_coverage_spec),              ## no specific vaccine for non-elderly whilst that vaccination campaign is ongoing
    rep(0, specific_protection_delay),          ## time between initiation of specific vaccine campaign and people first being protected by that first dose
    rep(daily_doses, runtime - specific_protection_delay - time_to_coverage_spec - specific_vaccine_start - detection_time)) # specific vaccination of all other ages until end of runtime

## in this new formulation, we don't give any secondary doses in the disease-specific vaccination campaign
## this is a little bit of a fudge - required because when specific vaccination campaign starts, we need to change the vaccine efficacy profile so that 
## primary vaccination for elderly becomes disease specific profile - required because we need to give specific vaccine to all elderly eligible,
## not just those meeting the size of the bpsv stockpile (which is a coverage < than target disease specific vaccination coverage)
second_doses <- 
  c(rep(0, 2),                                  ## arbitrary two day delay between primary and secondary doses
    rep(0, detection_time),                     ## time between epidemic start and detection
    rep(0, bpsv_start),                         ## time between detection and initiation of BPSV campaign
    rep(0, bpsv_protection_delay),              ## time between initiation of BPSV campaign and people first being protected by that first dose
    rep(daily_doses, time_to_coverage_bpsv),    ## protection (if any) emerges in BPSV-vaccinated primary vaccinated folks
    rep(0, specific_vaccine_start - time_to_coverage_bpsv - bpsv_protection_delay - bpsv_start), # specific vaccine campaign starts specific_vaccine_start days after detection
    rep(0, time_to_coverage_spec),              ## no specific vaccine for non-elderly whilst that vaccination campaign is ongoing
    rep(0, specific_protection_delay),          ## time between initiation of specific vaccine campaign and people first being protected by that first dose
    rep(0, runtime - specific_protection_delay - time_to_coverage_spec - specific_vaccine_start - detection_time)) ## no secondary doses in disease-specific regime (primary doses are treated as the complete course)
second_doses <- second_doses[1:length(primary_doses)]
# plot(primary_doses, type = "l")
# lines(second_doses, type = "l", col = "red")

## Booster Doses (only for elderly, this is the disease specific vaccine)
booster_doses <- 
  c(rep(0, detection_time),
    rep(0, specific_vaccine_start),
    rep(0, specific_protection_delay),
    rep(daily_doses, runtime - specific_protection_delay - specific_vaccine_start - detection_time))

## Creating matrices of infection and disease efficacy for each age-group and vaccine received

## Major change compared to previous iterations: 
## During period of bpsv campaign, vaccine efficacy matrix has primary, secondary = bpsv and booster = specific for elderly population
## During period of specific campaign, matrix has primary = specific, secondary = bpsv and booster = specific for elderly population
## Required as some eligible elderly won't have received bpsv during initial campaign due to limited size of bpsv stockpile (coverage bpsv < coverage specific)
## No second doses are distributed during the specific campaign. 

## Ordering in efficacy vector
## 1 = first dose received
## 2 = second dose received 1; 3 = second dose received 2
## 4 = primary/secondary waned
## 5 = booster 1; 6 = booster 2
## 7 = booster waned
ve_bpsv <- list(infection = efficacy_infection_bpsv, disease = efficacy_disease_bpsv)
ve_spec <- list(infection = efficacy_infection_spec, disease = efficacy_disease_spec)

### Infection Efficacy
ve_i_elderly_bpsv_campaign <- c(ve_bpsv$infection, ve_bpsv$infection, ve_bpsv$infection, 0, ve_spec$infection, ve_spec$infection, 0) 
ve_i_elderly_spec_campaign <- c(ve_spec$infection, ve_bpsv$infection, ve_bpsv$infection, 0, ve_spec$infection, ve_spec$infection, 0) 
ve_i_non_elderly <- c(ve_spec$infection, ve_spec$infection, ve_spec$infection, 0, ve_spec$infection, ve_spec$infection, 0)

vaccine_efficacy_infection_bpsv_campaign <- matrix(c(rep(ve_i_non_elderly, 17 - length(priority_age_groups)), 
                                                     rep(ve_i_elderly_bpsv_campaign, length(priority_age_groups)) ), ncol = 17)
vaccine_efficacy_infection_spec_campaign <- matrix(c(rep(ve_i_non_elderly, 17 - length(priority_age_groups)), 
                                                     rep(ve_i_elderly_spec_campaign, length(priority_age_groups)) ), ncol = 17)
### Disease Efficacy
ve_d_elderly_bpsv_campaign <- c(ve_bpsv$disease, ve_bpsv$disease, ve_bpsv$disease, 0, ve_spec$disease, ve_spec$disease, 0) 
ve_d_elderly_spec_campaign <- c(ve_spec$disease, ve_bpsv$disease, ve_bpsv$disease, 0, ve_spec$disease, ve_spec$disease, 0) 
ve_d_non_elderly <- c(ve_spec$disease, ve_spec$disease, ve_spec$disease, 0, ve_spec$disease, ve_spec$disease, 0)

vaccine_efficacy_disease_bpsv_campaign <- matrix(c(rep(ve_d_non_elderly, 17 - length(priority_age_groups)), 
                                                   rep(ve_d_elderly_bpsv_campaign, length(priority_age_groups)) ), ncol = 17)
vaccine_efficacy_disease_spec_campaign <- matrix(c(rep(ve_d_non_elderly, 17 - length(priority_age_groups)), 
                                                   rep(ve_d_elderly_spec_campaign, length(priority_age_groups)) ), ncol = 17)

### Vaccine protection waning
dur_V <- matrix(data = c(c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_bpsv, length(priority_age_groups))),    ## duration of primary series protection (bpsv for elderly, specific vaccine for everyone else in this scenario)
                         c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_bpsv, length(priority_age_groups))),    ## duration of primary series protection (bpsv for elderly, specific vaccine for everyone else in this scenario)
                         c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_spec, length(priority_age_groups))),           ## duration of booster protection (specific for elderly, not used for everyone else)
                         c(rep(dur_spec, min(priority_age_groups) - 1), rep(dur_spec, length(priority_age_groups)))),         ## duration of booster protection (specific for elderly, not used for everyone else)
                nrow = 4, ncol = 17, byrow = TRUE)

## Eligibility for Boosters (Elderly Only)
vaccine_booster_follow_up_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(0, length(priority_age_groups))) # no booster follow ons
vaccine_booster_initial_coverage <- c(rep(0, min(priority_age_groups) - 1), rep(1, length(priority_age_groups))) ## note that we set coverage to be 100% as it's already bounded by whatever coverage the 1st/2nd dose coverage got to

## Running the Model 
mm <- squire::get_mixing_matrix(country = country)    
rel_infectiousness_vaccinated <- squire.page:::probs_booster$rel_infectiousness_vaccinated
rel_infectiousness_vaccinated[rel_infectiousness_vaccinated < 1] <- 1
mod_run <- run_booster(time_period = runtime,
                       population = standard_pop,                                                 
                       contact_matrix_set = mm,                                                   
                       R0 = Rt,     
                       tt_R0 = tt_Rt, 
                       hosp_bed_capacity = hosp_bed_capacity,                                     
                       ICU_bed_capacity = ICU_bed_capacity,                                       
                       prob_hosp = prob_hosp,
                       dur_IMild = TgVary_dur_IMild,
                       dur_ICase = TgVary_dur_ICase,
                       dur_R = dur_R,                                                        
                       seeding_cases = seeding_cases,
                       vaccine_coverage_mat = vaccine_coverage_mat,
                       
                       ## Old
                       # vaccine_efficacy_infection = vaccine_efficacy_infection,                   
                       # vaccine_efficacy_disease = vaccine_efficacy_disease,
                       
                       ## New
                       vaccine_efficacy_infection = list(vaccine_efficacy_infection_bpsv_campaign, 
                                                         vaccine_efficacy_infection_spec_campaign),
                       tt_vaccine_efficacy_infection = c(0, specific_vaccine_start - 10), # calculate this ahead of time to make sure it doesn't bleed into time when bpsv is being given out
                       vaccine_efficacy_disease = list(vaccine_efficacy_disease_bpsv_campaign, 
                                                       vaccine_efficacy_disease_spec_campaign),
                       tt_vaccine_efficacy_disease = c(0, specific_vaccine_start - 10), # calculate this ahead of time to make sure it doesn't bleed into time when bpsv is being given out
                       
                       rel_infectiousness_vaccinated = rel_infectiousness_vaccinated, 
                       dur_V = dur_V,                                              
                       primary_doses = primary_doses,  
                       second_doses = second_doses,
                       booster_doses = booster_doses,                                             
                       vaccine_booster_follow_up_coverage = vaccine_booster_follow_up_coverage,
                       vaccine_booster_initial_coverage = vaccine_booster_initial_coverage)

check <- nimue::format(mod_run, compartments = c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_booster_dose"),
                       reduce_age = FALSE) %>%
  filter(t > 1, compartment == "deaths" |  
           compartment == "vaccinated_first_dose" | compartment == "vaccinated_second_dose" | compartment == "vaccinated_booster_dose") %>%
  group_by(replicate, t) 


ggplot() +
  geom_line(data = subset(check, compartment == "vaccinated_booster_dose"),
            aes(x = t, y = value, col = compartment)) +
  geom_line(data = subset(check, compartment == "vaccinated_first_dose"),
            aes(x = t, y = value, col = compartment)) +
  geom_line(data = subset(check, compartment == "vaccinated_second_dose"),
            aes(x = t, y = value, col = compartment)) +
  facet_wrap(~age_group)

# one check will be to run this with bpsv_coverage and coverage set to the same,
# and compare it to the current/(I guess with this in place, old) framework to check they produce almost identical results

