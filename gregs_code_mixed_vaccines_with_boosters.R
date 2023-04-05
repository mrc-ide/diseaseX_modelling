devtools::install_github("mrc-ide/squire.page", ref = "three-compartment-primary") #get the version with three primary series compartments instead of two 
library(squire.page) 
library(tidyverse) 
country <- "France"

#assumed scenario
#vaccine 1 stockplied and immediately rollout out to elderly populations
#vaccine 2 produced after a year with higher efficacy and rolled out first to the elderly then to the rest of the pop

duration_of_protection <- 1.5 * 365 #1.5 year duration have to assume this is identical between vaccines
ve_v1 <- list(infection = 0.5, disease = 0.5)
ve_v2 <- list(infection = 0.75, disease = 0.95)
coverage <- 0.8
time_to_coverage_v1 <- 300 #how long it will take to reach v1 coverage in the target_age_groups
start_time_v1_campaign <- 50 #starts 50 days in 
target_age_groups <- 14:17 #55+ 
start_time_v2_campaign <- 365 
v2_daily_doses <- 20000

#calculate how many doses of v1 is needed each day to reach the coverage 
elderly_pop_to_vaccinate <- sum(squire::get_population(country)$n[target_age_groups]) * coverage 
v1_daily_doses <- elderly_pop_to_vaccinate/time_to_coverage_v1 
#calculate how long it will take to boost the elderly pop (who are vaccinated) 
time_to_coverage_v2_eldery <- ceiling(elderly_pop_to_vaccinate/v2_daily_doses)

#convert these into parameters for the booster model 
vaccine_coverage_mat <- matrix(c( 
  rep(0, 17 - length(target_age_groups)), rep(coverage, length(target_age_groups)), 
  rep(0, 3), rep(coverage, 17 - 3) ), 
  ncol = 17, byrow = TRUE)

#format for ve in the booster model (first dose, second dose, second dose, waned, booster, booster , waned) 
ve_i_booster_elderly <- c(ve_v1$infection, ve_v1$infection, ve_v1$infection, 0, ve_v2$infection, ve_v2$infection, 0) 
ve_i_booster_non_elderly <- c(ve_v2$infection, ve_v2$infection, ve_v2$infection, 0, ve_v2$infection, ve_v2$infection, 0)

#setup as a matrix giving VE per age group

vaccine_efficacy_infection <- matrix(c( rep(ve_i_booster_non_elderly, 17 - length(target_age_groups)), 
                                        rep(ve_i_booster_elderly, length(target_age_groups)) ), ncol = 17)

#repeat for disease 
ve_i_booster_elderly <- c(ve_v1$disease, ve_v1$disease, ve_v1$disease, 0, ve_v2$disease, ve_v2$disease, 0) 
ve_i_booster_non_elderly <- c(ve_v2$disease, ve_v2$disease, ve_v2$disease, 0, ve_v2$disease, ve_v2$disease, 0)

vaccine_efficacy_disease <- matrix(c( rep(ve_i_booster_non_elderly, 17 - length(target_age_groups)), 
                                      rep(ve_i_booster_elderly, length(target_age_groups)) ), ncol = 17)

#since we're ignoring first doses set a one day delay between first dose and second dose (in primary series)

second_dose_delay <- 1 
dur_V <- rep(duration_of_protection/2, 4) #spend half the duration in each compartment giving erlang-2 waning
primary_doses <- c( 
  0, v1_daily_doses, 0, v2_daily_doses ) 
tt_primary_doses <- c( 
  0, start_time_v1_campaign, start_time_v1_campaign + time_to_coverage_v1, start_time_v2_campaign + time_to_coverage_v2_eldery ) 
booster_doses <- c( 
  0, v2_daily_doses, 0 ) 
tt_booster_doses <- c( 
  0, start_time_v2_campaign, start_time_v2_campaign + time_to_coverage_v2_eldery )

#other notes: 
#all other parameters should be as in nimue::run 
#the booster doesn't have compartments to simulate a delay in the development of 
#protection. this can be approximated with the parameters protection_delay_rate 
#and protection_delay_shape which approximate an gamma delay in development of protection


multi_vaccine_sim <- run_booster( 
  country = country, 
  time_period = 1500, 
  vaccine_coverage_mat = vaccine_coverage_mat, 
  vaccine_efficacy_infection = vaccine_efficacy_infection, 
  vaccine_efficacy_disease = vaccine_efficacy_disease, 
  second_dose_delay = second_dose_delay, 
  dur_V = dur_V, 
  primary_doses = primary_doses, 
  tt_primary_doses = tt_primary_doses, 
  booster_doses = booster_doses, 
  tt_booster_doses = tt_booster_doses )
  #Rest of the parameter here see ?run_booster for details )

# # This setup doesn't model into endemicity, i.e. further doses 
# 
#but this can be implemented by just adding more doses to the booster, which should apply to the entire population: 
# time_to_coverage_v2_non_elderly <- (sum(squire::get_population(country)$n[4:17]) * coverage - elderly_pop_to_vaccinate)/v2_daily_doses 
# booster_doses <- c( 
    # 0, v2_daily_doses, 0, v2_daily_doses # ) 
# tt_booster_doses <- c( 
  # 0, start_time_v2_campaign, start_time_v2_campaign + time_to_coverage_v2_eldery, start_time_v2_campaign + time_to_coverage_v2_eldery + time_to_coverage_v2_non_elderly 
# )


#plots to demonstrate vaccinations
nimue_format(multi_vaccine_sim, c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_booster_dose")) %>%
  select(compartment, y, t) %>%
  mutate(compartment = str_remove(compartment, "vaccinated_") %>%
           str_replace("_", " ") %>%
           str_to_title()) %>%
  ggplot(aes(x = t, y = y, colour = compartment)) +
  geom_line() +
  labs(x = "Time", y = "Vaccinated", colour = "Dose")

#now to convert this into our multivaccines
nimue_format(multi_vaccine_sim, c("vaccinated_second_dose", "vaccinated_booster_dose"), reduce_age = FALSE) %>%
  mutate(
    vaccine = if_else(
      compartment == "vaccinated_second_dose" & age_group %in% sort(unique(age_group))[target_age_groups],
      "V1",
      "V2"
    )
  ) %>%
  group_by(t, vaccine) %>%
  summarise(
    y = sum(y),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = t, y = y, colour = vaccine)) +
  geom_line() +
  labs(x = "Time", y = "Vaccinated", colour = "Vaccine:")

#Alternatively we could model a uniform distribution by just using primary as V1
#and boosters as V2, however we'd need everyone to get V1 before they can get V2

