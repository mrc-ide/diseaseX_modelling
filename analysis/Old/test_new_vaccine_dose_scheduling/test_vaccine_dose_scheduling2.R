library(dplyr); library(ggplot2); library(squire.page)

probs_booster <- squire.page:::probs_booster
durs_booster <- squire.page:::durs_booster
vaccine_pars_booster <- squire.page:::vaccine_pars_booster

primary_campaign_start <- 10
primary_series_protection_delay <- 5
second_dose_delay <- 1
primary_series_vaccination_campaign_duration <- 100
booster_campaign_start <- 40
booster_protection_delay <- 10
booster_campaign_duration <- 250
doses_per_day <- 5 * 10^4
runtime <- 700

# vaccine_coverage_mat <- matrix(rep(0.95, 17), ncol = 17, byrow = TRUE)
vaccine_coverage_mat <- matrix(c(rep(0, 16), rep(0.95, 1)), ncol = 17, byrow = TRUE)

vaccine_coverage_mat <- matrix(c(rep(0, 16), rep(0.95, 1),
                                 rep(0, 5), rep(0.95, 12)), ncol = 17, byrow = TRUE)


vaccine_efficacy_disease <- vaccine_pars_booster$vaccine_efficacy_disease
vaccine_efficacy_disease[vaccine_efficacy_disease<1] <- 0

## Q1 when I have vaccine_efficacy_disease set to 1, how is it possible that my numbers of vaccinated_first_dose folk decline over time?
## Q2 when I vary number of 0s in vaccine_coverage_mat, how comes I get weird behaviour where boosters doesn't reach the top?
#### ==> let's try and replicate this with Greg's baseline model

## note currently assumes primary series vaccination begins at 0, you'd need another timepoint in for primary doses if you wanted it to start it later
x <- run_booster(# Misc Params
  country = "Argentina",
  time_period = runtime,
  seed = 100,
  R0 = 4,
  prob_hosp = probs_booster$prob_hosp * 3,
  hosp_bed_capacity = 10^8,
  ICU_bed_capacity = 10^8,
  dur_R = 365000000,
  dur_V = rep(365000000, 4),
  vaccine_coverage_mat = vaccine_coverage_mat,
  vaccine_booster_initial_coverage = c(rep(0, 16), 1),
  vaccine_efficacy_disease = vaccine_efficacy_disease,

  # New vaccine dose params
  primary_doses = c(rep(0, primary_campaign_start + primary_series_protection_delay),                              ## primary doses injected but folks not yet moved into pV_1 compartment
                    rep(doses_per_day, primary_series_vaccination_campaign_duration),                                       ## protection developing - folks being moved into pV_1
                    rep(0, runtime - primary_series_vaccination_campaign_duration - primary_series_protection_delay - primary_campaign_start)), ## no more primary doses delivered
  second_doses = c(rep(0, primary_campaign_start + primary_series_protection_delay),                                                        ## primary doses injected but folks not yet moved into pV_1 compartment
                   rep(0, second_dose_delay),                                                                      ## delay between receiving primary dose and second dose
                   rep(doses_per_day, primary_series_vaccination_campaign_duration),                                        ## moving folks into the second dose primary series compartment - fV_1
                   rep(0, runtime - primary_series_vaccination_campaign_duration - second_dose_delay - primary_series_protection_delay - primary_campaign_start)),  ## no more secondary doses delivered
  booster_doses = c(rep(0, booster_campaign_start + booster_protection_delay),                                     ## booster campaign begins later
                    rep(doses_per_day, booster_campaign_duration),
                    rep(0, runtime - booster_campaign_duration - booster_protection_delay - booster_campaign_start)))

check <- nimue::format(x, compartments = c("vaccinated_first_dose", "vaccinated_second_dose", "vaccinated_booster_dose"),
                       reduce_age = FALSE) %>%
  filter(t > 1,
         compartment == "vaccinated_first_dose" | compartment == "vaccinated_second_dose" | compartment == "vaccinated_booster_dose") %>%
  group_by(replicate, t)

pop_df <- data.frame(age_group = sort(unique(check$age_group)), population = x$parameters$population)
check <- check %>%
  left_join(pop_df, by = "age_group") %>%
  mutate(prop = value / population)

ggplot() +
  geom_line(data = check[check$age_group == "80+", ], aes(x = t, y = value, col = compartment)) +
  facet_wrap(~age_group)

# ggplot() +
#   geom_line(data = check, aes(x = t, y = prop, col = compartment)) +
#   facet_wrap(~age_group)


# vaccine_coverage_mat <- matrix(c(rep(0, 13), rep(0.95, 4),
#                                  rep(0, 5), rep(0.95, 12)), ncol = 17, byrow = TRUE)
# vaccine_coverage_mat <- matrix(c(rep(0, 5), rep(0.95, 12)), ncol = 17, byrow = TRUE)

