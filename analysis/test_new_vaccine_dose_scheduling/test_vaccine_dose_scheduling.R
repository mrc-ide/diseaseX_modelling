library(dplyr); library(ggplot2); library(squire.page)

probs_booster <- squire.page:::probs_booster
durs_booster <- squire.page:::durs_booster
vaccine_pars_booster <- squire.page:::vaccine_pars_booster


primary_campaign_start <- 10
primary_series_protection_delay <- 5
second_dose_delay <- 50
primary_series_vaccination_campaign_duration <- 50
booster_campaign_start <- 150
booster_protection_delay <- 10
booster_campaign_duration <- 50
doses_per_day <- 5 * 10^5

## note currently assumes primary series vaccination begins at 0, you'd need another timepoint in for primary doses if you wanted it to start it later
x <- run_booster(# Misc Params
                 country = "Argentina",
                 time_period = 365,
                 seed = 100,
                 hosp_bed_capacity = 10^8,
                 ICU_bed_capacity = 10^8,
                 dur_R = 365000000,
                 dur_V = rep(365000000, 4),
                 vaccine_coverage_mat = vaccine_pars_booster$vaccine_coverage_mat,
                 vaccine_booster_initial_coverage = rep(0.99, 17),

                 # New vaccine dose params
                 primary_doses = c(rep(0, primary_campaign_start + primary_series_protection_delay),                              ## primary doses injected but folks not yet moved into pV_1 compartment
                                   rep(doses_per_day, primary_series_vaccination_campaign_duration),                                       ## protection developing - folks being moved into pV_1
                                   rep(0, 365 - primary_series_vaccination_campaign_duration - primary_series_protection_delay - primary_campaign_start)), ## no more primary doses delivered
                 second_doses = c(rep(0, primary_campaign_start + primary_series_protection_delay),                                                        ## primary doses injected but folks not yet moved into pV_1 compartment
                                  rep(0, second_dose_delay),                                                                      ## delay between receiving primary dose and second dose
                                  rep(doses_per_day, primary_series_vaccination_campaign_duration),                                        ## moving folks into the second dose primary series compartment - fV_1
                                  rep(0, 365 - primary_series_vaccination_campaign_duration - second_dose_delay - primary_series_protection_delay - primary_campaign_start)),  ## no more secondary doses delivered
                 booster_doses = c(rep(0, booster_campaign_start + booster_protection_delay),                                     ## booster campaign begins later
                                   rep(doses_per_day, booster_campaign_duration),
                                   rep(0, 365 - booster_campaign_duration - booster_protection_delay - booster_campaign_start)))

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
  geom_line(data = check, aes(x = t, y = prop, col = compartment)) +
  facet_wrap(~age_group)

ggplot() +
  geom_line(data = check[check$age_group == "80+", ], aes(x = t, y = prop, col = compartment)) +
  facet_wrap(~age_group) +
  scale_x_continuous(breaks = seq(0, 370, 10))
