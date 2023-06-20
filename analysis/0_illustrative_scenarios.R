# Notes:
## 1) Double check IFR calc (fine to just alter prob_hosp?)
## 2) Are we fine assuming unlimited healthcare capacity?

# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))

### Defining Central Scenarios

#### Demographic Parameters
target_pop <- 1e6
country <- "Argentina"

#### Healthcare Parameters
hosp_bed_capacity <- 100000000                                         
ICU_bed_capacity <- 100000000    

#### Epidemiological an NPI Parameters
Rt <- c(1.5, 2, 3)
tt_Rt <- 1:9 # 9 different NPI scenarios
Tg <- c(7, 14)
IFR <- c(0.5, 1.5)

## Vaccine-Related Parameters
vaccine_scenario <- c("specific_only", "both_vaccines") # which scenario to explore
bpsv_start <- c(30, 50)                        # BPSV distribution start
specific_vaccine_start <- c(100, 200, 365)     # specific vaccine distribution start
efficacy_infection_bpsv <- 0.35                # vaccine efficacy against infection - BPSV
efficacy_disease_bpsv <- 0.8                   # vaccine efficacy against disease - BPSV
efficacy_infection_spec <- 0.55                # vaccine efficacy against infection - specific vaccine
efficacy_disease_spec <- 0.9                   # vaccine efficacy against disease - specific vaccine
dur_R <- 1000 * 365                            # duration of infection-induced immunity
dur_V_vec <- rep(1000 * 365, 4)                # duration of vaccine-induced immunity for both vaccines
second_dose_delay <- 7                         # controls how many days after "1st dose" people receive second dose; see here: https://github.com/mrc-ide/squire.page/blob/main/inst/odin/nimue_booster.R#L427-L430
dur_vacc_delay <- 7                            # mean duration from vaccination to protection
coverage <- 0.8                                # proportion of the population vaccinated
vaccination_rate <- 0.035                       # vaccination rate per week as percentage of population
min_age_group_index_priority <- 13             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
min_age_group_index_non_priority <- 4          # index of the youngest age group that *receives* vaccines (4 = 15+)


# Expand_grid creates bpsv_start * specific_only scenarios, but they're all the same as no BPSV, so keep only first set
baseline_scenarios <- expand_grid(population_size = 1e6,
                                  country = "Argentina",
                                  hosp_bed_capacity = hosp_bed_capacity,                                         
                                  ICU_bed_capacity = ICU_bed_capacity,
                                  Rt = Rt,
                                  tt_Rt = tt_Rt,  ## placeholder for now, used to specify NPI scenario later
                                  Tg = Tg,
                                  IFR = IFR,
                                  vaccine_scenario = vaccine_scenario,
                                  bpsv_start = bpsv_start,
                                  specific_vaccine_start = specific_vaccine_start,
                                  efficacy_infection_bpsv = efficacy_infection_bpsv,
                                  efficacy_disease_bpsv = seq(0.35, 0.9, 0.05), #efficacy_disease_bpsv,
                                  efficacy_infection_spec = efficacy_infection_spec, 
                                  efficacy_disease_spec = efficacy_disease_spec,
                                  dur_R = dur_R, 
                                  dur_V = 1, # placeholder to be filled in below
                                  second_dose_delay = second_dose_delay, 
                                  dur_vacc_delay = dur_vacc_delay,
                                  coverage = coverage,
                                  vaccination_rate = vaccination_rate, 
                                  min_age_group_index_priority = min_age_group_index_priority,
                                  min_age_group_index_non_priority = min_age_group_index_non_priority,
                                  runtime = 2*365,
                                  seeding_cases = 2,
                                  NPI_scenario = "temp") %>% 
  filter(vaccine_scenario == "both_vaccines" | (bpsv_start == bpsv_start[1] & vaccine_scenario == "specific_only"))

## Coercing dur_V into correct format
baseline_scenarios <- baseline_scenarios %>%
  mutate(dur_V = list(dur_V_vec))

## Generating NPI scenarios
### 1 = NPIs reducing Rt < 1 until BPSV campaign is done, followed by minimum mandate until specific campaign is done, then full release
### 2 = NPIs reducing Rt < 1 until BPSV campaign is done, followed by full release
### 3 = NPIs reducing Rt < 1 until spec vaccine campaign is done, followed by full release
### 4 = NPIs minimal mandate until BPSV campaign is done, followed by full release
### 5 = NPIs minimal mandate until spec vaccine campaign is done, followed by full release
### 6 = NPIs reducing Rt < 1 for fixed amount of time, followed by minimum mandate until BPSV campaign is done, then full release 
### 7 = NPIs reducing Rt < 1 for fixed amount of time, followed by minimum mandate until spec campaign is done, then full release 
### 8 = NPIs for a fixed calendar amount of time e.g. 30 days followed by full release
### 9 = No NPIs
#### NOTE - REMEMBER THAT SPECIFIC VACCINE DEVELOPMENT MUST BE *AFTER* THE BPSV CAMPAIGN IS DONE
#### NOTE - MIGHT NEED TO TWEAK THESE SOMEWHAT DEPENDING ON WHETHER EMERGENCE OR SECONDARY COUNTRY
detection_time <- 14
minimal_mandate_reduction <- 0.25
lockdown_Rt <- 0.9
fixed_lockdown_time <- 14
standard_pop <- generate_standard_pop(country = unique(baseline_scenarios$country), population_size = unique(baseline_scenarios$population_size))
daily_doses <- unique(baseline_scenarios$vaccination_rate) * unique(baseline_scenarios$population_size) / 7    # rate of vaccination with primary series
priority_age_groups <- unique(baseline_scenarios$min_age_group_index_priority):17            
elderly_pop_to_vaccinate <- sum(standard_pop[priority_age_groups]) * coverage # 60+s receive primary (BNPCV) and booster (diseaseX-specific); under 60s receive just primary (diseaseX-specific)
time_to_coverage_bpsv <- ceiling(elderly_pop_to_vaccinate/daily_doses) + 1
time_to_coverage_spec <- time_to_coverage_bpsv ## assumed same rate for now

tic()
scenarios <- baseline_scenarios %>%
  rowwise() %>%
  mutate(dur_V = list(dur_V_vec)) %>%
  mutate(NPI_scenario_int = tt_Rt) %>% 
  mutate(NPI_scenario = case_when(tt_Rt == 1 ~ "LockdownBPSVFinish_minMandateSpecFinish_fullRelease",
                                  tt_Rt == 2 ~ "LockdownBPSVFinish_fullRelease",
                                  tt_Rt == 3 ~ "LockdownSpecFinish_fullRelease",
                                  tt_Rt == 4 ~ "minMandateBPSVFinish_fullRelease",
                                  tt_Rt == 5 ~ "minMandateSpecFinish_fullRelease",
                                  tt_Rt == 6 ~ "fixedLockdownTime_minMandateBPSVFinish_fullRelease",  ## watch out for ordering of Rts and tt_Rts here
                                  tt_Rt == 7 ~ "fixedLockdownTime_minMandateSpecFinish_fullRelease",  ## watch out for ordering of Rts and tt_Rts here
                                  tt_Rt == 8 ~ "fixedLockdownTime_fullRelease",
                                  tt_Rt == 9 ~ "Nothing")) %>%
  mutate(Rt = case_when(tt_Rt == 1 ~ list(c(Rt, lockdown_Rt, Rt * (1 - minimal_mandate_reduction), Rt)),
                        tt_Rt == 2 ~ list(c(Rt, lockdown_Rt, Rt)),
                        tt_Rt == 3 ~ list(c(Rt, lockdown_Rt, Rt)),
                        tt_Rt == 4 ~ list(c(Rt, Rt * (1 - minimal_mandate_reduction), Rt)),
                        tt_Rt == 5 ~ list(c(Rt, Rt * (1 - minimal_mandate_reduction), Rt)),
                        tt_Rt == 6 ~ list(c(Rt, lockdown_Rt, Rt * (1 - minimal_mandate_reduction), Rt)), 
                        tt_Rt == 7 ~ list(c(Rt, lockdown_Rt, Rt * (1 - minimal_mandate_reduction), Rt)),  
                        tt_Rt == 8 ~ list(c(Rt, lockdown_Rt, Rt)),
                        tt_Rt == 9 ~ list(Rt))) %>%
  mutate(tt_Rt = case_when(tt_Rt == 1 ~ list(c(0, detection_time, detection_time + time_to_coverage_bpsv, detection_time + specific_vaccine_start + time_to_coverage_spec)),
                           tt_Rt == 2 ~ list(c(0, detection_time, detection_time + time_to_coverage_bpsv)),
                           tt_Rt == 3 ~ list(c(0, detection_time, detection_time + specific_vaccine_start + time_to_coverage_spec)),
                           tt_Rt == 4 ~ list(c(0, detection_time, detection_time + time_to_coverage_bpsv)),
                           tt_Rt == 5 ~ list(c(0, detection_time, detection_time + specific_vaccine_start + time_to_coverage_spec)),
                           tt_Rt == 6 ~ list(c(0, detection_time, detection_time + fixed_lockdown_time, detection_time + time_to_coverage_bpsv)), # watch out for this one - poss for 2nd element to be > 3rd element - still unsure what to do with this 
                           tt_Rt == 7 ~ list(c(0, detection_time, detection_time + fixed_lockdown_time, detection_time + specific_vaccine_start + time_to_coverage_spec)),  
                           tt_Rt == 8 ~ list(c(0, detection_time, detection_time + fixed_lockdown_time)),
                           tt_Rt == 9 ~ list(0))) # No NPIs
toc()

tic()
plan(multisession, workers = 55) # multicore does nothing on windows as multicore isn't supported
system.time({out <- future_pmap(scenarios, run_sars_x, .progress = TRUE)})
toc()
#329 seconds to run 11664 simulations with 40 cores - effectively 0.028 seconds per iteration 
#271 seconds to run 11664 simulations with 55 cores - effectively 0.023 seconds per iteration 

#### code here to take "out" and create summary dataframe
### need to add NPI scenario in here
cl <- makeCluster(5)
clusterEvalQ(cl, {
  library(data.table)
})
tic()
data <- parLapply(cl, out, function(x) {
  y <- data.frame(x$summary_metrics, 
                  country = x$model_arguments$country,
                  population_size = x$model_arguments$population_size,
                  hosp_bed_capacity = x$model_arguments$hosp_bed_capacity,
                  ICU_bed_capacity = x$model_arguments$ICU_bed_capacity,
                  R0 = x$model_arguments$Rt[1],
                  Tg = x$model_arguments$Tg,
                  IFR = x$model_arguments$IFR,
                  vaccine_scenario = x$model_arguments$vaccine_scenario,
                  bpsv_start = ifelse(x$model_arguments$vaccine_scenario == "specific_only", NA, x$model_arguments$bpsv_start),
                  specific_vaccine_start = x$model_arguments$specific_vaccine_start,
                  efficacy_infection_bpsv = x$model_arguments$efficacy_infection_bpsv,
                  efficacy_disease_bpsv = x$model_arguments$efficacy_disease_bpsv,
                  efficacy_infection_spec = x$model_arguments$efficacy_infection_spec,
                  efficacy_disease_spec = x$model_arguments$efficacy_disease_spec,
                  dur_R = x$model_arguments$dur_R,
                  dur_V = x$model_arguments$dur_V[1],
                  second_dose_delay = x$model_arguments$second_dose_delay,
                  dur_vacc_delay = x$model_arguments$dur_vacc_delay,
                  coverage = x$model_arguments$coverage,
                  vaccination_rate = x$model_arguments$vaccination_rate,
                  min_age_group_index_priority = x$model_arguments$min_age_group_index_priority,
                  min_age_group_index_non_priority = x$model_arguments$min_age_group_index_non_priority,
                  runtime = x$model_arguments$runtime,
                  seeding_cases = x$model_arguments$seeding_cases,
                  NPI_scenario = x$model_arguments$NPI_scenario,
                  NPI_scenario_int = x$model_arguments$NPI_scenario_int) 
})
toc()
stopCluster(cl) 

combined_data <- rbindlist(data)
table(combined_data$vaccine_scenario, useNA = "ifany")

scenarios_to_plot <- combined_data

two_vax <- scenarios_to_plot %>%
  filter(vaccine_scenario == "both_vaccines") %>%
  rename(deaths_bpsv = deaths,
         time_under_NPIs_bpsv = time_under_NPIs,
         composite_NPI_bpsv = composite_NPI)

one_vax <- scenarios_to_plot %>%
  filter(vaccine_scenario == "specific_only") %>%
  select(R0, Tg, IFR, NPI_scenario, specific_vaccine_start, deaths, time_under_NPIs, composite_NPI, efficacy_disease_bpsv) %>% ## IMPORTANT - YOU NEED TO INCLUDE THE VARIABLE YOU'RE DOING THE SENSITIVITY ANALYSIS OVER
  rename(deaths_spec = deaths,
         time_under_NPIs_spec = time_under_NPIs,
         composite_NPI_spec = composite_NPI)

joined <- two_vax %>%
  left_join(one_vax, by = c("R0", "Tg", "IFR", "NPI_scenario", "specific_vaccine_start", "efficacy_disease_bpsv")) %>% ## IMPORTANT - YOU NEED TO INCLUDE THE VARIABLE YOU'RE DOING THE SENSITIVITY ANALYSIS OVER
  mutate(deaths_saved = deaths_spec - deaths_bpsv)

plot(joined$deaths_bpsv, joined$deaths_spec, xlim = c(0, 17500), ylim = c(0, 17500))
lines(1:17500, 1:17500)
plot(joined$efficacy_disease_bpsv, joined$deaths_saved)

test_plot <- joined %>%
  filter(R0 == 2,
         Tg == 7,
         IFR == 1.5) %>%
  mutate(specific_vaccine_start = factor(specific_vaccine_start),
         bpsv_start = factor(bpsv_start),
         NPI_scenario_int = factor(NPI_scenario_int))

ggplot(test_plot, aes(x = efficacy_disease_bpsv, y = deaths_saved, colour = specific_vaccine_start)) +
  geom_path() +
  facet_grid(bpsv_start ~ NPI_scenario_int)

check <- test_plot[test_plot$NPI_scenario_int == 9 & 
                     test_plot$specific_vaccine_start == 100 &
                     #test_plot$efficacy_disease_bpsv == 0.5 &
                     test_plot$bpsv_start == 50, ]
ggplot(check, aes(x = efficacy_disease_bpsv, y = deaths_saved, colour = specific_vaccine_start)) +
  geom_path() +
  facet_grid(bpsv_start ~ NPI_scenario_int)

ggplot(check, aes(x = efficacy_disease_bpsv, y = deaths_bpsv, colour = specific_vaccine_start)) +
  geom_path() 
ggplot(check, aes(x = efficacy_disease_bpsv, y = deaths_spec, colour = specific_vaccine_start)) +
  geom_path() 

check2 <- scenarios_to_plot %>%
  filter(R0 == 2,
         Tg == 7,
         IFR == 1.5,
         NPI_scenario_int == 9,
         specific_vaccine_start == 100,
         is.na(bpsv_start),
         vaccine_scenario == "specific_only")

check_og <- baseline_scenarios %>%
  filter(Rt == 2,
         Tg == 7,
         IFR == 1.5,
         tt_Rt == 9,
         specific_vaccine_start == 100,
         bpsv_start == 30,
         vaccine_scenario == "specific_only")

index <- which(baseline_scenarios$Rt == 2 & 
      baseline_scenarios$Tg == 7 &
      baseline_scenarios$IFR == 1.5 & 
      baseline_scenarios$tt_Rt == 9 &
      baseline_scenarios$specific_vaccine_start == 100 &
      baseline_scenarios$bpsv_start == 30 &
      baseline_scenarios$vaccine_scenario == "specific_only")

test_run <- scenarios[index, ]
tic()
plan(multisession, workers = 5) 
system.time({out <- future_pmap(test_run, run_sars_x, .progress = TRUE)})
toc()

data <- lapply(out, function(x) {
  y <- data.frame(x$summary_metrics,
                  R0 = x$model_arguments$Rt[1],
                  Tg = x$model_arguments$Tg,
                  IFR = x$model_arguments$IFR,
                  vaccine_scenario = x$model_arguments$vaccine_scenario,
                  bpsv_start = x$model_arguments$bpsv_start,
                  specific_vaccine_start = x$model_arguments$specific_vaccine_start,
                  efficacy_infection_bpsv = x$model_arguments$efficacy_infection_bpsv,
                  efficacy_disease_bpsv = x$model_arguments$efficacy_disease_bpsv,
                  efficacy_infection_spec = x$model_arguments$efficacy_infection_spec,
                  efficacy_disease_spec = x$model_arguments$efficacy_disease_spec,
                  NPI_scenario_int = x$model_arguments$NPI_scenario_int)
})
x <- rbindlist(data)

check3 <- check2 %>%
  select(deaths, R0, Tg, IFR, bpsv_start, specific_vaccine_start, efficacy_disease_bpsv) %>%
  mutate(bpsv_start = 30)

z <- x %>%
  left_join(check3, by = c("R0", "Tg", "IFR", "bpsv_start", "specific_vaccine_start", "efficacy_disease_bpsv")) %>%
  relocate(deaths.y)

plot(z$deaths.x, z$deaths.y)

check2$deaths[order(check2$deaths)]
x$deaths[order(x$deaths)]


## Checking weirdness:
index <- which(baseline_scenarios$Rt == 2 & 
                 baseline_scenarios$Tg == 7 &
                 baseline_scenarios$IFR == 1.5 & 
                 baseline_scenarios$tt_Rt == 9 &
                 baseline_scenarios$specific_vaccine_start == 100 &
                 baseline_scenarios$bpsv_start == 30 &
                 baseline_scenarios$vaccine_scenario == "specific_only")
test_run <- scenarios[index, ]
tic()
plan(multisession, workers = 5) 
system.time({out <- future_pmap(test_run, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
toc()
data <- lapply(out, function(x) {
  y <- data.frame(x$summary_metrics,
                  R0 = x$model_arguments$Rt[1],
                  Tg = x$model_arguments$Tg,
                  IFR = x$model_arguments$IFR,
                  vaccine_scenario = x$model_arguments$vaccine_scenario,
                  bpsv_start = x$model_arguments$bpsv_start,
                  specific_vaccine_start = x$model_arguments$specific_vaccine_start,
                  efficacy_infection_bpsv = x$model_arguments$efficacy_infection_bpsv,
                  efficacy_disease_bpsv = x$model_arguments$efficacy_disease_bpsv,
                  efficacy_infection_spec = x$model_arguments$efficacy_infection_spec,
                  efficacy_disease_spec = x$model_arguments$efficacy_disease_spec,
                  NPI_scenario_int = x$model_arguments$NPI_scenario_int)
})
x <- rbindlist(data)
deaths1 <- x$deaths

[1] 7374.526 7411.961 7371.021 7328.408 7543.190 7578.105 7543.190 7374.526 7583.370 7415.128 7328.408 7583.370
[1] 7374.526 7411.961 7371.021 7328.408 7543.190 7578.105 7543.190 7374.526 7583.370 7415.128 7328.408 7583.370

index <- which(baseline_scenarios$Rt == 2 & 
                 baseline_scenarios$Tg == 7 &
                 baseline_scenarios$IFR == 1.5 & 
                 baseline_scenarios$tt_Rt == 9 &
                 baseline_scenarios$specific_vaccine_start == 100 &
                 baseline_scenarios$bpsv_start == 30 &
                 baseline_scenarios$vaccine_scenario == "specific_only")
test_run <- scenarios[index, ]
# test_run$bpsv_start <- 30
tic()
plan(multisession, workers = 5) 
system.time({out <- future_pmap(test_run, run_sars_x, .progress = TRUE)})
toc()
data <- lapply(out, function(x) {
  y <- data.frame(x$summary_metrics,
                  R0 = x$model_arguments$Rt[1],
                  Tg = x$model_arguments$Tg,
                  IFR = x$model_arguments$IFR,
                  vaccine_scenario = x$model_arguments$vaccine_scenario,
                  bpsv_start = x$model_arguments$bpsv_start,
                  specific_vaccine_start = x$model_arguments$specific_vaccine_start,
                  efficacy_infection_bpsv = x$model_arguments$efficacy_infection_bpsv,
                  efficacy_disease_bpsv = x$model_arguments$efficacy_disease_bpsv,
                  efficacy_infection_spec = x$model_arguments$efficacy_infection_spec,
                  efficacy_disease_spec = x$model_arguments$efficacy_disease_spec,
                  NPI_scenario_int = x$model_arguments$NPI_scenario_int)
})
x <- rbindlist(data)
deaths2 <- x$deaths

deaths1
deaths2
plot(deaths1, deaths2)

## Plotting out the NPI scenarios
Rt_check <- Rt[3]
specific_vaccine_start_check <- specific_vaccine_start[1]

Rt_S1 <- c(Rt_check, lockdown_Rt, Rt_check * (1 - minimal_mandate_reduction), Rt_check)
tt_Rt_S1 <- c(0, detection_time, detection_time + time_to_coverage_bpsv, detection_time + specific_vaccine_start_check + time_to_coverage_spec)

Rt_S2 <- c(Rt_check, lockdown_Rt, Rt_check)
tt_Rt_S2 <- c(0, detection_time, detection_time + time_to_coverage_bpsv)

Rt_S3 <- c(Rt_check, lockdown_Rt, Rt_check)
tt_Rt_S3 <- c(0, detection_time, detection_time + specific_vaccine_start_check + time_to_coverage_spec)

Rt_S4 <- c(Rt_check, Rt_check * (1 - minimal_mandate_reduction), Rt_check)
tt_Rt_S4 <- c(0, detection_time, detection_time + time_to_coverage_bpsv)

Rt_S5 <- c(Rt_check, Rt_check * (1 - minimal_mandate_reduction), Rt_check)
tt_Rt_S5 <- c(0, detection_time, detection_time + specific_vaccine_start_check + time_to_coverage_spec)

Rt_S6 <- c(Rt_check, lockdown_Rt, Rt_check * (1 - minimal_mandate_reduction), Rt_check) ## unsure about this one as depending on fixed lockdown time, potentially shorter than time_to_coverage
tt_Rt_S6 <- c(0, detection_time, detection_time + fixed_lockdown_time, detection_time + time_to_coverage_bpsv) 

Rt_S7 <- c(Rt_check, lockdown_Rt, Rt_check * (1 - minimal_mandate_reduction), Rt_check) ## less worried about this one because specific vaccine will take min 100 days
tt_Rt_S7 <- c(0, detection_time, detection_time + fixed_lockdown_time, detection_time + specific_vaccine_start_check + time_to_coverage_spec) 

Rt_S8 <- c(Rt_check, lockdown_Rt, Rt_check)
tt_Rt_S8 <- c(0, detection_time, detection_time + fixed_lockdown_time)

Rt_S9 <- Rt_check
tt_Rt_S9 <- 0

NPI_df <- data.frame(Rt = c(Rt_S1, Rt_S2, Rt_S3, Rt_S4, Rt_S5, Rt_S6, Rt_S7, Rt_S8, c(Rt_S9, Rt_S9)),
                     tt_Rt = c(tt_Rt_S1, tt_Rt_S2, tt_Rt_S3, tt_Rt_S4, tt_Rt_S5, tt_Rt_S6, tt_Rt_S7, tt_Rt_S8, c(tt_Rt_S9, 10)),
                     scenario = c(rep("1", length(Rt_S1)), 
                                  rep("2", length(Rt_S2)),
                                  rep("3", length(Rt_S3)),
                                  rep("4", length(Rt_S4)),
                                  rep("5", length(Rt_S5)),
                                  rep("6", length(Rt_S6)),
                                  rep("7", length(Rt_S7)),
                                  rep("8", length(Rt_S8)),
                                  rep("9", length(Rt_S9) + 1)))
NPI_df$scenario <- paste0("Scenario ", NPI_df$scenario)
NPI_df <- NPI_df %>%
  group_by(scenario) %>%
  mutate(next_time = lead(tt_Rt),
         next_value = lead(Rt)) %>%
  mutate(next_time = ifelse(is.na(next_time), unique(baseline_scenarios$runtime), next_time),
         next_value = ifelse(is.na(next_value), Rt_check, next_value))
overplot_factor <- 1
example_NPI <- subset(NPI_df, scenario == "Scenario 1")
example_NPI$next_time[length(example_NPI$next_time)] <- 150
a <- ggplot(example_NPI) +
  geom_segment(aes(x = tt_Rt - overplot_factor, xend = next_time + overplot_factor, y = Rt, yend = Rt), size = 2) +
  geom_segment(aes(x = next_time, xend = next_time, y = Rt, yend = next_value), size = 2) +
  geom_hline(aes(yintercept = 1)) +
  geom_hline(aes(yintercept = lockdown_Rt), linetype = "dashed") +
  geom_hline(aes(yintercept = Rt_check * (1 - minimal_mandate_reduction)), linetype = "dashed") +
  geom_vline(aes(xintercept = detection_time)) +
  geom_vline(aes(xintercept = detection_time + time_to_coverage_bpsv)) +
  geom_vline(aes(xintercept = detection_time + specific_vaccine_start_check + time_to_coverage_spec)) +
  theme_bw() +
  geom_segment(aes(x = detection_time, y = 4, xend = detection_time, yend = Rt_check + 0.7), arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = detection_time, y = 4.1, label = "Detection") +
  geom_segment(aes(x = detection_time + time_to_coverage_bpsv, y = 4, xend = detection_time + time_to_coverage_bpsv, yend = Rt_check + 0.7), arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = detection_time + time_to_coverage_bpsv, y = 4.2, label = "BPSV Vaccination\nCompleted", hjust = 0) +
  geom_segment(aes(x = detection_time + specific_vaccine_start_check + time_to_coverage_spec, y = 4, xend = detection_time + specific_vaccine_start_check + time_to_coverage_spec, yend = Rt_check + 0.7), arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = detection_time + specific_vaccine_start_check + time_to_coverage_spec, y = 4.2, label = "Spec. Vaccination\nCompleted", hjust = 0.5) +
  geom_segment(aes(x = -20, y = Rt_check * (1 - minimal_mandate_reduction), xend = -10, yend = Rt_check * (1 - minimal_mandate_reduction)), arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = -20, y = Rt_check * (1 - minimal_mandate_reduction), label = "Min.\nMandate", hjust = 1) +
  geom_segment(aes(x = -20, y = lockdown_Rt, xend = -10, yend = 0.9), arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = -20, y = lockdown_Rt, label = "Lockdown", hjust = 1) +
  geom_segment(aes(x = -20, y = Rt_check, xend = -10, yend = Rt_check), arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = -20, y = Rt_check, label = "R0", hjust = 1) +
  theme(plot.margin = margin(2.5, 1, 2.5, 2.5, "cm")) +
  coord_cartesian(clip = 'off', xlim = c(0, 150), ylim = c(0.5, Rt_check + 0.5)) +
  scale_y_continuous(position = "right") +
  labs(x = "Time (Days)")

b <- ggplot(NPI_df, aes(x = tt_Rt - overplot_factor, colour = scenario)) +
  geom_hline(aes(yintercept = 1), linewidth = 0.2) +
  geom_hline(aes(yintercept = lockdown_Rt), linetype = "dashed", linewidth = 0.2) +
  geom_hline(aes(yintercept = Rt_check * (1 - minimal_mandate_reduction)), linetype = "dashed", linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time + time_to_coverage_bpsv), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time + specific_vaccine_start_check + time_to_coverage_spec), linewidth = 0.2) +
  geom_segment(aes(xend = next_time + overplot_factor, y = Rt, yend = Rt), size = 1) +
  geom_segment(aes(x = next_time, xend = next_time, y = Rt, yend = next_value), size = 1) +
  theme_bw() +
  facet_wrap(~scenario) +
  coord_cartesian(xlim = c(0, 150), ylim = c(0.5, Rt_check + 0.5)) +
  theme(legend.position = "none")

cowplot::plot_grid(a, b, rel_widths = c(1, 1.2)) # 11.3 x 5.7 dimensions is good

