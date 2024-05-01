# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/branching_process_spatial_vaccination.R"))

### Epi parameters
# SC1_generation_time <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
# SC1_infection_to_onset <- function(n) { rgamma(n, shape = 0.1, rate = 1) } ## Ask Azra for values (negligible assumed currently)
# SC1_prop_asymptomatic <- 0
# SC1_prob_hosp <- 0.95 # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
# SC1_hospitalisation_delay <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/

SC1_generation_time <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
SC1_infection_to_onset <- function(n) { rgamma(n, shape = 0.1, rate = 1) } ## Ask Azra for values (negligible assumed currently)
SC1_prop_asymptomatic <- 0
SC1_prob_hosp <- 1 # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
SC1_hospitalisation_delay <- function(n) { 12 } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/

### Spatial kernel parameters
mu <- 1
size <- 5
spatial_kernel <- function(n) { rnbinom(n, size = 5, mu = mu) }
spatial_kernel <- function(n) { 0 }
vaccination_radius <- 1000 * mu

### Vaccine related parameters
vaccine_coverage <- 1
vaccine_efficacy_infection <- 0.5
vaccine_efficacy_transmission <- 0.5
vaccine_efficacy_disease <- 0.5
vaccine_logistical_delay <- 2
vaccine_protection_delay <- 7

### Other parameters
pop <- 10^10
check_final_size <- 2000
initial_immune <- 0
seeding_cases <- 3

set.seed(10)
test <- spatial_bp_geog_vacc(mn_offspring = 6,
                             generation_time = SC1_generation_time,
                             spatial_kernel = spatial_kernel,
                             t0 = 0, tf = Inf,
                             check_final_size = check_final_size,
                             seeding_cases = seeding_cases,
                             prop_asymptomatic = SC1_prop_asymptomatic,
                             prob_hosp = SC1_prob_hosp,
                             hospitalisation_delay = SC1_hospitalisation_delay,
                             detection_threshold = 10,
                             vaccine_campaign_radius = vaccination_radius,
                             vaccine_coverage = vaccine_coverage,
                             vaccine_efficacy_infection = vaccine_efficacy_infection,
                             vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                             vaccine_efficacy_disease = vaccine_efficacy_disease,
                             vaccine_logistical_delay = vaccine_logistical_delay,
                             vaccine_protection_delay = vaccine_protection_delay)
sum(!is.na(test$time_infection))

x <- test %>%
  select(ancestor, n_offspring, n_offspring_new, n_offspring_new_new, hospitalised, time_infection,
         vaccinated, time_vaccinated, vaccinated_before_infection, time_protected, protected_before_infection)
table(x$vaccinated)

hist(x$n_offspring)
hist(x$n_offspring_new)
hist(x$n_offspring_new_new)

y <- x %>%
  filter(!is.na(vaccinated) & vaccinated == 1)

mean(y$n_offspring, na.rm = TRUE)
mean(y$n_offspring_new, na.rm = TRUE)
mean(y$n_offspring_new_new, na.rm = TRUE)

### Spatial kernel parameters
mu <- 25
size <- 5
spatial_kernel <- function(n) { rnbinom(n, size = 5, mu = mu) }
# spatial_kernel <- function(n) { 0 }
vaccination_radius <- 35 * mu

### Vaccine related parameters
vaccine_coverage <- 0.8
vaccine_efficacy_infection <- 0.35
vaccine_efficacy_transmission <- 0.35
vaccine_efficacy_disease <- 0.95
vaccine_logistical_delay <- 2
vaccine_protection_delay <- 7

### Other parameters
pop <- 10^10
check_final_size <- 2000
initial_immune <- 0
seeding_cases <- 3

### Sensitivity analysis parameters
R0_scan <- c(1, 1.25, 1.5, 1.75, 2)
surveillance_scan <- c(1, 10, 100)
iterations <- 10

storage <- array(data = NA, dim = c(iterations, length(R0_scan), length(surveillance_scan)))
for (i in 1:length(R0_scan)) {
  for (j in 1:length(surveillance_scan)) {
    for (k in 1:iterations) {
      
      temp <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                   generation_time = SC1_generation_time,
                                   spatial_kernel = spatial_kernel,
                                   t0 = 0, tf = Inf,
                                   check_final_size = check_final_size,
                                   seeding_cases = seeding_cases,
                                   prop_asymptomatic = SC1_prop_asymptomatic,
                                   prob_hosp = SC1_prob_hosp,
                                   hospitalisation_delay = SC1_hospitalisation_delay,
                                   detection_threshold = surveillance_scan[j],
                                   vaccine_campaign_radius = vaccination_radius,
                                   vaccine_coverage = vaccine_coverage,
                                   vaccine_efficacy_infection = vaccine_efficacy_infection,
                                   vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                   vaccine_efficacy_disease = vaccine_efficacy_disease,
                                   vaccine_logistical_delay = vaccine_logistical_delay,
                                   vaccine_protection_delay = vaccine_protection_delay)

      storage[k, i, j] <- sum(!is.na(temp$time_infection))
  
      print(paste0("i = ", i, ", j = ", j, ", k = ", k))
    }
  }
}

x <- reshape2::melt(storage)
colnames(x) <- c("iteration", "R0", "surveillance", "outbreak_size")
test <- x %>%
  mutate(contained = ifelse(outbreak_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(R0, surveillance) %>%
  summarise(proportion_contained = sum(contained) / iterations)
ggplot(test, aes(x = R0, y = proportion_contained, col = factor(surveillance))) +
  geom_line(linewidth = 1) +
  scale_x_continuous(labels = R0_scan) +
  scale_colour_manual(labels = surveillance_scan,
                      values = c("#FEC9F1", "#E899DC", "#B279A7"),
                      name = "Surveillance\nThreshold\nTrigger") +
  theme_bw()
  

# 
# ## note that vax campaign can never be triggered until all infections are generated from index case.
# 
# mn_offspring = 6
# generation_time = SC1_generation_time
# spatial_kernel = spatial_kernel
# t0 = 0
# tf = Inf
# check_final_size = check_final_size
# seeding_cases = seeding_cases
# prop_asymptomatic = SC1_prop_asymptomatic
# prob_hosp = SC1_prob_hosp
# hospitalisation_delay = SC1_hospitalisation_delay
# detection_threshold = 5
# vaccine_campaign_radius = vaccination_radius
# vaccine_coverage = vaccine_coverage
# vaccine_efficacy_infection = vaccine_efficacy_infection
# vaccine_efficacy_transmission = vaccine_efficacy_transmission
# vaccine_efficacy_disease = vaccine_efficacy_disease
# vaccine_logistical_delay = vaccine_logistical_delay
# vaccine_protection_delay = vaccine_protection_delay