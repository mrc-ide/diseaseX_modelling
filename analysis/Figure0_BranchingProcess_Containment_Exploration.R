# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/branching_process.R"))

## Currently missing:
### Hospitalisations as a trigger for implementing the vaccine 
### Vaccine efficacy waning

## Multiple simulations
generation_time <- function(n) { rgamma(n, shape = 13.5, rate = 2) } # 6.75 day generation time Gamam distributed (as per Walker et al, Science, 2020)
infection_to_onset <- function(n) { rgamma(n, shape = 6.75, rate = 2) } # ~50% of transmission presymptomatic (per SARS-CoV-2, see https://bmjopen.bmj.com/content/11/6/e041240)
hist(generation_time(10^5))
hist(infection_to_onset(10^5))
median(generation_time(10^5) / infection_to_onset(10^5))
iterations <- 100
R0_scan <- c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5)
check_final_size <- 1000
vaccine_start <- 25
vaccine_coverage <- 0.8
storage_nothing <- matrix(nrow = iterations, ncol = length(R0_scan))
storage_monoclonal <- matrix(nrow = iterations, ncol = length(R0_scan))
storage_vacc <- matrix(nrow = iterations, ncol = length(R0_scan))

for (i in 1:length(R0_scan)) {
  for (j in 1:iterations) {
    
    vacc <- chain_sim_susc_ring_vacc(offspring = "pois",
                                     mn_offspring = R0_scan[i],
                                     generation_time = generation_time,
                                     t0 = 0, tf = Inf, pop = 10^7, check_final_size = check_final_size, initial_immune = 0,
                                     seeding_cases = 5, prop_asymptomatic = 0.2,
                                     infection_to_onset = infection_to_onset,
                                     vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                     vaccine_efficacy_infection = 0.75,
                                     vaccine_efficacy_transmission = 0.55,
                                     vaccine_logistical_delay = 2,
                                     vaccine_protection_delay = 5)
    storage_vacc[j, i] <- nrow(vacc)
    
    monoclonal <- chain_sim_susc_ring_vacc(offspring = "pois",
                                           mn_offspring = R0_scan[i],
                                           generation_time = generation_time,
                                           t0 = 0, tf = Inf, pop = 10^7, check_final_size = check_final_size, initial_immune = 0,
                                           seeding_cases = 5, prop_asymptomatic = 0.2,
                                           infection_to_onset = infection_to_onset,
                                           vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                           vaccine_efficacy_infection = 0.75,
                                           vaccine_efficacy_transmission = 0.55,
                                           vaccine_logistical_delay = 2,
                                           vaccine_protection_delay = 0)
    storage_monoclonal[j, i] <- nrow(monoclonal)
    
    no_vacc <- chain_sim_susc_ring_vacc(offspring = "pois",
                                        mn_offspring = R0_scan[i],
                                        generation_time = generation_time,
                                        t0 = 0, tf = Inf, pop = 10^7, check_final_size = check_final_size, initial_immune = 0,
                                        seeding_cases = 5, prop_asymptomatic = 0,
                                        infection_to_onset = infection_to_onset,
                                        vaccine_start = 1000, vaccine_coverage = vaccine_coverage,
                                        vaccine_efficacy_infection = 0,
                                        vaccine_efficacy_transmission = 0,
                                        vaccine_logistical_delay = 100,
                                        vaccine_protection_delay = 100)
    storage_nothing[j, i] <- nrow(no_vacc)
    if (j %% 25 == 0) {
      print(j)
    }
  }
  print(i)
}

## Quick check plotting
plot(R0_scan, apply(storage_vacc, 2, median), type = "l", xlab = "R0", ylab = paste0("Final Epidemic Size (Capped at ", check_final_size),
     ylim = c(0, check_final_size))
lines(R0_scan, apply(storage_nothing, 2, median), col = "blue")
lines(R0_scan, apply(storage_monoclonal, 2, median), col = "red")

## Creating overall dataframe with all the results
no_vacc <- data.frame(iteration = 1:iterations, scenario = "no_vaccination", storage_nothing)
colnames(no_vacc) <- c("iteration", "scenario", paste0("R0=", R0_scan))
vacc <- data.frame(iteration = 1:iterations, scenario = "vaccination", storage_vacc)
colnames(vacc) <- c("iteration", "scenario", paste0("R0=", R0_scan))
monoclonal <- data.frame(iteration = 1:iterations, scenario = "monoclonal", storage_monoclonal)
colnames(monoclonal) <- c("iteration", "scenario", paste0("R0=", R0_scan))

overall_bp_df <- rbind(no_vacc, vacc, monoclonal) %>%
  pivot_longer(cols = starts_with("R0"), names_to = "R0", values_to = "Epidemic Size") %>%
  mutate(actualR0 = as.numeric(gsub("R0=", "", R0)))
saveRDS(object = overall_bp_df, file = "outputs/branching_process_outputs/branchingProcess_R0_Scan.rds")
containment_df <- overall_bp_df %>%
  mutate(contained = ifelse(`Epidemic Size` < (0.9 * check_final_size), 1, 0)) %>%
  group_by(actualR0, R0, scenario) %>%
  summarise(proportion_contained = sum(contained) / iterations)
# ggplot(overall_bp_df) +
#   geom_boxplot(aes(x = actualR0, y = `Epidemic Size`, group = actualR0, col = scenario)) +
#   facet_wrap(scenario ~ .)

## Percentage of Outbreaks Contained Plot
containment_plot <- ggplot(containment_df, aes(x = actualR0, y = 100 * proportion_contained, col = scenario)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_colour_manual(values = c("#6DB937", "#3E9BDB", "#CA2E6B"), 
                      labels = c("Monoclonal", "No Vaccination", "BPSV"),
                      name = "Intervention") +
  labs(x = "R0", y = "% Outbreaks (Self-)Contained")

## Added Value Plot
df_m_v <- containment_df %>% filter(scenario %in% c('monoclonal', 'vaccination'))
df_nv <- containment_df %>% filter(scenario == 'no_vaccination')
df_merged <- left_join(df_m_v, df_nv, by = c('R0', "actualR0"))
df_final <- df_merged %>%
  select(-scenario.y) %>%
  rename(proportion_contained_no_vaccination = "proportion_contained.y",
         proportion_contained_intervention = "proportion_contained.x",
         scenario = "scenario.x") %>%
  mutate(proportion_contained_diff = proportion_contained_intervention - proportion_contained_no_vaccination)

containment_diff_plot <- ggplot(df_final, aes(x = actualR0, y = 100 * proportion_contained_diff, col = scenario)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_colour_manual(values = c("#6DB937", "#CA2E6B"), 
                      labels = c("Monoclonal", "BPSV"),
                      name = "Intervention") +
  labs(x = "R0", y = "Increased % of Outbreaks Contained")

raw_containment_plot <- cowplot::plot_grid(containment_plot + theme(legend.position = "none"), 
                                               containment_diff_plot + theme(legend.position = "none"),
                                               nrow = 1, ncol = 2, rel_widths = c(1.5, 1), 
                                               labels = c("A", "B"))
legend <- cowplot::get_legend(containment_plot)
overall_containment_plot <- plot_grid(raw_containment_plot, legend, rel_widths = c(3, .5))
## approx 10 wide * 3.75 tall seems good

### STILL NEED TO ADD IN THE DETECTION TIME NOT BASED ON CALENDAR TIME AND INSTEAD NUMBER OF CASES

## Sensitivity Analysis - R0 vs Vaccine Efficacy
R0_seq <- R0_scan[-c(1, 2)]
vaccine_efficacy_seq <- seq(0, 1, 0.2)
iterations <- 25
storage_R0_efficacy_sensitivity <- array(data = NA, dim = c(length(R0_seq), 
                                                            length(vaccine_efficacy_infection_seq),
                                                            iterations))
for (i in 1:length(R0_seq)) {
  for (j in 1:length(vaccine_efficacy_infection_seq)) {
    for (k in 1:iterations) {
      bp_out <- chain_sim_susc_ring_vacc(offspring = "pois",
                                         mn_offspring = R0_seq[i],
                                         generation_time = generation_time,
                                         t0 = 0, tf = Inf, pop = 10^7, 
                                         check_final_size = check_final_size, initial_immune = 0,
                                         seeding_cases = 5, prop_asymptomatic = 0.2,
                                         infection_to_onset = infection_to_onset,
                                         vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                         vaccine_efficacy_infection = vaccine_efficacy_seq[j],
                                         vaccine_efficacy_transmission = vaccine_efficacy_seq[j],
                                         vaccine_logistical_delay = 2,
                                         vaccine_protection_delay = 7)
      storage_R0_efficacy_sensitivity[i, j, k] <- nrow(bp_out)
    }
  }
  print(i)
}

storage_R0_efficacy_df <- reshape2::melt(storage_R0_efficacy_sensitivity)
storage_R0_efficacy_df$R0 <- rep(R0_seq, length(vaccine_efficacy_infection_seq) * iterations)
storage_R0_efficacy_df$vaccine_efficacy <- rep(vaccine_efficacy_infection_seq, each = length(R0_seq), times = iterations)
storage_R0_efficacy_df$iteration <- rep(1:iterations, each = length(R0_seq) * length(vaccine_efficacy_infection_seq))
storage_R0_efficacy_df <- storage_R0_efficacy_df %>%
  select(-(Var1:Var3)) %>%
  rename(final_size = value) %>%
  mutate(contained = ifelse(final_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(R0, vaccine_efficacy) %>%
  summarise(proportion_contained = sum(contained) / n())
saveRDS(object = storage_R0_efficacy_df, file = "outputs/branching_process_outputs/branchingProcess_R0_Efficacy_scan.rds")

R0_vaccine_efficacy_plot <- ggplot(storage_R0_efficacy_df, aes(x = R0, y = vaccine_efficacy, fill = proportion_contained)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", name = "Proportion\nContained",
                      limits = c(0, 1)) +
  labs(x = "R0",
       y = "Vaccine Efficacy") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  coord_cartesian(expand = FALSE)

## Sensitivity Analysis - R0 vs Protection Delay
R0_seq <- R0_scan[-c(1, 2)]
vaccine_protection_delay_seq <- seq(0, 14, 2)
iterations <- 25
storage_R0_delay_sensitivity <- array(data = NA, dim = c(length(R0_seq), 
                                                            length(vaccine_protection_delay_seq),
                                                            iterations))
for (i in 1:length(R0_seq)) {
  for (j in 1:length(vaccine_protection_delay_seq)) {
    for (k in 1:iterations) {
      bp_out <- chain_sim_susc_ring_vacc(offspring = "pois",
                                         mn_offspring = R0_seq[i],
                                         generation_time = generation_time,
                                         t0 = 0, tf = Inf, pop = 10^7, 
                                         check_final_size = check_final_size, initial_immune = 0,
                                         seeding_cases = 5, prop_asymptomatic = 0.2,
                                         infection_to_onset = infection_to_onset,
                                         vaccine_start = vaccine_start, vaccine_coverage = vaccine_coverage,
                                         vaccine_efficacy_infection = 0.75,
                                         vaccine_efficacy_transmission = 0.55,
                                         vaccine_logistical_delay = 2,
                                         vaccine_protection_delay = vaccine_protection_delay_seq[j])
      storage_R0_delay_sensitivity[i, j, k] <- nrow(bp_out)
    }
  }
  print(i)
}

storage_R0_delay_df <- reshape2::melt(storage_R0_delay_sensitivity)
storage_R0_delay_df$R0 <- rep(R0_seq, length(vaccine_protection_delay_seq) * iterations)
storage_R0_delay_df$protection_delay <- rep(vaccine_protection_delay_seq, each = length(R0_seq), times = iterations)
storage_R0_delay_df$iteration <- rep(1:iterations, each = length(R0_seq) * length(vaccine_protection_delay_seq))
storage_R0_delay_df <- storage_R0_delay_df %>%
  select(-(Var1:Var3)) %>%
  rename(final_size = value) %>%
  mutate(contained = ifelse(final_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(R0, protection_delay) %>%
  summarise(proportion_contained = sum(contained) / n())
saveRDS(object = storage_R0_delay_df, file = "outputs/branching_process_outputs/branchingProcess_R0_Delay_scan.rds")

R0_protection_delay_plot <- ggplot(storage_R0_delay_df, aes(x = R0, y = protection_delay, fill = proportion_contained)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", name = "Proportion\nContained",
                      limits = c(0, 1)) +
  labs(x = "R0",
       y = "Vaccine Protection Delay") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  coord_cartesian(expand = FALSE)

raw_heatmaps <- plot_grid(R0_vaccine_efficacy_plot + theme(legend.position = "none"), 
                          R0_protection_delay_plot + theme(legend.position = "none"), 
                          rel_widths = c(1, 1), labels = c("C", "D"))
heatmaps_legend <- cowplot::get_legend(R0_vaccine_efficacy_plot)
heatmaps <- plot_grid(raw_heatmaps, heatmaps_legend, rel_widths = c(3, .5))

overall_plot <- cowplot::plot_grid(overall_containment_plot,
                                   heatmaps, nrow = 2, rel_heights = c(1, 1))




