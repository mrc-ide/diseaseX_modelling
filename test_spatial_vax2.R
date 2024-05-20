# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/branching_process_spatial_vaccination.R"))

### SC1 parameters
SC1_generation_time <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
SC1_infection_to_onset <- function(n) { rgamma(n, shape = 0.1, rate = 1) } ## Ask Azra for values (negligible assumed currently)
SC1_prop_asymptomatic <- 0
SC1_prob_hosp <- 0.9 # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/
SC1_hospitalisation_delay <- function(n) { 12 } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/

### SC2 parameters
SC2_generation_time <- function(n) { rgamma(n, shape = 13.5, rate = 2) } # 6.75 day generation time Gamam distributed (as per Walker et al, Science, 2020)
SC2_infection_to_onset <- function(n) { rgamma(n, shape = 13.5/3, rate = 2) } # ~35% of transmission presymptomatic (per SARS-CoV-2, slightly lower than but roughly aligned with: https://bmjopen.bmj.com/content/11/6/e041240)
SC2_prop_asymptomatic <- 0.15
SC2_prob_hosp <- 0.05
SC2_hospitalisation_delay <- function(n) { rgamma(n, shape = 24, rate = 2) } # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169193/

### Spatial kernel parameters
mu <- 10
size <- 4
spatial_kernel <- function(n) { rnbinom(n, size = eval(size), mu = eval(mu)) }
spatial_ratio_scan <- c(1, 10, 25, 50, 100)

### Vaccine related parameters
vaccine_coverage <- 0.8
vaccine_efficacy_infection <- 0.35
vaccine_efficacy_transmission <- 0.35
vaccine_efficacy_disease <- 0.95
vaccine_logistical_delay <- 2
vaccine_protection_delay <- 7

### Other parameters
pop <- 10^10
check_final_size <- 4000
initial_immune <- 0
seeding_cases <- 3

### Sensitivity analysis parameters
R0_scan <- c(0.75, 1, 1.25, 1.5, 1.75, 2)
surveillance_scan <- c(1, 10, 25, 50, 100)
iterations <- 100

library(future)
library(future.apply)
plan(multisession, workers = 2, future.seed=TRUE)
SC1_storage <- array(data = NA, dim = c(iterations, length(R0_scan), length(surveillance_scan), length(spatial_ratio_scan)))
SC2_storage <- array(data = NA, dim = c(iterations, length(R0_scan), length(surveillance_scan), length(spatial_ratio_scan)))
### need to re-run this with no vaccination in place
for (i in 1:length(R0_scan)) {
  for (j in 1:length(surveillance_scan)) {
    for (k in 1:length(spatial_ratio_scan)) {
      
      # Setup parallel processing for the iterations
      results <- future_lapply(1:iterations, function(l) {
        vaccination_radius <- spatial_ratio_scan[k] * mu
        
        # SARS-CoV-1 Pathogen Archetype
        SC1_temp <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
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
        SC1_count <- sum(!is.na(SC1_temp$time_infection))
        
        # SARS-CoV-2 Pathogen Archetype
        SC2_temp <- spatial_bp_geog_vacc(mn_offspring = R0_scan[i],
                                         generation_time = SC2_generation_time,
                                         spatial_kernel = spatial_kernel,
                                         t0 = 0, tf = Inf,
                                         check_final_size = check_final_size,
                                         seeding_cases = seeding_cases,
                                         prop_asymptomatic = SC2_prop_asymptomatic,
                                         prob_hosp = SC2_prob_hosp,
                                         hospitalisation_delay = SC2_hospitalisation_delay,
                                         detection_threshold = surveillance_scan[j],
                                         vaccine_campaign_radius = vaccination_radius,
                                         vaccine_coverage = vaccine_coverage,
                                         vaccine_efficacy_infection = vaccine_efficacy_infection,
                                         vaccine_efficacy_transmission = vaccine_efficacy_transmission,
                                         vaccine_efficacy_disease = vaccine_efficacy_disease,
                                         vaccine_logistical_delay = vaccine_logistical_delay,
                                         vaccine_protection_delay = vaccine_protection_delay)
        SC2_count <- sum(!is.na(SC2_temp$time_infection))
        
        list(SC1_count = SC1_count, SC2_count = SC2_count)
      })
      
      # Extract results and store them in the respective storage arrays
      for (l in 1:iterations) {
        SC1_storage[l, i, j, k] <- results[[l]]$SC1_count
        SC2_storage[l, i, j, k] <- results[[l]]$SC2_count
      }
    }
  }
  print(paste0("i = ", i, ", j = ", j, ", k = ", k))
}

SC1_reshaped <- reshape2::melt(SC1_storage)
colnames(SC1_reshaped) <- c("iteration", "R0", "surveillance", "spatial_ratio", "outbreak_size")
SC1_reshaped$pathogen <- "SARS-CoV-1"
saveRDS(SC1_reshaped, "outputs/Figure1_branchingProcess_Containment/spatial_vax_sens_SC1.rds")

SC2_reshaped <- reshape2::melt(SC2_storage)
colnames(SC2_reshaped) <- c("iteration", "R0", "surveillance", "spatial_ratio", "outbreak_size")
SC2_reshaped$pathogen <- "SARS-CoV-2"
saveRDS(SC2_reshaped, "outputs/Figure1_branchingProcess_Containment/spatial_vax_sens_SC2.rds")


SC1_reshaped <- readRDS("outputs/Figure1_branchingProcess_Containment/spatial_vax_sens_SC1.rds")
SC2_reshaped <- readRDS("outputs/Figure1_branchingProcess_Containment/spatial_vax_sens_SC2.rds")

overall <- rbind(SC1_reshaped, SC2_reshaped) %>%
  mutate(contained = ifelse(outbreak_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(R0, surveillance, spatial_ratio, pathogen) %>%
  summarise(proportion_contained = sum(contained) / iterations)

ggplot(overall, aes(x = R0, y = proportion_contained, col = factor(surveillance))) +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_x_continuous(breaks = 1:length(R0_scan), labels = R0_scan) +
  facet_grid(pathogen ~ spatial_ratio,
             labeller = as_labeller(c(`1`= "Radius 1x Average Inf. Dist",
                                      `2`= "Radius 10x Average Inf. Dist",
                                      `3`= "Radius 25x Average Inf. Dist",
                                      `4`= "Radius 50x Average Inf. Dist",
                                      `5`= "Radius 100x Average Inf. Dist",
                                      `SARS-CoV-1` = "SARS-CoV-1",
                                      `SARS-CoV-2` = "SARS-CoV-2"))) +
  scale_colour_manual(labels = surveillance_scan,
                      values = c("#FEC9F1", "#E899DC", "#D387AB", "#B279A7", "#948D9B"),
                      name = "Surveillance\nThreshold\nTrigger")

### subsetting
subset <- overall %>%
  filter(spatial_ratio == 3)
colours <- c("#E4AFC2", "#D37598", "#C33F6F", "#B31A51", "#810331") # "black" for no vaccination
colours <- c("#E3AFD6", "#D474B5", "#B52F8D", "#9C1068", "#820366") # "black" for no vaccination
colours <- c("#E3AFCB", "#D474A4", "#B52F7B", "#9C105A", "#6B0045") # "black" for no vaccination

fig1_c <- ggplot(subset, aes(x = R0, y = 100 * proportion_contained, col = factor(surveillance))) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = 1:length(R0_scan), labels = R0_scan) +
  theme_bw() +
  facet_grid(.~pathogen) + 
  scale_colour_manual(labels = surveillance_scan,
                      values = colours,
                      name = "Surveillance\nThreshold\nTrigger") +
  labs(x = "R0", y = "% Outbreaks Contained") +
  theme(strip.background = element_rect(fill = "white"))

R0_surveillance <- overall %>%
  filter(pathogen == "SARS-CoV-2" & spatial_ratio == 5 & R0 > 2)
R0_surveillance_plot <- ggplot(R0_surveillance, aes(x = R0, y = surveillance , fill = 100 * proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 100), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  scale_y_continuous(breaks = 1:5, labels = surveillance_scan) +
  scale_x_continuous(breaks = 1:6, labels = R0_scan) +
  labs(x = "R0",
       y = "Surveillance Threshold") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

R0_spatial <- overall %>%
  filter(pathogen == "SARS-CoV-2" & surveillance == 2 & R0 > 2)
R0_spatial_plot <- ggplot(R0_spatial, aes(x = R0, y = spatial_ratio, fill = 100 * proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 100), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  scale_y_continuous(breaks = 1:5, labels = spatial_ratio_scan) +
  scale_x_continuous(breaks = 1:6, labels = R0_scan) +
  labs(x = "R0",
       y = "Surveillance Threshold") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE)

heatmap_legend <- cowplot::get_legend(R0_spatial_plot)
heatmaps <- cowplot::plot_grid(R0_surveillance_plot + theme(legend.position = "none"), 
                               R0_spatial_plot + theme(legend.position = "none"), 
                               heatmap_legend,
                               ncol = 3, rel_widths = c(1, 1, 0.3))
cowplot::plot_grid(fig1_c, heatmaps, nrow = 2)


SC1_no_vaccination <- data.frame(expand_grid(R0_scan, iterations = 1:iterations), surveillance = 0, spatial_ratio = 0, vaccine_efficacy = 0, outbreak_size = NA_real_)

