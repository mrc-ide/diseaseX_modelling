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
spatial_kernel <- function(n) { rnbinom(n, size = 4, mu = 10) }
spatial_ratio_scan <- c(1, 10, 25, 50, 75, 100)

### Vaccine related parameters
vaccine_coverage <- 0.8
vaccine_efficacy_infection_scan <- c(0.3, 0.35, seq(0.4, 0.9, 0.1))
vaccine_efficacy_transmission_scan <- c(0.3, 0.35, seq(0.4, 0.9, 0.1))
vaccine_efficacy_disease <- 0.95
vaccine_logistical_delay <- 2
vaccine_protection_delay <- 7

### Other parameters
pop <- 10^10
check_final_size <- 10000
initial_immune <- 0
seeding_cases <- 3

### Sensitivity analysis parameters
R0_scan <- c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5)
surveillance_scan <- c(1, 10, 25, 50, 75, 100)
iterations <- 100

library(parallel)
no_cores <- 40
cl <- makeCluster(no_cores)
clusterExport(cl, list("mu", "R0_scan", "SC1_generation_time", "spatial_kernel",
                       "check_final_size", "seeding_cases", "SC1_prop_asymptomatic",
                       "SC1_prob_hosp", "SC1_hospitalisation_delay", "surveillance_scan",
                       "vaccine_coverage", "vaccine_efficacy_infection_scan", "vaccine_efficacy_transmission_scan",
                       "vaccine_efficacy_disease", "vaccine_logistical_delay", "vaccine_protection_delay",
                       "spatial_ratio_scan", "spatial_bp_geog_vacc", "spatial_calc",
                       "SC2_generation_time", "SC2_prop_asymptomatic", "SC2_prob_hosp", "SC2_hospitalisation_delay"))

SC1_storage <- array(data = NA, dim = c(iterations, length(R0_scan), length(surveillance_scan), length(spatial_ratio_scan), length(vaccine_efficacy_infection_scan)))
SC2_storage <- array(data = NA, dim = c(iterations, length(R0_scan), length(surveillance_scan), length(spatial_ratio_scan), length(vaccine_efficacy_infection_scan)))
for (i in 4:length(R0_scan)) {
  for (j in 1:length(surveillance_scan)) {
    for (k in 1:length(spatial_ratio_scan)) {
      for (l in 1:length(vaccine_efficacy_infection_scan)) {
        
        clusterExport(cl, list("i", "j", "k", "l"))
        
        # Setup parallel processing for the iterations
        results <- parLapply(cl, 1:iterations, function(m) {
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
                                           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[l],
                                           vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[l],
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
                                           vaccine_efficacy_infection = vaccine_efficacy_infection_scan[l],
                                           vaccine_efficacy_transmission = vaccine_efficacy_transmission_scan[l],
                                           vaccine_efficacy_disease = vaccine_efficacy_disease,
                                           vaccine_logistical_delay = vaccine_logistical_delay,
                                           vaccine_protection_delay = vaccine_protection_delay)
          SC2_count <- sum(!is.na(SC2_temp$time_infection))
          
          list(SC1_count = SC1_count, SC2_count = SC2_count)
        })
        
        # Extract results and store them in the respective storage arrays
        for (m in 1:iterations) {
          SC1_storage[m, i, j, k, l] <- results[[m]]$SC1_count
          SC2_storage[m, i, j, k, l] <- results[[m]]$SC2_count
        }
        print(paste0("i = ", i, ", j = ", j, ", k = ", k, ", l = " , l))
      }
    }
  }
}
stopCluster(cl)
saveRDS(SC1_storage, "outputs/Figure1_branchingProcess_Containment/temp_spatial_vax_sens_SC1_temp.rds")
saveRDS(SC2_storage, "outputs/Figure1_branchingProcess_Containment/temp_spatial_vax_sens_SC2_temp.rds")

SC1_no_vaccination <- data.frame(expand_grid(R0_scan, iterations = 1:iterations), surveillance = 0, spatial_ratio = 0, vaccine_efficacy = 0, outbreak_size = NA_real_, vaccine = "no_vaccine")
SC2_no_vaccination <- data.frame(expand_grid(R0_scan, iterations = 1:iterations), surveillance = 0, spatial_ratio = 0, vaccine_efficacy = 0, outbreak_size = NA_real_, vaccine = "no_vaccine")
counter <- 1
for (i in 1:length(R0_scan)) {
  for (j in 1:iterations) {
    
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
                                     detection_threshold = 1000,
                                     vaccine_campaign_radius = 1,
                                     vaccine_coverage = 0,
                                     vaccine_efficacy_infection = 0,
                                     vaccine_efficacy_transmission = 0,
                                     vaccine_efficacy_disease = 0,
                                     vaccine_logistical_delay = 100,
                                     vaccine_protection_delay = 100)
    SC1_count <- sum(!is.na(SC1_temp$time_infection))
    SC1_no_vaccination[counter, "outbreak_size"] <- SC1_count
    
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
                                     detection_threshold = 1000,
                                     vaccine_campaign_radius = 1,
                                     vaccine_coverage = 0,
                                     vaccine_efficacy_infection = 0,
                                     vaccine_efficacy_transmission = 0,
                                     vaccine_efficacy_disease = 0,
                                     vaccine_logistical_delay = 100,
                                     vaccine_protection_delay = 100)
    SC2_count <- sum(!is.na(SC2_temp$time_infection))
    SC2_no_vaccination[counter, "outbreak_size"] <- SC1_count
    counter <- counter + 1
  }
  print(i)
}

SC1_reshaped <- reshape2::melt(SC1_storage)
colnames(SC1_reshaped) <- c("iteration", "R0", "surveillance", "spatial_ratio", "vaccine_efficacy", "outbreak_size")
SC1_reshaped$vaccine <- "yes_vaccine"
SC1_no_vaccination <- SC1_no_vaccination %>%
  select(iterations, R0_scan, everything()) %>%
  rename(iteration = iterations, R0 = R0_scan)
SC1_no_vaccination$vaccine <- "no_vaccine"
SC1_reshaped2 <- rbind(SC1_reshaped, SC1_no_vaccination)
SC1_reshaped2$pathogen <- "SARS-CoV-1"
saveRDS(SC1_reshaped2, "outputs/Figure1_branchingProcess_Containment/spatial_vax_sens_SC1.rds")

SC2_reshaped <- reshape2::melt(SC2_storage)
colnames(SC2_reshaped) <- c("iteration", "R0", "surveillance", "spatial_ratio", "vaccine_efficacy", "outbreak_size")
SC2_reshaped$vaccine <- "yes_vaccine"
SC2_no_vaccination <- SC2_no_vaccination %>%
  select(iterations, R0_scan, everything()) %>%
  rename(iteration = iterations, R0 = R0_scan)
SC2_no_vaccination$vaccine <- "no_vaccine"
SC2_reshaped2 <- rbind(SC2_reshaped, SC2_no_vaccination)
SC2_reshaped2$pathogen <- "SARS-CoV-2"
saveRDS(SC2_reshaped2, "outputs/Figure1_branchingProcess_Containment/spatial_vax_sens_SC2.rds")

SC1_reshaped <- readRDS("outputs/Figure1_branchingProcess_Containment/spatial_vax_sens_SC1.rds")
SC2_reshaped <- readRDS("outputs/Figure1_branchingProcess_Containment/spatial_vax_sens_SC2.rds")

overall <- rbind(SC1_reshaped, SC2_reshaped) %>%
  mutate(contained = ifelse(outbreak_size < (0.9 * check_final_size), 1, 0)) %>%
  group_by(R0, surveillance, spatial_ratio, vaccine_efficacy, vaccine, pathogen) %>%
  summarise(proportion_contained = sum(contained) / iterations)

### Plotting the output
vaccine_efficacy_index <- which(vaccine_efficacy_infection_scan == 0.35)
spatial_ratio_index <- which(spatial_ratio_scan == 50)
surveillance_scan_index <- which(surveillance_scan == 10)
surveillance_scan_exclude <- which(surveillance_scan == 25)

overall2 <- overall %>%
  filter(vaccine_efficacy == vaccine_efficacy_index | vaccine_efficacy == 0 ,
         spatial_ratio == spatial_ratio_index | spatial_ratio == 0,
         surveillance != surveillance_scan_exclude) %>%
  mutate(R0 = ifelse(vaccine == "no_vaccine", R0, R0_scan[R0]))

fig1c <- ggplot(overall2, aes(x = R0, y = 100 * proportion_contained, col = factor(surveillance))) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_grid(. ~ pathogen,
             labeller = as_labeller(c(`SARS-CoV-1` = "SARS-CoV-1",
                                      `SARS-CoV-2` = "SARS-CoV-2"))) +
  scale_colour_manual(labels = c("No\nVaccine", paste0(surveillance_scan, " Hosp.")),
                      values = c("#474747", c("#E3AFCB", "#D474A4", "#B52F7B", "#9C105A", "#6B0045")),
                      name = "Surveillance\nThreshold\nTrigger") +
  labs(y = "% Outbreaks Contained") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))
ggsave(filename = "figures/Figure_1_BranchingProcess/Fig1HI_SpatialVaxContainmentPlot.pdf", 
       plot = fig1c, 
       width = 7.7, height = 3.1)

R0_surveillance <- overall %>%
  mutate(R0 = ifelse(vaccine == "no_vaccine", R0, R0_scan[R0])) %>%
  filter(pathogen == "SARS-CoV-2" & 
           spatial_ratio == spatial_ratio_index & 
           vaccine_efficacy == vaccine_efficacy_index & 
           R0 > 1 & 
           vaccine == "yes_vaccine")
R0_surveillance_plot <- ggplot(R0_surveillance, aes(x = R0, y = surveillance , fill = 100 * proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 100), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  scale_y_continuous(breaks = 1:length(unique(R0_surveillance$surveillance)), labels = surveillance_scan) +
  scale_x_continuous(breaks = c(1.5, 2.0, 2.5), labels = c(1.5, 2, 2.5)) +
  labs(x = "R0",
       y = "Surveillance Threshold") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE) +
  theme(legend.position = "none")
ggsave(filename = "figures/Figure_1_BranchingProcess/Figure1L_SurvThresh_SpatialVax.pdf", 
       plot = R0_surveillance_plot, width = 2.5, height = 2.34)


R0_spatial <- overall %>%
  mutate(R0 = ifelse(vaccine == "no_vaccine", R0, R0_scan[R0])) %>%
  filter(pathogen == "SARS-CoV-2" & 
           surveillance == surveillance_scan_index & 
           vaccine_efficacy == vaccine_efficacy_index & 
           R0 > 1 & 
           vaccine == "yes_vaccine")
R0_spatial_plot <- ggplot(R0_spatial, aes(x = R0, y = spatial_ratio, fill = 100 * proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 100), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  scale_y_continuous(breaks = 1:length(unique(R0_spatial$spatial_ratio)), labels = spatial_ratio_scan) +
  scale_x_continuous(breaks = c(1.5, 2.0, 2.5), labels = c(1.5, 2, 2.5)) +
  labs(x = "R0",
       y = "Ratio Spatial Vax Radius") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE) +
  theme(legend.position = "none")
ggsave(filename = "figures/Figure_1_BranchingProcess/Figure1J_SpatRatio_SpatialVax.pdf", 
       plot = R0_spatial_plot, width = 2.5, height = 2.34)



vaccine_efficacy_central_index <- which(vaccine_efficacy_infection_scan == 0.35)
R0_efficacy <- overall %>%
  mutate(R0 = ifelse(vaccine == "no_vaccine", R0, R0_scan[R0])) %>%
  filter(pathogen == "SARS-CoV-2" & 
           surveillance == surveillance_scan_index & 
           spatial_ratio == spatial_ratio_index & 
           R0 > 1 & 
           vaccine == "yes_vaccine" &
           vaccine_efficacy != vaccine_efficacy_central_index) %>%
  mutate(vaccine_efficacy = ifelse(vaccine_efficacy > 2, vaccine_efficacy - 1, vaccine_efficacy))
R0_efficacy_plot <- ggplot(R0_efficacy, aes(x = R0, y = vaccine_efficacy, fill = 100 * proportion_contained)) +
  geom_tile(colour = "black") +
  scale_fill_viridis_c(option = "rocket", limits = c(0, 100), begin = 0.175, end = 1, name = "Proportion\nContained",
                       direction = -1) +
  scale_y_continuous(breaks = c(2, 4, 6), 
                     labels = c(40, 60, 80)) +
  scale_x_continuous(breaks = c(1.5, 2.0, 2.5), labels = c(1.5, 2, 2.5)) +
  labs(x = "R0",
       y = "Vaccine Efficacy (%)") +
  theme(axis.text = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 0.5)) +  # Add black border
  coord_cartesian(expand = FALSE) +
  theme(legend.position = "none")
ggsave(filename = "figures/Figure_1_BranchingProcess/Figure1K_Efficacy_SpatialVax.pdf", 
       plot = R0_efficacy_plot, width = 2.5, height = 2.34)


heatmap_legend <- cowplot::get_legend(R0_spatial_plot)
heatmaps <- cowplot::plot_grid(R0_spatial_plot + theme(legend.position = "none"), 
                               R0_efficacy_plot + theme(legend.position = "none"),
                               R0_surveillance_plot + theme(legend.position = "none"), 
                               heatmap_legend,
                               ncol = 4, rel_widths = c(1, 1, 1, 0.3))
cowplot::plot_grid(fig1c, heatmaps, nrow = 2)

ggsave(filename = "figures/Figure_1_BranchingProcess/Fig1H_SpatialVaxContainmentPlot.pdf", 
       plot = fig1c, 
       width = 8, height = 3.1)

legend <- R0_spatial_plot + theme(legend.position = "bottom")
ggsave(file = "figures/Figure_1_BranchingProcess/LegendSpatialVaxHeatmap.pdf", plot = legend, width = 2.4, height = 2.4)


# ### subsetting
# subset <- overall %>%
#   filter(spatial_ratio == 3)
# colours <- c("#E4AFC2", "#D37598", "#C33F6F", "#B31A51", "#810331") # "black" for no vaccination
# colours <- c("#E3AFD6", "#D474B5", "#B52F8D", "#9C1068", "#820366") # "black" for no vaccination
# colours <- c("#E3AFCB", "#D474A4", "#B52F7B", "#9C105A", "#6B0045") # "black" for no vaccination
# 
# fig1_c <- ggplot(subset, aes(x = R0, y = 100 * proportion_contained, col = factor(surveillance))) +
#   geom_line() +
#   geom_point() +
#   scale_x_continuous(breaks = 1:length(R0_scan), labels = R0_scan) +
#   theme_bw() +
#   facet_grid(.~pathogen) + 
#   scale_colour_manual(labels = surveillance_scan,
#                       values = colours,
#                       name = "Surveillance\nThreshold\nTrigger") +
#   labs(x = "R0", y = "% Outbreaks Contained") +
#   theme(strip.background = element_rect(fill = "white"))
