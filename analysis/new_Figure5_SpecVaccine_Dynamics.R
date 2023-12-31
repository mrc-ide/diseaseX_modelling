# Load required libraries
library(sf); library(ggplot2); library(dplyr); library(rnaturalearth)
library(ggspatial); library(rgdal)
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/helper_functions.R"))

# Loading in bp based detection and calculating detection times for the the different R0 values
bp_df_long <- readRDS("outputs/Figure1_bp_detection_times.rds")
prob_hosp <- squire.page.sarsX:::probs_booster$prob_hosp
arg_pop <- squire::get_population("Argentina")
IHR <- sum(prob_hosp * arg_pop$n / sum(arg_pop$n)) 
num_hosp <- 1:20
detection_hosp <- round(num_hosp / IHR, digits = 0)
num_hosp <- c(1, 5, 10, 20)
infection_thresholds <- detection_hosp[num_hosp]
bp_df_mean_subset <- bp_df_long %>%
  filter(!is.infinite(value),
         detection %in% infection_thresholds) %>%
  group_by(R0, detection, metric) %>%
  summarise(mean = mean(value)) %>%
  filter(detection == 209) %>%
  filter(metric == "Daily Incidence")
R0_detection_time_pairs <- bp_df_mean_subset %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  ungroup() %>%
  select(R0, detection_time, metric)

# NPI Relevant Parameters
default <- define_default_params()
lockdown_Rt <- default$lockdown_Rt                                # Rt achieved under lockdown
minimal_mandate_reduction <- default$minimal_mandate_reduction    # Fold-reduction in R0 achieved under minimal mandate restrictions

## Figure 5 - Varying Disease-Specific Vaccine Development Time (and associated access)
raw_vacc_delay_scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5), specific_vaccine_start = 100 + seq(0, 720, 5))
raw_vacc_delay_scenarios2 <- expand_grid(raw_vacc_delay_scenarios, detection_threshold = unique(bp_df_mean_subset$detection)) %>%
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean) 
NPIs_vacc_delay <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, 
                                         minimal_mandate_reduction = minimal_mandate_reduction, 
                                         NPI_scenarios = c(4, 7, 8), 
                                         scenarios = raw_vacc_delay_scenarios2) 
vacc_delay_scenarios <- raw_vacc_delay_scenarios2 %>%
  full_join(NPIs_vacc_delay, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",
                                    "specific_vaccine_start", "vaccination_rate_bpsv", "vaccination_rate_spec", 
                                    "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all")
final_vacc_delay_scenarios <- vacc_delay_scenarios %>%
  semi_join(R0_detection_time_pairs, by = c("R0", "detection_time", "metric")) %>%
  mutate(main_varied = "specific_development_time")

vars_for_index <- c(variable_columns(vacc_delay_scenarios))
final_vacc_delay_scenarios2 <- final_vacc_delay_scenarios %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
cores <- parallel::detectCores() - 5
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(final_vacc_delay_scenarios2, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/Figure4_VaccineRollout_Characteristics/NEW_Figure_5_specific_vaccination_dev_access_delay_scenarios.rds")
} else {
  model_outputs <- readRDS("outputs/Figure4_VaccineRollout_Characteristics/NEW_Figure_5_specific_vaccination_dev_access_delay_scenarios.rds")
}

## Downstream here I need to create columns for different access timings based on an assumed development time
## E.g. if I set specific dev time to 250, then 300 is a delay to access of 50 days
##      anything below 250 is ignored

## Joining back in the detection metrics
detection_df <- final_vacc_delay_scenarios2 %>%
  select(scenario_index, specific_vaccine_start, detection_threshold, all_of(vars_for_index)) %>%
  filter(vaccine_scenario == "specific_only") %>%
  ungroup() %>%
  select(-vaccine_scenario) %>%
  select(scenario_index, specific_vaccine_start, NPI_int, detection_time, detection_threshold, metric)
model_outputs2 <- model_outputs %>%
  left_join(detection_df, by = c("scenario_index", "specific_vaccine_start", "detection_time", "NPI_int")) %>%
  mutate(detection_threshold_hosp = round(detection_threshold * IHR)) %>%
  mutate(detection_timing = case_when(detection_threshold_hosp == 1 ~ "Early",
                                      detection_threshold_hosp == 5 ~ "Intermediate",
                                      detection_threshold_hosp == 10 ~ "Late",
                                      detection_threshold_hosp == 20 ~ "Very Late")) %>%
  filter(metric == "Daily Incidence") 
labeller_lookup <- c(`4` = "Minimal NPIs", `7` = "Intermediate NPIs", `8` = "Stringent NPIs")


#### Figure Delay to Access

# Get high-quality natural earth data
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(region_un != "Antarctica") %>%
  filter(region_un != "Seven seas (open ocean)") %>%
  filter(continent != "Seven seas (open ocean)")

data <- data.frame(country = world$geounit, value = runif(n = nrow(world))) ## add in BPSV impact here
world_robinson <- st_transform(world, crs = "ESRI:54030")
merged_data <- merge(world_robinson, data, by.x = 'geounit', by.y = 'country')

# Plotting the output 
ggplot() +
  geom_sf(data = merged_data, aes(fill = continent), color = NA, border = NA) +
  coord_sf(crs = st_crs("ESRI:54030")) +
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none") +
  theme_void() +
  scale_fill_brewer(type = "qual", palette = 2, direction = -1) +
  guides(fill = "none")


# Figure 4C - Empirical Delays to Access for COVID-19 
owid_vacc_data <- read.csv("data/owid-covid-data.csv") %>%
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>%
  filter(!is.na(total_vaccinations_per_hundred)) %>%
  filter(date > as.Date("01/11/2020", format = "%d/%m/%Y"))
time_to_one <- owid_vacc_data %>%
  filter(continent != "") %>%
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>%
  group_by(continent, location) %>%
  filter(!is.na(total_vaccinations_per_hundred) & total_vaccinations_per_hundred >= 1) %>%
  summarise(time = min(date)) %>%
  ungroup() %>%
  mutate(delay = as.numeric(time - min(time)))
time_to_one$continent <- factor(time_to_one$continent, levels = rev(levels(factor(time_to_one$continent))))
time_to_one_continent_df <- time_to_one %>%
  group_by(continent) %>%
  summarise(median_delay = median(delay))

access_delay_empirical_boxplot <- ggplot(time_to_one, aes(x = continent, y = delay)) + 
  geom_boxplot(aes(col = continent), position = position_dodge(0.8), outlier.shape = NA, linewidth = 0.75) +
  geom_jitter(aes(fill = continent), position = position_jitterdodge(1.75), size = 2, pch = 21) +
  theme_bw() +
  labs(x = "", y = "Delay to 1% Population Vaccinated (Days)") +
  coord_flip() +
  scale_y_continuous(position = "right", limits = c(0, 300)) + 
  scale_x_discrete(position = "top") + 
  theme(legend.position = "none") +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_fill_brewer(type = "qual", palette = 2) 




delay_access_plot <- ggplot(data = subset(model_outputs2, detection_timing == "Intermediate" & NPI_int %in% c(4, 7, 8) & specific_vaccine_start <= 400)) +
  geom_line(aes(x = specific_vaccine_start - min(model_outputs2$specific_vaccine_start), 
                y = bpsv_deaths_averted * 1000 / unique(model_outputs2$population_size), col = factor(R0)),
            linewidth = 1) +
  scale_colour_manual(values = c("#E6D6CD", "#C9ADA7", "#9A8C98", "#4A4E69", "#22223B")) +
  new_scale("colour") +
  geom_vline(data = time_to_one_continent_df, aes(xintercept = median_delay, col = continent), linewidth = 1) +
  scale_color_brewer(type = "qual", palette = 2) +
  labs(y = "Deaths Averted by Stockpiled BPSV", x = "Delay to Access Disease-Specific Vaccine (Days)") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill="#F5F5F5")) +
  scale_y_continuous(position = "right") +
  facet_wrap(NPI_int ~ ., nrow = 3, labeller = labeller(NPI_int = labeller_lookup), strip.position = "top")

delay_access_combined_plot <- cowplot::plot_grid(access_delay_empirical_boxplot, delay_access_plot, 
                                                 nrow = 2, align = "v", axis = "lr", rel_heights = c(1, 2.5),
                                                 labels = c("C", "D"))
ggsave(filename = "figures/Figure_4_VaccineRollout_Characteristics/Figure_4CD_VaccineAccess_Delay_plot.pdf",
       width = 7, height = 8)

figure4_overall <- cowplot::plot_grid(rate_and_stockpile_plot, delay_access_combined_plot,
                                      nrow = 1, rel_widths = c(0.66, 1))
ggsave(filename = "figures/Figure_4_VaccineRollout_Characteristics/Figure_4_Overall.pdf",
       plot = figure4_overall, 
       width = 9.3, height = 7.5)
# 9.26 * 8.14


##### Scrap ##### 
## Barplots for development time
# development_times <- unique(model_outputs2$specific_vaccine_start)
# 
# delay_100day_development <- unique(model_outputs2$specific_vaccine_start)[sapply(time_to_one_continent_df$median_delay, function(x) {
#   differences <- abs(development_times - (x + 100))
#   return(which.min(differences))
# })]
# time_to_one_continent_df$delay_100 <- delay_100day_development
# 
# delay_220day_development <- unique(model_outputs2$specific_vaccine_start)[sapply(time_to_one_continent_df$median_delay, function(x) {
#   differences <- abs(development_times - (x + 220))
#   return(which.min(differences))
# })]
# time_to_one_continent_df$delay_220 <- delay_220day_development
# time_to_one_continent_df2 <- time_to_one_continent_df %>%
#   pivot_longer(cols = delay_100:delay_220, names_to = "scenario", values_to = "specific_vaccine_start")
# 
# model_outputs_delay <- model_outputs2 %>%
#   filter(specific_vaccine_start %in% unique(c(delay_100day_development, delay_220day_development))) %>%
#   left_join(time_to_one_continent_df2, by = "specific_vaccine_start")
# 
# ggplot(data = subset(model_outputs_delay, R0 == 2)) +
#   geom_bar(aes(x = continent,
#                y = bpsv_deaths_averted * 1000 / unique(model_outputs2$population_size), fill = factor(R0)), stat = "identity") +
#   scale_colour_manual(values = c("#F2E9E4",
#                                  "#C9ADA7",
#                                  "#9A8C98",
#                                  "#4A4E69",
#                                  "#22223B")) +
#   facet_grid(NPI_int ~ scenario, labeller = labeller(NPI_int = labeller_lookup))
# 
# ggplot(data = subset(model_outputs_delay, R0 == 2)) +
#   geom_point(aes(x = composite_NPI_bpsv,
#                y = bpsv_deaths_averted * 1000 / unique(model_outputs2$population_size), 
#                col = scenario), stat = "identity") +
#   facet_grid(NPI_int ~ ., labeller = labeller(NPI_int = labeller_lookup))
# 
# model_outputs_delay2 <- model_outputs_delay %>%
#   select(continent, scenario, NPI_int, bpsv_deaths_averted, R0) %>%
#   filter(continent %in% c("Europe", "Africa")) %>%
#   pivot_wider(names_from = c(continent, scenario), values_from = bpsv_deaths_averted) %>%
#   mutate(afr_eur_100 = Africa_delay_100 - Europe_delay_100,
#          afr_eur_220 = Africa_delay_220 - Europe_delay_220) %>%
#   pivot_longer(cols = afr_eur_100:afr_eur_220, names_to = "scenario", values_to = "difference")
# 
# ggplot(data = model_outputs_delay2) +
#   geom_bar(aes(x = R0,
#                y = difference * 1000 / unique(model_outputs2$population_size), 
#                fill = factor(R0)), stat = "identity") +
#   scale_colour_manual(values = c("#F2E9E4",
#                                  "#C9ADA7",
#                                  "#9A8C98",
#                                  "#4A4E69",
#                                  "#22223B")) +
#   facet_grid(NPI_int ~ scenario, labeller = labeller(NPI_int = labeller_lookup))

















