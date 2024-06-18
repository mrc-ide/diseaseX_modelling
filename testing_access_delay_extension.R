# Load required libraries
library(sf); library(ggplot2); library(dplyr); library(rnaturalearth)
library(ggspatial) #; library(rgdal)
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/helper_functions.R"))

# Loading in bp based detection and calculating detection times for the the different R0 values
bp_df_long <- readRDS("outputs/Figure1_branchingProcess_Containment/Figure1_bp_detection_times.rds")
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

#### Figure Delay to Access

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

### Figure 4D Delay to Access Specific Development Time
raw_vacc_delay_scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5), specific_vaccine_start = 100 + seq(0, 720, 5), runtime = 1000)
raw_vacc_delay_scenarios2 <- expand_grid(raw_vacc_delay_scenarios, detection_threshold = unique(bp_df_mean_subset$detection)) %>%
  left_join(bp_df_mean_subset, by = c("R0" = "R0", "detection_threshold" = "detection")) %>%
  mutate(detection_time = round(mean, digits = 0)) %>%
  select(-mean)
NPIs_vacc_delay <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt,
                                         minimal_mandate_reduction = minimal_mandate_reduction,
                                         NPI_scenarios = 1:9, #c(4, 7, 8),
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
cores <- parallel::detectCores() - 2
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(final_vacc_delay_scenarios2, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/Figure5_DiseaseSpecific_Dev_Access/updated_NEW_Figure_5_vax_implementation_Delayscenarios.rds")
} else {
  model_outputs <- readRDS("outputs/Figure5_DiseaseSpecific_Dev_Access/updated_NEW_Figure_5_vax_implementation_Delayscenarios.rds")
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

## Access Delay Summary
delay_df <- time_to_one_continent_df %>%
  mutate(delay = 5 * round(median_delay / 5))

dev_start <- 250
delay_plotting <- model_outputs2 %>%
  mutate(delay = specific_vaccine_start - dev_start) %>%
  filter(delay <= dev_start & delay >= 0) %>%
  # filter(NPI_int == 7) %>%
  left_join(delay_df, by = "delay") %>%
  filter(!is.na(continent)) %>%
  mutate(NPI_int = factor(NPI_int))

ggplot(model_outputs2, aes(x = specific_vaccine_start, y = bpsv_deaths_averted * 1000 / unique(model_outputs2$population_size), 
                          col = factor(NPI_int))) +
          geom_line() +
  facet_grid(.~R0)

ggplot(data = delay_plotting) +
  geom_bar(aes(x = factor(continent), y = bpsv_deaths_averted * 1000 / unique(model_outputs2$population_size), 
               fill = continent), position = "dodge", stat = "identity", width = 0.5) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Dark2")) +
  guides(fill = guide_legend(title = "Continent", reverse = TRUE)) +
  facet_wrap(R0 ~ NPI_int, nrow = 3, scales = "free_y")

ggplot(data = delay_plotting) +
  geom_point(aes(x = time_under_NPIs_bpsv , 
               y = bpsv_deaths_averted * 1000 / unique(model_outputs2$population_size), 
               colour = NPI_int)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Dark2")) +
  guides(fill = guide_legend(title = "Continent", reverse = TRUE)) +
  facet_grid(R0 ~ continent)


ggplot(data = delay_plotting) +
  geom_point(aes(x = composite_NPI_bpsv , 
                 y = deaths_bpsv * 1000 / unique(model_outputs2$population_size), 
                 colour = NPI_int)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Dark2")) +
  guides(fill = guide_legend(title = "Continent", reverse = TRUE)) +
  facet_grid(R0 ~ continent)

delay_plotting2 <- model_outputs2 %>%
  mutate(delay = specific_vaccine_start - dev_start) %>%
  filter(delay <= dev_start & delay >= 0) %>%
  # filter(NPI_int == 7) %>%
  left_join(delay_df, by = "delay") %>%
  filter(!is.na(continent)) %>%
  mutate(NPI_int = factor(NPI_int))

ggplot(data = subset(model_outputs2, specific_vaccine_start < 400)) +
  geom_path(aes(x = as.numeric(specific_vaccine_start), 
                 y = deaths_bpsv * 1000 / unique(model_outputs2$population_size))) +
  geom_line(aes(x = as.numeric(specific_vaccine_start),
                y = deaths_spec * 1000 / unique(model_outputs2$population_size))) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Dark2")) +
  # guides(fill = guide_legend(title = "Continent", reverse = TRUE)) +
  facet_grid(NPI_int ~ R0)

## find the difference between BPSV + spec deaths vs just spec deaths
## work out how many more NPI days 

# %>%
  # theme_bw() +
  # lims(y = c(0, 6)) +
  # theme(legend.position = "none") +
  # labs(x = "Basic Reproduction Number (R0)", y = "BPSV Deaths Averted Per 1,000 Population") 

