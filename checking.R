# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))

# NPI Relevant Parameters
lockdown_Rt <- 0.9                   # Rt achieved under lockdown
minimal_mandate_reduction <- 0.25    # Fold-reduction in R0 achieved under minimal mandate restrictions

### BPSV Delay Between Vaccination Receipt & Protection
raw_bpsv_delay_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3, 3.5),                   # Basic reproduction number
                                             IFR = 1,                                       # IFR
                                             population_size = 10^10,
                                             Tg = 5.5,                                      # Tg
                                             detection_time = 14,                           # detection time
                                             bpsv_start = 7,                                # BPSV distribution start (time after detection time)
                                             bpsv_protection_delay = seq(10, 150, 10),        # delay between receipt of BPSV dose and protection
                                             specific_vaccine_start = c(200, 265),          # specific vaccine distribution start (time after detection time)
                                             specific_protection_delay = 7,                 # delay between receipt of specific dose and protection
                                             efficacy_infection_bpsv = 0.35,                # vaccine efficacy against infection - BPSV
                                             efficacy_disease_bpsv = 0.75,                  # vaccine efficacy against disease - BPSV
                                             efficacy_infection_spec = 0.55,                # vaccine efficacy against infection - specific vaccine
                                             efficacy_disease_spec = 0.9,                   # vaccine efficacy against disease - specific vaccine
                                             dur_R = 365000000,                                # duration of infection-induced immunity
                                             dur_bpsv = 365000000,                          # duration of BPSV vaccine immunity
                                             dur_spec = 365000000,                          # duration of disease-specific vaccine immunity
                                             coverage = 0.8,                                # proportion of the population vaccinated
                                             vaccination_rate = 0.035,                      # vaccination rate per week as percentage of population
                                             min_age_group_index_priority = 13,             # index of the youngest age group given priority w.r.t vaccines (13 = 60+)
                                             min_age_group_index_non_priority = 4)          # index of the youngest age group that *receives* vaccines (4 = 15+)
NPIs_bpsv_delay <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                         NPI_scenarios = c(4, 7, 8, 9), scenarios = raw_bpsv_delay_scenarios)
bpsv_delay_scenarios <- raw_bpsv_delay_scenarios %>%
  full_join(NPIs_bpsv_delay, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                    "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")

## Creating index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vars_for_index <- c(variable_columns(bpsv_delay_scenarios), "NPI_int")
bpsv_delay_scenarios <- bpsv_delay_scenarios %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
cores <- parallel::detectCores() - 3
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(bpsv_delay_scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/bpsv_delay_check.rds")
} else {
  model_outputs <- readRDS("outputs/bpsv_delay_check.rds")
}

## Plotting the output
colour_func <- scales::hue_pal()(max(model_outputs$NPI_int))
NPI_colours <- c("#C64191", "#F0803C", "#0D84A9", "#41A07D")
population_size <- unique(model_outputs$population_size)
runtime <- unique(model_outputs$runtime)
NPI_to_include <- c(4, 7, 8, 9)

NPI_df <- NPIs_bpsv_delay %>%
  filter(R0 == 2.5, specific_vaccine_start == 200) %>%
  select(R0, detection_time, bpsv_start, specific_vaccine_start, time_to_coverage_bpsv, time_to_coverage_spec, NPI_int, Rt, tt_Rt) %>%
  rowwise() %>%
  mutate(scenario_info = list(tibble(Rt = Rt, tt_Rt = tt_Rt))) %>%
  select(-Rt, -tt_Rt) %>%
  unnest(cols = c(scenario_info)) %>%
  mutate(scenario = paste0("Scenario ", NPI_int)) %>%
  group_by(scenario) %>%
  mutate(next_time = lead(tt_Rt),
         next_value = lead(Rt)) %>%
  mutate(next_time = ifelse(is.na(next_time), runtime, next_time),
         next_value = ifelse(is.na(next_value), R0, next_value))
overplot_factor <- 1

NPI_plot <- ggplot(NPI_df, aes(x = tt_Rt - overplot_factor, colour = scenario)) +
  geom_hline(aes(yintercept = 1), linewidth = 0.2) +
  geom_hline(aes(yintercept = lockdown_Rt), linetype = "dashed", linewidth = 0.2) +
  geom_hline(aes(yintercept = R0 * (1 - minimal_mandate_reduction)), linetype = "dashed", linewidth = 0.2) +
  geom_vline(aes(xintercept = 0), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time + bpsv_start + time_to_coverage_bpsv), linewidth = 0.2) +
  geom_vline(aes(xintercept = detection_time + specific_vaccine_start + time_to_coverage_spec), linewidth = 0.2) +
  geom_segment(aes(xend = next_time + overplot_factor, y = Rt, yend = Rt), size = 1) +
  geom_segment(aes(x = next_time, xend = next_time, y = Rt, yend = next_value), size = 1) +
  theme_bw() +
  scale_colour_manual(values = NPI_colours) +
  scale_x_continuous(breaks = c(0, unique(NPI_df$detection_time),
                                unique(NPI_df$detection_time) + unique(NPI_df$bpsv_start) + unique(NPI_df$time_to_coverage_bpsv),
                                unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start) + unique(NPI_df$time_to_coverage_spec)),
                     labels = c("", "", "BPSV\nFinish", "Spec\nFinish")) +
  scale_y_continuous(breaks = c(0, 1, unique(NPI_df$R0)),
                     labels = c("", "1", "R0")) +
  facet_wrap(scenario~., nrow = 3,
             labeller = as_labeller(c(`Scenario 4`='Minimal', `Scenario 7`='Moderate', `Scenario 8`='Stringent',
                                      `Scenario 9`='Nothing'))) +
  coord_cartesian(xlim = c(0, unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start) + unique(NPI_df$time_to_coverage_spec) + 10),
                  ylim = c(0.5, unique(NPI_df$R0) + 0.5)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill="#F5F5F5"))

## Delay to protection plot
delay_protect_plotting <- model_outputs %>%
  filter(NPI_int %in% NPI_to_include,
         specific_vaccine_start == specific_vaccine_start_fixed) %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "specific_vaccine_start", "bpsv_protection_delay")))) %>%
  group_by(R0, specific_vaccine_start, bpsv_protection_delay , NPI_int) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted) * 1000 / population_size,
            max_deaths_averted = max(bpsv_deaths_averted) * 1000 / population_size,
            central_deaths_averted = bpsv_deaths_averted * 1000 / population_size,
            perc_deaths_averted = 100 * bpsv_deaths_averted / deaths_spec,
            total_deaths_spec = deaths_spec * 1000 / population_size,
            total_deaths_bpsv = deaths_bpsv * 1000 / population_size,
            time_under_NPIs_bpsv = time_under_NPIs_bpsv,
            composite_NPI_bpsv = composite_NPI_bpsv)

delay_protect_plot <- ggplot(subset(delay_protect_plotting, R0 == 2.5)) +
  geom_line(aes(x = bpsv_protection_delay, y = central_deaths_averted, col = factor(NPI_int)), size = 1) +
  geom_point(aes(x = bpsv_protection_delay, y = central_deaths_averted, fill = factor(NPI_int)), 
             size = 2, pch = 21, col = "black") +
  scale_colour_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                      n_colours = 3))[2],
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3))[2],
                                 rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                      n_colours = 3))[2],
                                 rev(generate_palette(NPI_colours[4], modification = "go_lighter", 
                                                      n_colours = 3))[2]))  +
  scale_fill_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                    n_colours = 3))[2],
                               rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                    n_colours = 3))[2],
                               rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                    n_colours = 3))[2],
                               rev(generate_palette(NPI_colours[4], modification = "go_lighter", 
                                                    n_colours = 3))[2]))  +
  theme_bw() +
  lims(y = c(0, max(subset(delay_protect_plotting, R0 == 2.5)$central_deaths_averted))) +
  labs(x = "Protection Delay (Days)", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none")

ribbon_plotting_delay <- delay_protect_plotting %>%
  filter(R0 != 1.5 & R0 != 3) %>%
  group_by(bpsv_protection_delay, R0) %>%
  summarise(lower = ifelse(min(central_deaths_averted) - 0.2 < 0, 0, min(central_deaths_averted) - 0.2),
            upper = max(central_deaths_averted) + 0.2)

delay_protect_plot2 <- ggplot(subset(delay_protect_plotting, R0 != 1.5 & R0 != 3)) +
  geom_ribbon(data = ribbon_plotting_delay, aes(x = bpsv_protection_delay, ymin = lower, ymax = upper, group = R0), 
              alpha = 0.1, colour = "black", linetype = "dashed") +
  geom_line(aes(x = bpsv_protection_delay, y = central_deaths_averted, 
                col = interaction(factor(R0), factor(NPI_int))), size = 1) +
  geom_point(aes(x = bpsv_protection_delay, y = central_deaths_averted, 
                 fill = interaction(factor(R0), factor(NPI_int))), size = 2, pch = 21, col = "black") +
  scale_colour_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                      n_colours = 3)),
                                 rev(generate_palette(NPI_colours[4], modification = "go_lighter", 
                                                      n_colours = 3)))) +
  scale_fill_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                    n_colours = 3)),
                               rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                    n_colours = 3)),
                               rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                    n_colours = 3)),
                               rev(generate_palette(NPI_colours[4], modification = "go_lighter", 
                                                    n_colours = 3)))) +
  theme_bw() +
  lims(y = c(0, max(subset(delay_protect_plotting, R0 != 1.5)$central_deaths_averted) + 0.2)) +
  labs(x = "Protection Delay (Days)", y = "Deaths Averted By BPSV Per 1000") +
  guides(fill = guide_legend("NPI\nScenario"), colour = "none") +
  theme(legend.position = "none")

### scenarios for close inspection
scen_1 <- bpsv_delay_scenarios %>%
  filter(NPI_int == 9 & bpsv_protection_delay == 20 & R0 == 3.5 & specific_vaccine_start == 200)
scen_1_both <- scen_1[scen_1$vaccine_scenario == "both_vaccines", ]

s1_both <- run_sars_x(population_size = scen_1_both$population_size,
                      country = scen_1_both$country,
                      hosp_bed_capacity = scen_1_both$hosp_bed_capacity,
                      ICU_bed_capacity = scen_1_both$ICU_bed_capacity,
                      Rt = scen_1_both$Rt, 
                      tt_Rt = scen_1_both$tt_Rt, 
                      Tg = scen_1_both$Tg,
                      IFR = scen_1_both$IFR,
                      vaccine_scenario = "both_vaccines", 
                      detection_time = scen_1_both$detection_time,
                      bpsv_start = scen_1_both$bpsv_start,
                      bpsv_protection_delay = scen_1_both$bpsv_protection_delay,
                      specific_vaccine_start = scen_1_both$specific_vaccine_start,             
                      specific_protection_delay = scen_1_both$specific_protection_delay,            
                      efficacy_infection_bpsv = scen_1_both$efficacy_infection_bpsv,           
                      efficacy_disease_bpsv = scen_1_both$efficacy_disease_bpsv,              
                      efficacy_infection_spec = scen_1_both$efficacy_infection_spec,           
                      efficacy_disease_spec = scen_1_both$efficacy_disease_spec,              
                      dur_R = scen_1_both$dur_R,                       
                      dur_bpsv = scen_1_both$dur_bpsv,                    
                      dur_spec = scen_1_both$dur_spec,                   
                      coverage = scen_1_both$coverage,                          
                      vaccination_rate = scen_1_both$vaccination_rate,                  
                      min_age_group_index_priority = scen_1_both$min_age_group_index_priority,        
                      min_age_group_index_non_priority = scen_1_both$min_age_group_index_non_priority,     
                      runtime = 1000, 
                      seeding_cases = scen_1_both$seeding_cases,
                      output = "full", 
                      NPI_int = 0, 
                      scenario_index = 0, 
                      varied = "")

check <- nimue::format(s1_both$model_output, compartments = "D", summaries = "infections") %>%
  filter(t > 1, compartment != "D")
plot(check$t, check$value, xlim = c(0, 200))



