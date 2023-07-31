# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))

# NPI Relevant Parameters
lockdown_Rt <- 0.9                   # Rt achieved under lockdown
minimal_mandate_reduction <- 0.25    # Fold-reduction in R0 achieved under minimal mandate restrictions

# Generate parameter combinations for model running
raw_R0_spec_start_scenarios <- create_scenarios(R0 = c(1.5, 2, 2.5, 3, 3.5),                   
                                                IFR = 1,                          
                                                population_size = 10^10,
                                                Tg = 5.5,                                      
                                                detection_time = 14,                           
                                                bpsv_start = 7,                                
                                                bpsv_protection_delay = 7,                     
                                                specific_vaccine_start = seq(100, 730, 10),
                                                specific_protection_delay = 7,                 
                                                efficacy_infection_bpsv = 0.35,                
                                                efficacy_disease_bpsv = 0.75,    
                                                efficacy_infection_spec = 0.55,                
                                                efficacy_disease_spec = 0.9,                   
                                                dur_R = 365000000,                             
                                                dur_bpsv = 365000000,                          
                                                dur_spec = 365000000,                          
                                                coverage = 0.8,                                
                                                vaccination_rate = 0.035,                      
                                                min_age_group_index_priority = 13,             
                                                min_age_group_index_non_priority = 4,
                                                runtime = 1000)          
NPIs_R0_spec <- default_NPI_scenarios(lockdown_Rt = lockdown_Rt, minimal_mandate_reduction = minimal_mandate_reduction, 
                                      NPI_scenarios = c(4, 7, 8), scenarios = raw_R0_spec_start_scenarios)
R0_spec_start_scenarios <- raw_R0_spec_start_scenarios %>%
  full_join(NPIs_R0_spec, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence NPI scenario timing
                                  "specific_vaccine_start", "vaccination_rate", "coverage", "min_age_group_index_priority"), multiple = "all")

## Creating index for output (important as it orders dataframe so that pairs of identical scenarios save for BPSV Y/N are next to each other)
vars_for_index <- c(variable_columns(R0_spec_start_scenarios), "NPI_int")
scenarios <- R0_spec_start_scenarios %>%
  group_by(vaccine_scenario) %>%
  arrange_at(vars_for_index) %>%
  mutate(scenario_index = 1:n())

## Generating/Loading Outputs as Required
cores <- parallel::detectCores() - 3
fresh_run <- FALSE
if (fresh_run) {
  plan(multisession, workers = cores) # multicore does nothing on windows as multicore isn't supported
  system.time({out <- future_pmap(scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))})
  model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = cores)
  saveRDS(model_outputs, "outputs/R0_specific_start_bivariate.rds")
} else {
  model_outputs <- readRDS("outputs/R0_specific_start_bivariate.rds")
}

## Processing raw model outputs
proc_model_outputs <- model_outputs %>%
  filter(map_lgl(varied, ~ setequal(., c("R0", "specific_vaccine_start")))) %>%
  group_by(R0, specific_vaccine_start, efficacy_disease_bpsv, NPI_int) %>%
  summarise(min_deaths_averted = min(bpsv_deaths_averted) * 1000 / population_size,
            max_deaths_averted = max(bpsv_deaths_averted) * 1000 / population_size,
            central_deaths_averted = bpsv_deaths_averted * 1000 / population_size,
            perc_deaths_averted = 100 * bpsv_deaths_averted / deaths_spec,
            total_deaths_spec = deaths_spec * 1000 / population_size,
            total_deaths_bpsv = deaths_bpsv * 1000 / population_size,
            time_under_NPIs_bpsv = time_under_NPIs_bpsv,
            composite_NPI_bpsv = composite_NPI_bpsv)

## Plotting the output
colour_func <- scales::hue_pal()(max(model_outputs$NPI_int))
NPI_colours <- colour_func[unique(model_outputs$NPI_int)]
population_size <- unique(model_outputs$population_size)
runtime <- unique(model_outputs$runtime)
NPI_colours <- c("#C64191", "#F0803C", "#0D84A9")

## NPI Plots
NPI_df <- NPIs_R0_spec %>%
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
                     labels = c("", "Detect", "BPSV\nFinish", "Spec\nFinish")) +
  scale_y_continuous(breaks = c(0, 1, unique(NPI_df$R0)),
                     labels = c("", "1", "R0")) +
  facet_wrap(scenario~., nrow = 3,
             labeller = as_labeller(c(`Scenario 4`='Minimal', `Scenario 7`='Moderate', `Scenario 8`='Stringent'))) +
  labs(x = "Time Since Spillover") +
  coord_cartesian(xlim = c(0, unique(NPI_df$detection_time) + unique(NPI_df$specific_vaccine_start) + unique(NPI_df$time_to_coverage_spec) + 10),
                  ylim = c(0.5, unique(NPI_df$R0) + 0.5)) +
  theme(legend.position = "none",
        strip.background = element_rect(fill="#F5F5F5"))

## 1st Option for R0 vs Specific Vaccine Time Plot
pal <- scales::viridis_pal(alpha = 1, begin = 0.15, end = 1, direction = -1, option = "magma")
x <- ggplot(subset(proc_model_outputs, R0 %in% c(1.5, 2, 2.5, 3, 3.5))) +
  geom_line(aes(x = specific_vaccine_start, y = central_deaths_averted, 
                col = interaction(factor(R0), factor(NPI_int))), size = 1) +
  scale_colour_manual(values = c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                                                      n_colours = 5, view_palette = TRUE)),
                                 rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                                                      n_colours = 5, view_palette = TRUE)),
                                 rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                                                      n_colours = 5, view_palette = TRUE))))  +
  scale_y_continuous(position = "right") +
  scale_x_continuous(breaks = seq(100, 700, 100)) +
  facet_wrap(. ~ NPI_int, nrow = 3) +
  geom_vline(aes(xintercept = 100), linetype = "dashed") +
  geom_vline(aes(xintercept = 365), linetype = "dashed") +
  theme_bw() +
  labs(x = "Time to Disease-Specific\nVaccine", y = "Additional Deaths Averted By BPSV Per 1000") +
  guides(colour = guide_legend("R0")) +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(colour = "white"))

first_option <- cowplot::plot_grid(NPI_plot, x, rel_widths = c(1, 2), labels = c('A','B'), axis = "h", align = "bt")
ggsave(filename = "figures/R0_specStart_V1.pdf", plot = first_option, width = 11, height = 8)

## 2nd Option for R0 vs Specific Vaccine Time Plot
cols <- c(rev(generate_palette(NPI_colours[1], modification = "go_lighter", 
                               n_colours = 3, view_palette = TRUE)),
          rev(generate_palette(NPI_colours[2], modification = "go_lighter", 
                               n_colours = 3, view_palette = TRUE)),
          rev(generate_palette(NPI_colours[3], modification = "go_lighter", 
                               n_colours = 3, view_palette = TRUE)))
x <- ggplot(subset(proc_model_outputs, R0 %in% c(1.5, 2.5, 3.5))) +
  geom_line(aes(x = specific_vaccine_start, y = central_deaths_averted, 
                col = interaction(factor(NPI_int), factor(R0))), size = 1) +
  scale_colour_manual(values = c(cols[1], cols[4], cols[7],
                                 cols[2], cols[5], cols[8],
                                 cols[3], cols[6], cols[9]))  +
  scale_y_continuous(position = "right") +
  theme_bw() +
  labs(x = "Time to Disease-Specific\nVaccine", y = "Additional Deaths Averted By BPSV Per 1000") +
  guides(colour = guide_legend("R0")) +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_blank())

second_option <- cowplot::plot_grid(NPI_plot, x, rel_widths = c(1, 2), labels = c('A','B'), 
                                    axis = "h", align = "bt")
ggsave(filename = "figures/R0_specStart_V2.pdf", plot = second_option, width = 9.5, height = 7.5)
