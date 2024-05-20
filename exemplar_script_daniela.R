# Load required libraries
source(here::here("main.R"))

# Load required functions
source(here::here("functions/run_sars_x.R"))
source(here::here("functions/helper_functions.R"))

## Generating default parameters
default <- define_default_params()

# Generate parameter combinations for model running
scenarios <- create_scenarios(R0 = c(1.5, 2.5, 3.5),                      # Basic reproduction number
                              specific_vaccine_start = c(100, 250, 365))  # Specific vaccine distribution start (time after detection time)
# 3 x R0 * 3 x disease specific vaccine start * 2 scenarios (both vaccines and just disease-specific)

## Generating default NPI scenarios (i.e. Rt and tt_Rt for model parameter combinations) and joining to parameter combos
NPIs <- default_NPI_scenarios(lockdown_Rt = default$lockdown_Rt, 
                              minimal_mandate_reduction = default$minimal_mandate_reduction, 
                              NPI_scenarios = 1:9, 
                              scenarios = scenarios)

## Scenarios ready to run
scenarios_NPIs <- scenarios %>%
  full_join(NPIs, by = c("R0", "country", "population_size", "detection_time", "bpsv_start",    # joining by all columns which influence timing of NPI scenarios
                         "specific_vaccine_start", "vaccination_rate_bpsv", "vaccination_rate_spec",
                         "coverage_bpsv", "coverage_spec", "min_age_group_index_priority"), multiple = "all")

## Create index for each model run
final_scenarios <- scenarios_NPIs %>%
  group_by(vaccine_scenario) %>%
  mutate(scenario_index = 1:n())

## Running the model and summarising the output
single_run <- final_scenarios[1, ]
example_single_run <- run_sars_x(population_size = single_run$population_size,
                                 country = single_run$country,
                                 # etc 
                                 # etc
                                 )

## Running all scenarios
plan(multisession, workers = 4) # multicore does nothing on windows as multicore isn't supported
out <- future_pmap(final_scenarios, run_sars_x, .progress = TRUE, .options = furrr_options(seed = 123))
out[[1]]

## Summarising model outputs 
model_outputs <- format_multirun_output(output_list = out, parallel = TRUE, cores = 2)
