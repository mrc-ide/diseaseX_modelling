# Format output
format_out <- function(out, scenarios, target_pop = 50e6){
  # Combine_inputs and outputs
  out1 <- bind_cols(scenarios, bind_rows(out))
  # Isolate counterfactual (Coverage == 0)
  if ("coverage" %in% colnames(out1)) {
    outcf <- filter(out1, coverage == 0) %>%
      select(-coverage) %>%
      rename(output_cf = output,
             output_age_cf = output_age) %>%
      unique()
  } else {
    outcf <- filter(out1, (coverage_children == 0 & coverage_middle == 0 & coverage_old == 0)) %>%
      select(-coverage_children, -coverage_middle, -coverage_old) %>%
      rename(output_cf = output,
             output_age_cf = output_age) %>%
      unique()
  }
  
  # Combine runs and counterfactual and estimate summaries
  summaries <- left_join(out1, outcf)
  
  m <- ncol(summaries)+1
  n <- ncol(summaries)+14
  
  summarise_2021 <- summarise_outputs_age(summaries, p = 2021)
  colnames(summarise_2021)[m:n] <- paste0(colnames(summarise_2021)[m:n], "_2021")
  
  summarise_2022 <- summarise_outputs_age(summaries, p = 2022)
  colnames(summarise_2022)[m:n] <- paste0(colnames(summarise_2022)[m:n], "_2022")
  
  summarise_all_t <- summarise_outputs_age(summaries, p = "all_t")
  colnames(summarise_all_t)[m:n] <- paste0(colnames(summarise_all_t)[m:n], "_all_t")
  
  summaries <- left_join(summaries, select(summarise_2021, -contains("output"))) %>%
    left_join(select(summarise_2022, -contains("output"))) %>%
    left_join(select(summarise_all_t, -contains("output")))
}

# summarise by age and time period
summarise_by_age <- function(x, t_start, t_end, period){
  filter(x, t >= t_start, t < t_end) %>%
    group_by(age_group, compartment) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop_last") %>%
    mutate(period = factor(period))
} 

# Summarise outputs by age over different time period
summarise_outputs_age <- function(x, p) {
  mutate(x, 
         infections = round(map_dbl(output_age, pull_total, outcome = "infections", time_period = p), 2),
         hospitalisations = round(map_dbl(output_age, pull_total, outcome = "hospitalisations", time_period = p), 2),
         deaths = round(map_dbl(output_age, pull_total, outcome = "deaths", time_period = p), 2),
         yll = round(map_dbl(output_age, summarise_yll, time_period = p), 2),
         infections_cf = round(map_dbl(output_age_cf, pull_total, outcome = "infections", time_period = p), 2),
         hospitalisations_cf = round(map_dbl(output_age_cf, pull_total, outcome = "hospitalisations", time_period = p), 2),
         deaths_cf = round(map_dbl(output_age_cf, pull_total, outcome = "deaths", time_period = p), 2),
         yll_cf = round(map_dbl(output_age_cf, summarise_yll, time_period = p), 2),
         infections_averted = infections_cf - infections,
         hospitalisations_averted = hospitalisations_cf - hospitalisations,
         deaths_averted = deaths_cf - deaths,
         deaths_averted_prop = deaths_averted / deaths_cf,
         years_life_saved = yll_cf - yll,
         vaccine_n = round(map_dbl(output_age, pull_total, outcome = "vaccines", time_period = p)))
}



# Estimate total years of life lost
summarise_yll <- function(x, lifespan = 86.6, time_period){
  filter(x, compartment == "deaths", period == time_period) %>%
    mutate(mid_age = (((as.integer(age_group) - 1) * 5) + 2.5),
           yll = pmax(0, (lifespan - mid_age) * value)) %>%
    pull(yll) %>%
    sum()
}


# Pull sum totals
pull_total <- function(x, outcome, time_period){
  filter(x, compartment == outcome, period == time_period) %>%
    pull(value) %>%
    sum()
}
