time_to_nth_infection <- function(tdf, n) {
  
  # Round time_infection to the nearest whole number (day)
  daily_summary <- tdf %>%
    mutate(infection_day = round(time_infection)) %>%            # round to nearest day
    group_by(infection_day) %>%
    summarise(incidence = n(), .groups = "drop") %>%              # count infections per day
    complete(infection_day = seq(min(infection_day, na.rm = TRUE), max(infection_day, na.rm = TRUE), by = 1), # include missing days
             fill = list(incidence = 0)) %>%
    arrange(infection_day) %>%
    mutate(cum_infections = cumsum(incidence)) %>%                   # cumulative sum of infections
    filter(!is.na(infection_day))
    
  # Check if the total infections are less than n
  if (max(daily_summary$cum_infections) < n) {
    warning("The total number of infections in the simulation is less than n.")
    return(NA)
  }
  
  # Find the first day where cumulative infections reach or exceed n
  nth_day <- daily_summary$infection_day[which(daily_summary$cum_infections >= n)[1]]
  return(list(nth_day, 
              daily_summary))
}
