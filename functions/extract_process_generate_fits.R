grab_fit <- function(iso3c, excess_mortality, booster = FALSE){
  
  if (excess_mortality) {
    path <- paste0("https://github.com/mrc-ide/covid-vaccine-impact-orderly/raw/main/data/excess_mortality/model_fits/", iso3c, ".Rds")
  } else {
    path <- paste0("https://github.com/mrc-ide/covid-vaccine-impact-orderly/raw/main/data/reported_deaths/model_fits/", iso3c, ".Rds")
  }
  if (booster) {
    if (excess_mortality) {
      path <- paste0("https://github.com/GBarnsley/booster_model_fits/blob/main/", iso3c, ".Rds?raw=true")
    } else {
      path <- paste0("https://github.com/mrc-ide/nimue_global_fits/raw/main/reported_deaths/", iso3c, ".Rds")
    }
  }
  
  download.file(path, "temp.Rds", mode = "wb", quiet = TRUE)
  fit <- readRDS("temp.Rds")
  unlink("temp.Rds")
  fit
}

simple_Rt <- function (model_out) {
  date_0 <- model_out$inputs$start_date
  iso3c <- squire::get_population(model_out$parameters$country)$iso3c[1]
  return(
    lapply(seq_along(model_out$samples),
           function(y) {
             Rt <- model_out$samples[[y]]$R0
             tt <- list(change = seq_along(Rt), dates = date_0 +
                          model_out$samples[[y]]$tt_R0)
             df <- data.frame(
               Rt = Rt,
               date = date_0 + model_out$samples[[y]]$tt_R0
             ) %>%
               dplyr::mutate(t = as.numeric(.data$date - min(.data$date))) %>%
               dplyr::mutate(iso3c = iso3c, rep = y)
             return(df)
           }))
}

quick_format <- function(x, var_select, date_0) {
  
  d <- nimue:::odin_index(x$model)
  
  do.call(rbind,lapply(var_select, function(i){
    do.call(rbind, lapply(seq_len(dim(x$output)[3]), function(y) {
      df <- data.frame(y = rowSums(x$output[,d[[i]],y]), compartment = i)
      df$t <- seq_len(nrow(df)) - nrow(df)
      df$replicate <- y
      df$date <- df$t + date_0
      return(df)
    }))
  }))
  
}

get_deaths_infections_hosps_time <- function(out){
  value <- quick_format(out, c("D", "infections_cumu"), out$inputs$start_date)
  value$date <- as.Date(rownames(value))
  value <- value %>%
    group_by(replicate, compartment) %>%
    arrange(date) %>%
    filter(y > 0) %>%
    transmute(
      y = c(0, diff(y)),
      date = date,
      replicate = replicate
    ) %>%
    ungroup() %>%
    pivot_wider(names_from = compartment, values_from = y) %>%
    rename(deaths = D, infections = infections_cumu)
  value
}

# Using overall Rt and seeding cases as inputs to re-running squire.page to recreate the deaths curves
seed_infections <- function(squire_model, country, seeding_cases){
  init <- squire.page:::assign_infections(do.call(squire_model$parameter_func, list(country = country)), seeding_cases)
  init_vars <- str_subset(names(init), "_0")
  names(init_vars) <- init_vars
  map(init_vars, function(var, init){
    init[[var]][1:17, 1:6] #have to cut it so that it works with the legacy code from nimue, these dimensions get readded
  }, init = init)
}