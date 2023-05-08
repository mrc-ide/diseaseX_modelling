# Load required libraries
library(odin); library(ggpubr)
source(here::here("main.R"))


# Simple ODE Model to Simulate Waning Efficacy Over Time
### This matches the model currently used in squire.page: 
###  https://github.com/mrc-ide/squire.page/blob/main/inst/odin/nimue_booster.R#L20-L30
calc_eff_gen <- odin({
  initial(C1) <- 1
  initial(C2) <- 0
  initial(C3) <- 0
  deriv(C1) <- -w_1*C1
  deriv(C2) <- w_1*C1 - w_2*C2
  deriv(C3) <- w_2*C2
  output(ve_d) <- (C1 * ved + C2 * ved_2 + C3 * ved_3)
  output(ve_i) <- (C1 * vei + C2 * vei_2 + C3 * vei_3)
  w_1 <- user()
  w_2 <- user()
  ved <- user()
  ved_2 <- user()
  ved_3 <- user()
  vei <- user()
  vei_2 <- user()
  vei_3 <- user()
})
set.seed(1000101)
n_samples <- 50
simulate_time <- 2*365

t_plot <- seq(0, simulate_time, length.out = 100)
calc_eff <- calc_eff_gen$new(  user = list(
  w_1 = 0.011889010,
  w_2 = 0.001492792,
  ved = 0.99,
  ved_2 = 0.95,
  ved_3 = 0.654333312,
  vei = 0.996074975,
  vei_2 = 0.522821603,
  vei_3 = 0.279789215))

### Note that as far as my current understanding goes, when someone receives first 
### dose, they're automatically protected (i.e. there's no delay between receiving vaccination)
### and getting protection) - the same applies when they're given second dose, booster etc.
mod_value <- calc_eff$run(t = t_plot)
second_dose_timing <- 90
first_dose_efficacy_infection <- 0.65
first_dose_efficacy_hospitalisation <- 0.80
output <- tibble(
  t = rep(c(0:(second_dose_timing-1), t_plot + second_dose_timing), 2),
  value = c(c(rep(first_dose_efficacy_hospitalisation, second_dose_timing), mod_value[, "ve_d"]),
            c(rep(first_dose_efficacy_infection, second_dose_timing), mod_value[, "ve_i"])),
  endpoint = c(rep("Hospitalisation", length(t_plot) + second_dose_timing), 
               rep("Infection", length(t_plot) + second_dose_timing)))

ggplot(data = output, aes(x = t, y = value, color = endpoint)) +
  geom_line() +
  labs(y = "Vaccine Efficacy", x = "Days Since First Dose") +
  theme_bw() +
  ylim(c(0, 1))

# ### Note: much of this code is taken from: 
# ###       https://github.com/mrc-ide/global-lmic-reports-orderly/blob/minimal_fitting/src/parameters_vaccines/script.R
# simulate_ab <- function(t, initial_ab, h_s, h_l, t_s) {
#   pi1 <- -log(2)/h_s
#   pi2 <- -log(2)/h_l
#   initial_ab * (
#     (exp(pi1 * t + pi2 * t_s) + exp(pi1 * t_s + pi2 * t))/
#       (exp(pi1 * t_s) + exp(pi2 * t_s))
#   )
# }
# ab_to_ve <- function(ab, n50, k){
#   1/(1 + exp(-k * (log10(ab) - n50)))
# }
# ab_params <- list(
#   ni50 = log10(0.052),
#   ns50 = log10(0.01),
#   k = 2.8,
#   hl_s = 33,
#   hl_l = 580,
#   period_s = 86
# )
# t_check <- seq(1:500)
# ab <- simulate_ab(t = t_check, initial_ab = 0.1, h_s = ab_params$hl_s, h_l = ab_params$hl_l, t_s = ab_params$period_s)
# plot(ab_to_ve(ab, ab_params$ni50, ab_params$k)) # protection against infection
# plot(ab_to_ve(ab, ab_params$ns50, ab_params$k)) # protection against death
