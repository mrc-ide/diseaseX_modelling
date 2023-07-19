# if (!require("pak")) {
#   install.packages("pak")
#   pak::pkg_install("epiverse-trace/bpmodels")
# }
# Load required libraries
library(bpmodels); library(dplyr)

# Set seed and run branching process with susceptible depletion up until day 50
set.seed(457)
chain_sim_eg <- chain_sim_susc(offspring = "pois", 
                               mn_offspring = 2,
                               t0 = 0, 
                               tf = 50,
                               serial = function(n) {
                                 rgamma(n, shape = 5, rate = 1)},
                               pop = 10^5,
                               initial_immune = 0)

# Get number of people infected by each infector
indiv_num_infected <- chain_sim_eg %>%
  filter(!is.na(ancestor)) %>%
  group_by(ancestor) %>%
  summarise(num_infected = n())
mean(indiv_num_infected$num_infected) ## double checking approximately matches the mean above

# Get overall incidence
incidence <- chain_sim_eg %>%
  mutate(daily = round(time, digits = 0)) %>%
  group_by(daily) %>%
  summarise(incidence = n()) %>%
  tidyr::complete(daily = min(daily):max(daily), fill = list(incidence = 0))
plot(incidence$daily, incidence$incidence)

# Get number of individuals in each generation
num_gen <- chain_sim_eg %>%
  group_by(generation) %>%
  summarise(count = n())

# Get total number of infections arising from the individuals in the selected generation
selected_generation <- 3
id_select_gen <- chain_sim_eg %>%
  filter(generation == selected_generation) %>%
  pull(id)
individual_transmission_branch_counts <- list()
for (i in 1:length(id_select_gen)) {
  x <- chain_sim_eg %>%
    filter(id %in% traverse_tree(id_select_gen[i], df, c())) %>%
    mutate(daily = round(time, digits = 0)) %>%
    group_by(daily) %>%
    summarise(incidence = n()) %>%
    tidyr::complete(daily = 0:max(daily), fill = list(incidence = 0)) %>%
    mutate(cumulative = cumsum(incidence))
  individual_transmission_branch_counts[[i]] <- x
  print(i)
}

sum(individual_transmission_branch_counts[[1]]$incidence) + 
  sum(individual_transmission_branch_counts[[2]]$incidence) + 
  sum(individual_transmission_branch_counts[[3]]$incidence) + 
  sum(individual_transmission_branch_counts[[4]]$incidence) ## total = 2596 (total number of people in sim, minus the 4 
                                                            ##         individuals in starting gen, minus individuals in
                                                            ##         previous gen - 2 in this case)

plot(individual_transmission_branch_counts[[1]]$incidence, type = "l")
lines(individual_transmission_branch_counts[[2]]$incidence, col = "red")
lines(individual_transmission_branch_counts[[3]]$incidence, col = "blue")
lines(individual_transmission_branch_counts[[4]]$incidence, col = "purple")

plot(individual_transmission_branch_counts[[1]]$cumulative, type = "l")
lines(individual_transmission_branch_counts[[2]]$cumulative, col = "red")
lines(individual_transmission_branch_counts[[3]]$cumulative, col = "blue")
lines(individual_transmission_branch_counts[[4]]$cumulative, col = "purple")

