# if (!require("pak")) {
#   install.packages("pak")
#   pak::pkg_install("epiverse-trace/bpmodels")
# }
# Load required libraries
library(bpmodels); library(dplyr)

# Sourcing functions
source("functions/run_sars_x.R")

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
  summarise(num_infected = n()) %>%
  rename(id = ancestor)
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

# Get total number of infections arising from infections at time = t (ft removing all infectors for whom all infectees have happened)
num_initial_infections <- 8
gen_select <- chain_sim_eg %>%
  group_by(generation) %>%
  summarise(size = n()) %>%
  filter(size > num_initial_infections) %>%
  filter(generation == min(generation)) %>%
  pull(generation)

id_select_gen <- chain_sim_eg %>%
  filter(generation %in% c(gen_select - 1, gen_select)) %>% 
  left_join(indiv_num_infected, by = "id") %>%
  mutate(num_infected = ifelse(is.na(num_infected), 0, num_infected)) %>%
  arrange(time)

potential_ids <- id_select_gen$id
counter <- 1
selected_ids <- c()
i <- 1
while(length(selected_ids) < num_initial_infections) {
  temp <- potential_ids[i]
  if (i == 1) {
    num_infected <- id_select_gen$num_infected[id_select_gen$id %in% temp]
    infected_ids <- id_select_gen$id[id_select_gen$ancestor %in% temp]
    selected_ids <- c(selected_ids, temp)
  } else {
    ancestor <- id_select_gen$ancestor[id_select_gen$id %in% temp]
    ancestor_num_infected <- id_select_gen$num_infected[id_select_gen$id %in% ancestor]
    ancestor_infected <- id_select_gen$id[id_select_gen$ancestor == ancestor]
    if (sum(ancestor_infected %in% c(selected_ids, temp)) == length(ancestor_infected)) {
      selected_ids <- selected_ids[which(selected_ids != ancestor)]
      selected_ids <- c(selected_ids, temp)
    } else {
      selected_ids <- c(selected_ids, temp)
    }  
  }
  i <- i + 1
}

# for the selected ids, getting their number infected (minus those also in selected_ids)
selected_id_df <- id_select_gen %>%
  filter(id %in% selected_ids) 
num_anc_occurence <- data.frame(ancestor = as.numeric(names(table(selected_id_df$ancestor))), 
                                ancestor_num_infected = indiv_num_infected$num_infected[indiv_num_infected$id %in% as.numeric(names(table(selected_id_df$ancestor)))], 
                                ancestor_count = as.vector(unname(table(selected_id_df$ancestor))))

num_id_occurence <- data.frame(id = as.numeric(names(table(selected_id_df$id))))
num_id_occurence$id_offspring_present <- NA
for (i in 1:length(num_id_occurence$id)) {
  num_id_occurence$id_offspring_present[i] <- sum(selected_id_df$ancestor %in% num_id_occurence$id[i])
}

selected_id_df2 <- selected_id_df %>%
  left_join(num_id_occurence, by = "id") %>%
  mutate(diff = ifelse(id_offspring_present == 0, 0, num_infected - id_offspring_present))

# if some (not all) of infector's offspring have also been selected,
# work out which haven't been selected so we can calculate branch size for them and then
# sum them together later on.
selected_id_df3 <- tibble(selected_id_df2)
selected_id_df3$offspring_not_already_in_df <- NA_real_
for (i in 1:length(selected_id_df3$id)) {
  if (selected_id_df3$diff[i] == 0) {
    selected_id_df3$offspring_not_already_in_df[i] <- NA_real_
  } else {
    temp_infector <- selected_id_df3$id[i]
    temp_infectees <- chain_sim_eg$id[chain_sim_eg$ancestor == temp_infector & !is.na(chain_sim_eg$ancestor)]
    not_present <- selected_id_df3 %>%
      filter(ancestor == temp_infector)
    to_consider <- temp_infectees[!(temp_infectees %in% not_present$id)]
    selected_id_df3$offspring_not_already_in_df[i] <- list(to_consider)
  }
}

individual_transmission_branch_counts <- list()
for (i in 1:length(selected_id_df3$id)) {
  check <- sum(is.na(selected_id_df3$offspring_not_already_in_df[[i]]))
  if (check != 0) {
    cluster_ids <- traverse_tree(selected_id_df3$id[i], chain_sim_eg, c())
    x <- chain_sim_eg %>%
      filter(id %in% cluster_ids) %>%
      mutate(daily = round(time, digits = 0)) %>%
      group_by(daily) %>%
      summarise(incidence = n()) %>%
      tidyr::complete(daily = 0:max(daily), fill = list(incidence = 0)) %>%
      mutate(cumulative = cumsum(incidence))
    individual_transmission_branch_counts[[i]] <- x
  } else {
    ids <- selected_id_df3$offspring_not_already_in_df[[i]]
    for (j in 1:length(ids)) {
      if (j == 1) {
        children <- chain_sim_eg %>% 
          filter(ancestor == ids[j]) %>% 
          pull(id)
        if(is.integer(children) & length(children) == 0) {
          x <- data.frame(daily = 0, incidence = 0, cumulative = 0, chain_branch = j)
        } else {
          cluster_ids <- traverse_tree(ids[j], chain_sim_eg, c())
          x <- chain_sim_eg %>%
            filter(id %in% cluster_ids) %>%
            mutate(daily = round(time, digits = 0)) %>%
            group_by(daily) %>%
            summarise(incidence = n()) %>%
            tidyr::complete(daily = 0:max(daily), fill = list(incidence = 0)) %>%
            mutate(cumulative = cumsum(incidence))  %>%
            mutate(chain_branch = j)
        }
      } else {
        children <- chain_sim_eg %>% 
          filter(ancestor == ids[j]) %>% 
          pull(id)
        cluster_ids <- traverse_tree(ids[j], chain_sim_eg, c())
        if(is.integer(children) & length(children) == 0) {
          y <- data.frame(daily = 0, incidence = 0, cumulative = 0, chain_branch = j)
        } else {
          y <- chain_sim_eg %>%
            filter(id %in% cluster_ids) %>%
            mutate(daily = round(time, digits = 0)) %>%
            group_by(daily) %>%
            summarise(incidence = n()) %>%
            tidyr::complete(daily = 0:max(daily), fill = list(incidence = 0)) %>%
            mutate(cumulative = cumsum(incidence)) %>%
            mutate(chain_branch = j)
          x <- rbind(x, y)
        }
      }
      print(j)
    }
    x <- x %>%
      group_by(daily) %>%
      summarise(incidence = sum(incidence),
                cumulative = sum(cumulative))
  }
  individual_transmission_branch_counts[[i]] <- x
  print(i)
}

colour_func <- scales::hue_pal()(length(individual_transmission_branch_counts))
for (i in 1:length(individual_transmission_branch_counts)) {
  if (i == 1) {
    plot(individual_transmission_branch_counts[[i]]$incidence, type = "l", col = colour_func[i])
  } else {
    lines(individual_transmission_branch_counts[[i]]$incidence, col = colour_func[i])
  }
}
## need to double check I'm counting initial individuals
## need to double check what I do re time (perhaps rejig??) - I'm doing
## a complete call above when defining these so perhaps want to rejig time so that it starts at time of that first seeding case 
## (and is equal to 0)

calc_time_to_cluster_size(individual_transmission_branch_counts, 10)

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

# Get total number of infections arising from the individuals in the selected generation
# selected_generation <- 3
# id_select_gen <- chain_sim_eg %>%
#   filter(generation == selected_generation) %>%
#   pull(id)
# individual_transmission_branch_counts <- list()
# for (i in 1:length(id_select_gen)) {
#   cluster_ids <- traverse_tree(id_select_gen[i], chain_sim_eg, c())
#   x <- chain_sim_eg %>%
#     filter(id %in% cluster_ids) %>%
#     mutate(daily = round(time, digits = 0)) %>%
#     group_by(daily) %>%
#     summarise(incidence = n()) %>%
#     tidyr::complete(daily = 0:max(daily), fill = list(incidence = 0)) %>%
#     mutate(cumulative = cumsum(incidence))
#   individual_transmission_branch_counts[[i]] <- x
#   print(i)
# }

