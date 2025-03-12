offspring = c("pois") 
mn_offspring = 3
disp_offspring = 1
generation_time = function(n) {return(rep(10, n))}
prop_asymptomatic = 0.1
infection_to_onset = function(n) {return(rep(5, n))}
vaccine_start = 0
vaccine_coverage = 0.2
vaccine_efficacy_infection = 0.75
vaccine_efficacy_transmission = 0.2
vaccine_logistical_delay = 1
vaccine_protection_delay = 1
prob_quarantine = 0
onset_to_quarantine = function(n) {return (rep(1, n))}
quarantine_efficacy = 0.5
t0 = 0
tf = Inf
population = 10^6
check_final_size = 10000
initial_immune = 0
seeding_cases = 5
seed = 200

x <- ring_vax_bp_sim(offspring = offspring,
                     mn_offspring = mn_offspring,
                     disp_offspring = disp_offspring, 
                     generation_time = generation_time,                   
                     prop_asymptomatic = prop_asymptomatic,                 
                     infection_to_onset = infection_to_onset,                
                     vaccine_start = vaccine_start,                     
                     vaccine_coverage = vaccine_coverage,                  
                     vaccine_efficacy_infection = vaccine_efficacy_infection,        
                     vaccine_efficacy_transmission = vaccine_efficacy_transmission,     
                     vaccine_logistical_delay = vaccine_logistical_delay,          
                     vaccine_protection_delay = vaccine_protection_delay,          
                     prob_quarantine = prob_quarantine,                   
                     onset_to_quarantine = onset_to_quarantine,               
                     quarantine_efficacy = quarantine_efficacy,               
                     t0 = 0,                            
                     tf = Inf,                          
                     population = population,                        
                     check_final_size = check_final_size,                  
                     initial_immune = 0,                
                     seeding_cases = seeding_cases,                     
                     seed = seed)

head(x)

y <- time_to_nth_infection(x, 100)
