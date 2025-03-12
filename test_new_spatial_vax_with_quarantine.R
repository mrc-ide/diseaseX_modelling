

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

x <- spatial_vax_bp_sim(offspring = c("pois"),   
                        mn_offspring = 3,                      
                        disp_offspring =1,                    
                        spatial_kernel = function(n) { return(rnbinom(n, size = 4, mu = 10)) },                    
                        generation_time = function(n) {return(rep(10, n))},                   
                        prop_asymptomatic = 0.1,                 
                        infection_to_onset = function(n) {return(rep(5, n))},                
                        prob_hosp = 0.5,                         
                        hospitalisation_delay = function(n) {return(rep(10, n))},             
                        vaccine_campaign_radius = 100,           
                        vaccine_coverage = 0.2,                  
                        vaccine_efficacy_infection = 0.75,        
                        vaccine_efficacy_transmission = 0.2,     
                        vaccine_logistical_delay = 1,          
                        vaccine_protection_delay = 1,          
                        vaccine_efficacy_disease = 0.5,          
                        detection_threshold = 1,               
                        t0 = 0, 
                        tf = Inf, 
                        population = 10^6,
                        check_final_size = 10000,
                        initial_immune = 0,
                        seeding_cases = 5,
                        seed = 200)

head(x)

y <- time_to_nth_infection(x, 100)
