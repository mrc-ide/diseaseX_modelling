library(dplyr)


# 1. Calculate life expectancy discounted 

discounted_le <- function (x, d= 0.03){
  out <- 0
  x<- round(x,0)
  for (i in 1:x) {out <- out + 1/(1+d)^i}
  return(out)
}

life_expectancy <- read.csv("data/life_expectancy.csv")
life_exp_discounted <- as.matrix(life_expectancy[,2:5])
life_exp_discounted <- data.frame(apply(X = life_exp_discounted,FUN = discounted_le,MARGIN = c(1,2)))
life_exp_discounted$age_group <- life_expectancy$age_group

write.csv(life_exp_discounted, "./data/life_expectancy_discounted.csv")

working_le <- read.csv("./data/working_life_expectacy.csv") 
x<- sapply(X = working_le$life_expectancy ,FUN = discounted_le)
working_le$disc_le <- x

write.csv(working_le, "./data/working_life_expectacy.csv")


#2. Calculate weighted discounted life expectancy per income group 

source("./R/functions_multivaccine_secondary.R")

w_le <- function (income_group, le_discounted ){ 
  
  rep_country <- get_representative_country(income_group = income_group)
  pop <- squire::get_population(country = rep_country)$n
  
  le_vector <- le_discounted %>% 
    select({{income_group}}) 
  
  w_le <- sum(le_vector*pop)/sum(pop)
  
  return (w_le)
}

vsly <- read.csv("./data/vsl.csv") %>% mutate (w_le = 0)

for (i in 1:nrow(vsly)) { vsly$w_le [i] <- w_le(vsly$income_group[i], life_exp_discounted)}

vsly <- vsly %>% mutate( vsly= vsl/w_le)

write.csv(vsly,"./data/vsly.csv")

