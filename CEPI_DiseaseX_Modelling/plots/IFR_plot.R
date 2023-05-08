library(readxl)
library(tidyr)
library (dplyr)
library(ggplot2)
library(forcats)
library(wesanderson)

order = c("0 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49",
          "50 to 54", "55 to 59","60 to 64",  "65 to 69" ,"70 to 74" ,"75 to 79", "80+")

order2 <- c ("Severe disease","Hospitalisation", "Death no severe disease","Death severe disease")

IFR_Distributions <- read_excel("data/IFR Distributions.xlsx")

to_plot <- IFR_Distributions %>% select(c("Ages","Prob_hosp", "Prob_severe","Prob_death_ns", "Prob_death_s")) %>%
  rename("Hospitalisation" = Prob_hosp,
         "Severe disease" = Prob_severe,
         "Death severe disease" = Prob_death_s,
         "Death no severe disease" =Prob_death_ns) %>% na.omit() %>%  mutate_at(vars(2:5), as.numeric) %>%
  pivot_longer(- Ages) 


to_plot$Ages <- factor(to_plot$Ages, levels= order)
to_plot$name <- factor(to_plot$name, levels=order2)



ggplot (to_plot) + geom_col(aes(x=Ages, y=value, fill=name), position="dodge", width=1) + theme_classic(base_size = 13) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position= "top") + labs(x = "Age group", y= "Probability", fill ="")+
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4))

