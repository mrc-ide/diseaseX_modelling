
df_main <- readRDS("output/4_secondary_mv_scenarios.rds") 

#df_main <- out_combine

# Adding cumulative values to the output df. 

cum <- function(df){
  
  out <- df %>% group_by (compartment) %>% mutate (cum= cumsum(value))
  return (out)
}
days_above_capacity <- function (out, income_group,target_pop)
{
  country  <- get_representative_country(income_group)
  capacity <- squire::get_healthcare_capacity(country=country)$hosp_beds *target_pop/1000
  
  r <- out %>% filter (compartment == "hospital_occupancy") %>% summarise (days_above = length(value[value>capacity]))
  
  return(as.numeric(r$days_above))
  
}

df_main <- df_main %>%
  rowwise() %>%
  mutate(output = list(cum(output)),
         output_cf= list(cum(output_cf)),
         day_above_capacity = days_above_capacity(output, income_group,target_pop),
         day_above_capacity_cf = days_above_capacity(output_cf, income_group,target_pop))


target_pop <- 1e6


# Time series plots

income_to_plot = "LMIC"

# limit to an illustrative scenario
df1 <- df_main %>%
  filter(income_group == income_to_plot & ifr_scaling == 1) %>%
  filter(timing1a == 50 | timing1a == 110 | timing1a == 170) %>%
  unnest(cols = output) %>%
  filter(t>1) %>%
  mutate(importation_time_label = case_when ((timing1a == 50) ~ "50 days",
                                             (timing1a == 110) ~ "110 days",
                                             (timing1a == 170) ~ "170 days"
  )) %>%
  mutate(importation_time_label = factor(importation_time_label, levels=c("50 days", "110 days", "170 days"))) %>%
  mutate(vaccine_start_label = case_when( (vaccine_1_start < runtime) ~ "Broadly protective sarbecovirus vaccine",
                                          (vaccine_1_start == runtime & vaccine_2_start < runtime) ~ "SARS-X vaccine",
                                          (vaccine_1_start == runtime & vaccine_2_start == runtime) ~ "No vaccine",
  )) %>%
  mutate(importation_time_label = factor(importation_time_label)) %>%
  mutate(vaccine_start_label = factor(vaccine_start_label)) %>%
  mutate(vaccine_start_label = fct_relevel(vaccine_start_label, "Broadly protective sarbecovirus vaccine", 
                                           "SARS-X vaccine", "No vaccine")) %>%
  filter(Rt1a==R0) ##temporarily in - remove once re-run batch

df_deaths <- df1 %>%
  filter(compartment == "deaths") %>%
  filter(R0==1.5) %>%
  filter((vaccine_start_label=="No vaccine" & two_vaccines==0) | vaccine_start_label == "SARS-X vaccine" 
         | vaccine_start_label == "Broadly protective sarbecovirus vaccine")

df_rec<- df_deaths %>% filter(t==100)


g1a <- ggplot(data = df_deaths) +
  geom_line(aes(x = t, y = (value),col=vaccine_start_label),linewidth=1.1) +
  geom_rect(aes(xmin=timing1b, xmax=timing2, ymin=0,ymax=Inf) ,alpha=0.8, fill="#d1e0e0") +
  facet_grid(~importation_time_label, scales="free") +
  #  facet_wrap(vars(detection_time), labeller = as_labeller(detection_time_label)) +
  labs(x = "Time (days)", y = "Deaths per million per day", col = "Vaccine Strategy", linetype = "") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = cols2) +
  scale_linetype_manual(values = c("dashed", "dotted")) 
g1a
ggsave("plots/Secondary_Time_Series_Deaths_R0_1.5.png",g1a, height=3, width=10)

df_deaths <- df1 %>%
  filter(compartment == "deaths") %>%
  filter(R0==2.5) %>%
  filter((vaccine_start_label=="No vaccine" & two_vaccines==0) | vaccine_start_label == "SARS-X vaccine" 
         | vaccine_start_label == "Broadly protective sarbecovirus vaccine")

g1b <- ggplot(data = df_deaths) +
  geom_rect(data = df_rec %>% filter(vaccine_start_label == "SARS-X vaccine"),aes(xmin=min(timing1b,timing2), xmax=timing2, ymin=90,ymax=93), fill="#bbdf27", alpha=0.6) +
  geom_rect(data = df_rec %>% filter(vaccine_start_label == "No vaccine") ,aes(xmin=min(timing1b,timing2), xmax=timing2, ymin=87,ymax=90),alpha=0.6,fill="grey60") +
  geom_rect(data = df_rec %>% filter(vaccine_start_label == "Broadly protective sarbecovirus vaccine"),aes(xmin=min(timing1b,timing2), xmax=timing2, ymin=84,ymax=87), fill="#25848e",alpha=0.6) +
  geom_line(aes(x = t, y = (value),col=vaccine_start_label),linewidth=1.1) +
  facet_grid(~importation_time_label, scales="free") +
  #  facet_wrap(vars(detection_time), labeller = as_labeller(detection_time_label)) +
  labs(x = "Time (days)", y = "Deaths per million per day", col = "Vaccine Strategy", linetype = "") +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = cols2) +
  scale_linetype_manual(values = c("dashed", "dotted"))

g1b
ggsave("plots/Secondary_Time_Series_Deaths_R0_2.5.png",g1b, height=3, width=10)

df_deaths <- df1 %>%
  filter(compartment == "deaths") %>%
  filter(R0==4) %>%
  filter((vaccine_start_label=="No vaccine" & two_vaccines==0) | vaccine_start_label == "SARS-X vaccine" 
         | vaccine_start_label == "Broadly protective sarbecovirus vaccine")

g1c <- ggplot() +
  geom_line(data = df_deaths, aes(x = t, y = (value),col=vaccine_start_label),linewidth=1.2) +
  facet_wrap(vars(importation_time_label)) +
  #  facet_wrap(vars(detection_time), labeller = as_labeller(detection_time_label)) +
  labs(x = "Time (days)", y = "Deaths per million per day", col = "Vaccine Strategy", linetype = "") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = cols2) +
  scale_linetype_manual(values = c("dashed", "dotted")) 
g1c
ggsave("plots/Secondary_Time_Series_Deaths_R0_4.png",g1c, height=3, width=10)

# Vaccinations plot 

df_v1 <- df1 %>%
  filter(timing1a==50) %>%
  filter(compartment == "V1") %>%
  #  filter(t<=250) %>%
  filter(R0 == 2.5) %>%
  group_by(vaccine_start_label)


df_v2 <- df1 %>%
  filter(timing1a==50) %>%
  filter(compartment == "V2") %>%
  #  filter(t<=250) %>%
  filter(R0 == 2.5) %>%
  group_by(vaccine_start_label)


g2 <- ggplot() +
  geom_line(data = df_v1, aes(x = t, y = value/target_pop*100, linetype = "Sarbecovirus", col=vaccine_start_label),linewidth=1.2) +
  geom_line(data = df_v2, aes(x = t, y = value/target_pop*100, linetype = "SARS-X", col=vaccine_start_label),linewidth=1.2) +
  scale_linetype_manual(name = "Vaccine", values = c("solid", "dashed")) +
  scale_y_continuous(limits=c(0,100), breaks = c(0,20,40,60,80,100) ) +
  scale_color_manual(name = "Vaccine Strategy", values = cols2a) + 
  labs(x = "Time (days)", y = "Vaccine coverage (%)", col = "Vaccine Strategy", linetype = "") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) 
g2
ggsave("plots/Secondary_Vaccine_Rollout_emergence.png",g2, height=3, width=6)



# Hospitalisation plot 


rep_country <- get_representative_country(income_group = income_to_plot)
capacity <- squire::get_healthcare_capacity(country=rep_country)$hosp_beds


df_hosp <- df1 %>%
  filter(compartment == "hospital_occupancy") %>%
  filter(R0==2.5) %>%
  filter((vaccine_start_label=="No vaccine" & two_vaccines==0) | vaccine_start_label == "SARS-X vaccine" 
         | vaccine_start_label == "Broadly protective sarbecovirus vaccine")



g3b <- ggplot(data = df_hosp) +
  geom_rect(data = df_rec %>% filter(vaccine_start_label == "SARS-X vaccine"),aes(xmin=min(timing1b,timing2), xmax=timing2, ymin=4500,ymax=4400), fill="#bbdf27", alpha=0.6) +
  geom_rect(data = df_rec %>% filter(vaccine_start_label == "No vaccine") ,aes(xmin=min(timing1b,timing2), xmax=timing2, ymin=4400,ymax=4300),alpha=0.6,fill="grey60") +
  geom_rect(data = df_rec %>% filter(vaccine_start_label == "Broadly protective sarbecovirus vaccine"),aes(xmin=min(timing1b,timing2), xmax=timing2, ymin=4300,ymax=4200), fill="#25848e",alpha=0.6) +
  geom_line(aes(x = t, y = (value),col=vaccine_start_label),linewidth=1.1) +
  geom_hline(yintercept = capacity*1000)+
  facet_grid(~importation_time_label, scales="free") +
  #  facet_wrap(vars(detection_time), labeller = as_labeller(detection_time_label)) +
  labs(x = "Time (days)", y = "Hospital occupancy per million per day", col = "Vaccine Strategy", linetype = "") +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = cols2) +
  scale_linetype_manual(values = c("dashed", "dotted")) 
g3b
ggsave("plots/Secondary_Time_Series_Hosp_R0_2.5.png",g1b, height=3, width=10)


days_above_capacity <- df_hosp %>% group_by (vaccine_start_label,importation_time_label) %>% 
  select (vaccine_start_label,importation_time_label,day_above_capacity,day_above_capacity_cf) %>% distinct()%>%
  flextable()

save_as_docx("Secondary" = days_above_capacity,
             path = "tables/Secondary_Hospital_Capacity_Table.docx")

##############################################################
# Summary plots

# remove time series for faster processing
df2 <- select(df_main,-contains("output")) %>%
  filter(timing1a == 50 | timing1a == 110 | timing1a == 170) %>%
  mutate(importation_time_label = case_when ((timing1a == 50) ~ "50 days",
                                             (timing1a == 110) ~ "110 days",
                                             (timing1a == 170) ~ "170 days"
  )) %>%
  mutate(importation_time_label = factor(importation_time_label, levels=c("50 days", "110 days", "170 days"))) %>%
  mutate(vaccine_start_label = case_when( (vaccine_1_start == runtime & vaccine_2_start == runtime) ~ "No vaccine",
                                          (vaccine_1_start == runtime & vaccine_2_start < runtime) ~ "SARS-X vaccine",
                                          (vaccine_1_start < runtime) ~ "Broadly protective sarbecovirus vaccine"
  )) %>%
  mutate(importation_time_label = factor(importation_time_label)) %>%
  filter(Rt1a==R0) ##temporarily in - remove once re-run batch


# select only vaccine scenarios

df2a <-filter(df2, vaccine_2_start<runtime) 

s0 <- ggplot(df2a, aes(x = importation_time_label, y = hospitalisations_averted)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = round(stat(y), 0), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Day of importation", y = "Hospitalisations averted") +
  theme_classic(base_size = 13) +
  scale_fill_manual("Vaccination Strategy", values = cols2)
s0
ggsave("plots/Secondary_HospitalisationsAverted.png",s1, height=3, width=6)



s1 <- ggplot(df2a, aes(x = importation_time_label, y = infections_averted)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = round(stat(y), 0), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Day of importation", y = "Infections averted") +
  theme_classic() +
  scale_fill_manual("Vaccination Strategy", values = cols2)
s1
ggsave("plots/Secondary_InfectionsAverted.png",s1, height=3, width=6)


s2 <- ggplot(df2a, aes(x = importation_time_label, y = deaths_averted)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = round(stat(y), 0), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Day of importation", y = "Deaths averted") +
  theme_classic(base_size = 13) +
  scale_fill_manual("Vaccination Strategy", values = cols2)
s2
ggsave("plots/Secondary_DeathsAverted.png",s2, height=3, width=6)


s3 <- ggplot(df2a, aes(x = importation_time_label, y = years_life_saved)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = round(stat(y), 0), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Day of importation", y = "Life years saved") +
  theme_classic() +
  scale_fill_manual("Vaccination Strategy", values = cols2)
s3
ggsave("plots/Secondary_YLL.png",s3, height=3, width=6)

s4 <- ggplot(df2a, aes(x = importation_time_label, y = vsl_sum/1e9)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = paste0("$", round(stat(y), 2), "B"), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Day of importation", y = expression("Financial value of lives saved")) +
  theme_classic() +
  scale_fill_manual("Vaccination Strategy", values = cols2) +
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0.05)) # set the upper limit to 10 and add some padding
s4
ggsave("plots/Secondary_VSL.png",s4, height=3, width=6)


s5 <- ggplot(df2a, aes(x = importation_time_label, y = prod_loss/1e6)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = paste0("$", round(stat(y), 2), "M"), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Detection time", y = expression("Productivity loss due to mortality (USD)")) +
  theme_classic() +
  scale_fill_manual("Vaccination Strategy", values = cols2) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0.05)) # set the upper limit to 10 and add some padding
s5



layout <- "
AA
BB
CD
"
combined <- g1b + 
  g3b + 
  s0 + 
  s2 + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect", design = layout)

combined
ggsave("plots/Secondary_combined.png", combined, height=16, width=12)


##############################################################
# Economic tables 

#econ <- df2 %>% filter (vaccine_start_label != "No vaccine") %>% filter( ifr_scaling ==1 & R0==2.5) %>% 
#  group_by(income_group, importation_time_label, vaccine_start_label)

#write.csv (econ, file= "./table_secondary.csv")

# Economic tables 


df_summary <- df2 %>%
  filter (vaccine_start_label != "No vaccine") %>%
  mutate(vaccine_start_label = case_when( (vaccine_1_start < runtime) ~ "Sarbecovirus_vaccine",
                                          (vaccine_1_start == runtime & vaccine_2_start < runtime) ~ "SARS_3_vaccine",
                                          (vaccine_1_start == runtime & vaccine_2_start == runtime) ~ "No vaccine"
  ))%>%
  mutate(total_econ = vsl_sum + prod_loss + coi_sum) %>%
  group_by(income_group, importation_time_label, vaccine_start_label) %>%
  summarize(
    median_total_econ = median(total_econ)/1e6,
    lower_total_econ = quantile(total_econ, 0.25)/1e6,
    upper_total_econ = quantile(total_econ, 0.75)/1e6,
    median_vsly = median(vsly_sum)/1e6,
    lower_vsly = quantile(vsly_sum, 0.25)/1e6,
    upper_vsly = quantile(vsly_sum, 0.75)/1e6,
    median_prod_loss = median(prod_loss)/1e6,
    lower_prod_loss = quantile(prod_loss, 0.25)/1e6,
    upper_prod_loss = quantile(prod_loss, 0.75)/1e6,
    median_coi = median(coi_sum)/1e6,
    lower_coi = quantile(coi_sum, 0.25)/1e6,
    upper_coi = quantile(coi_sum, 0.75)/1e6,
    median_vsl = median(vsl_sum)/1e6,
    lower_vsl = quantile(vsl_sum, 0.25)/1e6,
    upper_vsl = quantile(vsl_sum, 0.75)/1e6
  ) %>%
  ungroup()

table <- df_summary %>%
  mutate(median_total_econ = round(median_total_econ,digits=1), lower_total_econ = round(lower_total_econ,digits=1), upper_total_econ = round(upper_total_econ,digits=1)) %>%
  mutate(total_econ = paste0(median_total_econ, " (",lower_total_econ,"-",upper_total_econ,")")) %>%
  mutate(median_vsly = round(median_vsly,digits=1), lower_vsly = round(lower_vsly,digits=1), upper_vsly = round(upper_vsly,digits=1)) %>%
  mutate(vsly = paste0(median_vsly, " (",lower_vsly,"-",upper_vsly,")")) %>%
  mutate(median_prod_loss = round(median_prod_loss,digits=1), lower_prod_loss = round(lower_prod_loss,digits=1), upper_prod_loss = round(upper_prod_loss,digits=1)) %>%
  mutate(prod_loss = paste0(median_prod_loss, " (",lower_prod_loss,"-",upper_prod_loss,")")) %>%
  mutate(median_coi = round(median_coi,digits=1), lower_coi = round(lower_coi,digits=1), upper_coi = round(upper_coi,digits=1)) %>%
  mutate(coi = paste0(median_coi, " (",lower_coi,"-",upper_coi,")")) 

tb.econ <- table %>%
  select(income_group, importation_time_label, vaccine_start_label,
         total_econ, vsly, prod_loss, coi) %>%
  pivot_wider(names_from = vaccine_start_label, values_from = c(total_econ, vsly, prod_loss, coi)) %>%
  flextable()         

save_as_docx("Secondary" = tb.econ,
             path = "tables/Secondary_Table.docx")


