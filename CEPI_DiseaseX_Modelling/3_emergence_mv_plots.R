
source("./R/functions_multivaccine.R")

df_main <- readRDS("output/3_emergence_mv_scenarios.rds") 
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


# Time series plots for LMIC

# limit to an illustrative scenario

income_to_plot = "HIC"


df1 <- df_main %>%
  filter(income_group == income_to_plot & ifr_scaling == 1) %>%
  unnest(cols = output) %>%
  filter(t>1) %>%
  mutate(detection_time_label = case_when ((timing1 == 30) ~ "30 days",
                                           (timing1 == 60) ~ "60 days",
                                           (timing1 == 90) ~ "90 days"
  )) %>%
  mutate(vaccine_start_label = case_when( (vaccine_1_start < runtime) ~ "BPSV + SARS-X",
                                          (vaccine_1_start == runtime & vaccine_2_start < runtime) ~ "SARS-X ",
                                          (vaccine_1_start == runtime & vaccine_2_start == runtime) ~ "No vaccine",
  )) %>%
  mutate(detection_time_label = factor(detection_time_label)) %>%
  mutate(vaccine_start_label = factor(vaccine_start_label)) %>%
  mutate(vaccine_start_label = fct_relevel(vaccine_start_label, "BPSV + SARS-X", 
                                           "SARS-X ", "No vaccine"))


df_deaths <- df1 %>%
  filter(compartment == "deaths") %>%
  filter(R0==1.5) %>%
  filter((vaccine_start_label=="No vaccine" & two_vaccines==0) | vaccine_start_label == "SARS-X " 
         | vaccine_start_label == "BPSV + SARS-X")

df_deaths_tot <- df_deaths %>%
  group_by(vaccine_start_label, detection_time_label) %>%
  summarise(total = sum(value))

g1a <- ggplot() +
  geom_line(data = df_deaths, aes(x = t, y = (value),col=vaccine_start_label),linewidth=1.1) +
  facet_wrap(vars(detection_time_label)) +
  #  facet_wrap(vars(detection_time), labeller = as_labeller(detection_time_label)) +
  labs(x = "Time (days)", y = "Deaths per million per day", col = "Vaccine Strategy", linetype = "") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = cols1) +
  scale_linetype_manual(values = c("dashed", "dotted")) 
g1a
ggsave("plots/Time_Series_Deaths_R0_1.5.png",g1a, height=3, width=10)

df_deaths <- df1 %>%
  filter(compartment == "deaths") %>%
  filter(R0==2.5) %>%
  filter((vaccine_start_label=="No vaccine" & two_vaccines==0) | vaccine_start_label == "SARS-X " 
         | vaccine_start_label == "BPSV + SARS-X")

df_rec<- df_deaths %>% filter(t==100)


g1b <- ggplot() +
  geom_rect(data = df_rec %>% filter(vaccine_start_label == "SARS-X "),aes(xmin=timing1, xmax=timing2, ymin=90,ymax=93), fill="#440154ff", alpha=0.8) +
  geom_rect(data = df_rec %>% filter(vaccine_start_label == "No vaccine") ,aes(xmin=timing1, xmax=timing2, ymin=87,ymax=90),alpha=0.8,fill="grey40") +
  geom_rect(data = df_rec %>% filter(vaccine_start_label == "BPSV + SARS-X"),aes(xmin=timing1, xmax=timing2, ymin=84,ymax=87), fill="#08519c",alpha=0.8) +
  #geom_rect (data= df_rec, aes(xmin=timing1, xmax=timing2, ymin=0,ymax=53,fill=vaccine_start_label ), alpha=0.1)+ 
  geom_line(data = df_deaths,aes(x = t, y = (value),col=vaccine_start_label),linewidth=1.2) +
  facet_wrap(vars(detection_time_label)) +
  labs(x = "Time (days)", y = "Deaths per million per day", col = "Vaccine Strategy", linetype = "") +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = cols1) +
  scale_fill_manual(values=cols1) +
scale_linetype_manual(values = c("dashed", "dotted")) 
g1b
ggsave("plots/Time_Series_Deaths_R0_2.5.png",g1b, height=3, width=10)


g1b_cum <- ggplot() +
   geom_line(data = df_deaths,aes(x = t, y = (cum),col=vaccine_start_label),linewidth=1.2) +
  facet_wrap(vars(detection_time_label)) +
  labs(x = "Time (days)", y = "Cumulative deaths per million", col = "Vaccine Strategy", linetype = "") +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = cols1) +
  scale_fill_manual(values=cols1) +
  scale_linetype_manual(values = c("dashed", "dotted")) 
g1b_cum
ggsave("plots/Time_Series_Deaths_cum_R0_2.5.png",g1b, height=3, width=10)


df_deaths <- df1 %>%
  filter(compartment == "deaths") %>%
  filter(R0==4) %>%
  filter((vaccine_start_label=="No vaccine" & two_vaccines==0) | vaccine_start_label == "SARS-X " 
         | vaccine_start_label == "BPSV + SARS-X")

g1c <- ggplot() +
  geom_line(data = df_deaths, aes(x = t, y = (value),col=vaccine_start_label),linewidth=1.2) +
  facet_wrap(vars(detection_time_label)) +
  #  facet_wrap(vars(detection_time), labeller = as_labeller(detection_time_label)) +
  labs(x = "Time (days)", y = "Deaths per million per day", col = "Vaccine Strategy", linetype = "") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = cols1) +
  scale_linetype_manual(values = c("dashed", "dotted")) 
g1c
ggsave("plots/Time_Series_Deaths_R0_4.png",g1c, height=3, width=10)


df_v1 <- df1 %>%
  filter(timing1==60) %>%
  filter(compartment == "V1") %>%
  #  filter(t<=250) %>%
  filter(R0 == 2.5) %>%
  group_by(vaccine_start_label)


df_v2 <- df1 %>%
  filter(timing1==60) %>%
  filter(compartment == "V2") %>%
  #  filter(t<=250) %>%
  filter(R0 == 2.5) %>%
  group_by(vaccine_start_label)


g2 <- ggplot() +
  geom_line(data = df_v1, aes(x = t, y = value/target_pop, linetype = "BPSV + SARS-X", col=vaccine_start_label),linewidth=1.2) +
  geom_line(data = df_v2, aes(x = t, y = value/target_pop, linetype = "SARS-X", col=vaccine_start_label),linewidth=1.2) +
  scale_linetype_manual(name = "Vaccine", values = c("solid", "dashed")) +
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1.0) ) +
  scale_color_manual(name = "Strategy", values = cols1a) + 
  labs(x = "Time", y = "Vaccine coverage", col = "Strategy", linetype = "", title = income_to_plot) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) 

ggsave(paste0("./plots/Vaccine_Rollout_emergence_",income_to_plot,".png"),g2, height=3, width=6)


g2_HIC <- g2

layout <- "
AB
CD
"
combined <- g2_LIC +
  g2_LMIC +
  g2_UMIC +
  g2_HIC +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", design = layout)

combined

ggsave("plots/Emergence_vaccRates.png", combined, height=13, width=16)



# Hospitalisation plot 

rep_country <- get_representative_country(income_group = income_to_plot)
capacity <- squire::get_healthcare_capacity(country=rep_country)$hosp_beds

df_hosp <- df1 %>%
  filter(compartment == "hospital_occupancy") %>%
  filter(R0==2.5) %>%
  filter((vaccine_start_label=="No vaccine" & two_vaccines==0) | vaccine_start_label == "SARS-X " 
         | vaccine_start_label == "BPSV + SARS-X")

df_rec<- df_hosp %>% filter(t==100)

g3b <- ggplot(data = df_hosp) +
  geom_rect(data = df_rec %>% filter(vaccine_start_label == "SARS-X "),aes(xmin=timing1, xmax=timing2, ymin=5000,ymax=4900), fill="#440154ff", alpha=0.8) +
  geom_rect(data = df_rec %>% filter(vaccine_start_label == "No vaccine") ,aes(xmin=timing1, xmax=timing2, ymin=4900,ymax=4800),alpha=0.8,fill="grey20") +
  geom_rect(data = df_rec %>% filter(vaccine_start_label == "BPSV + SARS-X"),aes(xmin=timing1, xmax=timing2, ymin=4800,ymax=4700), fill="#08519c",alpha=0.8) +
  #geom_rect (data= df_rec, aes(xmin=timing1, xmax=timing2, ymin=0,ymax=3300,fill=vaccine_start_label ), alpha=0.1)+ 
  geom_line(aes(x = t, y = (value),col=vaccine_start_label),linewidth=1.1) +
  geom_hline(yintercept = capacity*1000)+
  facet_grid(~detection_time_label, scales="free") +
  #  facet_wrap(vars(detection_time), labeller = as_labeller(detection_time_label)) +
  labs(x = "Time (days)", y = "Hospital ocuppancy per million per day", col = "Vaccine Strategy", linetype = "") +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank(),
      legend.position = "none") +
  scale_color_manual(values = cols1) +
  scale_fill_manual(values=cols1) +
scale_linetype_manual(values = c("dashed", "dotted")) 
g3b

ggsave("plots/Emergence_Time_Series_Hosp_R0_2.5.png",g1b, height=3, width=10)

# Table days above capacity

days_above_capacity <- df_hosp %>% group_by (vaccine_start_label,detection_time_label) %>% 
  select (vaccine_start_label,detection_time_label,day_above_capacity,day_above_capacity_cf) %>% distinct() %>%
  flextable()

save_as_docx("Emergence" = days_above_capacity,
             path = "tables/Emergence_Hospital_Capacity_Table.docx")



df_hosp <- df1 %>%
  filter(compartment == "hospitalisations") %>%
  filter(R0==2.5) %>%
  filter((vaccine_start_label=="No vaccine" & two_vaccines==0) | vaccine_start_label == "SARS-X " 
         | vaccine_start_label == "BPSV + SARS-X")

g3b_cum <- ggplot(data = df_hosp) +
    geom_line(aes(x = t, y = (cum),col=vaccine_start_label),linewidth=1.1) +
   facet_grid(~detection_time_label, scales="free") +
  #  facet_wrap(vars(detection_time), labeller = as_labeller(detection_time_label)) +
  labs(x = "Time (days)", y = "Cumulative hospitalisations per million", col = "Vaccine Strategy", linetype = "") +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = cols1) +
  scale_fill_manual(values=cols1) +
  scale_linetype_manual(values = c("dashed", "dotted")) 
g3b_cum


##############################################################
# Summary plots

# remove time series for faster processing
df2 <- select(df_main,-contains("output")) %>%
  mutate(detection_time_label = case_when ((timing1 == 30) ~ "30 days",
                                           (timing1 == 60) ~ "60 days",
                                           (timing1 == 90) ~ "90 days"
  )) %>%
  mutate(vaccine_start_label = case_when( (vaccine_1_start == runtime & vaccine_2_start == runtime) ~ "No vaccine",
                                          (vaccine_1_start == runtime & vaccine_2_start < runtime) ~ "SARS-X ",
                                          (vaccine_1_start < runtime) ~ "BPSV + SARS-X"
  )) %>%
  mutate(detection_time_label = factor(detection_time_label)) 

# select only vaccine scenarios

df2a <-filter(df2, vaccine_2_start<runtime) 

s0 <- ggplot(df2a, aes(x = detection_time_label, y = hospitalisations_averted)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = round(stat(y), 0), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Detection time", y = " Hospitalisations averted") +
  theme_classic(base_size = 13) +
  scale_fill_manual("Vaccination Strategy", values = cols1)
s0
ggsave("plots/Emergence_HospitalisationsAverted.png",s0, height=3, width=6)

s1 <- ggplot(df2a, aes(x = detection_time_label, y = infections_averted)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = round(stat(y), 0), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Detection time", y = "Infections averted") +
  theme_classic(base_size = 13) +
  scale_fill_manual("Vaccination Strategy", values = cols1)
s1
ggsave("plots/Emergence_InfectionsAverted.png",s1, height=3, width=6)


s2 <- ggplot(df2a, aes(x = detection_time_label, y = deaths_averted)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = round(stat(y), 0), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Detection time", y = "Deaths averted") +
  theme_classic(base_size = 13) +
  scale_fill_manual("Vaccination Strategy", values = cols1)
s2
ggsave("plots/Emergence_DeathsAverted.png",s2, height=3, width=6)


s3 <- ggplot(df2a, aes(x = detection_time_label, y = years_life_saved)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = round(stat(y), 0), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Detection time", y = "Life years saved") +
  theme_classic(base_size = 13) +
  scale_fill_manual("Vaccination Strategy", values = cols1)
s3
ggsave("plots/Emergence_YLL.png",s3, height=3, width=6)

s4 <- ggplot(df2a, aes(x = detection_time_label, y = vsl_sum/1e9)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = paste0("$", round(stat(y), 2), "B"), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Detection time", y = expression("Financial value of lives saved (USD)")) +
  theme_classic(base_size = 13) +
  scale_fill_manual("Vaccination Strategy", values = cols1) +
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0.05)) # set the upper limit to 10 and add some padding
s4
ggsave("plots/Emergence_VSL.png",s4, height=3, width=6)


s5 <- ggplot(df2a, aes(x = detection_time_label, y = prod_loss/1e6)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = paste0("$", round(stat(y), 2), "M"), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Detection time", y = expression("Productivity loss due to mortality (USD)")) +
  theme_classic(base_size = 13) +
  scale_fill_manual("Vaccination Strategy", values = cols1) +
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


ggsave(paste0("plots/Emergence_combined_" ,income_to_plot,".png"), combined, height=15, width=11)


#Appendix 


layout <- "
AAA
BBB
CDE
"
combined <- g1b_cum + 
  g3b_cum + 
  s3 +
  s4 + 
  s5 +
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect", design = layout)

combined


ggsave(paste0("plots/Emergence_appendix_" ,income_to_plot,".png"), combined, height=16, width=16)



##############################################################
# Economic tables 

#econ <- df2 %>% filter (vaccine_start_label != "No vaccine") %>% filter( ifr_scaling ==1 & R0==2.5) %>% 
#  group_by(income_group, detection_time_label, vaccine_start_label)

#write.csv (econ, file= "./table_emergence.csv")

df_summary <- df2 %>%
  filter (vaccine_start_label != "No vaccine") %>%
  mutate(vaccine_start_label = case_when( (vaccine_1_start < runtime) ~ "Sarbecovirus_vaccine",
                                          (vaccine_1_start == runtime & vaccine_2_start < runtime) ~ "SARS_3_vaccine",
                                          (vaccine_1_start == runtime & vaccine_2_start == runtime) ~ "No vaccine"
                                           ))%>%
  mutate(total_econ = vsl_sum + prod_loss + coi_2) %>%
  group_by(income_group, detection_time_label, vaccine_start_label) %>%
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
  mutate(median_vsly = round(median_vsly,digits=1), lower_vsly = round(lower_vsly,digits=0), upper_vsly = round(upper_vsly,digits=1)) %>%
  mutate(vsly = paste0(median_vsly, " (",lower_vsly,"-",upper_vsly,")")) %>%
  mutate(median_prod_loss = round(median_prod_loss,digits=1), lower_prod_loss = round(lower_prod_loss,digits=1), upper_prod_loss = round(upper_prod_loss,digits=1)) %>%
  mutate(prod_loss = paste0(median_prod_loss, " (",lower_prod_loss,"-",upper_prod_loss,")")) %>%
  mutate(median_coi = round(median_coi,digits=1), lower_coi = round(lower_coi,digits=1), upper_coi = round(upper_coi,digits=1)) %>%
  mutate(coi = paste0(median_coi, " (",lower_coi,"-",upper_coi,")")) 
  
tb.econ <- table %>%
  select(income_group, detection_time_label, vaccine_start_label,
         total_econ, vsly, prod_loss, coi) %>%
  pivot_wider(names_from = vaccine_start_label, values_from = c(total_econ, vsly, prod_loss, coi)) %>%
  flextable()         

save_as_docx("Emergence" = tb.econ,
             path = "tables/Emergence_Table1.docx")


