
df_main <- readRDS("output/3_emergence_mv_scenarios.rds") 
#df_main <- out_combine

# Time series plots for UMIC

# limit to an illustrative scenario
df1 <- df_main %>%
  filter(income_group == "UMIC" & ifr_scaling == 1) %>%
  unnest(cols = output) %>%
  filter(t>1) %>%
  mutate(detection_time_label = case_when ((timing1 == 30) ~ "30 days",
                                           (timing1 == 60) ~ "60 days",
                                           (timing1 == 90) ~ "90 days"
                                           )) %>%
  mutate(vaccine_start_label = case_when( (vaccine_1_start < runtime) ~ "Sarbecovirus vaccine",
                                          (vaccine_1_start == runtime & vaccine_2_start < runtime) ~ "SARS-3 vaccine",
                                          (vaccine_1_start == runtime & vaccine_2_start == runtime) ~ "No vaccine",
                                           )) %>%
  mutate(detection_time_label = factor(detection_time_label)) %>%
  mutate(vaccine_start_label = factor(vaccine_start_label)) %>%
  mutate(vaccine_start_label = fct_relevel(vaccine_start_label, "Sarbecovirus vaccine", 
                                           "SARS-3 vaccine", "No vaccine"))


df_deaths <- df1 %>%
            filter(compartment == "deaths") %>%
            filter(R0==1.5) %>%
            filter((vaccine_start_label=="No vaccine" & two_vaccines==0) | vaccine_start_label == "SARS-3 vaccine" 
                   | vaccine_start_label == "Sarbecovirus vaccine")

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
  filter((vaccine_start_label=="No vaccine" & two_vaccines==0) | vaccine_start_label == "SARS-3 vaccine" 
         | vaccine_start_label == "Sarbecovirus vaccine")

g1b <- ggplot() +
  geom_line(data = df_deaths, aes(x = t, y = (value),col=vaccine_start_label),linewidth=1.2) +
  facet_wrap(vars(detection_time_label)) +
  #  facet_wrap(vars(detection_time), labeller = as_labeller(detection_time_label)) +
  labs(x = "Time (days)", y = "Deaths per million per day", col = "Vaccine Strategy", linetype = "") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = cols1) +
  scale_linetype_manual(values = c("dashed", "dotted")) 
g1b
ggsave("plots/Time_Series_Deaths_R0_2.5.png",g1b, height=3, width=10)

df_deaths <- df1 %>%
  filter(compartment == "deaths") %>%
  filter(R0==4) %>%
  filter((vaccine_start_label=="No vaccine" & two_vaccines==0) | vaccine_start_label == "SARS-3 vaccine" 
         | vaccine_start_label == "Sarbecovirus vaccine")

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
  geom_line(data = df_v1, aes(x = t, y = value/target_pop, linetype = "Sarbecovirus", col=vaccine_start_label),linewidth=1.2) +
  geom_line(data = df_v2, aes(x = t, y = value/target_pop, linetype = "SARS-3", col=vaccine_start_label),linewidth=1.2) +
  scale_linetype_manual(name = "Vaccine", values = c("solid", "dashed")) +
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1.0) ) +
  scale_color_manual(name = "Strategy", values = cols1a) + 
  labs(x = "Time", y = "Vaccine coverage", col = "Strategy", linetype = "") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) 
g2
ggsave("plots/Vaccine_Rollout_emergence.png",g2, height=3, width=6)

##############################################################
# Summary plots

# remove time series for faster processing
df2 <- select(df_main,-contains("output")) %>%
  mutate(detection_time_label = case_when ((timing1 == 30) ~ "30 days",
                                           (timing1 == 60) ~ "60 days",
                                           (timing1 == 90) ~ "90 days"
  )) %>%
  mutate(vaccine_start_label = case_when( (vaccine_1_start == runtime & vaccine_2_start == runtime) ~ "No vaccine",
                                          (vaccine_1_start == runtime & vaccine_2_start < runtime) ~ "SARS-3 vaccine",
                                          (vaccine_1_start < runtime) ~ "Sarbecovirus vaccine"
  )) %>%
  mutate(detection_time_label = factor(detection_time_label)) 

# select only vaccine scenarios

df2a <-filter(df2, vaccine_2_start<runtime) 

s1 <- ggplot(df2a, aes(x = detection_time_label, y = infections_averted)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = round(stat(y), 0), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Detection time", y = "Infections averted") +
  theme_classic() +
  scale_fill_manual("Vaccination Strategy", values = cols1)
s1
ggsave("plots/Emergence_InfectionsAverted.png",s1, height=3, width=6)


s2 <- ggplot(df2a, aes(x = detection_time_label, y = deaths_averted)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = round(stat(y), 0), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Detection time", y = "Deaths averted") +
  theme_classic() +
  scale_fill_manual("Vaccination Strategy", values = cols1)
s2
ggsave("plots/Emergence_DeathsAverted.png",s2, height=3, width=6)


s3 <- ggplot(df2a, aes(x = detection_time_label, y = years_life_saved)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = round(stat(y), 0), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Detection time", y = "Life years saved") +
  theme_classic() +
  scale_fill_manual("Vaccination Strategy", values = cols1)
s3
ggsave("plots/Emergence_YLL.png",s3, height=3, width=6)

s4 <- ggplot(df2a, aes(x = detection_time_label, y = vsl_sum/1e9)) +
  geom_boxplot(aes(fill = vaccine_start_label)) +
  stat_summary(aes(label = paste0("$", round(stat(y), 2), "B"), fill = vaccine_start_label),
               fun = median, geom = "text", size = 3.5, vjust = -0.5, color = "white",
               position = position_dodge(width = 0.75), show.legend = FALSE) +
  labs(x = "Detection time", y = expression("Financial value of lives saved")) +
  theme_classic() +
  scale_fill_manual("Vaccination Strategy", values = cols1) +
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0.05)) # set the upper limit to 10 and add some padding
s4
ggsave("plots/Emergence_VSL.png",s4, height=3, width=6)



layout <- "
AA
BC
DE
"
combined <- g1b + 
  s1 + 
  s2 + 
  s3 + 
  s4 +
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect", design = layout)

combined
ggsave("plots/Emergence_combined.png", combined, height=16, width=12)

