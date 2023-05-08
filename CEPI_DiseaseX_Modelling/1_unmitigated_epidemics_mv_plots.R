
df_main <- readRDS("output/1_unmitigated_mv_scenarios.rds") 
#df_main <- out_combine

# Time series plots for HIC 

# limit to an illustrative scenario
df1 <- df_main %>%
  select(-output_cf) %>%
  filter(income_group == "HIC" & ifr_scaling == 1 ) %>%
  unnest(cols = output) %>%
  filter(t>1) %>%
  mutate(vaccine_start_label = case_when( (vaccine_1_start == 365*2 & vaccine_2_start == 365*2) ~ "No vaccine",
                                          (vaccine_1_start == 365*2 & vaccine_2_start == 100) ~ "SARS-3 vaccine",
                                          (vaccine_1_start < 365*2) ~ paste0("Sarbecovirus vaccine day ",vaccine_1_start)
  )
  ) %>%
  mutate(vaccine_start_label = factor(vaccine_start_label)) %>%
  mutate(vaccine_start_label = fct_relevel(vaccine_start_label, "Sarbecovirus vaccine day 10", "Sarbecovirus vaccine day 20", "Sarbecovirus vaccine day 30", "SARS-3 vaccine", "no vaccine"))

df_deaths <- df1 %>%
            filter(compartment == "deaths" & t<=250)

R0_labs <- c('1.5' = "R0 = 1.5",'2.5'="R0 = 2.5",'4'="R0 =4 ")
g1 <- ggplot() +
 # geom_segment(data = df_deaths, aes(x = 10, y = 240, xend = 10, yend = 240), col = "#08519c") +
#  geom_segment(data = df_deaths, aes(x = 20, y = 210, xend = 20, yend = 240), col = "#440154ff") +
#  geom_segment(data = df_deaths, aes(x = 50, y = 210, xend = 50, yend = 240), col = "#9ecae1") +
#  geom_segment(data = df_deaths, aes(x = 100, y = 210, xend = 100, yend = 240), col = "#993366") +
  geom_line(data = df_deaths, aes(x = t, y = (value),col=vaccine_start_label),linewidth=1.2) +
  facet_wrap(vars(R0), labeller = as_labeller(R0_labs)) +
  labs(x = "Time (days)", y = "Deaths per million per day", col = "Day of vaccine introduction", linetype = "") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = cols1) +
  scale_linetype_manual(values = c("dashed")) 
g1
ggsave("plots/Unmitigated_DeathsOverTime_mv.png",g1, height=3, width=10)


df_v1 <- df1 %>%
  filter(compartment == "V1") %>%
 # filter(t<=250) %>%
  filter(R0 == 2.5) %>%
  group_by(vaccine_start_label) 


df_v2 <- df1 %>%
  filter(compartment == "V2") %>%
  # filter(t<=250) %>%
  filter(R0 == 2.5) %>%
  group_by(vaccine_start_label) 
 
g2 <- ggplot() +
  geom_line(data = df_v1, aes(x = t, y = value/sum(pop[]), linetype = "Sarbecovirus", col=vaccine_start_label),linewidth=1.2) +
  geom_line(data = df_v2, aes(x = t, y = value/sum(pop[]), linetype = "SARS-3", col=vaccine_start_label),linewidth=1.2) +
  scale_linetype_manual(name = "Vaccine", values = c("solid", "dashed")) +
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1.0) ) +
  labs(x = "Time", y = "Vaccine coverage", col = "Strategy", linetype = "") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_manual(values = cols1)
g2
ggsave("plots/Vaccine_Rollout_mv.png",g2, height=3, width=6)


layout <- "
A
B
"
combined <- (g2+ theme(legend.position = "none")) + g1 +
  plot_annotation(tag_levels = "A") + 
  plot_layout(design = layout)
combined
ggsave("plots/Figure1_mv.png",combined, height=3, width=9)

##############################################################
# Summary plots

# remove time series for faster processing
df2 <- select(df_main,-contains("output")) %>%
  mutate(vaccine_start_label = case_when( (vaccine_1_start == 365*2 & vaccine_2_start == 365*2) ~ "No vaccine",
                                          (vaccine_1_start == 365*2 & vaccine_2_start == 100) ~ "SARS-3 vaccine",
                                          (vaccine_1_start < 365*2) ~ paste0("Sarbecovirus vaccine day ",vaccine_1_start)
                                          )) %>%
  mutate(vaccine_start_label = factor(vaccine_start_label)) %>%
#  mutate(vaccine_start_label = fct_relevel(vaccine_start_label, "10", "20", "50", "100", "no vaccine")) %>%
  mutate(income_group = fct_relevel(income_group, "HIC", "UMIC", "LMIC", "LIC"))


s1 <- ggplot() +
  geom_boxplot(data = df2, aes(x = income_group, y = deaths, fill = vaccine_start_label)) +
  labs(x = "World Bank Income Setting", y = "Deaths") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_fill_manual( "Day of vaccine introduction",values = cols1) 
  
s1
ggsave("plots/Unmitigated_Deaths_mv.png",s1, height=3, width=6)

df2a <-filter(df2, vaccine_start_label != "No vaccine")

s2 <- ggplot() +
  geom_boxplot(data = df2a, aes(x = income_group, y = deaths_averted, fill = vaccine_start_label)) +
  labs(x = "World Bank Income Setting", y = "Deaths Averted") +
  theme_classic() +
  scale_fill_manual( "Day of vaccine introduction",values = cols1) 
s2
ggsave("plots/Unmitigated_DeathAverted_mv.png",s2, height=3, width=6)


s3 <- ggplot() +
  geom_boxplot(data = df2a, aes(x = income_group, y = years_life_saved, fill = vaccine_start_label)) +
  labs(x = "World Bank Income Setting", y = "Years of Life Saved") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_fill_manual( "Day of vaccine introduction",values = cols1) 
s3
ggsave("plots/Unmitigated_YLL_mv.png",s3, height=3, width=6)

