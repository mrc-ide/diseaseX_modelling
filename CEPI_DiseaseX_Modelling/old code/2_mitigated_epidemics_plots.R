
df_main <- readRDS("output/2_mitigated_scenarios.rds") 
#df_main <- out_combine

# Time series plots for HIC 

# limit to an illustrative scenario
df1 <- df_main %>%
  filter(income_group == "HIC" & ifr_scaling == 1 ) %>%
  unnest(cols = output) %>%
  filter(t>1) %>%
  mutate(vaccine_start_label = ifelse(vaccine_start == 365*2, "no vaccine", paste0(vaccine_start))) %>%
  mutate(vaccine_start_label = factor(vaccine_start_label)) %>%
  mutate(vaccine_start_label = fct_relevel(vaccine_start_label, "10", "20", "50", "100", "no vaccine"))

df_deaths <- df1 %>%
            filter(compartment == "deaths" & t<=365*2)
R0_labs <- c('1.5' = "R0 = 1.5",'2.5'="R0 = 2.5",'4'="R0 =4 ")
g1 <- ggplot() +
  geom_line(data = df_deaths, aes(x = t, y = (value),col=vaccine_start_label),linewidth=1.2) +
  facet_wrap(vars(R0), labeller = as_labeller(R0_labs)) +
  labs(x = "Time (days)", y = "Deaths per million per day", col = "Day of vaccine introduction", linetype = "") +
  theme_classic() +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = cols1) +
  scale_linetype_manual(values = c("dashed", "dotted")) 
g1
ggsave("plots/Mitigated_DeathsOverTime.png",g1, height=3, width=10)


df_vaccines <- df1 %>%
  filter(compartment == "vaccines" & t<=250) %>%
  filter(R0 == 2.5) %>%
  group_by(vaccine_start) %>%
  mutate(cum_vaccines = cumsum(value)) 

g2 <- ggplot() +
   geom_line(data = df_vaccines, aes(x = t, y = cum_vaccines/target_pop,col=vaccine_start_label),linewidth=1.2) +
  labs(x = "Time", y = "Vaccine coverage", col = "Day of vaccine introduction", linetype = "") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_manual(values = cols1) +
  scale_linetype_manual(values = c("dashed", "dotted"))
g2

plot <- (g2 | g1 )
ggsave("plots/Figure2.png",plot, height=6, width=18)

##############################################################
# Summary plots

# remove time series for faster processing
df2 <- select(df_main,-contains("output")) %>%
  mutate(vaccine_start_label = ifelse(vaccine_start == 365*2, "no vaccine", paste0(vaccine_start))) %>%
  mutate(vaccine_start_label = factor(vaccine_start_label)) %>%
  mutate(vaccine_start_label = fct_relevel(vaccine_start_label, "10", "20", "50", "100", "no vaccine")) %>%
  mutate(income_group = fct_relevel(income_group, "HIC", "UMIC", "LMIC", "LIC")) %>%
  mutate(prop_mitigation_avoided = (contact_reduction.y-contact_reduction.x)/contact_reduction.y) %>%
  mutate(R0 = factor(R0))



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
ggsave("plots/Mitigated_Deaths.png",s1, height=3, width=6)

s4 <- ggplot() + 
  geom_boxplot(data = df2, aes(x = vaccine_start_label, y = contact_reduction.x, fill = vaccine_start_label)) +
  labs(x = "", y = "Weighted mitigation days") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_fill_manual( "Day of vaccine introduction",values = cols1) 
s4

s5 <- ggplot(data = df2, aes(x = R0, y=contact_reduction.x, fill=vaccine_start_label)) + 
  geom_col(position = "dodge") +
  labs(x = "R0", y = "Contact-weighted NPI days") +
  theme_classic() +
  scale_fill_manual( "Day of vaccine introduction",values = cols1) 
s5
ggsave("plots/Mitigated_NPI_Impact.png",s5, height=3, width=6)


df2a <-filter(df2, vaccine_start<365*2)
s2 <- ggplot() +
  geom_boxplot(data = df2a, aes(x = income_group, y = deaths_averted, fill = vaccine_start_label)) +
  labs(x = "World Bank Income Setting", y = "Deaths Averted") +
  theme_classic() +
  scale_fill_manual( "Day of vaccine introduction",values = cols1) 
s2
ggsave("plots/Mitigated_DeathsAverted.png",s2, height=3, width=6)


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
ggsave("plots/Mitigated_YLL.png",s3, height=3, width=6)


s5 <- ggplot(data = df2a, aes(x = vaccine_start_label, y=prop_mitigation_avoided, fill=vaccine_start_label)) + 
  geom_col(width=0.5) +
  labs(x = "", y = "Reduction in NPIs") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_fill_manual( "Day of vaccine introduction",values = cols1) 
s5

