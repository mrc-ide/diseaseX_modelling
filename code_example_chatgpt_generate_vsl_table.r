library(dplyr)
library(tidyr)
library(officer)

# example dataset
df <- data.frame(
  strategy = rep(c("Vaccine A", "Vaccine B"), each = 36),
  R0 = rep(rep(c(1.5, 2, 2.5), each = 12), 2),
  IFR = rep(c(0.01, 0.02, 0.03), each = 24),
  income_setting = rep(c("Low", "Medium", "High"), each = 12, times = 2),
  date_of_detection = rep(c("2022-01-01", "2022-02-01", "2022-03-01"), each = 8, times = 2),
  VSLY = rnorm(72, mean = 100, sd = 10)
)

# calculate median and upper/lower range of VSLY by strategy, R0, IFR, income_setting, and date_of_detection
df_summary <- df %>%
  group_by(strategy, R0, IFR, income_setting, date_of_detection) %>%
  summarize(
    median_VSLY = median(VSLY),
    lower_VSLY = quantile(VSLY, 0.25),
    upper_VSLY = quantile(VSLY, 0.75)
  ) %>%
  ungroup()

# pivot the table to get columns for each strategy
df_summary_wide <- df_summary %>%
  pivot_wider(
    names_from = strategy,
    values_from = c(median_VSLY, lower_VSLY, upper_VSLY),
    names_prefix = "VSLY_"
  )

# reorder columns for readability
df_summary_wide <- df_summary_wide %>%
  select(income_setting, date_of_detection, starts_with("Vaccine"))

# create a new Word document
doc <- read_docx()

# add a table to the document
doc <- doc %>%
  body_add_table(df_summary_wide, style = "Table Grid")

# save the document
print(doc, target = "summary_table.docx")
