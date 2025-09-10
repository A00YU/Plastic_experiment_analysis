# load library
library(survival)
library(coxme)
library(survminer)
library(lme4)
library(glmmTMB)
library(emmeans)
library(glm2)
library(tidyr)
library(dplyr)
library(ggplot2)

####### load data#########
# Load the master data sheet
weekly_data_r2 <- read.csv("/Users/aoyu/Desktop/Snail_Data/Plastic Tank Experiment Data Sheet Round 2.csv", skip = 2)

all_data_r2 <- weekly_data_r2 %>%
  mutate(
    Treatment = recode(Color,
                       "G" = "Control", "Y" = "Seasoned Macro - Low", "O" = "Seasoned Macro - High", "W-1" = "Virgin Macro - Low", "W-3" = "Virgin Macro - High"),
    Treatment = factor(Treatment, levels = c("Control", "Seasoned Macro - Low", "Seasoned Macro - High",
                                             "Virgin Macro - Low", "Virgin Macro - High")),
    # Unique physical jar (color x jar number)
    jar_id = interaction(Color, Jar_num, drop = TRUE)
  )

# table(all_data_r2$Color, all_data_r2$Week_num)
