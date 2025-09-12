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
library(mgcv)

####### load data#########
source("/Users/aoyu/Desktop/Snail_Data/R1_data_cleaning.R") # cleaned master data sheet called all_data

####### measurement error correction due to orientation #########
all_data_growth <- all_data %>%
  select(Color, Jar_num, Week_num, 
         starts_with("Adult_death_"), starts_with("Snail_length_"), 
         starts_with("snail_matching_"),starts_with("Unidentifed_"),
         starts_with("snail_mating_"), starts_with("snail_orientation_"), Treatment)

# table(all_data$Color,all_data$Week_num)
# replace NA with 0 for non-mating snails
all_data_growth <- all_data_growth %>%
  mutate(
    snail_mating_1 = ifelse(is.na(snail_mating_1), 0, snail_mating_1),
    snail_mating_2 = ifelse(is.na(snail_mating_2), 0, snail_mating_2),
    snail_mating_3 = ifelse(is.na(snail_mating_3), 0, snail_mating_3),
    snail_orientation_1 = ifelse(is.na(snail_orientation_1), 0, snail_orientation_1),
    snail_orientation_2 = ifelse(is.na(snail_orientation_2), 0, snail_orientation_2),
    snail_orientation_3 = ifelse(is.na(snail_orientation_3), 0, snail_orientation_3)
  )
# if mating, snail length is inaccurate, thus delete
all_data_growth <- all_data_growth %>%
  mutate(
    Snail_length_1 = ifelse(snail_mating_1 == 1, NA, as.numeric(Snail_length_1)),
    Snail_length_2 = ifelse(snail_mating_2 == 1, NA, as.numeric(Snail_length_2)),
    Snail_length_3 = ifelse(snail_mating_3 == 1, NA, as.numeric(Snail_length_3))
  )

# count number of mating snails
num_mating <- sum(all_data_growth$snail_mating_1 == 1, na.rm = TRUE) + sum(all_data_growth$snail_mating_2 == 1, na.rm = TRUE) + sum(all_data_growth$snail_mating_3 == 1, na.rm = TRUE)
num_mating # 212 snails

# if not mating but simply orientation is wrong, correct it with a correction ratio: correct/wrong = 1.030 = 1.02997
# round the product into 3 decimal places
all_data_growth <- all_data_growth %>%
  mutate(
    Snail_length_1 = round(ifelse(snail_orientation_1 == 1 & snail_mating_1 != 1, Snail_length_1 * 1.03, Snail_length_1), 3),
    Snail_length_2 = round(ifelse(snail_orientation_2 == 1 & snail_mating_2 != 1, Snail_length_2 * 1.03, Snail_length_2), 3),
    Snail_length_3 = round(ifelse(snail_orientation_3 == 1 & snail_mating_3 != 1, Snail_length_3 * 1.03, Snail_length_3), 3)
  )

# count number of snails with initial wrong orientation
num_orientation <- sum(all_data_growth$snail_orientation_1 == 1 & all_data_growth$snail_mating_1 != 1, na.rm = TRUE) + 
  sum(all_data_growth$snail_orientation_2 == 1 & all_data_growth$snail_mating_2 != 1, na.rm = TRUE) + 
  sum(all_data_growth$snail_orientation_3 == 1 & all_data_growth$snail_mating_3 != 1, na.rm = TRUE)
num_orientation # 10 snails

# count number of snails that have a length measurement
num_length <- sum(!is.na(all_data_growth$Snail_length_1)) + sum(!is.na(all_data_growth$Snail_length_2)) + sum(!is.na(all_data_growth$Snail_length_3))
num_length # 6000 snails have length measurement

# wrong orientation rate
num_orientation / num_length # 0.001666667

####### calculate growth method 1: jar level average length if no snails death that week #########
all_data_growth <- all_data_growth %>%
  mutate(snail_length_sum = Snail_length_1 + Snail_length_2 + Snail_length_3,
         snail_count = (!is.na(Snail_length_1)) + (!is.na(Snail_length_2)) + (!is.na(Snail_length_3))) %>%
  arrange(Jar_num, Color, Week_num) %>%
  group_by(Jar_num, Color) %>%
  mutate(prev_snail_length_sum = lag(snail_length_sum),
         prev_snail_count = lag(snail_count),
         weekly_growth = if_else(snail_count == prev_snail_count,
                                 (snail_length_sum - prev_snail_length_sum) / snail_count,
                                 NA_real_)) %>%
  ungroup()

na_growth_jar_death <- count(all_data_growth %>%
  filter(Adult_death_1 == 1 & Adult_death_2 == 1 & Adult_death_3 == 1)) # 34 dead jars (valid NAs)

na_growth <- count(all_data_growth %>%
  filter(is.na(weekly_growth))) - na_growth_jar_death  # (1181-34)/2418 = 47.4% NA growth due to snail death in previous/current week, mating, or missing data

negative_growth<- count(all_data_growth %>%
  filter(weekly_growth < 0)) # 33/2418 = 1.4% negative growth
# hist(negative_growth$weekly_growth)
# table(negative_growth$Color,negative_growth$Week_num)

all_data_growth <- all_data_growth %>%
  filter(!is.na(weekly_growth) & weekly_growth >= 0) 

# GAMM with Autoregressive Error Structure
gamm_model_jar_avgerage <- gamm(weekly_growth ~ Treatment + s(Week_num, bs = "cs", k = 5), correlation = corAR1(), data = all_data_growth)
summary(gamm_model_jar_avgerage$gam)


# plot weekly growth
ggplot(all_data_growth, aes(x = factor(Week_num), y = weekly_growth, color = Treatment)) +
  geom_boxplot() +
  labs(title = "Weekly Snail Growth by Treatment Group",
       x = "Week Number",
       y = "Average Weekly Growth per Snail",
       color = "Treatment Group") +
  scale_color_manual(values = c("green", "yellow", "orange", "darkgray")) +
  theme_minimal() +
  facet_wrap(~ Treatment)

#check model Independence of Errors assumption
gamm_resid <- resid(gamm_model_jar_avgerage$lme, type = "normalized")
acf(gamm_resid, main = "ACF of Residuals") # residuals are approximately independent, sufficient model adjustment

# check Homoscedasticity assumption
plot(gamm_model_jar_avgerage$gam$fitted.values, residuals(gamm_model_jar_avgerage$gam), 
     xlab = "Predicted Values", ylab = "Residuals")
abline(h = 0, col = "red")

#check model normality of residual assumption
qqnorm(gamm_resid, main = "Q-Q Plot of Residuals")



####### (re-run correction chunk before this chunk) calculate growth method 2: growth data extraction by mannual matching #########
# clean up the mannual matching data
all_data_growth$snail_matching_1 <- sub(".*_", "", all_data_growth$snail_matching_1)
all_data_growth$snail_matching_2 <- sub(".*_", "", all_data_growth$snail_matching_2)
all_data_growth$snail_matching_3 <- sub(".*_", "", all_data_growth$snail_matching_3)
# Convert Unidentified snails NAs with 0
all_data_growth$Unidentifed_1 <- as.numeric(all_data_growth$Unidentifed_1)
all_data_growth$Unidentifed_1 <- ifelse(is.na(all_data_growth$Unidentifed_1), 0, all_data_growth$Unidentifed_1)
all_data_growth$Unidentifed_2 <- as.numeric(all_data_growth$Unidentifed_2)
all_data_growth$Unidentifed_2 <- ifelse(is.na(all_data_growth$Unidentifed_2), 0, all_data_growth$Unidentifed_2)
all_data_growth$Unidentifed_3 <- as.numeric(all_data_growth$Unidentifed_3)
all_data_growth$Unidentifed_3 <- ifelse(is.na(all_data_growth$Unidentifed_3), 0, all_data_growth$Unidentifed_3)

# step 2: check the ratio of mannually identifiable to unidentifiable
# 92 snails unidentifiable
un_id <- sum(all_data_growth$Unidentifed_1==1) + sum(all_data_growth$Unidentifed_2==1) + sum(all_data_growth$Unidentifed_3==1) 
# 7161 snails are 0
id_0 <- sum(all_data_growth$Unidentifed_1==0) + sum(all_data_growth$Unidentifed_2==0) + sum(all_data_growth$Unidentifed_3==0)
# 558 snails from week 1 are not applicable dispite 0
id_0_week1 <- sum(all_data_growth$Week_num==1)*3
# snails that are dead are not applicable distite 0
id_0_dead <- sum(all_data_growth$Adult_death_1==1) + sum(all_data_growth$Adult_death_2==1) + sum(all_data_growth$Adult_death_3==1)
# 6603 snails are identifiable
id <- id_0 - id_0_week1
# identification rate: 
id/(id + un_id) # 98.6%

# creat a function to get growth data
calculate_growth <- function(data) {
  data <- data %>%
    arrange(Jar_num, Color, Week_num) %>%
    group_by(Jar_num, Color) %>%
    mutate(
      growth_1 = ifelse(!is.na(snail_matching_1), Snail_length_1 - case_when(
        snail_matching_1 == 1 ~ lag(Snail_length_1),
        snail_matching_1 == 2 ~ lag(Snail_length_2),
        snail_matching_1 == 3 ~ lag(Snail_length_3),
        TRUE ~ NA_real_
      ), NA),
      growth_2 = ifelse(!is.na(snail_matching_2), Snail_length_2 - case_when(
        snail_matching_2 == 1 ~ lag(Snail_length_1),
        snail_matching_2 == 2 ~ lag(Snail_length_2),
        snail_matching_2 == 3 ~ lag(Snail_length_3),
        TRUE ~ NA_real_
      ), NA),
      growth_3 = ifelse(!is.na(snail_matching_3), Snail_length_3 - case_when(
        snail_matching_3 == 1 ~ lag(Snail_length_1),
        snail_matching_3 == 2 ~ lag(Snail_length_2),
        snail_matching_3 == 3 ~ lag(Snail_length_3),
        TRUE ~ NA_real_
      ), NA)
    ) %>%
    ungroup()
  
  return(data)
}

all_data_growth <- calculate_growth(all_data_growth)  


# Transform to long format
all_data_growth <- all_data_growth %>%
  pivot_longer(
    cols = starts_with("growth_"),
    names_to = "snail_id",
    names_prefix = "growth_",
    values_to = "match_growth"
  ) %>% 
  filter(!is.na(match_growth) & match_growth >= 0)

summary(all_data_growth$match_growth) 
# filtered out: 2078 NAs due to not identifiable (92) + week 1 (558) + no length:dead/mating/missing for this or previous week(1428)
# filtered out: 1151 negative growths

# na_growth <- all_data_growth %>% filter(is.na(match_growth))
# table(na_growth$Color, na_growth$Week_num)
# 
# negative_growth <- all_data_growth %>% filter(match_growth < 0)
# table(negative_growth$Color, negative_growth$Week_num)

# GAMM analysis will only use rows that have a match_growth value
gamm_model <- gamm(match_growth ~ Treatment + s(Week_num, bs = "cs", k = 5), 
                   random = list(Jar_num = ~1), 
                   correlation = corAR1(), 
                   data = all_data_growth)

# You can then inspect the model:
summary(gamm_model$gam) # 4025 sample size

# plot weekly growth
ggplot(all_data_growth, aes(x = factor(Week_num), y = match_growth, color = Treatment)) +
  geom_boxplot() +
  labs(title = "Weekly Snail Growth by Treatment Group",
       x = "Week Number",
       y = "Average Weekly Growth per Snail",
       color = "Treatment Group") +
  scale_color_manual(values = c("green", "yellow", "orange", "darkgray")) +
  theme_minimal() +
  facet_wrap(~ Treatment)

#check model Independence of Errors assumption
gamm_resid <- resid(gamm_model$lme, type = "normalized")
acf(gamm_resid, main = "ACF of Residuals") # residuals are approximately independent, sufficient model adjustment

# check Homoscedasticity assumption: not perfect
plot(gamm_model$gam$fitted.values, residuals(gamm_model$gam), 
     xlab = "Predicted Values", ylab = "Residuals")
abline(h = 0, col = "red")

#check model normality of residual assumption: not good
qqnorm(gamm_resid, main = "Q-Q Plot of Residuals")

