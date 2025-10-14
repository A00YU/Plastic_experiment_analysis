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
source("/Users/aoyu/Desktop/Snail_Data/R2_data_cleaning.R") # cleaned master data sheet called all_data_r2

####### measurement error correction due to orientation #########
all_data_growth <- all_data_r2 %>%
  select(Color, Jar_num, Week_num, starts_with("Adult_death_"), starts_with("Snail_length_"), starts_with("snail_mating_"), starts_with("snail_orientation_"), Treatment)

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
num_mating # 130 snails

# if not mating but simply orientation is wrong, correct it with a correction ratio: correct/wrong = 1.030 = 1.02997
# round the product into 3 decimal places
all_data_growth <- all_data_growth %>%
  mutate(
    Snail_length_1 = round(ifelse(snail_orientation_1 == 1 & snail_mating_1 != 1, Snail_length_1 * 1.03, Snail_length_1), 3),
    Snail_length_2 = round(ifelse(snail_orientation_2 == 1 & snail_mating_2 != 1, Snail_length_2 * 1.03, Snail_length_2), 3),
    Snail_length_3 = round(ifelse(snail_orientation_3 == 1 & snail_mating_3 != 1, Snail_length_3 * 1.03, Snail_length_3), 3)
  )

# count number of snails with initial wrong orientation (not including mating as "wrong" orientation)
num_orientation <- sum(all_data_growth$snail_orientation_1 == 1 & all_data_growth$snail_mating_1 != 1, na.rm = TRUE) + 
  sum(all_data_growth$snail_orientation_2 == 1 & all_data_growth$snail_mating_2 != 1, na.rm = TRUE) + 
  sum(all_data_growth$snail_orientation_3 == 1 & all_data_growth$snail_mating_3 != 1, na.rm = TRUE)
num_orientation # 1279 snails

# count number of snails that have a length measurement
num_length <- sum(!is.na(all_data_growth$Snail_length_1)) + sum(!is.na(all_data_growth$Snail_length_2)) + sum(!is.na(all_data_growth$Snail_length_3))
num_length # 6489 snails have length measurement

# wrong orientation rate
num_orientation / num_length # 19.7%

# to get the snail_length_avg for each jar at each week (including the ones with snail death)
snail_length_avg_df <- all_data_growth %>%
  mutate(
    Snail_length_1 = ifelse(Adult_death_3 == 1, 0, as.numeric(Snail_length_1)),
    Snail_length_2 = ifelse(Adult_death_2 == 1, 0, as.numeric(Snail_length_2)),
    Snail_length_3 = ifelse(Adult_death_1 == 1, 0, as.numeric(Snail_length_3))) %>% 
  arrange(Jar_num, Color, Week_num) %>%
  mutate(snail_length_sum = Snail_length_1 + Snail_length_2 + Snail_length_3,
         snail_count = rowSums(!is.na(cbind(Adult_death_1, Adult_death_2, Adult_death_3)) & cbind(Adult_death_1, Adult_death_2, Adult_death_3) != 1)) %>%
  arrange(Jar_num, Color, Week_num) %>%
  group_by(Jar_num, Color) %>%
  mutate(prev_snail_length_sum = lag(snail_length_sum),
         prev_snail_count = lag(snail_count),
         weekly_growth = if_else(snail_count == prev_snail_count,
                                 (snail_length_sum - prev_snail_length_sum) / snail_count,
                                 NA_real_),
         prev_snail_length_avg = prev_snail_length_sum/prev_snail_count) %>%
  ungroup() 

####### calculate growth method 1: jar level average length if no snails death that week #########
all_data_growth <- all_data_growth %>%
  mutate(snail_length_sum = Snail_length_1 + Snail_length_2 + Snail_length_3,
         snail_count = rowSums(!is.na(cbind(Adult_death_1, Adult_death_2, Adult_death_3)) & cbind(Adult_death_1, Adult_death_2, Adult_death_3) != 1)) %>%
  arrange(Jar_num, Color, Week_num) %>%
  group_by(Jar_num, Color) %>%
  mutate(prev_snail_length_sum = lag(snail_length_sum),
         prev_snail_count = lag(snail_count),
         weekly_growth = if_else(snail_count == prev_snail_count,
                                 (snail_length_sum - prev_snail_length_sum) / snail_count,
                                 NA_real_)) %>%
  ungroup()

# join the prev_snail_length_avg from snail_length_avg_df to all_data_growth
all_data_growth <- all_data_growth %>%
  left_join(snail_length_avg_df %>% select(Jar_num, Color, Week_num, prev_snail_length_avg), by = c("Jar_num", "Color", "Week_num"))

# Create a new column for combined groups
all_data_growth$Group <- factor(with(all_data_growth, 
                                     ifelse(Color %in% c("Y", "O"), "Seasoned Macro",
                                            ifelse(Color %in% c("W-1", "W-3"), "Virgin Macro","Control"))))

# Verify the new grouping
table(all_data_growth$Group, all_data_growth$Week_num)

# data counts for each senario
na_growth_jar_death <- count(all_data_growth %>%
                               filter(Adult_death_1 == 1 & Adult_death_2 == 1 & Adult_death_3 == 1)) # 131 dead jars (valid NAs)

na_growth <- count(all_data_growth %>%
                     filter(is.na(weekly_growth))) - na_growth_jar_death  # 1192/2704 = 44% NA growth due to snail death in previous/current week, mating, or missing data

negative_growth<- all_data_growth %>%
                          filter(weekly_growth < 0)

negative_growth <-count(all_data_growth %>%
                          filter(weekly_growth < 0)) # 243/2704 = 9% negative growth
# hist(negative_growth$weekly_growth)
# table(negative_growth$Color,negative_growth$Week_num)

# save a copy of all_data_growth before filtering NAs and negative growths for later test
all_data_growth_length_covariate <- all_data_growth

all_data_growth <- all_data_growth %>%
  mutate(weekly_growth = ifelse(weekly_growth < 0, 0, weekly_growth)) %>% # include negative growth as 0 growth
  filter(!is.na(weekly_growth) & weekly_growth >= 0) # 1381 data points for growth analysis

# Calculate the IQR
Q1 <- quantile(all_data_growth$weekly_growth, 0.25, na.rm = TRUE)
Q3 <- quantile(all_data_growth$weekly_growth, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1

# Define the lower and upper bounds for outliers
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Identify outliers
outliers <- all_data_growth %>%
  filter(weekly_growth < lower_bound | weekly_growth > upper_bound)

# GAMM with Autoregressive Error Structure
gamm_model_jar_avgerage <- gamm(weekly_growth ~ Group + s(prev_snail_length_avg, bs = "cs", k = 5), correlation = corAR1(), data = all_data_growth)
summary(gamm_model_jar_avgerage$gam)

gamm_model_jar_avgerage1 <- gamm(weekly_growth ~ Group + s(Week_num, bs = "cs", k = 5), correlation = corAR1(), data = all_data_growth)
summary(gamm_model_jar_avgerage1$gam)
AIC(gamm_model_jar_avgerage$lme, gamm_model_jar_avgerage1$lme) # model with prev length is better

gamm_model_jar_avgerage_abundance <- gamm(weekly_growth ~ Treatment + s(prev_snail_length_avg, bs = "cs", k = 5), correlation = corAR1(), data = all_data_growth)
summary(gamm_model_jar_avgerage_abundance$gam)

gamm_model_jar_avgerage_abundance1 <- gamm(weekly_growth ~ Treatment + s(Week_num, bs = "cs", k = 5), correlation = corAR1(), data = all_data_growth)

AIC(gamm_model_jar_avgerage_abundance$lme, gamm_model_jar_avgerage_abundance1$lme) # model with prev length is better

# plot weekly growth
ggplot(all_data_growth, aes(x = factor(Week_num), y = weekly_growth, color = Group)) +
  geom_boxplot() +
  labs(title = "Weekly Snail Growth by Treatment Group",
       x = "Week Number",
       y = "Average Weekly Growth per Snail",
       color = "Treatment Group") +
  scale_color_manual(values = c("#6baf78", "#a5855f", "#8a8a8a")) +
  theme_minimal() +
  facet_wrap(~ Group)

#check model Independence of Errors assumption
gamm_resid <- resid(gamm_model_jar_avgerage$lme, type = "normalized")
acf(gamm_resid, main = "ACF of Residuals") # residuals are approximately independent, sufficient model adjustment

# check Homoscedasticity assumption
plot(gamm_model_jar_avgerage$gam$fitted.values, residuals(gamm_model_jar_avgerage$gam), 
     xlab = "Predicted Values", ylab = "Residuals")
abline(h = 0, col = "red")

#check model normality of residual assumption
qqnorm(gamm_resid, main = "Q-Q Plot of Residuals")

################check for potential errors or anomalies ##################
test<- all_data_growth %>% filter(Group == "Seasoned Macro" | Group == "Virgin Macro") # hist comparison
test<- all_data_growth %>% filter(Group == "Virgin Macro") # single hist
ggplot(all_data_growth, aes(x = weekly_growth, color = Group)) +
  geom_histogram() +
  theme_minimal()  +
  facet_wrap(~ Group) 

ggplot(test, aes(x = weekly_growth, color = Group)) +
  geom_histogram() +
  theme_minimal()  +
  facet_wrap(~ Week_num) 
       