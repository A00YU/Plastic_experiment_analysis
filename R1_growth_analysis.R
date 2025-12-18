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
library(broom)
library(modelsummary)

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

# specifying the 17 jars (19 snails) + 6 jars to eliminate due to missing/duplicate picture or within treatment misplacement
eliminate_jars <- data.frame(
  Week_num = c(1, 1, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 9, 8, 8, 9, 9, 11, 11), 
  Color = c("W", "W", "O", "O", "O", "Y", "O", "O", "O", "O", "Y", "Y", "O", "W", "Y", "Y", "G", "G", "G", "G", "G", "G", "G"),  
  Jar_num = c(34, 55, 16, 18, 32, 11, 3, 7, 12, 40, 35, 57, 29, 19, 39, 70, 20, 6, 67, 33, 45, 8, 11) 
)

# Filter the data frame to exclude the specified jars in the specific week, color, and jar number combinations
all_data_growth <- all_data_growth[!with(all_data_growth, paste(Week_num, Color, Jar_num) %in% paste(eliminate_jars$Week_num, eliminate_jars$Color, eliminate_jars$Jar_num)), ]

# count number of snails that have a length measurement
num_length <- sum(!is.na(all_data_growth$Snail_length_1)) + sum(!is.na(all_data_growth$Snail_length_2)) + sum(!is.na(all_data_growth$Snail_length_3))
num_length # 5958 snails have length measurement

# wrong orientation rate
num_orientation / num_length # 0.001678416

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
         snail_count = (!is.na(Snail_length_1)) + (!is.na(Snail_length_2)) + (!is.na(Snail_length_3))) %>%
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

na_growth_jar_death <- count(all_data_growth %>%
  filter(Adult_death_1 == 1 & Adult_death_2 == 1 & Adult_death_3 == 1)) # 34 dead jars (valid NAs)

na_growth <- count(all_data_growth %>%
  filter(is.na(weekly_growth))) - na_growth_jar_death  # (1181-34)/2418 = 47.4% NA growth due to snail death in previous/current week, mating, or missing data

negative_growth<- count(all_data_growth %>%
  filter(weekly_growth < 0)) # 33/2418 = 1.4% negative growth

negative_growth<- all_data_growth %>%
                          filter(weekly_growth < 0)
# hist(negative_growth$weekly_growth)
table(negative_growth$Color,negative_growth$Week_num)

# save a copy of all_data_growth before filtering NAs and negative growths for later test
all_data_growth_length_covariate <- all_data_growth

all_data_growth <- all_data_growth %>%
  mutate(weekly_growth = ifelse(weekly_growth < 0, 0, weekly_growth)) %>% 
  filter(!is.na(weekly_growth) & weekly_growth >= 0) 

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
gamm_model_jar_avgerage <- gamm(weekly_growth ~ Treatment + s(prev_snail_length_avg, bs = "cs", k = 5), correlation = corAR1(), data = all_data_growth)
summary(gamm_model_jar_avgerage$gam) # 1234 sample size (unit in cm)

# Saves model results as a CSV
R1_growth_model <- data.frame(
  term = c("Virgin Micro", "Biofouled Macro", "Biofouled Macro + Virgin Micro"),
  estimate = c(0.0040116, 0.0093312, 0.0051165),
  standard.error = c(0.0011426, 0.0020414,0.0017913),
  p.value = c(0.000463,5.34e-06, 0.004358))

write.csv(R1_growth_model,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R1_growth_model.csv")

# plot weekly growth
R1_growth_plot <- ggplot(all_data_growth, aes(x = factor(Week_num), y = weekly_growth, color = Treatment)) +
  geom_boxplot() +
  labs(x = "Weeks",
       y = "Mean Shell Length Growth per Snail (cm)",
       color = "Treatment") +
  scale_color_manual(values = c("Control" = "#6baf78", 
                                "Micro" = "#5a9fd6", 
                                "Macro" = "#a5855f", 
                                "Macro+Micro" = "#9d7ca5"),
                     labels = c("Control", 
                                "Virgin Micro", 
                                "Biofouled Macro", 
                                "Biofouled Macro + Virgin Micro")) +
  theme_minimal() +
  theme(strip.text.x = element_blank())+
  facet_wrap(~ Treatment)

# Save the plot to a PNG file
ggsave(filename = "R1_growth_plot.png", plot = R1_growth_plot, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)

#check model Independence of Errors assumption
gamm_resid <- resid(gamm_model_jar_avgerage$lme, type = "normalized")
acf(gamm_resid, main = "ACF of Residuals") # residuals are approximately independent, sufficient model adjustment

# check Homoscedasticity assumption
plot(gamm_model_jar_avgerage$gam$fitted.values, residuals(gamm_model_jar_avgerage$gam), 
     xlab = "Predicted Values", ylab = "Residuals")
abline(h = 0, col = "red")

#check model normality of residual assumption
qqnorm(gamm_resid, main = "Q-Q Plot of Residuals")

source("/Users/aoyu/Desktop/Snail_Data/R1_data_cleaning.R") # cleaned master data sheet called all_data
# final size and whether diff between treatment:
randomization_check <- all_data %>% 
  filter(Week_num == 13) %>%
  select(Color, Jar_num, Snail_length_1, Snail_length_2, Snail_length_3) %>% 
  filter(Color == "G"|Color == "Y" |Color ==  "W"|Color ==  "O") 

table(randomization_check$Color, randomization_check$Jar_num)

final_size_long <- randomization_check %>%
  pivot_longer(cols = starts_with("Snail_length_"),
               names_to = "Snail_ID",
               values_to = "Snail_length") %>% 
  mutate(Snail_length = as.numeric(Snail_length))

# ANOVA to check if initial snail size differs by treatment
anova_result <- aov(Snail_length ~ Color, data = final_size_long)
summary(anova_result) # significant difference numerically (could be due to snail picture angle), but acceptable biologically
TukeyHSD(anova_result) # post-hoc test

# Boxplot to visualize the distribution of initial snail sizes by treatment
ggplot(final_size_long, aes(x = Color, y = Snail_length, color = Color)) +
  geom_boxplot() +
  labs(title = "Distribution of Initial Snail Sizes by Treatment",
       x = "Treatment Group",
       y = "Snail Length") +
  theme_minimal() 
