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
library(ggsignif)
library(mgcv)
library(car)

####### load data#########
data <- read.csv("/Users/aoyu/Desktop/Snail_Data/Seasoned plastic dimensions & algae measurements.csv")

####### chlorophyll a analysis #########
chlor.a <- data %>% 
  filter(test=="chlor.a") %>% 
  select(batch_number, test, sample_number, chlor_absolute_ug_per_unit,pheophytin_absolute_ug_per_unit,Chlor_pheophytin_ratio) %>%
  mutate(Chlor_pheophytin_ratio = as.numeric(Chlor_pheophytin_ratio),
         batch_method = ifelse(batch_number %in% c(1,2), "lab", "field"))

summary(chlor.a$Chlor_pheophytin_ratio) # >1, Healthy, active algae with chlorophyll a dominates; little degradation
hist(chlor.a$Chlor_pheophytin_ratio)

wilcox.test(Chlor_pheophytin_ratio ~ batch_method, data = chlor.a) # p<0.05, significantly differ by method

ggplot(chlor.a, aes(x = batch_method, y = Chlor_pheophytin_ratio, fill = batch_method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("lab" = "green", "field" = "lightgreen")) +
  labs(title = "Chlorophyll a. to pheophytin a. ratio",
       x = "seasoning method",
       y = "Chlorophyll a. to pheophytin a. ratio") +
  theme_minimal()

####### DM & AFDM analysis #########
AFDM <- data %>% 
  filter(test=="AFDM") %>% 
  select(batch_number, test, sample_number, DM_ug_per_unit, AFDM_ug_per_unit, AFDM_DM_ratio) %>% 
  mutate(AFDM_DM_ratio = as.numeric(AFDM_DM_ratio),
         DM_ug_per_unit = as.numeric(DM_ug_per_unit),
         batch_method = ifelse(batch_number %in% c(1,2), "lab", "field"))

summary(AFDM$AFDM_DM_ratio) # < 0.3, low to moderate organic contents with some algae + microbes and high inorganic such as sediments + mineral
hist(AFDM$AFDM_DM_ratio)

summary(AFDM$DM_ug_per_unit)
hist(AFDM$DM_ug_per_unit)

wilcox.test(AFDM_DM_ratio ~ batch_method, data = AFDM) # p<0.05, significantly differ by method
wilcox.test(DM_ug_per_unit ~ batch_method, data = AFDM) # p<0.05, significantly differ by method

ggplot(AFDM, aes(x = batch_method, y = AFDM_DM_ratio, fill = batch_method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("lab" = "green", "field" = "lightgreen")) +
  labs(title = "Ash free dry mass to dry mass ratio",
       x = "seasoning method",
       y = "Ash free dry mass to dry mass ratio") +
  theme_minimal() # Field less *organic matter* load then lab (0.5 times the abundance)

ggplot(AFDM, aes(x = batch_method, y = DM_ug_per_unit, fill = batch_method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("lab" = "green", "field" = "lightgreen")) +
  labs(title = "Dry mass",
       x = "seasoning method",
       y = "Dry mass (ug)") +
  theme_minimal() # Field much more periphyton load then lab (23 times the abundance)

####### autotrophic index: autotrophic to heterotrophic ratio analysis #########
chlor_absolute <- data %>% 
  select(batch_number, sample_number, chlor_absolute_ug_per_unit) %>%
  filter(!is.na(chlor_absolute_ug_per_unit), batch_number %in% c(1,2,3,4)) %>% 
  mutate(chlor_absolute_ug_per_unit = as.numeric(chlor_absolute_ug_per_unit),
         batch_method = ifelse(batch_number %in% c(1,2), "lab", "field"))

AFDM_absolute <- data %>% 
  select(batch_number, sample_number, AFDM_ug_per_unit) %>%
  filter(!is.na(AFDM_ug_per_unit), batch_number %in% c(1,2,3,4)) %>% 
  mutate(AFDM_ug_per_unit = as.numeric(AFDM_ug_per_unit),
         batch_method = ifelse(batch_number %in% c(1,2), "lab", "field"))

AI <- left_join(chlor_absolute, AFDM_absolute, by = c("batch_number", "sample_number", "batch_method")) %>% 
  mutate(AI =  AFDM_ug_per_unit / chlor_absolute_ug_per_unit)
hist(AI$AI)

wilcox.test(AI ~ batch_method, data = AI)

ggplot(AI, aes(x = batch_method, y = AI, fill = batch_method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("lab" = "green", "field" = "lightgreen")) +
  labs(title = "Autotrophic Index ",
       x = "seasoning method",
       y = "Autotrophic Index") +
  theme_minimal()

####### algae and lettuce CN analysis #########
CN_data <- read.csv("/Users/aoyu/Desktop/Snail_Data/periphyton_lettuce_CN_data.csv")
CN_data <- CN_data[!is.na(CN_data$average_CN), ]  # Remove rows with NA in average_CN

# Subset data 
lettuce_data <- subset(CN_data, food_type == "lettuce") # 18 samples
algae_data <- subset(CN_data, food_type == "algae") # 23 samples

# lettuce:
# Calculate IQR
Q1 <- quantile(lettuce_data$average_CN, 0.25)
Q3 <- quantile(lettuce_data$average_CN, 0.75)
IQR <- Q3 - Q1

# Define outlier thresholds
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Identify outliers: no outliers
outliers <- lettuce_data$average_CN < lower_bound | lettuce_data$average_CN > upper_bound

# algae:
# Calculate IQR
Q1 <- quantile(algae_data$average_CN, 0.25)
Q3 <- quantile(algae_data$average_CN, 0.75)
IQR <- Q3 - Q1

# Define outlier thresholds
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Identify outliers: 2 outliers
outliers <- algae_data$average_CN < lower_bound | algae_data$average_CN > upper_bound

# Subset data end up with 18 lettuce data and 21 algae data

# normality
qqnorm(lettuce_data$average_CN, main = "QQ Plot for CN - Lettuce") # normal
qqnorm(algae_data$average_CN, main = "QQ Plot for CN - Algae") # not normal


# equal variance check
leveneTest(average_CN ~ food_type, data = CN_data) #  p>0.05 the variances are equal across the groups 

# gives mean of both group:
# t.test(average_CN ~ food_type, data = CN_data, var.equal = FALSE) # not a valid test given non-normality

# nonparametric test:
wilcox.test(average_CN ~ food_type, data = CN_data) # p>0.05 not significantly differ


# plot the CN data
ggplot(CN_data, aes(x = food_type, y = average_CN, fill = food_type)) +
  geom_boxplot() +
  scale_fill_manual(values = c("algae" = "green", "lettuce" = "lightgreen")) +
  labs(title = "Carbon to Nitrogen Ratio by Food Type",
       x = "Food Type",
       y = "Average C:N Ratio") +
  theme_minimal()
