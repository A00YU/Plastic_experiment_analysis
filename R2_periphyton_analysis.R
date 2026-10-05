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
data <- read.csv("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/Seasoned plastic dimensions & algae measurements.csv")

####### chlorophyll a analysis #########
chlor.a <- data %>% 
  filter(test=="chlor.a") %>% 
  select(batch_number, test, sample_number, chlor_absolute_ug_per_unit,pheophytin_absolute_ug_per_unit,Chlor_pheophytin_ratio) %>%
  mutate(Chlor_pheophytin_ratio = as.numeric(Chlor_pheophytin_ratio),
         batch_method = ifelse(batch_number %in% c(1,2), "lab", "field"))

summary(chlor.a$Chlor_pheophytin_ratio) # >1, Healthy, active algae with chlorophyll a dominates; little degradation
hist(chlor.a$Chlor_pheophytin_ratio)

wilcox.test(Chlor_pheophytin_ratio ~ batch_method, data = chlor.a) # p<0.05, significantly differ by method

chlorophyll_pheophytin_ratio <- ggplot(chlor.a, aes(x = batch_method, y = Chlor_pheophytin_ratio, fill = batch_method)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6, color = "darkblue") +
  scale_fill_manual(values = c("lab" = "#4c78a8", "field" = "#4ead9f"),
                    labels = c("Field", "Lab")) +
  scale_x_discrete(labels = c(lab = "Lab", field = "Field"))+
  labs(x = "Biofouling method",
       y = "Chlorophyll a. to pheophytin a. ratio",
       fill = "Biofouling method") +
  theme_minimal()

chlorophyll_pheophytin_ratio

# ggsave(filename = "Chlorophyll_pheophytin_ratio.png", plot = chlorophyll_pheophytin_ratio, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/FinalOutputs", width = 6, height = 6, dpi = 300)
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

AFDM_DM_ratio <- ggplot(AFDM, aes(x = batch_method, y = AFDM_DM_ratio, fill = batch_method)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6, color = "darkblue") +
  scale_fill_manual(values = c("lab" = "#4c78a8", "field" = "#4ead9f"),
                    labels = c("Field", "Lab")) +
  scale_x_discrete(labels = c(lab = "Lab", field = "Field"))+
  labs(x = "Biofouling method",
       y = "Ash free dry mass to dry mass ratio",
       fill = "Biofouling method") +
  theme_minimal()

AFDM_DM_ratio

# ggsave(filename = "AFDM_DM_ratio.png", plot = AFDM_DM_ratio, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/FinalOutputs", width = 6, height = 6, dpi = 300)

DM_ug_per_unit <- ggplot(AFDM, aes(x = batch_method, y = DM_ug_per_unit, fill = batch_method)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6, color = "darkblue") +
  scale_fill_manual(values = c("lab" = "#4c78a8", "field" = "#4ead9f"),
                    labels = c("Field", "Lab")) +
  scale_x_discrete(labels = c(lab = "Lab", field = "Field"))+
  labs(x = "Biofouling method",
       y = "Periphyton dry mass (ug)",
       fill = "Biofouling method") +
  theme_minimal()

DM_ug_per_unit

# ggsave(filename = "DM_ug_per_unit.png", plot = DM_ug_per_unit, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/FinalOutputs", width = 6, height = 6, dpi = 300)

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

# cannot randomly match two results to calculate AI given not from exatly the same sample
# AI <- left_join(chlor_absolute, AFDM_absolute, by = c("batch_number", "sample_number", "batch_method")) %>% 
#   mutate(AI =  AFDM_ug_per_unit / chlor_absolute_ug_per_unit)

# instead, calculate mean of each batch and then calculate AI
chlor_mean <- chlor_absolute %>%
  group_by(batch_number, batch_method) %>%
  summarize(mean_chlor = mean(chlor_absolute_ug_per_unit, na.rm = TRUE))

AFDM_mean <- AFDM_absolute %>%
  group_by(batch_number, batch_method) %>%
  summarize(mean_AFDM = mean(AFDM_ug_per_unit, na.rm = TRUE))

AI <- left_join(chlor_mean, AFDM_mean, by = c("batch_number", "batch_method")) %>%
  mutate(AI =  mean_AFDM / mean_chlor)

summary(AI$AI) 
wilcox.test(AI ~ batch_method, data = AI)

AI <- ggplot(AI, aes(x = batch_method, y = AI, fill = batch_method)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6, color = "darkblue") +
  scale_fill_manual(values = c("lab" = "#4c78a8", "field" = "#4ead9f"),
                    labels = c("Field", "Lab")) +
  scale_x_discrete(labels = c(lab = "Lab", field = "Field"))+
  labs(x = "Biofouling method",
       y = "Autotrophic index",
       fill = "Biofouling method") +
  theme_minimal()

AI

# ggsave(filename = "AI.png", plot = AI, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/FinalOutputs", width = 6, height = 6, dpi = 300)

####### algae and lettuce CN analysis #########
CN_data <- read.csv("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/periphyton_lettuce_CN_data.csv")
CN_data <- CN_data[!is.na(CN_data$average_CN), ]  # Remove rows with NA in average_CN

# Subset data 
lettuce_data <- subset(CN_data, food_type == "lettuce") # 18 samples
algae_data <- subset(CN_data, food_type == "algae") # 23 samples
algae_data$batch_method <- c(rep("field", 12), rep("lab", 11))

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
CN_data %>%
  filter(food_type %in% c("algae", "lettuce")) %>%
  group_by(food_type) %>%
  summarise(
    n = sum(!is.na(average_CN)),
    mean_value = mean(average_CN, na.rm = TRUE),
    sd_value = sd(average_CN, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_value = sd_value / sqrt(n),
    t_critical = qt(0.975, df = pmax(n - 1, 1)),
    lower95 = if_else(n > 1, mean_value - t_critical * se_value, NA_real_),
    upper95 = if_else(n > 1, mean_value + t_critical * se_value, NA_real_)
  ) %>%
  select(food_type, n, mean_value, lower95, upper95)


# nonparametric test:
wilcox.test(average_CN ~ food_type, data = CN_data) # p>0.05 not significantly differ
wilcox.test(average_CN ~ batch_method, data = algae_data) # p<0.005  significantly differ

# plot the CN data on food type
CN_foodtype <- ggplot(CN_data, aes(x = food_type, y = average_CN, fill = food_type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6, color = "darkblue") +
  scale_fill_manual(values = c("algae" = "#a5855f", "lettuce" = "#2CA25F"),
                    labels = c("Periphyton", "Lettuce")) +
  scale_x_discrete(labels = c(algae = "Periphyton", lettuce = "Lettuce"))+
  labs(x = "Food type",
       y = "Mean C:N ratio",
       fill = "Food type") +
  theme_minimal()

CN_foodtype

# ggsave(filename = "CN_foodtype.png", plot = CN_foodtype, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/FinalOutputs", width = 6, height = 6, dpi = 300)

# plot the CN data on algae batch method
CN_lab_field <- ggplot(algae_data, aes(x = batch_method, y = average_CN, fill = batch_method)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6, color = "darkblue") +
  scale_fill_manual(values = c("lab" = "#4c78a8", "field" = "#4ead9f"),
                    labels = c("Field", "Lab")) +
  scale_x_discrete(labels = c(lab = "Lab", field = "Field"))+
  labs(x = "Biofouling method",
       y = "Mean C:N ratio",
       fill = "Biofouling method") +
  theme_minimal()

CN_lab_field

# ggsave(filename = "CN_lab_field.png", plot = CN_lab_field, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/FinalOutputs", width = 6, height = 6, dpi = 300)
