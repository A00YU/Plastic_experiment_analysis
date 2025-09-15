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

####### load data#########
source("/Users/aoyu/Desktop/Snail_Data/R1_data_cleaning.R") # cleaned master data sheet called all_data

####### Offspring hatching analysis #########
all_data_hatch <- all_data %>% 
  filter(Week_num != 1) %>% # Remove week 1 given no eggs yet
  select(Jar_num, Color, Week_num, Total_egg_mass_num, Enrolled_egg_num, Unhatched_egg_num, Hatch_check_1, Hatch_check_2, Treatment,Adult_death_1, Adult_death_2, Adult_death_3)%>%  
  mutate(Hatch_check_1 = as.numeric(Hatch_check_1), Hatch_check_2 = as.numeric(Hatch_check_2))

# table(all_data_hatch$Treatment,all_data_hatch$Week_num)

#replace all NA with 0
all_data_hatch$Enrolled_egg_num[is.na(all_data_hatch$Enrolled_egg_num)] <- 0
all_data_hatch$Unhatched_egg_num[is.na(all_data_hatch$Unhatched_egg_num)] <- 0
all_data_hatch$Hatch_check_1[is.na(all_data_hatch$Hatch_check_1)] <- 0
all_data_hatch$Hatch_check_2[is.na(all_data_hatch$Hatch_check_2)] <- 0

# Check for missing or problematic values (aim to enroll around 10 eggs & we see median = 10)
summary(all_data_hatch$Enrolled_egg_num) # enrolled num
summary(all_data_hatch$Unhatched_egg_num) # enrolled and unhatched (sometimes a few snails are already hatched before enrolled)
summary(all_data_hatch$Hatch_check_1) # hatched num at w1 check point
summary(all_data_hatch$Hatch_check_2) # hatched num at w3 check point (2nd check point)

# when Enrolled_egg_num = 0, Two different senarios should be considered:
# 1) 0 egg mass to enroll = a NA hatch success, which is when Total_egg_mass_num = 0 & Enrolled_egg_num = 0, should be excluded from the study
# 2) have some egg mass but 0 valid egg to enroll = 0 hatch rate, which is Total_egg_mass_num > 0 & Enrolled_egg_num = 0, should be included in the study 

#Among those who COULD hatch, meaning with a viable structure:
#check for measurement error:
#negative values are spoted, where some of the "un-hatcheables" are hatched.
count(all_data_hatch, Unhatched_egg_num - Hatch_check_1 - Hatch_check_2 < 0) # 11/(2221+11) mistakes = error rate 0.5%; ~1.00% error rate if assuming it is same likelihood to misjudge an egg to be alive or dead
# head(all_data_hatch %>% filter(Unhatched_egg_num - Hatch_check_1 - Hatch_check_2 < 0))

# make sure that the number of enrolled eggs is greater than or equal to the number of unhatched eggs
valid_0 <- all_data_hatch %>% 
  filter(Total_egg_mass_num > 0 & Enrolled_egg_num == 0)
all_data_hatch <- all_data_hatch %>% 
  filter(Enrolled_egg_num >= Unhatched_egg_num & Enrolled_egg_num > 0) 
all_data_hatch <- rbind(all_data_hatch, valid_0) # add back the valid 0 hatch rate data

# Perform the calculation for hatch success 
all_data_hatch$hatch <- all_data_hatch$Enrolled_egg_num - all_data_hatch$Unhatched_egg_num + all_data_hatch$Hatch_check_1 + all_data_hatch$Hatch_check_2
all_data_hatch$unhatch <- pmax(0, all_data_hatch$Unhatched_egg_num - all_data_hatch$Hatch_check_1 - all_data_hatch$Hatch_check_2)
all_data_hatch$hatch_success <- all_data_hatch$hatch / (all_data_hatch$hatch + all_data_hatch$unhatch) 

# Replace NaN values with 0 in the hatch_success column (due to 0/0 computation but still valid biological 0s)
all_data_hatch <- all_data_hatch %>%
  mutate(hatch_success = ifelse(is.nan(hatch_success), 0, hatch_success))

# Biomial: number of success out of total number of trials
# Calculate mean and variance of the dependent variable
mean_value <- mean(all_data_hatch$hatch_success)
variance_value <- var(all_data_hatch$hatch_success)

# Print the results
cat("Mean :", mean_value, "\n")
cat("Variance :", variance_value, "\n") # Variance < Mean

all_data_hatch$snails_num <- rowSums(all_data_hatch[, c("Adult_death_1", "Adult_death_2", "Adult_death_3")] == 0)

# binomial model with zero inflation
hatch_model_withWeek_num <- glmmTMB( cbind(hatch, unhatch) ~ Treatment + Week_num +  snails_num ,
                                     ziformula = ~Treatment + Week_num +  snails_num,
                                     family = binomial, data = all_data_hatch)

hatch_model <- glmmTMB( cbind(hatch, unhatch) ~ Treatment +  snails_num ,
                        ziformula = ~Treatment +  snails_num,
                        family = binomial, data = all_data_hatch)

anova(hatch_model_withWeek_num, hatch_model)
# smaller the AIC the better, so the model with parent age is better

# Summary of the model
summary(hatch_model_withWeek_num)

# # Post-hoc comparisons (emmeans)
emmeans_hatch <- emmeans(hatch_model_withWeek_num, pairwise ~ Treatment)
summary(emmeans_hatch)

######## Plot #########
# Summarize the data by Color and Week_num to calculate mean and 95% CI
fitness_summary <- all_data_hatch %>%
  group_by(Treatment, Week_num) %>%
  summarise(
    mean_fitness = mean(hatch_success, na.rm = TRUE),
    sd_fitness   = sd(hatch_success, na.rm = TRUE),
    n          = n(),
    se         = sd_fitness / sqrt(n),
    ci_lower   = mean_fitness - 1.96 * se,
    ci_upper   = mean_fitness + 1.96 * se
  ) %>%
  ungroup()

# Plot across weeks
ggplot(fitness_summary, aes(x = Week_num, y = mean_fitness, color = Treatment, group = Treatment)) +
  geom_line() +  # add line connecting dots for each Color group
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  # geom_smooth(method = "loess", se = FALSE) +
  # geom_smooth(method = "loess", se = TRUE, aes(fill = Treatment), alpha = 0.2) +  # with se
  scale_color_manual(values = c("Control" = "green", "Micro" = "yellow", "Macro" = "orange", "Macro+Micro" = "darkgray")) +
  scale_fill_manual(values = c("Control" = "green", "Micro" = "yellow", "Macro" = "orange", "Macro+Micro" = "darkgray")) +
  scale_x_continuous(breaks = 1:13, labels = 1:13)+
  labs(
    x = "Week Number",
    y = "Mean Hatching Success Rate",
    title = "Offsping Fitness Over Time by Treatment"
  ) +
  theme_minimal()

# plots with summary statistics over time masks timely variation:

# ggplot(all_data_hatch, aes(x = factor(Color, levels = c("G", "Y", "O", "W")), 
#                            y = hatch_success, fill = Color)) +
#   geom_boxplot() +
#   scale_x_discrete(labels = c("G" = "Negative Control",
#                               "Y" = "Microplastic",
#                               "O" = "Macroplastic",
#                               "W" = "Microplastic & Macroplastic")) +
#   scale_fill_manual(values = c("G" = "green",
#                                "Y" = "yellow",
#                                "O" = "orange",
#                                "W" = "darkgray")) +
#   theme_minimal() +
#   labs(title = "Offspring Hatching Success by Treatment Group",
#        x = "Treatment Group",
#        y = "Hatching Success Rate") +
#   geom_signif(comparisons = list(c("G", "Y")),
#               annotations = c("***"),
#               y_position = c(0.95),
#               tip_length = 0,
#               bracket.size = 0) +
#   scale_y_continuous(limits = c(0, 1))
