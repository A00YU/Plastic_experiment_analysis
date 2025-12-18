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
library(broom)
library(modelsummary)

####### load data#########
source("/Users/aoyu/Desktop/Snail_Data/R1_data_cleaning.R") # cleaned master data sheet called all_data

all_data_hatch <- all_data %>% 
  filter(Week_num != 1) %>% # Remove week 1 given no eggs yet 
  filter(Week_num != 13) %>% # & 13 with no hatching check
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
summary(all_data_hatch$Unhatched_egg_num) # enrolled and unhatched (sometimes a few snails are already hatched before enrolled over night)
summary(all_data_hatch$Hatch_check_1) # hatched num at w1 check point
summary(all_data_hatch$Hatch_check_2) # hatched num at w3 check point (2nd check point)

# when Enrolled_egg_num = 0, Two different senarios should be considered:
# 1) 0 egg mass to enroll = a NA hatch success, which is when Total_egg_mass_num = 0 & Enrolled_egg_num = 0, should be excluded from the study
# 2) have some egg mass but 0 valid egg to enroll = 0 hatch rate, which is Total_egg_mass_num > 0 & Enrolled_egg_num = 0, should be included in the study 

#Among those who COULD hatch, meaning with a viable structure:
#check for measurement error:
#negative values are spoted, where some of the "un-hatcheables" are hatched.
count(all_data_hatch, Unhatched_egg_num - Hatch_check_1 - Hatch_check_2 < 0) # 9/(2143+9) mistakes = error rate 0.4%; ~0.8% error rate if assuming it is same likelihood to misjudge an egg to be alive or dead
# head(all_data_hatch %>% filter(Unhatched_egg_num - Hatch_check_1 - Hatch_check_2 < 0))

# make sure that the number of enrolled eggs is greater than or equal to the number of unhatched eggs
no_lay<-all_data_hatch %>% filter(Total_egg_mass_num == 0)
valid_0 <- all_data_hatch %>% 
  filter(Total_egg_mass_num > 0 & Enrolled_egg_num == 0)
all_data_hatch <- all_data_hatch %>% 
  filter(Enrolled_egg_num >= Unhatched_egg_num & Enrolled_egg_num > 0) 
all_data_hatch <- rbind(all_data_hatch, valid_0) # add back the valid 0 hatch rate data

####### Offspring hatch success analysis #########
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

# estimate eggs per mass among treatments per weeks
egg_per_mass_estimates<- all_data_hatch %>% filter(Enrolled_egg_num != 0) %>% group_by(Treatment, Week_num) %>% summarise(egg_per_mass = median(Enrolled_egg_num)) #slight variation

# egg_per_mass_estimates<- all_data_hatch %>% filter(Enrolled_egg_num != 0) %>% group_by(Treatment) %>% summarise(egg_per_mass = median(Enrolled_egg_num)) # avg of 11 eggs per mass across treatments
ggplot(egg_per_mass_estimates, aes(x = Week_num, y = egg_per_mass, color = Treatment)) +
  geom_line() +
  geom_point()+
  scale_x_continuous(breaks = 1:12, labels = 1:12)

# save egg/mass as a variable for expected reproductive output model
write.csv(egg_per_mass_estimates,  file = "/Users/aoyu/Desktop/Snail_Data/FinalOutputs/R1_egg_per_mass_estimates.csv")

# add prevous week's snail size as covariate to control for
source("/Users/aoyu/Desktop/Snail_Data/R1_growth_analysis.R") # cleaned master data sheet called all_data_growth
covariate <- all_data_growth_length_covariate %>% select(Jar_num, Color, Week_num, prev_snail_length_avg)
all_data_hatch <- merge(all_data_hatch, covariate, by = c("Jar_num", "Color", "Week_num"), all.x = TRUE)
all_data_hatch <- all_data_hatch %>% filter(!is.na(prev_snail_length_avg)) # 116 missing prev_snail_length_avg

# binomial model with zero inflation
hatch_model <- glmmTMB( cbind(hatch, unhatch) ~ Treatment + prev_snail_length_avg +  snails_num,
                                     ziformula = ~Treatment + prev_snail_length_avg +  snails_num,
                                     family = binomial, data = all_data_hatch)

hatch_model_base <- glmmTMB( cbind(hatch, unhatch) ~ Treatment + prev_snail_length_avg,
                        ziformula = ~Treatment + prev_snail_length_avg,
                        family = binomial, data = all_data_hatch)

hatch_model_rm <- glmmTMB( cbind(hatch, unhatch) ~ Treatment + prev_snail_length_avg +  snails_num + (1|Jar_num),
                        ziformula = ~Treatment + prev_snail_length_avg +  snails_num + (1|Jar_num),
                        family = binomial, data = all_data_hatch)

AIC(hatch_model, hatch_model_base) # hatch_model with snail_num is better
AIC(hatch_model, hatch_model_rm) 

# Summary of the model
summary(hatch_model_rm)
# mannual data entry

# R1_hatch_model <- data.frame(
#   term = c("Virgin Micro", "Weathered Macro", "Weathered Macro + Virgin Micro", "Snail Shell Length Avg per Jar", "Number of Snails per Jar"),
#   estimate = c(-0.094980, -0.099884, 0.005612 , -2.734268, -0.086147),
#   standard.error = c(0.051304, 0.088186, 0.090227, 0.251424, 0.048805),
#   exp.estimate = c(exp(-0.094980), exp(-0.099884), exp(0.005612), exp(-2.734268), exp(-0.086147)),
#   p.value = c(0.0641,0.2574, 0.9504, "< 2e-16", 0.0775))
# 
# write.csv(R1_hatch_model,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R1_hatch_model.csv")

# # Post-hoc comparisons (emmeans)
emmeans_hatch <- emmeans(hatch_model_rm, pairwise ~ Treatment)
summary(emmeans_hatch)

# mannual data entry
# R1_hatch_pairwise <- data.frame(
#   term = c("Control - Virgin Micro","Control - Weathered Macro",
#            "Control - (Weathered Macro + Virgin Micro)",
#            "Virgin Micro - Weathered Macro",
#            "Virgin Micro - (Weathered Macro + Virgin Micro)", 
#            "Weathered Macro - (Weathered Macro + Virgin Micro)"),
#   estimate = c(0.09498, 0.09988, -0.00561 , 0.00490, -0.10059, 0.10550),
#   standard.error = c(0.0513, 0.0882, 0.0902, 0.0879, 0.0900, 0.1170),
#   exp.estimate = c(exp(0.09498), exp(0.09988), exp(-0.00561), exp(0.00490), exp(-0.10059), exp(0.10550)),
#   p.value = c(0.2494,0.6693, 0.9999, 0.9999, 0.6786, 0.8029))
# 
# write.csv(R1_hatch_pairwise,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R1_hatch_pairwise.csv")

######## Plot hatch success #########
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
R1_hatch_plot <- ggplot(fitness_summary, aes(x = Week_num, y = mean_fitness, color = Treatment, group = Treatment)) +
  geom_line() +  # add line connecting dots for each Color group
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  # geom_smooth(method = "loess", se = FALSE) +
  # geom_smooth(method = "loess", se = TRUE, aes(fill = Treatment), alpha = 0.2) +  # with se
  scale_color_manual(values = c("Control" = "#6baf78", 
                                "Micro" = "#5a9fd6", 
                                "Macro" = "#a5855f", 
                                "Macro+Micro" = "#9d7ca5"),
                     labels = c("Control", 
                                "Virgin Micro", 
                                "Biofouled Macro", 
                                "Biofouled Macro + Virgin Micro")) +
  scale_x_continuous(breaks = 1:13, labels = 1:13)+
  labs(
    x = "Weeks",
    y = "Mean Hatching Success Rate per Jar",) +
  theme_minimal()


# Save the plot to a PNG file
ggsave(filename = "R1_hatch_plot.png", plot = R1_hatch_plot, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/FinalOutputs", width = 9, height = 6, dpi = 300)

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

####### explore offspring hatch rate analysis #########
# Perform the calculation for hatch rate 
all_data_hatch$hatch_t1 <- all_data_hatch$Enrolled_egg_num - all_data_hatch$Unhatched_egg_num + all_data_hatch$Hatch_check_1 
all_data_hatch$hatch_t2 <- all_data_hatch$Hatch_check_2
all_data_hatch$unhatch_t1 <- pmax(0, all_data_hatch$Unhatched_egg_num - all_data_hatch$Hatch_check_1)
all_data_hatch$hatch_rate_t1 <- all_data_hatch$hatch_t1 / all_data_hatch$Enrolled_egg_num 
all_data_hatch$hatch_rate_t2 <- all_data_hatch$hatch_t2 / all_data_hatch$unhatch_t1

# Replace NaN values with 0 in the hatch_success column (due to 0/0 computation but still valid biological 0s)
# Replace the 9 rows of hatch rate bigger than 1 due to measurement error with 1
all_data_hatch <- all_data_hatch %>%
  mutate(hatch_rate_t1 = ifelse(is.nan(hatch_rate_t1), 0, hatch_rate_t1),
         hatch_rate_t1 = ifelse(hatch_rate_t1>1, 1, hatch_rate_t1),
         hatch_rate_t2 = ifelse(is.nan(hatch_rate_t2), 0, hatch_rate_t2),
         hatch_rate_t2 = ifelse(hatch_rate_t2>1, 1, hatch_rate_t2))

all_data_hatch$snails_num <- rowSums(all_data_hatch[, c("Adult_death_1", "Adult_death_2", "Adult_death_3")] == 0)
# Biomial: number of success out of total number of trials
# Calculate mean and variance of the dependent variable
mean(all_data_hatch$hatch_rate_t1)
var(all_data_hatch$hatch_rate_t1)
mean(all_data_hatch$hatch_rate_t2)
var(all_data_hatch$hatch_rate_t2)
# both Variance < Mean

# binomial regression model with zero inflation
hatchrate_t1 <- glmmTMB( cbind(hatch_t1, Enrolled_egg_num) ~ Treatment + prev_snail_length_avg +  snails_num ,
                                     ziformula = ~Treatment + prev_snail_length_avg +  snails_num,
                                     family = binomial, data = all_data_hatch)

hatchrate_t2 <- glmmTMB( cbind(hatch_t2, unhatch_t1) ~ Treatment + prev_snail_length_avg +  snails_num ,
                                      ziformula = ~Treatment + prev_snail_length_avg +  snails_num,
                                      family = binomial, data = all_data_hatch)

summary(hatchrate_t1)
summary(hatchrate_t2) # nothing interesting
