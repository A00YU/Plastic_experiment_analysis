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
source("/Users/aoyu/Desktop/Snail_Data/R2_data_cleaning.R") # cleaned master data sheet called all_data

all_data_hatch <- all_data_r2 %>% 
  filter(Week_num != 1) %>% # Remove week 1 given no eggs yet
  select(Jar_num, Color, Week_num, Total_egg_mass_num, Enrolled_egg_num, Unhatched_egg_num, Hatch_check_1, Hatch_check_2, Treatment,Adult_death_1, Adult_death_2, Adult_death_3)%>%  
  mutate(Enrolled_egg_num = as.numeric(Enrolled_egg_num), Unhatched_egg_num = as.numeric(Unhatched_egg_num),
         Hatch_check_1 = as.numeric(Hatch_check_1), Hatch_check_2 = as.numeric(Hatch_check_2))

all_data_hatch$snails_num <- rowSums(all_data_hatch[, c("Adult_death_1", "Adult_death_2", "Adult_death_3")] == 0)

table(all_data_hatch$Treatment,all_data_hatch$Week_num)

#replace all NA with 0
all_data_hatch$Enrolled_egg_num[is.na(all_data_hatch$Enrolled_egg_num)] <- 0
all_data_hatch$Unhatched_egg_num[is.na(all_data_hatch$Unhatched_egg_num)] <- 0
all_data_hatch$Hatch_check_1[is.na(all_data_hatch$Hatch_check_1)] <- 0
all_data_hatch$Hatch_check_2[is.na(all_data_hatch$Hatch_check_2)] <- 0

# Check for missing or problematic values (aim to enroll around 10 eggs & we see median = 11)
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
count(all_data_hatch, Unhatched_egg_num - Hatch_check_1 - Hatch_check_2 < 0) # 8/(2488+8) mistakes = error rate 0.3%; ~0.6% error rate if assuming it is same likelihood to misjudge an egg to be alive or dead

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

# Combine Y&O into SM and W-1&W-3 into VM
all_data_hatch$Group <- factor(with(all_data_hatch, 
                                    ifelse(Color %in% c("Y", "O"), "Seasoned Macro",
                                           ifelse(Color %in% c("W-1", "W-3"), "Virgin Macro","Control"))))

# Select week 1 to 12 to match 12 weeks of data in Round 1 
all_data_hatch <- all_data_hatch %>% filter(Week_num != 13)

# Verify the new grouping
table(all_data_hatch$Group, all_data_hatch$Week_num)

####### Offspring hatching analysis for each treatment group with varying abundance #########
# Biomial: number of success out of total number of trials
hist(all_data_hatch$hatch_success)

# Calculate mean and variance of the dependent variable
mean_value <- mean(all_data_hatch$hatch_success)
variance_value <- var(all_data_hatch$hatch_success)

# Print the results
cat("Mean :", mean_value, "\n")
cat("Variance :", variance_value, "\n") # Variance < Mean

# estimate eggs per mass among treatments per weeks
egg_per_mass_estimates<- all_data_hatch %>% filter(Enrolled_egg_num != 0) %>% group_by(Treatment, Week_num) %>% summarise(egg_per_mass = mean(Enrolled_egg_num)) #slight variation
egg_per_mass_estimates<- all_data_hatch %>%  group_by(Treatment, Week_num) %>% summarise(egg_per_mass = mean(Enrolled_egg_num)) #slight variation

# egg_per_mass_estimates<- all_data_hatch %>% filter(Enrolled_egg_num != 0) %>% group_by(Treatment) %>% summarise(egg_per_mass = median(Enrolled_egg_num)) # avg of 11 eggs per mass across treatments
ggplot(egg_per_mass_estimates, aes(x = Week_num, y = egg_per_mass, color = Treatment)) +
  geom_line() +
  geom_point()+
  scale_x_continuous(breaks = 1:12, labels = 1:12)

# save egg/mass as a variable for expected reproductive output model
write.csv(egg_per_mass_estimates,  file = "/Users/aoyu/Desktop/Snail_Data/FinalOutputs/R2_egg_per_mass_estimates.csv")

# add prevous week's snail size as covariate to control for
source("/Users/aoyu/Desktop/Snail_Data/R2_growth_analysis.R") # cleaned master data sheet called all_data_growth
covariate <- all_data_growth_length_covariate %>% select(Jar_num, Color, Week_num, prev_snail_length_avg)
all_data_hatch <- merge(all_data_hatch, covariate, by = c("Jar_num", "Color", "Week_num"), all.x = TRUE)
all_data_hatch <- all_data_hatch %>% filter(!is.na(prev_snail_length_avg)) # 976 missing prev_snail_length_avg

# binomial model with zero inflation
hatch_model_withWeek_num <- glmmTMB( cbind(hatch, unhatch) ~ Treatment + Week_num +  snails_num ,
                                     ziformula = ~Treatment + Week_num +  snails_num,
                                     family = binomial, data = all_data_hatch)
hatch_model_withWeek_num1 <- glmmTMB( cbind(hatch, unhatch) ~ Treatment + prev_snail_length_avg +  snails_num ,
                                     ziformula = ~Treatment + prev_snail_length_avg +  snails_num,
                                     family = binomial, data = all_data_hatch)

AIC(hatch_model_withWeek_num, hatch_model_withWeek_num1) # Week_num is better covariate than prev_snail_length_avg

hatch_model <- glmmTMB( cbind(hatch, unhatch) ~ Treatment +  snails_num ,
                        ziformula = ~Treatment +  snails_num,
                        family = binomial, data = all_data_hatch)

AIC(hatch_model_withWeek_num, hatch_model) # smaller the AIC the better, so the model with Week_num is better

# Summary of the model and calculation of odds for interpretation
summary(hatch_model_withWeek_num)

# R2_hatch_model_4treatments <- data.frame(
#   term = c("Seasoned Macro - Low", "Seasoned Macro - High", "Virgin Macro - Low", "Virgin Macro - High", "Experimental Week Number", "Number of Snails per Jar"),
#   estimate = c(1.01882, 1.37112, 0.06411 , 0.25767, -0.02329, -0.01715),
#   standard.error = c(0.05104, 0.05684, 0.05325, 0.05604, 0.00736, 0.03479),
#   exp.estimate = c(exp(1.01882), exp(1.37112), exp(0.06411), exp(0.25767), exp(-0.02329), exp(-0.01715)),
#   p.value = c("< 2e-16", "< 2e-16", 0.22865, 4.27e-06,0.00156, 0.62203))
# 
# write.csv(R2_hatch_model_4treatments,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_hatch_model_4treatments.csv")

# # Post-hoc comparisons (emmeans)
emmeans_hatch <- emmeans(hatch_model_withWeek_num, pairwise ~ Treatment)
summary(emmeans_hatch)

R2_hatch_pairwise_4treatments <- data.frame(
  term = c("Control - (Biofouled Macro - Low)",
           "Control - (Biofouled Macro - High)",
           "Control - (Virgin Macro - Low)",
           "Control - (Virgin Macro - High)",
           "(Biofouled Macro - Low) - (Biofouled Macro - High)", 
           "(Biofouled Macro - Low) - (Virgin Macro - Low)",
           "(Biofouled Macro - Low) - (Virgin Macro - High)",
           "(Biofouled Macro - High) - (Virgin Macro - Low) ",
           "(Biofouled Macro - High) - (Virgin Macro - High)",
           "(Virgin Macro - Low) - (Virgin Macro - High)"),
  estimate = c(-1.0278, -1.3670, -0.0669 , -0.2523, -0.3392, 0.9609,0.7755,1.3001,1.1147, -0.1854),
  standard.error = c(0.0513, 0.0571, 0.0537, 0.0562, 0.0605, 0.0582, 0.0602, 0.0633, 0.0651, 0.0624),
  exp.estimate = c(exp(-1.0278), exp(-1.3670), exp(-0.0669), exp(-0.2523), exp(-0.3392), exp(0.9609), exp(0.7755), exp(1.3001), exp(1.1147), exp(-0.1854)),
  p.value = c("< 0.0001","< 0.0001", 0.7245, 0.0001, "< 0.0001", "< 0.0001", "< 0.0001", "< 0.0001","< 0.0001", 0.0249))

write.csv(R2_hatch_pairwise_4treatments,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_hatch_pairwise_4treatments.csv")
######## Plot for each treatment with varying abundance #########
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
R2_hatch_plot_4treatments <- ggplot(fitness_summary, aes(x = Week_num, y = mean_fitness, color = Treatment, group = Treatment)) +
  geom_line() +  # add line connecting dots for each Color group
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  # geom_smooth(method = "loess", se = FALSE) +
  # geom_smooth(method = "loess", se = TRUE, aes(fill = Treatment), alpha = 0.2) +  # with se
  scale_color_manual(values = c("Control" = "#6baf78", "Seasoned Macro - Low" = "#d4c2a8", "Seasoned Macro - High" = "#a5855f", "Virgin Macro - Low" = "#bcbcbc", "Virgin Macro - High" = "#8a8a8a"),
                     labels = c("Control", 
                                "Biofouled Macro - Low", 
                                "Biofouled Macro - High",
                                "Virgin Macro - Low",
                                "Virgin Macro - High")) +
  scale_fill_manual(values = c("Control" = "#6baf78", "Seasoned Macro - Low" = "#d4c2a8", "Seasoned Macro - High" = "#a5855f", "Virgin Macro - Low" = "#bcbcbc", "Virgin Macro - High" = "#8a8a8a")) +
  scale_x_continuous(breaks = 1:13, labels = 1:13)+
  labs(
    x = "Weeks",
    y = "Mean Hatching Success Rate per Jar",
    color = "Treatment") +
  theme_minimal()

# Save the plot to a PNG file
ggsave(filename = "R2_hatch_plot_4treatments.png", plot = R2_hatch_plot_4treatments, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/FinalOutputs", width = 9, height = 6, dpi = 300)

####### Offspring hatching analysis for grouped treatment disregarding abundance #########
# Biomial: number of success out of total number of trials
hist(all_data_hatch$hatch_success)

# Calculate mean and variance of the dependent variable
mean_value <- mean(all_data_hatch$hatch_success)
variance_value <- var(all_data_hatch$hatch_success)

# Print the results
cat("Mean :", mean_value, "\n")
cat("Variance :", variance_value, "\n") # Variance < Mean

# binomial model with zero inflation
hatch_model_withWeek_num <- glmmTMB( cbind(hatch, unhatch) ~ Group + Week_num +  snails_num ,
                                     ziformula = ~ Group + Week_num +  snails_num,
                                     family = binomial, data = all_data_hatch)
hatch_model_withWeek_num1 <- glmmTMB( cbind(hatch, unhatch) ~ Group + prev_snail_length_avg +  snails_num ,
                                     ziformula = ~ Group + prev_snail_length_avg +  snails_num,
                                     family = binomial, data = all_data_hatch)
AIC(hatch_model_withWeek_num, hatch_model_withWeek_num1) # Week_num is better covariate than prev_snail_length_avg
hatch_model <- glmmTMB( cbind(hatch, unhatch) ~ Group +  snails_num ,
                        ziformula = ~ Group +  snails_num,
                        family = binomial, data = all_data_hatch)

anova(hatch_model_withWeek_num, hatch_model)# smaller the AIC the better, so the model with parent age is better

# Summary of the model
summary(hatch_model_withWeek_num)

# R2_hatch_model_grouped <- data.frame(
#   term = c("Seasoned Macro", "Virgin Macro", "Experimental Week Number", "Number of Snails per Jar"),
#   estimate = c(1.1741034, 0.1545806, -0.0206389 , -0.0002519),
#   standard.error = c(0.0445508, 0.0449325, 0.0073292, 0.0345753),
#   exp.estimate = c(exp(1.1741034), exp(0.1545806), exp(-0.0206389), exp(-0.0002519)),
#   p.value = c("< 2e-16", 0.000581, 0.004863, 0.994186))
# 
# write.csv(R2_hatch_model_grouped,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_hatch_model_grouped.csv")

# # Post-hoc comparisons (emmeans)
emmeans_hatch <- emmeans(hatch_model_withWeek_num1, pairwise ~ Group)
summary(emmeans_hatch)

R2_hatch_pairwise_grouped <- data.frame(
  term = c("Control - Biofouled Macro",
           "Control - Virgin Macro",
           "Biofouled Macro - Virgin Macro"),
  estimate = c(-1.15, -0.15, 1.00),
  standard.error = c(0.0447, 0.0452, 0.0431),
  exp.estimate = c(exp(-1.15), exp(-0.15), exp(1.00)),
  p.value = c("< 0.0001", 0.0026, "< 0.0001"))

write.csv(R2_hatch_pairwise_grouped,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_hatch_pairwise_grouped.csv")
######## Plot for grouped treatment disregarding abundance #########
# Summarize the data by Color and Week_num to calculate mean and 95% CI
fitness_summary <- all_data_hatch %>%
  group_by(Group, Week_num) %>%
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
R2_hatch_plot_grouped <- ggplot(fitness_summary, aes(x = Week_num, y = mean_fitness, color = Group, group = Group)) +
  geom_line() +  # add line connecting dots for each Color group
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  # geom_smooth(method = "loess", se = FALSE) +
  # geom_smooth(method = "loess", se = TRUE, aes(fill = Treatment), alpha = 0.2) +  # with se
  scale_color_manual(values = c("Control" = "#6baf78", "Seasoned Macro" = "#a5855f", "Virgin Macro" = "#8a8a8a"),
                     labels = c("Control", 
                                "Biofouled Macro", 
                                "Virgin Macro")) +
  scale_fill_manual(values = c("Control" = "#6baf78", "Seasoned Macro" = "#a5855f", "Virgin Macro" = "#8a8a8a")) +
  scale_x_continuous(breaks = 1:13, labels = 1:13)+
  labs(
    x = "Weeks",
    y = "Mean Hatching Success Rate",
    color = "Treatment") +
  theme_minimal()

# Save the plot to a PNG file
ggsave(filename = "R2_hatch_plot_grouped.png", plot = R2_hatch_plot_grouped, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)
