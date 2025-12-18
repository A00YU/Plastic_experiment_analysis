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
library(car)

####### load data#########
source("/Users/aoyu/Desktop/Snail_Data/R2_data_cleaning.R") # cleaned master data sheet called all_data_r2
all_data_repro <- all_data_r2[all_data_r2$Week_num != 1, ]  # Remove week 1 data
all_data_repro <- all_data_repro %>% select(Jar_num, Color, Week_num, Total_egg_mass_num, Adult_death_1, Adult_death_2, Adult_death_3, Treatment)
all_data_repro$total_egg_mass_num <- as.numeric(all_data_repro$Total_egg_mass_num)
# Combine Y&O into SM and W-1&W-3 into VM
all_data_repro$Group <- factor(with(all_data_repro, 
                                    ifelse(Color %in% c("Y", "O"), "Seasoned Macro",
                                           ifelse(Color %in% c("W-1", "W-3"), "Virgin Macro","Control"))))

# Verify the new grouping
table(all_data_repro$Group, all_data_repro$Week_num)

# distribution check
# hist(all_data_repro$total_egg_mass_num, breaks = 10, col = "lightblue", border = "black",
#      main = "Histogram", xlab = "egg mass num per jar", ylab = "counts")

# Calculate Spearman's correlation bewteen snail number and egg mass number
# Add a new column 'snails_num' 
all_data_repro$snails_num <- rowSums(all_data_repro[, c("Adult_death_1", "Adult_death_2", "Adult_death_3")] == 0)
# Eliminate rows where snails_num equals 0
all_data_repro <- all_data_repro[!all_data_repro$snails_num == 0, ]

spearman_cor_egg <- cor.test(all_data_repro$snails_num, all_data_repro$total_egg_mass_num, 
                             method = "pearson")
# moderate positive correlation between snail number and egg mass number, negligible correlation between snail number and egg mass number

####### Reproductive output analysis by grouped treatment disregard varying plastic abundance #########
# Decide on nbinom or piosson:
# Calculate mean and variance of the dependent variable
mean_value <- mean(all_data_repro$total_egg_mass_num)
variance_value <- var(all_data_repro$total_egg_mass_num)

# Print the results
cat("Mean of total_egg_mass_num:", mean_value, "\n")
cat("Variance of total_egg_mass_num:", variance_value, "\n") # Variance >> Mean, making it better suited for over dispersed data

# add prevous week's snail size as covariate to control for
source("/Users/aoyu/Desktop/Snail_Data/R2_growth_analysis.R") # cleaned master data sheet called all_data_growth
covariate <- all_data_growth_length_covariate %>% select(Jar_num, Color, Week_num, prev_snail_length_avg)
all_data_repro <- merge(all_data_repro, covariate, by = c("Jar_num", "Color", "Week_num"), all.x = TRUE)
all_data_repro <- all_data_repro %>% filter(!is.na(prev_snail_length_avg)) # 964 missing prev_snail_length_avg

# Negative Binomial Model for egg mass counts (glmmTMB package)
reproduction_model_1 <- glmmTMB(total_egg_mass_num ~ Treatment + Week_num , family = nbinom2, data = all_data_repro) # (1|jar) account for jar-specific variability over time
summary(reproduction_model_1)

reproduction_model <- glmmTMB(total_egg_mass_num ~ Treatment + prev_snail_length_avg , family = nbinom2, data = all_data_repro) # (1|jar) account for jar-specific variability over time
summary(reproduction_model)

AIC(reproduction_model, reproduction_model_1) # reproduction_model is better

reproduction_model_snail_num_1 <- glmmTMB(total_egg_mass_num ~ Treatment + Week_num + snails_num , family = nbinom2, data = all_data_repro)
summary(reproduction_model_snail_num_1)

reproduction_model_snail_num <- glmmTMB(total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num , family = nbinom2, data = all_data_repro)
summary(reproduction_model_snail_num)

AIC(reproduction_model_snail_num, reproduction_model_snail_num_1) # prev_snail_length_avg as covariate is better
anova(reproduction_model, reproduction_model_snail_num) # snail_num does have a big impact on egg mass production

pairwise_comparisons_repro <- emmeans(reproduction_model_snail_num, pairwise ~ Treatment)
# # Summarize pairwise comparisons
summary(pairwise_comparisons_repro)

R2_fecundity_pairwise_4treatments <- data.frame(
  term = c("Control - (Biofouled Macro - Low)", 
           "Control - (Biofouled Macro - High)", 
           "Control - (Virgin Macro - Low)", 
           "Control - (Virgin Macro - High)",
           "(Biofouled Macro - Low) - (Biofouled Macro - High)",
           "(Biofouled Macro - Low) - (Virgin Macro - Low)",
           "(Biofouled Macro - Low) - (Virgin Macro - High)",
           "(Biofouled Macro - High) - (Virgin Macro - Low)",
           "(Biofouled Macro - High) - (Virgin Macro - High)",
           "(Virgin Macro - Low) - (Virgin Macro - High)"),
  estimate = c(-0.14148, -0.17988, -0.06859 , -0.13214, -0.03840, 0.07289, 0.00934, 0.11129, 0.04775, -0.06355),
  standard.error = c(0.0333, 0.0335, 0.0392, 0.0393, 0.0327, 0.0387, 0.0386, 0.0388, 0.0386, 0.0439),
  exp.estimate = c(exp(-0.14148), exp(-0.17988), exp(-0.06859), exp(-0.13214), exp(-0.03840), exp(0.07289), exp(0.00934), exp(0.11129), exp(0.04775), exp(-0.06355)),
  p.value = c(0.0002,"<0.0001", 0.4025, 0.0068, 0.7656, 0.3249, 0.9992, 0.0335, 0.7306, 0.5958))

write.csv(R2_fecundity_pairwise_4treatments,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_fecundity_pairwise_4treatments.csv")

summary(reproduction_model_snail_num)

R2_fecundity_model_4treatments <- data.frame(
  term = c("Biofouled Macro - Low", "Biofouled Macro - High", "Virgin Macro - Low", "Virgin Macro - High", "Snail Shell Length Avg per Jar", "Number of Snails per Jar"),
  estimate = c(0.14148, 0.17988, 0.06859 , 0.13214, 0.89513, 0.35056),
  standard.error = c(0.03330, 0.03351, 0.03917, 0.03925, 0.18065, 0.02098),
  exp.estimate = c(exp(0.14148), exp(0.17988), exp(0.06859), exp(0.13214), exp(0.89513), exp(0.35056)),
  p.value = c(2.15e-05,7.98e-08, 0.079915, 0.000761, 7.23e-07, "< 2E-16"))

write.csv(R2_fecundity_model_4treatments,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_fecundity_model_4treatments.csv")
####### Visualization of reproductive output by individual treatment with varying plastic abundance #########
# Summarize the data by Color and Week_num to calculate mean and 95% CI
repro_summary <- all_data_repro %>%
  group_by(Treatment, Week_num) %>%
  summarise(
    mean_repro = mean(total_egg_mass_num, na.rm = TRUE),
    sd_repro   = sd(total_egg_mass_num, na.rm = TRUE),
    n          = n(),
    se         = sd_repro / sqrt(n),
    ci_lower   = mean_repro - 1.96 * se,
    ci_upper   = mean_repro + 1.96 * se
  ) %>%
  ungroup()

# Plot the data with ggplot2 using specified colors, connecting points with a smooth curve for each treatment
R2_fecundity_plot_4treatments <- ggplot(repro_summary, aes(x = Week_num, y = mean_repro, color = Treatment, group = Treatment)) +
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
    y = "Mean Total Number of Egg Mass per Jar") +
  theme_minimal()

# Save the plot to a PNG file
ggsave(filename = "R2_fecundity_plot_4treatments.png", plot = R2_fecundity_plot_4treatments, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)

####### Reproductive output analysis by grouped treatment disregard varying plastic abundance #########
# Negative Binomial Model for egg mass counts (glmmTMB package)
reproduction_model_group <- glmmTMB(total_egg_mass_num ~ Group + prev_snail_length_avg , family = nbinom2, data = all_data_repro) # (1|jar) account for jar-specific variability over time
summary(reproduction_model_group)

reproduction_model_snail_num_group <- glmmTMB(total_egg_mass_num ~ Group + prev_snail_length_avg + snails_num , family = nbinom2, data = all_data_repro)
summary(reproduction_model_snail_num_group)

reproduction_model_snail_num_group_rm <- glmmTMB(total_egg_mass_num ~ Group + prev_snail_length_avg + snails_num + (1|Jar_num), family = nbinom2, data = all_data_repro)
anova(reproduction_model_snail_num_group, reproduction_model_group) # the number of snails doesn't have an impact on egg mass production, even if the correlation coefficient is low.
anova(reproduction_model_snail_num_group, reproduction_model_snail_num_group_rm) 

pairwise_comparisons_repro_group <- emmeans(reproduction_model_snail_num_group, pairwise ~ Group)
# # Summarize pairwise comparisons
summary(pairwise_comparisons_repro_group)

R2_fecundity_pairwise_grouped <- data.frame(
  term = c("Control - Biofouled Macro", 
           "Control - Virgin Macro", 
           "Biofouled Macro - Virgin Macro"),
  estimate = c(-0.16, -0.10, 0.06),
  standard.error = c(0.0292, 0.0325, 0.0274),
  exp.estimate = c(exp(-0.16), exp(-0.10), exp(0.06)),
  p.value = c("<0.0001", 0.0058, 0.0724))

write.csv(R2_fecundity_pairwise_grouped,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_fecundity_pairwise_grouped.csv")

summary(reproduction_model_snail_num_group)

R2_fecundity_model_grouped <- data.frame(
  term = c("Biofouled Macro", "Virgin Macro", "Snail Shell Length Avg per Jar", "Number of Snails per Jar"),
  estimate = c(0.16027, 0.10024, 0.91637, 0.35342),
  standard.error = c(0.02916, 0.03253, 0.18042, 0.02093),
  exp.estimate = c(exp(0.16027), exp(0.10024), exp(0.91637), exp(0.35342)),
  p.value = c(3.89e-08,0.00206, 3.79e-07, "< 2e-16"))

write.csv(R2_fecundity_model_grouped,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_fecundity_model_grouped.csv")
####### Visualization of reproductive output by individual treatment with varying plastic abundance #########
# Summarize the data by Color and Week_num to calculate mean and 95% CI
repro_summary_group <- all_data_repro %>%
  group_by(Group, Week_num) %>%
  summarise(
    mean_repro = mean(total_egg_mass_num, na.rm = TRUE),
    sd_repro   = sd(total_egg_mass_num, na.rm = TRUE),
    n          = n(),
    se         = sd_repro / sqrt(n),
    ci_lower   = mean_repro - 1.96 * se,
    ci_upper   = mean_repro + 1.96 * se
  ) %>%
  ungroup()

# Plot the data with ggplot2 using specified colors, connecting points with a smooth curve for each treatment
R2_fecundity_plot_grouped <- ggplot(repro_summary_group, aes(x = Week_num, y = mean_repro, color = Group, group = Group)) +
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
    y = "Mean Total Number of Egg Mass per Jar",
    color = "Treatment"
  ) +
  theme_minimal()

# Save the plot to a PNG file
ggsave(filename = "R2_fecundity_plot_grouped.png", plot = R2_fecundity_plot_grouped, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)
