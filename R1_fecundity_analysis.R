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
library(broom)
library(modelsummary)

####### load data#########
source("/Users/aoyu/Desktop/Snail_Data/R1_data_cleaning.R") # cleaned master data sheet called all_data
all_data_repro <- all_data[all_data$Week_num != 1, ]  # Remove week 1 data
all_data_repro <- all_data_repro %>% select(Jar_num, Color, Week_num, Total_egg_mass_num, Enrolled_egg_num, Adult_death_1, Adult_death_2, Adult_death_3, Treatment)
all_data_repro$total_egg_mass_num <- as.numeric(all_data_repro$Total_egg_mass_num)

####### Reproductive output analysis at the jar level #########
# distribution check
# hist(all_data_repro$total_egg_mass_num, breaks = 10, col = "lightblue", border = "black",
     # main = "Histogram", xlab = "egg mass num per jar", ylab = "counts")

# Calculate Spearman's correlation bewteen snail number and egg mass number
# Add a new column 'snails_num' 
all_data_repro$snails_num <- rowSums(all_data_repro[, c("Adult_death_1", "Adult_death_2", "Adult_death_3")] == 0)
# Eliminate rows where snails_num equals 0
all_data_repro <- all_data_repro[!all_data_repro$snails_num == 0, ]

spearman_cor_egg <- cor.test(all_data_repro$snails_num, all_data_repro$total_egg_mass_num, 
                             method = "pearson")
# little magnitude of positive correlation between snail number and egg mass number, negligible correlation between snail number and egg mass number

# Decide on nbinom or piosson:
# Calculate mean and variance of the dependent variable
mean_value <- mean(all_data_repro$total_egg_mass_num)
variance_value <- var(all_data_repro$total_egg_mass_num)

# Print the results
cat("Mean of total_egg_mass_num:", mean_value, "\n")
cat("Variance of total_egg_mass_num:", variance_value, "\n") # Variance >> Mean, making it better suited for over dispersed data

# add prevous week's snail size as covariate to control for
source("/Users/aoyu/Desktop/Snail_Data/R1_growth_analysis.R") # cleaned master data sheet called all_data_growth
covariate <- all_data_growth_length_covariate %>% select(Jar_num, Color, Week_num, prev_snail_length_avg)
all_data_repro <- merge(all_data_repro, covariate, by = c("Jar_num", "Color", "Week_num"), all.x = TRUE)
all_data_repro <- all_data_repro %>% filter(!is.na(prev_snail_length_avg)) # 124 missing prev_snail_length_avg

# what if use enrolled eggs to proximate eggs per mass to get total eggs 
# all_data_repro <- all_data_repro %>% mutate(Enrolled_egg_num = if_else(is.na(Enrolled_egg_num), 0, Enrolled_egg_num),
#                                             total_egg_mass_num = total_egg_mass_num * Enrolled_egg_num)

# check estimated eggs per mass differences among treatments
anova_egg_per_mass <-aov(Enrolled_egg_num ~ Treatment, data = all_data_repro)
summary(anova_egg_per_mass)
TukeyHSD(anova_egg_per_mass)

# Negative Binomial Model for egg mass counts (glmmTMB package)
reproduction_model <- glmmTMB(total_egg_mass_num ~ Treatment + prev_snail_length_avg + (1| Jar_num), family = nbinom2, data = all_data_repro) # (1|jar) account for jar-specific variability over time
summary(reproduction_model)

reproduction_model_snail_num <- glmmTMB(total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num + (1| Jar_num), family = nbinom2, data = all_data_repro)
summary(reproduction_model_snail_num)

anova(reproduction_model_snail_num, reproduction_model) # the number of snails does have an impact on egg mass production

reproduction_model_snail_num_fix <- glmmTMB(total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num, family = nbinom2, data = all_data_repro)
summary(reproduction_model_snail_num_fix)

anova(reproduction_model_snail_num, reproduction_model_snail_num_fix) # random effect model better

pairwise_comparisons_repro <- emmeans(reproduction_model_snail_num, pairwise ~ Treatment)
# # Summarize pairwise comparisons
summary(pairwise_comparisons_repro)
# mannual data entry instead
# R1_fecundity_pairwise <- data.frame(
#   term = c("Virgin Micro - (Weathered Macro + Virgin Micro)", 
#            "Control - (Weathered Macro + Virgin Micro)", 
#            "Weathered Macro - (Weathered Macro + Virgin Micro)", 
#            "Virgin Micro - Weathered Macro",
#            "Control - Weathered Macro", "Control - Virgin Micro"),
#   estimate = c(-0.1986, -0.1848, -0.1584 , -0.0402, -0.0264, 0.0138),
#   standard.error = c(0.0480, 0.0481, 0.0624, 0.0492, 0.0489, 0.0297),
#   exp.estimate = c(exp(-0.1986), exp(-0.1848), exp(-0.1584), exp(-0.0402), exp(-0.0264), exp(0.0138)),
#   p.value = c(0.0002,0.0007, 0.0543, 0.8467, 0.9492, 0.9671))

# write.csv(R1_fecundity_pairwise,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R1_fecundity_pairwise.csv")

summary(reproduction_model_snail_num)
R1_fecundity_model <- data.frame(
  term = c("Virgin Micro", "Biofouled Macro", "Biofouled Macro + Virgin Micro", "Snail Shell Length Avg per Jar", "Number of Snails per Jar"),
  estimate = c(-0.01562, 0.01476, 0.18878 , 1.38546, 0.17769),
  standard.error = c(0.02974, 0.04894, 0.04790, 0.13936, 0.02464),
  exp.estimate = c(exp(-0.01562), exp(0.01476), exp(0.18878), exp(1.38546), exp(0.17769)),
  p.value = c(0.599,0.763, 8.11e-05, "< 2e-16", 5.48e-13))

write.csv(R1_fecundity_model,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R1_fecundity_model.csv")

####### Visualization of reproductive output #########
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

# Plot across weeks
R1_fecundity_plot <- ggplot(repro_summary, aes(x = Week_num, y = mean_repro, color = Treatment, group = Treatment)) +
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
    y = "Mean Total Number of Egg Mass per Jar") +
  theme_minimal()

# Save the plot to a PNG file
ggsave(filename = "R1_fecundity_plot.png", plot = R1_fecundity_plot, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)
