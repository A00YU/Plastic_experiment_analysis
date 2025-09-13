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
hist(all_data_repro$total_egg_mass_num, breaks = 10, col = "lightblue", border = "black",
     main = "Histogram", xlab = "egg mass num per jar", ylab = "counts")

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

# Negative Binomial Model for egg mass counts (glmmTMB package)
reproduction_model <- glmmTMB(total_egg_mass_num ~ Treatment + Week_num , family = nbinom2, data = all_data_repro) # (1|jar) account for jar-specific variability over time
summary(reproduction_model)

reproduction_model_snail_num <- glmmTMB(total_egg_mass_num ~ Treatment + Week_num + snails_num , family = nbinom2, data = all_data_repro)
summary(reproduction_model_snail_num)

anova(reproduction_model_snail_num, reproduction_model) # the number of snails does have an impact on egg mass production, even if the correlation coefficient is low.

pairwise_comparisons_repro <- emmeans(reproduction_model_snail_num, pairwise ~ Treatment)
# # Summarize pairwise comparisons
summary(pairwise_comparisons_repro)
summary(reproduction_model_snail_num)

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
ggplot(repro_summary, aes(x = Week_num, y = mean_repro, color = Treatment, group = Treatment)) +
  geom_line() +  # add line connecting dots for each Color group
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  # geom_smooth(method = "loess", se = FALSE) +
  # geom_smooth(method = "loess", se = TRUE, aes(fill = Treatment), alpha = 0.2) +  # with se
  scale_color_manual(values = c("Control" = "green", "Seasoned Macro - Low" = "yellow", "Seasoned Macro - High" = "orange", "Virgin Macro - Low" = "lightgray", "Virgin Macro - High" = "darkgray")) +
  scale_fill_manual(values = c("Control" = "green", "Seasoned Macro - Low" = "yellow", "Seasoned Macro - High" = "orange", "Virgin Macro - Low" = "lightgray", "Virgin Macro - High" = "darkgray")) +
  labs(
    x = "Week Number",
    y = "Mean Total Egg Mass Number",
    title = "Fecundity Over Time by Treatment"
  ) +
  theme_minimal()

####### Reproductive output analysis by grouped treatment disregard varying plastic abundance #########
# Negative Binomial Model for egg mass counts (glmmTMB package)
reproduction_model_group <- glmmTMB(total_egg_mass_num ~ Group + Week_num , family = nbinom2, data = all_data_repro) # (1|jar) account for jar-specific variability over time
summary(reproduction_model_group)

reproduction_model_snail_num_group <- glmmTMB(total_egg_mass_num ~ Group + Week_num + snails_num , family = nbinom2, data = all_data_repro)
summary(reproduction_model_snail_num_group)

anova(reproduction_model_snail_num_group, reproduction_model_group) # the number of snails does have an impact on egg mass production, even if the correlation coefficient is low.

pairwise_comparisons_repro_group <- emmeans(reproduction_model_snail_num_group, pairwise ~ Group)
# # Summarize pairwise comparisons
summary(pairwise_comparisons_repro_group)
summary(reproduction_model_snail_num_group)

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
ggplot(repro_summary_group, aes(x = Week_num, y = mean_repro, color = Group, group = Group)) +
  geom_line() +  # add line connecting dots for each Color group
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  # geom_smooth(method = "loess", se = FALSE) +
  # geom_smooth(method = "loess", se = TRUE, aes(fill = Treatment), alpha = 0.2) +  # with se
  scale_color_manual(values = c("Control" = "green", "Seasoned Macro" = "orange", "Virgin Macro" = "darkgray")) +
  scale_fill_manual(values = c("Control" = "green", "Seasoned Macro" = "orange", "Virgin Macro" = "darkgray")) +
  labs(
    x = "Week Number",
    y = "Mean Total Egg Mass Number",
    title = "Fecundity Over Time by Treatment"
  ) +
  theme_minimal()
