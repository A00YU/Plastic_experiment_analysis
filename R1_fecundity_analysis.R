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
all_data_repro <- all_data[all_data$Week_num != 1, ]  # Remove week 1 baseline data given all 0
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

library(performance)
r2_nakagawa(reproduction_model_snail_num)
# Conditional R2: 0.090 fixed + random effects
# Marginal R2: 0.073 fixed effects
VarCorr(reproduction_model_snail_num)
icc(reproduction_model_snail_num) # ICC < 0.05 → random effect is negligible

summary(reproduction_model_snail_num)

pairwise_comparisons_repro <- emmeans(reproduction_model_snail_num, pairwise ~ Treatment)
summary(pairwise_comparisons_repro)

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

### Plot across weeks ----
R1_fecundity_plot <- ggplot(repro_summary, aes(x = Week_num, y = mean_repro, color = Treatment, group = Treatment)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), linewidth = 0.5, width = 0.15) +
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
  theme_minimal() +
  theme(legend.position = "top",
                          legend.direction = "horizontal",
        strip.text.x = element_blank(),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(color = "black", size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(color = "black", size = 11))

R1_fecundity_plot
# Save the plot to a PNG file
# ggsave(filename = "R1_fecundity_plot_toplegend.png", plot = R1_fecundity_plot, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)

### marginalized effect plot r1 df ----
# marginalized mean egg mass per jar over a range of covariates values 
emm_repro <- emmeans(reproduction_model_snail_num, ~ Treatment)

# 2. Extract the 95% Confidence Intervals
df_95 <- as.data.frame(summary(emm_repro, level = 0.95)) %>% rename( "lower95"="asymp.LCL",
                                                                     "upper95"="asymp.UCL")

# 3. Extract the 50% Confidence Intervals
df_50 <- as.data.frame(summary(emm_repro, level = 0.50)) %>% rename( "lower95"="asymp.LCL",
                                                                     "upper95"="asymp.UCL")

# 4. Merge into a single final data frame
plot_df_repro_r1 <- data.frame(
  Treatment  = df_95[[1]],
  mean_egg = df_95[[2]],
  lower95    = df_95[[5]],
  upper95    = df_95[[6]],
  lower50    = df_50[[5]],
  upper50    = df_50[[6]]
)

saveRDS(plot_df_repro_r1, "/Users/aoyu/Desktop/Snail_Data/Outputs/plot_df_repro_r1.rds")


