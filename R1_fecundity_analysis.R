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
source("/Users/aoyu/Desktop/Snail_Data/R1_data_cleaning.R") # cleaned master data sheet called all_data
all_data_repro <- all_data[all_data$Week_num != 1, ]  # Remove week 1 data
all_data_repro <- all_data_repro %>% select(Jar_num, Color, Week_num, Total_egg_mass_num, Adult_death_1, Adult_death_2, Adult_death_3, Treatment)
all_data_repro$total_egg_mass_num <- as.numeric(all_data_repro$Total_egg_mass_num)

####### Reproductive output analysis at the jar level #########
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
# little magnitude of positive correlation between snail number and egg mass number, negligible correlation between snail number and egg mass number

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

# Plot the data with ggplot2 using specified colors, connecting points with a smooth curve for each treatment
ggplot(repro_summary, aes(x = Week_num, y = mean_repro, color = Treatment, group = Treatment)) +
  geom_line() +  # add line connecting dots for each Color group
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  # geom_smooth(method = "loess", se = FALSE) +
  # geom_smooth(method = "loess", se = TRUE, aes(fill = Treatment), alpha = 0.2) +  # with se
  scale_color_manual(values = c("Control" = "green", "Micro" = "yellow", "Macro" = "orange", "Macro+Micro" = "darkgray")) +
  scale_fill_manual(values = c("Control" = "green", "Micro" = "yellow", "Macro" = "orange", "Macro+Micro" = "darkgray")) +
  labs(
    x = "Week Number",
    y = "Mean Total Egg Mass Number",
    title = "Fecundity Over Time by Treatment"
  ) +
  theme_minimal()

# plot and summary statistics across time masks true variation and could be missleading:
# install.packages("geomViolinDiscrete") # unable to install; " package ‘geomViolinDiscrete’ is not available for this version of R"
# # library("geomViolinDiscrete")
# library(ggsignif)

# # Plot the raw data using a violin plot
# ggplot(all_data_repro, aes(x = factor(Color, levels = c("G", "Y", "O", "W")), y = total_egg_mass_num, fill = Color)) +
#   geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75), color = "black") + 
#   # geom_jitter(width = 0.1, alpha = 0.4, color = "darkblue") +  # Add individual data points
#   theme_minimal() +  # Clean theme
#   labs(title = "Egg Mass Counts by Treatment Group",
#        x = "Treatment Group",
#        y = "Total Egg Mass Count",
#        fill = "Treatment Group") +
#   scale_x_discrete(labels = c("O" = "Macroplastic", "Y" = "Microplastic", "W" = "Microplastic & Macroplastic", "G" = "Negative Control")) +
#   scale_fill_manual(values = c("O" = "orange", "Y" = "yellow", "W" = "darkgray", "G" = "green")) + 
#   stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black")  + # Add mean points
#   geom_signif(comparisons = list(c("G", "Y"), c("G", "O"), c("G", "W") ), # Add manual pair-wise comparisons
#               annotations = c("***", "**", "*"),  # p-value stars here
#               y_position = c(29, 34, 39),         # adjust y positions 
#               tip_length = 0,
#               bracket.size = 0)
# 
# #box plot:
# ggplot(all_data_repro, aes(x = factor(Color, levels = c("G", "Y", "O", "W")), y = total_egg_mass_num, fill = Color)) +
#   geom_boxplot(outlier.colour = "black", outlier.shape = 16, outlier.size = 2, 
#                notch = FALSE, width = 0.7) +  # Box plot with specified aesthetics
#   # geom_jitter(width = 0.1, alpha = 0.4, color = "darkblue") +  # Uncomment to add individual data points
#   theme_minimal() +  # Clean theme
#   labs(title = "Egg Mass Counts by Treatment Group",
#        x = "Treatment Group",
#        y = "Total Egg Mass Count",
#        fill = "Treatment Group") +
#   scale_x_discrete(labels = c("O" = "Macroplastic", "Y" = "Microplastic", "W" = "Microplastic & Macroplastic", "G" = "Negative Control")) +
#   scale_fill_manual(values = c("O" = "orange", "Y" = "yellow", "W" = "darkgray", "G" = "green"))+
#   geom_signif(comparisons = list(c("G", "Y"), c("G", "O"), c("G", "W")), # Add manual pair-wise comparisons
#               annotations = c("***", "**", "*"),  # p-value stars here
#               y_position = c(28, 30, 32),         # adjust y positions 
#               tip_length = 0,
#               bracket.size = 0)
