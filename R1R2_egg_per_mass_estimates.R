#egg per mass estimates for R1&R2 by treatments

egg_r1 <- read.csv("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/R1_egg_per_mass_estimates.csv")
egg_r2 <- read.csv("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/R2_egg_per_mass_estimates.csv")

egg_per_mass_r1 <- egg_r1 %>% group_by(Treatment) %>% summarise(egg_per_mass=mean(egg_per_mass)) # around 11
egg_per_mass_r2 <- egg_r2 %>% group_by(Treatment) %>% summarise(egg_per_mass=mean(egg_per_mass)) # around 11-12

# Load required libraries
library(dplyr)
library(emmeans)

# 1. Run the two-way ANOVA model with an interaction effect
anova_model_r1 <- aov(egg_per_mass ~ Treatment , data = egg_r1)
anova_model_r2 <- aov(egg_per_mass ~ Treatment , data = egg_r2)

# Check the overall ANOVA results
summary(anova_model_r1)
summary(anova_model_r2)

# 2. Extract pairwise comparisons between Treatments WITHIN each round, identity link used for Gaussian
treatment_comparisons_r1 <- emmeans(anova_model_r1, pairwise ~ Treatment)
treatment_comparisons_r2 <- emmeans(anova_model_r2, pairwise ~ Treatment)

# 3. View the specific contrast results (with adjusted p-values)
egg_per_mass_contrasts_r1 <-treatment_comparisons_r1$contrasts
egg_per_mass_contrasts_r2 <-treatment_comparisons_r2$contrasts

write.csv(egg_per_mass_contrasts_r1,  file = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/R1_egg_per_mass_contrast.csv")
write.csv(egg_per_mass_contrasts_r2,  file = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/R2_egg_per_mass_contrast.csv")
