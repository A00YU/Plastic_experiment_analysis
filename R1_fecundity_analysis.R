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
# library(modelsummary)

####### load data#########
source("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/R1_data_cleaning.R") # cleaned master data sheet called all_data
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
source("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/R1_growth_analysis.R") # cleaned master data sheet called all_data_growth
covariate <- all_data_growth_length_covariate %>% select(Jar_num, Color, Week_num, prev_snail_length_avg)
all_data_repro <- merge(all_data_repro, covariate, by = c("Jar_num", "Color", "Week_num"), all.x = TRUE)
all_data_repro <- all_data_repro %>% filter(!is.na(prev_snail_length_avg)) %>%  mutate(jar_id = interaction(Color, Jar_num, drop = TRUE))# 124 missing prev_snail_length_avg

# Negative Binomial Model for egg mass counts (glmmTMB package)
reproduction_model <- glmmTMB(total_egg_mass_num ~ Treatment + prev_snail_length_avg + (1| jar_id), family = nbinom2, data = all_data_repro) # (1|jar) account for jar-specific variability over time
summary(reproduction_model)

reproduction_model_snail_num <- glmmTMB(total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num + (1| jar_id), family = nbinom2, data = all_data_repro)
summary(reproduction_model_snail_num)

anova(reproduction_model_snail_num, reproduction_model) # the number of snails does have an impact on egg mass production

reproduction_model_snail_num_fix <- glmmTMB(total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num, family = nbinom2, data = all_data_repro)
summary(reproduction_model_snail_num_fix)

anova(reproduction_model_snail_num, reproduction_model_snail_num_fix) # random effect model better

summary(reproduction_model_snail_num)

# Fixed effects (conditional model) summary → CSV
glmm_model_summary <- summary(reproduction_model_snail_num)

coef_table <- as.data.frame(glmm_model_summary$coefficients$cond)
coef_table <- data.frame(term = rownames(coef_table), coef_table, row.names = NULL)

output_dir <- "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2"

# write.csv(
#   coef_table,
#   file.path(output_dir, "R1_fecundity_model.csv"),
#   row.names = FALSE
# )

# Random-effect variances and standard deviations (correct accessor for glmmTMB)
glmmTMB::VarCorr(reproduction_model_snail_num)
sigma(reproduction_model_snail_num)

# Confirm random-effect grouping names
names(glmmTMB::VarCorr(reproduction_model_snail_num)$cond)


pairwise_comparisons_repro <- emmeans(reproduction_model_snail_num, pairwise ~ Treatment, type = "response")
summary(pairwise_comparisons_repro)

fecundity_pairwise <- as.data.frame(
  summary(pairwise_comparisons_repro)$contrasts
)

# write.csv(
#   fecundity_pairwise,
#   file.path(output_dir, "R1_fecundity_pairwise.csv"),
#   row.names = FALSE
# )

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
# ggsave(filename = "R1_fecundity_plot_toplegend.png", plot = R1_fecundity_plot, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)

### marginalized effect plot r1 df ----
# marginalized mean egg mass per jar over a range of covariates values 
emm_repro <- emmeans(reproduction_model_snail_num, ~ Treatment, type = "response")

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

# saveRDS(plot_df_repro_r1, "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/plot_df_repro_r1.rds")




# ###explore: when does the plastic effect starts to be significant####
# library(dplyr)
# library(glmmTMB)
# library(emmeans)
# 
# all_data_repro <- all_data_repro %>%
#   filter(Treatment%in%c("Macro","Control")) %>% 
#   mutate(
#     Treatment = relevel(factor(Treatment), ref = "Control"),  # change "G" if control has another name
#     Jar_num = factor(Jar_num)
#   )
# 
# # week-specific test
# week2_df <- all_data_repro %>% filter(Week_num == 2)
# week3_df <- all_data_repro %>% filter(Week_num == 3)
# week4_df <- all_data_repro %>% filter(Week_num == 4)
# week5_df <- all_data_repro %>% filter(Week_num == 5)
# week6_df <- all_data_repro %>% filter(Week_num == 6)
# week7_df <- all_data_repro %>% filter(Week_num == 7)
# week8_df <- all_data_repro %>% filter(Week_num == 8)
# week9_df <- all_data_repro %>% filter(Week_num == 9)
# week10_df <- all_data_repro %>% filter(Week_num == 10)
# week11_df <- all_data_repro %>% filter(Week_num == 11)
# week12_df <- all_data_repro %>% filter(Week_num == 12)
# week13_df <- all_data_repro %>% filter(Week_num == 13)
# 
# week2_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week2_df
# )
# 
# week2_null <- glmmTMB(
#   total_egg_mass_num ~ prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week2_df
# )
# 
# anova(week2_null, week2_full) # treatment already is a sig factor of egg mass in week 2
# emmeans(week2_full, pairwise ~ Treatment, type = "response") # control/macro not diff.
# 
# week3_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week3_df
# )
# 
# week3_null <- glmmTMB(
#   total_egg_mass_num ~ prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week3_df
# )
# 
# anova(week3_null, week3_full) # treatment already is a sig factor of egg mass in week 2
# emmeans(week3_full, pairwise ~ Treatment, type = "response") # start to diff
# 
# week4_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week4_df
# )
# 
# week4_null <- glmmTMB(
#   total_egg_mass_num ~ prev_snail_length_avg + snails_num,
#   family = nbinom2,
#   data = week4_df
# )
# 
# anova(week4_null, week4_full) # NOT a sig factor of egg mass in week 4
# emmeans(week4_full, pairwise ~ Treatment, type = "response") # not sig
# 
# week5_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week5_df
# )
# 
# week5_null <- glmmTMB(
#   total_egg_mass_num ~ prev_snail_length_avg + snails_num,
#   family = nbinom2,
#   data = week5_df
# )
# 
# anova(week5_null, week5_full) # NOT a sig factor of egg mass in week 5
# emmeans(week5_full, pairwise ~ Treatment, type = "response") # not sig
# 
# week6_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week6_df
# )
# 
# week6_null <- glmmTMB(
#   total_egg_mass_num ~ prev_snail_length_avg + snails_num,
#   family = nbinom2,
#   data = week6_df
# )
# 
# anova(week6_null, week6_full) # NOT a sig factor
# emmeans(week6_full, pairwise ~ Treatment, type = "response") # not sig
# 
# week7_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week7_df
# )
# 
# week7_null <- glmmTMB(
#   total_egg_mass_num ~ prev_snail_length_avg + snails_num,
#   family = nbinom2,
#   data = week7_df
# )
# 
# anova(week7_null, week7_full) #  a sig factor of egg mass in week 5
# emmeans(week7_full, pairwise ~ Treatment, type = "response") # not sig
# 
# week8_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week8_df
# )
# 
# week8_null <- glmmTMB(
#   total_egg_mass_num ~ prev_snail_length_avg + snails_num,
#   family = nbinom2,
#   data = week8_df
# )
# 
# anova(week8_null, week8_full) # NOT a sig factor of egg mass in week 5
# emmeans(week8_full, pairwise ~ Treatment, type = "response") # not sig
# 
# week9_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week9_df
# )
# 
# week9_null <- glmmTMB(
#   total_egg_mass_num ~ prev_snail_length_avg + snails_num,
#   family = nbinom2,
#   data = week9_df
# )
# 
# anova(week9_null, week9_full) # NOT a sig factor 
# emmeans(week9_full, pairwise ~ Treatment, type = "response") # not sig
# 
# week10_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week10_df
# )
# 
# week10_null <- glmmTMB(
#   total_egg_mass_num ~ prev_snail_length_avg + snails_num,
#   family = nbinom2,
#   data = week10_df
# )
# 
# anova(week10_null, week10_full) # NOT a sig factor 
# emmeans(week10_full, pairwise ~ Treatment, type = "response") # not sig
# 
# week11_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week11_df
# )
# 
# week11_null <- glmmTMB(
#   total_egg_mass_num ~ prev_snail_length_avg + snails_num,
#   family = nbinom2,
#   data = week11_df
# )
# 
# anova(week11_null, week11_full) # NOT a sig factor 
# emmeans(week11_full, pairwise ~ Treatment, type = "response") # not sig
# 
# week12_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week12_df
# )
# 
# week12_null <- glmmTMB(
#   total_egg_mass_num ~ prev_snail_length_avg + snails_num,
#   family = nbinom2,
#   data = week12_df
# )
# 
# anova(week12_null, week12_full) # NOT a sig factor of egg mass in week 5
# emmeans(week12_full, pairwise ~ Treatment, type = "response") # not sig
# 
# week13_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num ,
#   family = nbinom2,
#   data = week13_df
# )
# 
# week13_null <- glmmTMB(
#   total_egg_mass_num ~ prev_snail_length_avg + snails_num,
#   family = nbinom2,
#   data = week13_df
# )
# 
# anova(week13_null, week13_full) # NOT a sig factor of egg mass in week 5
# emmeans(week13_full, pairwise ~ Treatment, type = "response") # not sig
# 
# # accumulative week test
# # week2_df is done above
# # accumulative observation-window test with repeated-measures random effect
# # Weeks 2-3 through Weeks 2-13
# 
# week2_3_df <- all_data_repro %>%
#   filter(Week_num %in% 2:3) %>%
#   mutate(
#     Treatment = relevel(factor(Treatment), ref = "Control"),
#     Week_num = factor(Week_num),
#     Jar_num = factor(Jar_num)
#   )
# 
# week2_3_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_3_df
# )
# 
# week2_3_null <- glmmTMB(
#   total_egg_mass_num ~ Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_3_df
# )
# 
# anova(week2_3_null, week2_3_full)
# emmeans(week2_3_full, ~ Treatment, type = "response") %>%
#   contrast(method = "trt.vs.ctrl", ref = "Control")
# 
# 
# week2_4_df <- all_data_repro %>%
#   filter(Week_num %in% 2:4) %>%
#   mutate(
#     Treatment = relevel(factor(Treatment), ref = "Control"),
#     Week_num = factor(Week_num),
#     Jar_num = factor(Jar_num)
#   )
# 
# week2_4_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_4_df
# )
# 
# week2_4_null <- glmmTMB(
#   total_egg_mass_num ~ Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_4_df
# )
# 
# anova(week2_4_null, week2_4_full)
# emmeans(week2_4_full, ~ Treatment, type = "response") %>%
#   contrast(method = "trt.vs.ctrl", ref = "Control")
# 
# 
# week2_5_df <- all_data_repro %>%
#   filter(Week_num %in% 2:5) %>%
#   mutate(
#     Treatment = relevel(factor(Treatment), ref = "Control"),
#     Week_num = factor(Week_num),
#     Jar_num = factor(Jar_num)
#   )
# 
# week2_5_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_5_df
# )
# 
# week2_5_null <- glmmTMB(
#   total_egg_mass_num ~ Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_5_df
# )
# 
# anova(week2_5_null, week2_5_full)
# emmeans(week2_5_full, ~ Treatment, type = "response") %>%
#   contrast(method = "trt.vs.ctrl", ref = "Control")
# 
# 
# week2_6_df <- all_data_repro %>%
#   filter(Week_num %in% 2:6) %>%
#   mutate(
#     Treatment = relevel(factor(Treatment), ref = "Control"),
#     Week_num = factor(Week_num),
#     Jar_num = factor(Jar_num)
#   )
# 
# week2_6_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_6_df
# )
# 
# week2_6_null <- glmmTMB(
#   total_egg_mass_num ~ Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_6_df
# )
# 
# anova(week2_6_null, week2_6_full)
# emmeans(week2_6_full, ~ Treatment, type = "response") %>%
#   contrast(method = "trt.vs.ctrl", ref = "Control")
# 
# 
# week2_7_df <- all_data_repro %>%
#   filter(Week_num %in% 2:7) %>%
#   mutate(
#     Treatment = relevel(factor(Treatment), ref = "Control"),
#     Week_num = factor(Week_num),
#     Jar_num = factor(Jar_num)
#   )
# 
# week2_7_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_7_df
# )
# 
# week2_7_null <- glmmTMB(
#   total_egg_mass_num ~ Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_7_df
# )
# 
# anova(week2_7_null, week2_7_full)
# emmeans(week2_7_full, ~ Treatment, type = "response") %>%
#   contrast(method = "trt.vs.ctrl", ref = "Control")
# 
# 
# week2_8_df <- all_data_repro %>%
#   filter(Week_num %in% 2:8) %>%
#   mutate(
#     Treatment = relevel(factor(Treatment), ref = "Control"),
#     Week_num = factor(Week_num),
#     Jar_num = factor(Jar_num)
#   )
# 
# week2_8_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_8_df
# )
# 
# week2_8_null <- glmmTMB(
#   total_egg_mass_num ~ Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_8_df
# )
# 
# anova(week2_8_null, week2_8_full)
# emmeans(week2_8_full, ~ Treatment, type = "response") %>%
#   contrast(method = "trt.vs.ctrl", ref = "Control")
# 
# 
# week2_9_df <- all_data_repro %>%
#   filter(Week_num %in% 2:9) %>%
#   mutate(
#     Treatment = relevel(factor(Treatment), ref = "Control"),
#     Week_num = factor(Week_num),
#     Jar_num = factor(Jar_num)
#   )
# 
# week2_9_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_9_df
# )
# 
# week2_9_null <- glmmTMB(
#   total_egg_mass_num ~ Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_9_df
# )
# 
# anova(week2_9_null, week2_9_full)
# emmeans(week2_9_full, ~ Treatment, type = "response") %>%
#   contrast(method = "trt.vs.ctrl", ref = "Control")
# 
# 
# week2_10_df <- all_data_repro %>%
#   filter(Week_num %in% 2:10) %>%
#   mutate(
#     Treatment = relevel(factor(Treatment), ref = "Control"),
#     Week_num = factor(Week_num),
#     Jar_num = factor(Jar_num)
#   )
# 
# week2_10_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_10_df
# )
# 
# week2_10_null <- glmmTMB(
#   total_egg_mass_num ~ Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_10_df
# )
# 
# anova(week2_10_null, week2_10_full)
# emmeans(week2_10_full, ~ Treatment, type = "response") %>%
#   contrast(method = "trt.vs.ctrl", ref = "Control")
# 
# 
# week2_11_df <- all_data_repro %>%
#   filter(Week_num %in% 2:11) %>%
#   mutate(
#     Treatment = relevel(factor(Treatment), ref = "Control"),
#     Week_num = factor(Week_num),
#     Jar_num = factor(Jar_num)
#   )
# 
# week2_11_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_11_df
# )
# 
# week2_11_null <- glmmTMB(
#   total_egg_mass_num ~ Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_11_df
# )
# 
# anova(week2_11_null, week2_11_full)
# emmeans(week2_11_full, ~ Treatment, type = "response") %>%
#   contrast(method = "trt.vs.ctrl", ref = "Control")
# 
# 
# week2_12_df <- all_data_repro %>%
#   filter(Week_num %in% 2:12) %>%
#   mutate(
#     Treatment = relevel(factor(Treatment), ref = "Control"),
#     Week_num = factor(Week_num),
#     Jar_num = factor(Jar_num)
#   )
# 
# week2_12_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_12_df
# )
# 
# week2_12_null <- glmmTMB(
#   total_egg_mass_num ~ Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_12_df
# )
# 
# anova(week2_12_null, week2_12_full)
# emmeans(week2_12_full, ~ Treatment, type = "response") %>%
#   contrast(method = "trt.vs.ctrl", ref = "Control")
# 
# 
# week2_13_df <- all_data_repro %>%
#   filter(Week_num %in% 2:13) %>%
#   mutate(
#     Treatment = relevel(factor(Treatment), ref = "Control"),
#     Week_num = factor(Week_num),
#     Jar_num = factor(Jar_num)
#   )
# 
# week2_13_full <- glmmTMB(
#   total_egg_mass_num ~ Treatment + Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_13_df
# )
# 
# week2_13_null <- glmmTMB(
#   total_egg_mass_num ~ Week_num + prev_snail_length_avg + snails_num + (1 | Jar_num),
#   family = nbinom2,
#   data = week2_13_df
# )
# 
# anova(week2_13_null, week2_13_full)
# emmeans(week2_13_full, ~ Treatment, type = "response") %>%
#   contrast(method = "trt.vs.ctrl", ref = "Control")

