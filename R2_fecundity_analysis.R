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

### load data####
source("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/R2_data_cleaning.R") # cleaned master data sheet called all_data_r2
all_data_repro <- all_data_r2[all_data_r2$Week_num != 1, ]  # Remove week 1 data
all_data_repro <- all_data_repro %>% select(Jar_num, Color, Week_num, Total_egg_mass_num, Adult_death_1, Adult_death_2, Adult_death_3, Treatment)
all_data_repro$total_egg_mass_num <- as.numeric(all_data_repro$Total_egg_mass_num)

plot_df_repro_r1 <- readRDS("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/plot_df_repro_r1.rds") # load plot_df_repro_r1 for joint plotting
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

### 4treatment: Reproductive output analysis by grouped treatment with varying plastic abundance ----
# Decide on nbinom or piosson:
# Calculate mean and variance of the dependent variable
mean_value <- mean(all_data_repro$total_egg_mass_num)
variance_value <- var(all_data_repro$total_egg_mass_num)

# Print the results
cat("Mean of total_egg_mass_num:", mean_value, "\n")
cat("Variance of total_egg_mass_num:", variance_value, "\n") # Variance >> Mean, making it better suited for over dispersed data

# add prevous week's snail size as covariate to control for
source("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/R2_growth_analysis.R") # cleaned master data sheet called all_data_growth
covariate <- all_data_growth_length_covariate_r2 %>% select(Jar_num, Color, Week_num, prev_snail_length_avg)
all_data_repro <- merge(all_data_repro, covariate, by = c("Jar_num", "Color", "Week_num"), all.x = TRUE)
all_data_repro <- all_data_repro %>% filter(!is.na(prev_snail_length_avg)) %>%  mutate(jar_id = interaction(Color, Jar_num, drop = TRUE))# 964 missing prev_snail_length_avg

# Negative Binomial Model for egg mass counts (glmmTMB package)
reproduction_model_re <- glmmTMB(total_egg_mass_num ~ Treatment + prev_snail_length_avg + (1|jar_id) , family = nbinom2, data = all_data_repro) # (1|jar) account for jar-specific variability over time
# summary(reproduction_model_1)

reproduction_model <- glmmTMB(total_egg_mass_num ~ Treatment + prev_snail_length_avg , family = nbinom2, data = all_data_repro) # (1|jar) account for jar-specific variability over time
# summary(reproduction_model)

reproduction_model_snail_num_re <- glmmTMB(total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num + (1|jar_id), family = nbinom2, data = all_data_repro)
# summary(reproduction_model_snail_num_1)

reproduction_model_snail_num <- glmmTMB(total_egg_mass_num ~ Treatment + prev_snail_length_avg + snails_num , family = nbinom2, data = all_data_repro)
# summary(reproduction_model_snail_num)

AIC(reproduction_model_re, reproduction_model, reproduction_model_snail_num, reproduction_model_snail_num_re) # AIC of reproduction_model_snail_num, reproduction_model_snail_num_re doesn't differ much and re make more statistical sense

# Fixed effects (conditional model) summary → CSV
glmm_model_summary <- summary(reproduction_model_snail_num_re)

coef_table <- as.data.frame(glmm_model_summary$coefficients$cond)
coef_table <- data.frame(term = rownames(coef_table), coef_table, row.names = NULL)

output_dir <- "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2"

# write.csv(
#   coef_table,
#   file.path(output_dir, "R2_fecundity_model.csv"),
#   row.names = FALSE
# )

# Random-effect variances and standard deviations (correct accessor for glmmTMB)
glmmTMB::VarCorr(reproduction_model_snail_num_re)
sigma(reproduction_model_snail_num_re)

# Confirm random-effect grouping names
names(glmmTMB::VarCorr(reproduction_model_snail_num_re)$cond)


pairwise_comparisons_repro <- emmeans(reproduction_model_snail_num_re, pairwise ~ Treatment, type = "response")
summary(pairwise_comparisons_repro)

fecundity_pairwise <- as.data.frame(
  summary(pairwise_comparisons_repro)$contrasts
)

# write.csv(
#   fecundity_pairwise,
#   file.path(output_dir, "R2_fecundity_pairwise.csv"),
#   row.names = FALSE
# )


### 4treatment: marginalized egg mass production ----

emm_repro <- emmeans(reproduction_model_snail_num_re, ~ Treatment, type = "response")

# 2. Extract the 95% Confidence Intervals
df_95 <- as.data.frame(summary(emm_repro, level = 0.95)) %>% rename( "lower95"="asymp.LCL",
                                                                     "upper95"="asymp.UCL")

# 3. Extract the 50% Confidence Intervals
df_50 <- as.data.frame(summary(emm_repro, level = 0.50)) %>% rename( "lower95"="asymp.LCL",
                                                                     "upper95"="asymp.UCL")

# 4. Merge into a single final data frame
plot_df_repro_r2 <- data.frame(
  Treatment  = df_95[[1]],
  mean_egg = df_95[[2]],
  lower95    = df_95[[5]],
  upper95    = df_95[[6]],
  lower50    = df_50[[5]],
  upper50    = df_50[[6]]
)

# # marginalized mean egg mass per jar over a range of covariates
# would need to read in R file R1_growth_anlysis.R to obtain the df marg_df_r1
plot_df_repro_r1r2 <- bind_rows(
  "Treatment Round 1" = plot_df_repro_r1, 
  "Treatment Round 2" = plot_df_repro_r2, 
  .id = "source"
)

#label pairwise comparison with compact letter display
letters_df <- data.frame(
  Treatment = c("Control", "Micro", "Macro", "Macro+Micro", 
                "Seasoned Macro - Low", "Seasoned Macro - High", 
                "Virgin Macro - Low", "Virgin Macro - High"),
  cld = c("a", "a", "ab", "b", "bc", "c", "ab", "bc") # Replace with your real letters
)

# Merge letters into your main plotting dataframe
plot_df_repro_r1r2 <- plot_df_repro_r1r2 %>%
  left_join(letters_df, by = "Treatment")

treatment_order <- c(
  "Control", 
  "Micro", 
  "Macro", 
  "Macro+Micro",
  "Seasoned Macro - Low", 
  "Seasoned Macro - High", 
  "Virgin Macro - Low", 
  "Virgin Macro - High"
)

# Apply order 
plot_df_repro_r1r2$Treatment <- factor(
  plot_df_repro_r1r2$Treatment, 
  levels = treatment_order
)

### plot marginalized effect ----
plot_repro_r1r2<- ggplot(plot_df_repro_r1r2, aes(x = Treatment, y = mean_egg, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower50, ymax = upper50),
    linewidth = 2.5,      # thicker than 95% CI
    width = 0,        # can be slightly wider
    # alpha = 0.5          # optional: slightly transparent
  ) +
  geom_errorbar(
    aes(ymin = lower95, ymax = upper95),
    linewidth = 1, width = 0.15
  ) +
  geom_text(aes(label = cld, y = upper95), 
            vjust = -0.4,           # Push letters above the error bar
            color = "black",        # Make letters black for readability
            fontface = "bold",
            size = 5) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Marginal mean egg mass count",
    color = "Treatment") + # Treatment effects marginalised over growth state
  scale_color_manual(values = c("Control" = "#6baf78", 
                                "Micro" = "#5a9fd6", 
                                "Macro" = "#a5855f", 
                                "Macro+Micro" = "#9d7ca5",
                                "Seasoned Macro - Low" = "#d4c2a8",
                                "Seasoned Macro - High" = "#a5855f",
                                "Virgin Macro - Low" = "#bcbcbc",
                                "Virgin Macro - High" = "#8a8a8a")) +
  theme(panel.spacing = unit(1.5, "lines"),
        strip.text.x = element_text(size = 15, face = "bold", color = "black"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(color = "black", size = 11),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  facet_wrap(~source, scales = "free_x", strip.position = "bottom") +
  scale_x_discrete(labels = c("Control"="Control", 
                              "Micro"="Virgin Micro", 
                              "Macro" = "Biofouled Macro", 
                              "Macro+Micro"="Biofouled Macro\n+ Virgin Micro",
                              "Seasoned Macro - Low" = "Biofouled Macro \n- Low",
                              "Seasoned Macro - High" ="Biofouled Macro \n- High",
                              "Virgin Macro - Low"="Virgin Macro \n- Low",
                              "Virgin Macro - High"="Virgin Macro \n- High")) 

plot_repro_r1r2


# Save the plot to a PNG file
ggsave(filename = "R1R2_fecundity_plot_marginalized_mean_seperated.png", plot = plot_repro_r1r2, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)

### Visualization of reproductive output by individual treatment with varying plastic abundance #########
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
  geom_line(linewidth = 1.5) +  # add line connecting dots for each Color group
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), linewidth = 0.5, width = 0.15) +
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
  theme_minimal()+ 
  theme(legend.position = "top",
        legend.direction = "horizontal",
        strip.text.x = element_blank(),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(color = "black", size = 11),
        legend.title = element_text(size = 10),
        legend.text = element_text(color = "black", size = 10))

R2_fecundity_plot_4treatments
# Save the plot to a PNG file
# ggsave(filename = "R2_fecundity_plot_4treatments_toplegend.png", plot = R2_fecundity_plot_4treatments, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)

### grouped: Reproductive output analysis by grouped treatment disregard varying plastic abundance #########
# Negative Binomial Model for egg mass counts (glmmTMB package)
reproduction_model_snail_num_re_group <- glmmTMB(total_egg_mass_num ~ Group + prev_snail_length_avg + snails_num + (1|jar_id), family = nbinom2, data = all_data_repro)
summary(reproduction_model_snail_num_re_group)

# Fixed effects (conditional model) summary → CSV
glmm_model_summary <- summary(reproduction_model_snail_num_re_group)

coef_table <- as.data.frame(glmm_model_summary$coefficients$cond)
coef_table <- data.frame(term = rownames(coef_table), coef_table, row.names = NULL)

output_dir <- "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2"

# write.csv(
#   coef_table,
#   file.path(output_dir, "R2_fecundity_model_grouped.csv"),
#   row.names = FALSE
# )

# Random-effect variances and standard deviations (correct accessor for glmmTMB)
glmmTMB::VarCorr(reproduction_model_snail_num_re_group)
sigma(reproduction_model_snail_num_re_group)

# Confirm random-effect grouping names
names(glmmTMB::VarCorr(reproduction_model_snail_num_re_group)$cond)


pairwise_comparisons_repro <- emmeans(reproduction_model_snail_num_re_group, pairwise ~ Group, type = "response")
summary(pairwise_comparisons_repro)

fecundity_pairwise <- as.data.frame(
  summary(pairwise_comparisons_repro)$contrasts
)

# write.csv(
#   fecundity_pairwise,
#   file.path(output_dir, "R2_fecundity_pairwise_group.csv"),
#   row.names = FALSE
# )


### Grouped: Visualization of reproductive output by individual treatment without varying plastic abundance #########
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
  geom_line(linewidth = 1.5) +  # add line connecting dots for each Color group
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), linewidth = 0.5, width = 0.15) +
  # geom_smooth(method = "loess", se = FALSE) +
  # geom_smooth(method = "loess", se = TRUE, aes(fill = Treatment), alpha = 0.2) +  # with se
  scale_color_manual(name = "Treatment",
                     values = c("Control" = "#6baf78", "Seasoned Macro" = "#a5855f", "Virgin Macro" = "#8a8a8a"),
                     labels = c("Control", 
                                "Biofouled Macro", 
                                "Virgin Macro")) +
  scale_fill_manual(values = c("Control" = "#6baf78", "Seasoned Macro" = "#a5855f", "Virgin Macro" = "#8a8a8a")) +
  scale_x_continuous(breaks = 1:13, labels = 1:13)+
  labs(
    x = "Weeks",
    y = "Mean Total Number of Egg Mass per Jar") +
  theme_minimal()+ 
  theme(legend.position = "top",
        legend.direction = "horizontal",
        strip.text.x = element_blank(),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(color = "black", size = 11),
        legend.title = element_text(size = 10),
        legend.text = element_text(color = "black", size = 10))

R2_fecundity_plot_grouped
# Save the plot to a PNG file
ggsave(filename = "R2_fecundity_plot_grouped_toplegend.png", plot = R2_fecundity_plot_grouped, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)


# ###when does the plastic effect starts to be significant####
# library(dplyr)
# library(glmmTMB)
# library(emmeans)
# 
# all_data_repro <- all_data_repro %>%
#   filter(Treatment%in%c("Seasoned Macro - High","Control")) %>% 
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
# anova(week2_null, week2_full) # not a sig factor of egg mass in week 2
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
# anova(week3_null, week3_full) # treatment is a sig factor of egg mass in week 3
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
# anova(week4_null, week4_full) #  a sig factor of egg mass in week 4
# emmeans(week4_full, pairwise ~ Treatment, type = "response") #  sig
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
# anova(week6_null, week6_full) # 
# emmeans(week6_full, pairwise ~ Treatment, type = "response") # 
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
# emmeans(week7_full, pairwise ~ Treatment, type = "response") #  sig
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
# anova(week8_null, week8_full) #  a sig factor of egg mass in week 
# emmeans(week8_full, pairwise ~ Treatment, type = "response") # sig
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
# anova(week9_null, week9_full) #  a sig factor 
# emmeans(week9_full, pairwise ~ Treatment, type = "response") #  sig
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
# anova(week11_null, week11_full) #  a sig factor 
# emmeans(week11_full, pairwise ~ Treatment, type = "response") #  sig
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
# anova(week12_null, week12_full) #  a sig factor of egg mass in week 5
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
