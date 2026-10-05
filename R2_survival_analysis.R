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

####### load data#########
source("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/R2_data_cleaning.R") # cleaned master data sheet called all_data_r2
source("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/R1_survival_analysis.R") 
####### build survival df #########
# Helper fxn to gather the weekly "Adult_death_k" and "Week_of_death_k" into long 
all_long_base <- all_data_r2 %>%
  select(Color, Jar_num, Week_num, starts_with("Adult_death_"), starts_with("Week_of_death_"))

# per-snail last week (based on Adult_death_k being recorded that week)
last_seen <- all_long_base %>%
  pivot_longer(starts_with("Adult_death_"),
               names_to = "snail_num", values_to = "death_flag",
               names_pattern = "Adult_death_(\\d+)") %>%
  group_by(Color, Jar_num, snail_num) %>%
  summarise(
    # Use last week where the snail column is present/recorded
    last_obs_week = suppressWarnings(max(Week_num[!is.na(death_flag)], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(last_obs_week = ifelse(is.finite(last_obs_week), last_obs_week, NA_real_))

# per-snail first recorded week of death
first_death_week <- all_long_base %>%
  pivot_longer(starts_with("Week_of_death_"),
               names_to = "snail_num", values_to = "week_of_death",
               names_pattern = "Week_of_death_(\\d+)") %>%
  group_by(Color, Jar_num, snail_num) %>%
  summarise(
    week_of_death = suppressWarnings(min(as.numeric(week_of_death), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(week_of_death = ifelse(is.finite(week_of_death), week_of_death, NA_real_))

# table(first_death_week$week_of_death, first_death_week$Color)

# join + add jar_id, treatment, and define outcome
per_snail <- last_seen %>%
  full_join(first_death_week, by = c("Color","Jar_num","snail_num")) %>%
  left_join(
    all_data_r2 %>% distinct(Color, Jar_num, jar_id, Treatment),
    by = c("Color","Jar_num")
  ) %>%
  mutate(
    survival_status = as.integer(!is.na(week_of_death)),
    week_survival   = ifelse(survival_status == 1, week_of_death, last_obs_week)
  ) %>%
  filter(is.finite(week_survival)) %>%
  mutate(
    Color    = factor(Color, levels = c("G","Y","O","W-1","W-3")),
    Treatment = factor(Treatment, levels = c("Control", "Seasoned Macro - Low", "Seasoned Macro - High",
                                             "Virgin Macro - Low", "Virgin Macro - High"))
  )

# Snails per jar (should be constant 3)
(per_snail %>% count(jar_id, name = "n_snails") %>% summarise(min=min(n_snails), median=median(n_snails), max=max(n_snails)))

# weeks survival
(per_snail %>% count(week_survival) %>% arrange(week_survival) %>% head(15))

# Create a new column for combined groups
per_snail$Group <- factor(with(per_snail, 
                               ifelse(Color %in% c("Y", "O"), "Seasoned Macro",
                                      ifelse(Color %in% c("W-1", "W-3"), "Virgin Macro","G"))))

# Verify the new grouping
table(per_snail$Group, per_snail$week_of_death)

#check if number of snail alive is a significant factor:
# Split follow-up at every recorded death/censoring time
cut_times <- sort(unique(per_snail$week_survival))
cut_times <- cut_times[cut_times < max(per_snail$week_survival)]

surv_density <- survival::survSplit(
  Surv(week_survival, survival_status) ~ .,
  data = per_snail,
  cut = cut_times,
  start = "tstart",
  id = "snail_id"
) %>%
  group_by(jar_id, tstart) %>%
  mutate(n_alive_start = n()) %>%
  ungroup()

# Deaths at an interval's end reduce counts only in later intervals.
# Both models retain jar-clustered standard errors.
fit_no_density <- coxph(
  Surv(tstart, week_survival, survival_status) ~
    Treatment + cluster(jar_id),
  data = surv_density
)

fit_with_density <- coxph(
  Surv(tstart, week_survival, survival_status) ~
    Treatment + n_alive_start + cluster(jar_id),
  data = surv_density
)

# Descriptive fit comparison: negative AIC change favors adding density
fit_comparison <- data.frame(
  model = c("Treatment only", "Treatment + living-snail count"),
  AIC = c(AIC(fit_no_density), AIC(fit_with_density))
)
fit_comparison$AIC_change_vs_original <-
  fit_comparison$AIC - fit_comparison$AIC[1]

print(fit_comparison)

# Cluster-robust Wald test and HR per additional living snail
b <- coef(fit_with_density)["n_alive_start"]
se <- sqrt(vcov(fit_with_density)["n_alive_start", "n_alive_start"])

density_result <- data.frame(
  HR = exp(b),
  lower_95 = exp(b - qnorm(0.975) * se),
  upper_95 = exp(b + qnorm(0.975) * se),
  robust_p = 2 * pnorm(abs(b / se), lower.tail = FALSE)
)

print(density_result) # HR 1.94-3.20; sig higher mortality hazard per additional living snail


#### Build weekly survival intervals with starting snail counts 
# Split each snail's follow-up into (week start, week end] intervals.
last_week <- max(per_snail$week_survival)
weekly_cuts <- seq_len(last_week)
weekly_cuts <- weekly_cuts[weekly_cuts < last_week]

per_snail_weekly <- survival::survSplit(
  Surv(week_survival, survival_status) ~ .,
  data = per_snail,
  cut = weekly_cuts,
  start = "week_start",
  id = "snail_id"
) %>%
  rename(week_end = week_survival) %>%
  group_by(jar_id, week_start) %>%
  mutate(n_alive_start = n()) %>%
  ungroup() %>%
  arrange(jar_id, snail_id, week_start)

###### Following analysis separate plastic abundance : "Treatment" #########
######## Cox survival models separate abundance : "Treatment" #########
# Cluster-robust Cox (jar-level clustering)
# Jar-clustered robust standard errors
fit_cr <- coxph(
  Surv(week_start, week_end, survival_status) ~
    Treatment + n_alive_start + cluster(jar_id),
  data = per_snail_weekly
)

# No jar adjustment
fit_fix <- coxph(
  Surv(week_start, week_end, survival_status) ~
    Treatment + n_alive_start,
  data = per_snail_weekly
)

# Jar random effect (shared frailty)
fit_re <- coxme(
  Surv(week_start, week_end, survival_status) ~
    Treatment + n_alive_start + (1 | jar_id),
  data = per_snail_weekly
)

(summary(fit_re))
(summary(fit_fix)) # no jar-level random effect

# boundary issues make LRT conservative; report AIC
(AIC(fit_fix))
(AIC(fit_re))
AIC(fit_cr)

anova(fit_fix, fit_re) # re model better

# # Saves re model results as a CSV
# R2_survival_model_4treatments <- data.frame()
# 
# write.csv(R2_survival_model_4treatments,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_survival_model_4treatments.csv")
####### pairwise comparisons HRs separate abundance : "Treatment" #########
# pairwise HRs with robust vcov for coxph (clustered)
emm  <- emmeans(fit_re, ~ Treatment, vcov. = vcov(fit_re))
pw   <- contrast(emm, "pairwise", adjust = "tukey")

pw_sum <- as.data.frame(summary(pw))
# 95% CI
pw_ci95 <- as.data.frame(confint(pw, adjust = "tukey", level = 0.95))
lower95 <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(pw_ci95))[1]
upper95 <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(pw_ci95))[1]

pw_ci95 <- pw_ci95 %>%
  dplyr::rename(
    LCL95 = !!lower95,
    UCL95 = !!upper95
  )

# 50% CI
pw_ci50 <- as.data.frame(confint(pw, adjust = "tukey", level = 0.50))
lower50 <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(pw_ci50))[1]
upper50 <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(pw_ci50))[1]

pw_ci50 <- pw_ci50 %>%
  dplyr::rename(
    LCL50 = !!lower50,
    UCL50 = !!upper50
  )

pairs_adj <- pw_ci95 %>%
  dplyr::select(contrast, estimate, SE, df, LCL95, UCL95) %>%
  dplyr::left_join(
    pw_ci50 %>% dplyr::select(contrast, LCL50, UCL50),
    by = "contrast"
  ) %>%
  dplyr::left_join(
    pw_sum %>% dplyr::select(contrast, p.value),
    by = "contrast"
  ) %>%
  dplyr::mutate(
    # emmeans reports log-HR as A − B → convert to B vs A
    HR        = exp(-estimate),
    
    HR.low95  = exp(-UCL95),
    HR.high95 = exp(-LCL95),
    
    HR.low50  = exp(-UCL50),
    HR.high50 = exp(-LCL50),
    
    A = sub(" - .*", "", contrast),
    B = sub(".* - ", "", contrast),
    label = sprintf("%s - %s", B, A),
    perc_reduction = 1 - HR,
    estimate = -estimate
  ) %>%
  dplyr::select(
    label, estimate, SE,
    HR,
    HR.low95, HR.high95,
    HR.low50, HR.high50,
    p.value
  )

# Pairwise HRs (Tukey-adjusted; robust vcov; B vs A)
(pairs_adj)
write.csv(pairs_adj,  file = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/R2_survival_pairwise_4treatments_corrected.csv")

pairs_adj_r2 <- pairs_adj[1:4,] %>% mutate(label = c("Biofouled Macro - Low", "Biofouled Macro - High", "Virgin Macro - Low", "Virgin Macro - High"))

# would need to read in R file R1_survival_aanlysis.R to obtain the df pairs_adj_r1 done in section "load data"

pairs_adj_r1r2 <- bind_rows(
  "Treatment Round 1" = pairs_adj_r1, 
  "Treatment Round 2" = pairs_adj_r2, 
  .id = "source"
) 
#label pairwise comparison with compact letter display
letters_df <- data.frame(
  label = c(
    "Biofouled Macro + Virgin Micro",
    "Biofouled Macro", 
    "Virgin Micro",
    "Virgin Macro - High",
    "Virgin Macro - Low", 
    "Biofouled Macro - High",
    "Biofouled Macro - Low"
  ),
  cld = c("ab", "a", "b", "a", "a", "a", "a") 
)

# Merge letters into your main plotting dataframe
pairs_adj_r1r2 <- pairs_adj_r1r2 %>%
  left_join(letters_df, by = "label")

treatment_order <- c("Virgin Micro",
  
  "Biofouled Macro", 
  "Biofouled Macro + Virgin Micro",
  "Biofouled Macro - Low",
  "Biofouled Macro - High",
  "Virgin Macro - Low",
  "Virgin Macro - High"
 
)

# Apply this order to your dataframe
pairs_adj_r1r2$label <- factor(
  pairs_adj_r1r2$label, 
  levels = rev(treatment_order)
)

hr_plot_r1r2 <- ggplot(
  pairs_adj_r1r2,
  aes(y = label, x = HR, color = label)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(xmin = HR.low50, xmax = HR.high50),
    orientation = "y",
    linewidth = 2.5,
    width = 0
  ) +
  geom_errorbar(
    aes(xmin = HR.low95, xmax = HR.high95),
    orientation = "y",
    linewidth = 1,
    width = 0.15
  ) +
  geom_text(
    aes(label = cld, x = HR.high95),
    hjust = -0.4,
    vjust = 0.5,
    color = "black",
    fontface = "bold",
    size = 5
  ) +
  theme_classic() +
  labs(
    x = "Hazard ratio",
    y = NULL,
    color = "Treatment"
  ) +
  scale_color_manual(values = c(
    "Control" = "#6baf78",
    "Virgin Micro" = "#5a9fd6",
    "Biofouled Macro" = "#a5855f",
    "Biofouled Macro + Virgin Micro" = "#9d7ca5",
    "Biofouled Macro - Low" = "#d4c2a8",
    "Biofouled Macro - High" = "#a5855f",
    "Virgin Macro - Low" = "#bcbcbc",
    "Virgin Macro - High" = "#8a8a8a"
  )) +
  theme(
    panel.spacing = unit(1.5, "lines"),
    strip.text.y = element_text(
      size = 15, face = "bold", color = "black"
    ),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(color = "black", size = 11),
    legend.position = "none"
  ) +
  facet_grid(
    rows = vars(source),
    scales = "free_y",
    switch = "y"
  ) +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  scale_y_discrete(labels = c(
    "Control" = "Control",
    "Virgin Micro" = "Virgin Micro",
    "Biofouled Macro" = "Biofouled Macro",
    "Biofouled Macro + Virgin Micro" = "Biofouled Macro\n+ Virgin Micro",
    "Biofouled Macro - Low" = "Biofouled Macro \n- Low",
    "Biofouled Macro - High" = "Biofouled Macro \n- High",
    "Virgin Macro - Low" = "Virgin Macro \n- Low",
    "Virgin Macro - High" = "Virgin Macro \n- High"
  ))

hr_plot_r1r2

ggsave(filename = "R1R2_survival_plot_HR_treatment.png", plot = hr_plot_r1r2, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)


####### Cox model diagnostics for "Treatment" #########
# porportional hazard assumption met?
ph_test <- cox.zph(fit_fix)
# Global and term-wise Schoenfeld PH tests (fixed-effects model)
(ph_test) # not significant with p > 0.05

####### plot survival curves for "Treatment"#########
# Surv Curves: KM (Kaplan-Meier estimator) by treatment, Cox-predicted (cluster-robust), stratified Cox (diagnostic)
sf_km <- survfit(Surv(week_survival, survival_status) ~ Treatment, data = per_snail)

km_plot <- ggsurvplot(
  sf_km, data = per_snail, conf.int = TRUE,conf.int.alpha = 0.15,
  risk.table = TRUE, risk.table.height = 0.22,
  ggtheme = theme_minimal(base_size = 10), 
  xlab = "Weeks", 
  ylab = "Survival Probability",
  legend.title = "Treatment",
  legend.labs = c("Control", "Biofouled Macro - Low", "Biofouled Macro - High",
                  "Virgin Macro - Low", "Virgin Macro - High"),
  palette = c("#6baf78", "#d4c2a8", "#a5855f", "#bcbcbc", "#8a8a8a"),
  xlim = c(0, 13),           
  break.time.by = 1)

km_plot$plot <- km_plot$plot + theme(strip.text.x = element_blank(),
                                     axis.title = element_text(size = 15, face = "bold"),
                                     axis.text = element_text(color = "black", size = 11),
                                     legend.title = element_text(size = 10),
                                     legend.text = element_text(color = "black", size = 10)) 

km_plot$plot <- km_plot$plot +  scale_y_continuous(limits = c(0.3, 1.0),
                                                   breaks = seq(0.3, 1, by = 0.1))

(km_plot$plot)
# Save the plot to a PNG file
ggsave(filename = "R2_survival_plot_4treatments_truncated.png", plot = km_plot$plot, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)

(km_plot$table)
ggsave(filename = "R2_survival_plot_4treatments_num_at_risk.png", plot = km_plot$table, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)

###### Following analysis group plastic abundance : "Group" #########
####### Cox survival models group abundance: "Group" #########
# Cluster-robust Cox (jar-level clustering)
fit_cr <- coxph(Surv(week_start, week_end, survival_status) ~ Group + n_alive_start + cluster(jar_id),
                data = per_snail_weekly, ties = "efron")

# Fixed-effects Cox (for PH tests and baseline comparison)
fit_fix <- coxph(Surv(week_start, week_end, survival_status) ~ Group + n_alive_start,
                 data = per_snail_weekly, ties = "efron")

# Random-effects Cox (frailty) — optional model comparison
fit_re  <- coxme(Surv(week_start, week_end, survival_status) ~ Group + n_alive_start + (1|jar_id),
                 data = per_snail_weekly)

(summary(fit_re))
# (summary(fit_fix)) # no jar-level random effect

# boundary issues make LRT conservative; report AIC
(AIC(fit_fix))
(AIC(fit_re))
(AIC(fit_cr))
anova(fit_fix, fit_re) # re significantly better

####### pairwise comparisons HRs group abundance: "Group" #########
# pairwise HRs with robust vcov for coxph (clustered)
emm  <- emmeans(fit_re, ~ Group, vcov. = vcov(fit_re))
pw   <- contrast(emm, "pairwise", adjust = "tukey")

pw_sum <- as.data.frame(summary(pw))
pw_ci  <- as.data.frame(confint(pw, adjust = "tukey"))

lower_name <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(pw_ci))[1]
upper_name <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(pw_ci))[1]

pairs_adj <- pw_ci %>%
  dplyr::select(contrast, estimate, SE, df, !!lower_name, !!upper_name) %>%
  dplyr::left_join(pw_sum %>% dplyr::select(contrast, p.value), by = "contrast") %>%
  dplyr::mutate(
    # emmeans on a Cox PH uses log-HR with the internal "A - B" direction
    HR      = exp(-estimate),                       # report B vs A
    HR.low  = exp(-.data[[upper_name]]),
    HR.high = exp(-.data[[lower_name]]),
    A = sub(" - .*", "", contrast),
    B = sub(".* - ", "", contrast),
    label = sprintf("%s - %s", B, A),
    perc_reduction = 1 - HR,
    estimate = -estimate
  ) %>%
  dplyr::select(label, estimate,
                SE, HR, HR.low, HR.high, p.value)

# Pairwise HRs (Tukey-adjusted; robust vcov; B vs A)
(pairs_adj)
write.csv(pairs_adj,  file = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/R2_survival_pairwise_grouped_corrected.csv")

####### Cox model diagnostics for "Group"#########
# porportional hazard assumption met?
ph_test <- cox.zph(fit_fix)
# Global and term-wise Schoenfeld PH tests (fixed-effects model)
(ph_test) # met assumption with p > 0.05

####### plot survival curves for "Group"#########
# Surv Curves: KM (Kaplan-Meier estimator) by treatment, Cox-predicted (cluster-robust), stratified Cox (diagnostic)
sf_km <- survfit(Surv(week_survival, survival_status) ~ Group, data = per_snail)

km_plot <- ggsurvplot(
  sf_km, data = per_snail, conf.int = TRUE,conf.int.alpha = 0.2,
  risk.table = TRUE, risk.table.height = 0.22,
  ggtheme = theme_minimal(base_size = 10), 
  xlab = "Weeks", 
  ylab = "Survival Probability",
  legend.title = "Treatment",
  legend.labs = c("Control", "Biofouled Macro", "Virgin Macro"),
  palette = c("#6baf78", "#a5855f", "#8a8a8a"),
  xlim = c(0, 13),           
  break.time.by = 1 
)

km_plot$plot <- km_plot$plot + theme(strip.text.x = element_blank(),
                                     axis.title = element_text(size = 15, face = "bold"),
                                     axis.text = element_text(color = "black", size = 11),
                                     legend.title = element_text(size = 11),
                                     legend.text = element_text(color = "black", size = 11)) 

km_plot$plot <- km_plot$plot +  scale_y_continuous(limits = c(0.3, 1.0),
                                                   breaks = seq(0.3, 1, by = 0.1))

# plots!
(km_plot$plot)# Save the plot to a PNG file
ggsave(filename = "R2_survival_plot_grouped_toplegend.png", plot = km_plot$plot, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)
(km_plot$table)
ggsave(filename = "R2_survival_plot_grouped_num_at_risk.png", plot = km_plot$table, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)


######## death by week and treatment plot #########

# Deaths table by week and treatment
death_table <- per_snail %>%
  mutate(week_of_death = ifelse(is.na(week_of_death), "alive", as.character(week_of_death))) %>% 
  group_by(week_of_death, Treatment) %>%
  summarise(deaths = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Treatment, values_from = deaths, values_fill = list(deaths = 0)) %>%
  mutate(week_of_death = factor(week_of_death, levels = c(as.character(2:13), "alive"))) %>%
  arrange(week_of_death) 

death_table

# Convert week_of_deadth to numeric for plotting (excluding "alive")
death_table <- death_table %>%
  filter(week_of_death != "alive") %>%
  mutate(week_of_death = as.numeric(week_of_death)+1)

# Reshape the data for plotting
death_long <- death_table %>%
  pivot_longer(cols = c("Control", "Seasoned Macro - Low", "Seasoned Macro - High",
                        "Virgin Macro - Low", "Virgin Macro - High"), names_to = "Treatment", values_to = "deaths")

# Plot for absolute death
ggplot(death_long, aes(x = week_of_death, y = deaths, color = Treatment)) +
  geom_line() +
  geom_point() +
  labs(title = "Number of Deaths per Week by Color",
       x = "Week of Death",
       y = "Number of Deaths") +
  scale_x_continuous(breaks = seq(0, 13, by = 1)) +
  theme_minimal()

# Calculate the total number of snails at the start for each color
total_snails <- per_snail %>%
  group_by(Treatment) %>%
  summarise(total = n(), .groups = 'drop')

# Merge with total snails to calculate relative deaths
death_long <- death_long %>%
  left_join(total_snails, by = "Treatment") %>%
  group_by(Treatment) %>%
  mutate(cumulative_death = cumsum(deaths)) 

death_long <- death_long %>%
  mutate(relative_deaths = deaths / total)

# Plot relative death
ggplot(death_long, aes(x = week_of_death, y = relative_deaths, color = Treatment)) +
  geom_line() +
  geom_point() +
  labs(title = "Relative Number of Deaths per Week",
       x = "Week of Death",
       y = "Relative Number of Deaths") +
  scale_x_continuous(breaks = seq(0, 13, by = 1)) +
  theme_minimal()
