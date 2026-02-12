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
source("/Users/aoyu/Desktop/Snail_Data/R2_data_cleaning.R") # cleaned master data sheet called all_data_r2
source("/Users/aoyu/Desktop/Snail_Data/R1_survival_analysis.R") 
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

###### Following analysis separate plastic abundance : "Treatment" #########
######## Cox survival models separate abundance : "Treatment" #########
# Cluster-robust Cox (jar-level clustering)
fit_cr <- coxph(Surv(week_survival, survival_status) ~ Treatment + cluster(jar_id),
                data = per_snail, ties = "efron")

# Fixed-effects Cox (for PH tests and baseline comparison)
fit_fix <- coxph(Surv(week_survival, survival_status) ~ Treatment,
                 data = per_snail, ties = "efron")

# Random-effects Cox (frailty) — optional model comparison
fit_re  <- coxme(Surv(week_survival, survival_status) ~ Treatment + (1|jar_id),
                 data = per_snail)

(summary(fit_cr))
(summary(fit_fix)) # no jar-level random effect

# boundary issues make LRT conservative; report AIC
(AIC(fit_fix))
(AIC(fit_re))
anova(fit_fix, fit_re) # fixed model better

# # Saves model results as a CSV
# R2_survival_model_4treatments <- data.frame(
#   term = c("Biofouled Macro - Low", "Biofouled Macro - High", "Virgin Macro - Low", "Virgin Macro - High"),
#   estimate = c(-0.10180, -0.04495, 0.01800, 0.06185),
#   standard.error = c(0.14710, 0.14666,0.16915, 0.16647),
#   HR = c(0.90321, 0.95605, 1.01816, 1.06381),
#   HR.low = c(0.6770, 0.7172, 0.7309, 0.7677),
#   HR.high = c(1.205, 1.274, 1.418, 1.474),
#   p.value = c(0.489,0.759, 0.915, 0.710 ))
# 
# write.csv(R2_survival_model_4treatments,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_survival_model_4treatments.csv")
####### pairwise comparisons HRs separate abundance : "Treatment" #########
# pairwise HRs with robust vcov for coxph (clustered)
Vrob <- sandwich::sandwich(fit_cr)   

emm  <- emmeans(fit_cr, ~ Treatment, vcov. = Vrob)
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
    perc_reduction = 1 - HR
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
write.csv(pairs_adj,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_survival_pairwise_4treatments.csv")

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
  # 95% CI (thin, background)
  geom_errorbarh(
    aes(xmin = HR.low95, xmax = HR.high95),
    width = 0.15,
    linewidth = 1,
    # alpha = 0.7
  ) +
  
  # 50% CI (thick, foreground)
  geom_errorbarh(
    aes(xmin = HR.low50, xmax = HR.high50),
    width = 0,
    linewidth = 2.5
  ) +
  
  # geom_vline(
  #   xintercept = 1,
  #   linetype = "dashed",
  #   color = "grey40"
  # ) +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4),
    limits = c(
      min(pairs_adj_r1r2$HR.low95,  na.rm = TRUE),
      max(pairs_adj_r1r2$HR.high95, na.rm = TRUE)
    )
  ) +
  facet_grid(
    rows = vars(source),
    scales = "free_y",
    switch = "y"
  ) +
  theme_classic() +
  labs(
    x = "Hazard ratio",
    y = NULL
  ) +
  theme(
    panel.spacing.y = unit(1.5, "lines"),
    strip.text.y = element_text(size = 15, face = "bold"),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(color = "black", size = 11),
    legend.position = "none"
  ) +
  scale_color_manual(values = c( 
                                "Virgin Micro" = "#5a9fd6", 
                                "Biofouled Macro" = "#a5855f", 
                                "Biofouled Macro + Virgin Micro" = "#9d7ca5",
                                "Biofouled Macro - Low" = "#d4c2a8",
                                "Biofouled Macro - High" = "#a5855f",
                                "Virgin Macro - Low" = "#bcbcbc",
                                "Virgin Macro - High" = "#8a8a8a")) +
  geom_text(
    aes(
      x = HR.high95,
      y = label,
      label = cld
    ),
    hjust = -1,        # push text to the right of the CI
    vjust = 0.5,
    color = "black",
    fontface = "bold",
    size = 5
  ) +
  scale_y_discrete(labels = c(
    "Virgin Micro" = "Virgin Micro",
    "Biofouled Macro" = "Biofouled Macro",
    "Biofouled Macro + Virgin Micro" = "Biofouled Macro\n+ Virgin Micro",
    "Biofouled Macro - Low" = "Biofouled Macro - Low",
    "Biofouled Macro - High" = "Biofouled Macro - High",
    "Virgin Macro - Low" = "Virgin Macro - Low",
    "Virgin Macro - High" = "Virgin Macro - High"
  ))

hr_plot_r1r2

ggsave(filename = "R1R2_survival_plot_HR_seperated.png", plot = hr_plot_r1r2, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)


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
ggsave(filename = "R2_survival_plot_4treatments_truncated.png", plot = km_plot$plot, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)
(km_plot$table)
ggsave(filename = "R2_survival_plot_4treatments_num_at_risk.png", plot = km_plot$table, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)

### Cox predicted curves by treatment from clustered model -----
#(survfit.coxph works with new data to get conditional curves by level)
# a model-based estimate that accounts for clustering within jars
levels_T <- levels(per_snail$Treatment)
pretty_labels <- levels_T
names(pretty_labels) <- levels_T

build_sf_row <- function(one_level) {
  nd_one <- data.frame(Treatment = factor(one_level, levels = levels_T))
  sf_one <- survfit(fit_cr, newdata = nd_one)
  ss_one <- survminer::surv_summary(sf_one, data = nd_one)
  ss_one$Treatment <- pretty_labels[[one_level]]
  ss_one
}
df_mod <- dplyr::bind_rows(lapply(levels_T, build_sf_row))

p_cox <- ggplot(df_mod, aes(x = time, y = surv, color = Treatment)) +
  geom_step() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Treatment),
              alpha = 0.20, colour = NA) +
  labs(x = "Weeks", y = "Survival probability",
       color = "Treatment", fill = "Treatment",
       title = "Cox-predicted survival (cluster-robust fit)") +
  theme_minimal(base_size = 10)

# check KM vs Cox by treatment
df_km <- survminer::surv_summary(sf_km, data = per_snail) %>%
  mutate(Treatment = sub("^Treatment=", "", as.character(strata))) %>%
  transmute(time, surv, lower, upper, Treatment)

cox_df <- df_mod %>% transmute(time, surv, lower, upper, Treatment)

km_df <- df_km; km_df$Source  <- "KM"
cx_df <- cox_df; cx_df$Source <- "Cox"
overlay_df <- bind_rows(km_df, cx_df)

p_overlay <- ggplot(overlay_df, aes(time, surv, color = Source, linetype = Source)) +
  geom_step() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Source),
              alpha = 0.12, colour = NA) +
  scale_color_manual(values = c("KM" = "black", "Cox" = "steelblue4")) +
  scale_fill_manual( values = c("KM" = "grey70", "Cox" = "steelblue2")) +
  scale_linetype_manual(values = c("KM" = "solid", "Cox" = "dashed")) +
  facet_wrap(~ Treatment, ncol = 2) +
  labs(x = "Weeks", y = "Survival probability", color = "", fill = "", linetype = "",
       title = "KM vs Cox (diagnostic only)") +
  theme_minimal(base_size = 10)

# stratified Cox (matches KM so just diagnostic)
sf_strat <- survfit(Surv(week_survival, survival_status) ~ strata(Treatment), data = per_snail)
df_strat <- survminer::surv_summary(sf_strat, data = per_snail) %>%
  mutate(Treatment = sub("^strata\\(Treatment\\)=?", "", strata))

p_strat <- ggplot() +
  geom_step(data = df_km,   aes(time, surv, color = Treatment), alpha = 0.7) +
  geom_step(data = df_strat, aes(time, surv, color = Treatment), linetype = "dotted") +
  labs(x = "Weeks", y = "Survival probability", color = "Treatment",
       title = "Stratified Cox (dotted) over KM (solid) — diagnostic") +
  theme_minimal(base_size = 10)

# plots!
# (p_cox)
# (p_overlay) # why macro+ micro not starts at week 1? 
# (p_strat) # why macro+ micro not starts at week 1? 


####### Omitted since assumption met:RMST (restricted mean survival time) for "Treatment"#########
# # RMST since PH not met (primary estimand), τ = 13 weeks
# tau <- 13  # common follow-up window; adjust as needed/wanted
# 
# # pairwise comparisons for RMST
# pair_levels <- combn(levels_T, 2, simplify = FALSE)
# pair_labels <- setNames(levels_T, levels_T)
# 
# robust_extract <- function(fit) {
#   # Extract RMST estimates for each arm and the unadjusted difference row
#   arm0 <- as.data.frame(fit$RMST.arm0$result)
#   arm1 <- as.data.frame(fit$RMST.arm1$result)
#   get_est <- function(df) {
#     idx <- if (length(rownames(df))) which(grepl("RMST", rownames(df), ignore.case = TRUE))[1] else NA_integer_
#     if (is.na(idx)) idx <- 1
#     est_col <- intersect(c("Est.", "Estimate", "est", "RMST"), names(df))
#     if (!length(est_col)) est_col <- names(df)[1]
#     as.numeric(df[idx, est_col[1]])
#   }
#   est0 <- get_est(arm0); est1 <- get_est(arm1)
#   
#   ud <- as.data.frame(fit$unadjusted.result)
#   r_ud <- if (length(rownames(ud))) which(grepl("^\\s*RMST", rownames(ud), ignore.case = TRUE))[1] else NA_integer_
#   if (is.na(r_ud)) r_ud <- 1
#   
#   est_col   <- intersect(c("Est.", "Estimate", "est", "Diff", "Difference"), names(ud)); if (!length(est_col)) est_col <- names(ud)[1]
#   lower_col <- names(ud)[grep("lower|lcl", names(ud), ignore.case = TRUE)][1]; if (is.na(lower_col)) lower_col <- names(ud)[max(2, ncol(ud)-2)]
#   upper_col <- names(ud)[grep("upper|ucl", names(ud), ignore.case = TRUE)][1]; if (is.na(upper_col)) upper_col <- names(ud)[max(3, ncol(ud)-1)]
#   p_col     <- names(ud)[grep("^p($|\\b)|p.value", names(ud), ignore.case = TRUE)][1]; if (is.na(p_col))     p_col     <- names(ud)[ncol(ud)]
#   
#   c(
#     RMST_ref   = est0,
#     RMST_comp  = est1,
#     dRMST      = as.numeric(ud[r_ud, est_col[1]]),
#     dRMST_low  = as.numeric(ud[r_ud, lower_col]),
#     dRMST_high = as.numeric(ud[r_ud, upper_col]),
#     p_value    = as.numeric(ud[r_ud, p_col])
#   )
# }
# 
# rmst_pairwise <- purrr::map_dfr(pair_levels, function(lv) {
#   lv1 <- lv[1]; lv2 <- lv[2]
#   subdat <- per_snail %>%
#     filter(Treatment %in% c(lv1, lv2)) %>%
#     mutate(grp = factor(Treatment, levels = c(lv1, lv2)))  # lv1 = 0, lv2 = 1
#   gnum <- as.integer(subdat$grp == lv2)
#   
#   fit <- survRM2::rmst2(time = subdat$week_survival, status = subdat$survival_status, arm = gnum, tau = tau)
#   vals <- robust_extract(fit)
#   
#   tibble::tibble(
#     comparison = sprintf("%s vs %s", pair_labels[[lv2]], pair_labels[[lv1]]),
#     RMST_ref   = vals["RMST_ref"],
#     RMST_comp  = vals["RMST_comp"],
#     dRMST      = vals["dRMST"],
#     dRMST_low  = vals["dRMST_low"],
#     dRMST_high = vals["dRMST_high"],
#     p_value    = vals["p_value"]
#   )
# })
# 
# rmst_pretty <- rmst_pairwise %>%
#   mutate(across(c(RMST_ref, RMST_comp, dRMST, dRMST_low, dRMST_high), ~round(as.numeric(.), 2))) %>%
#   mutate(`95% CI` = sprintf("[%.2f, %.2f]", dRMST_low, dRMST_high)) %>%
#   dplyr::select(comparison, RMST_ref, RMST_comp, dRMST, `95% CI`, p_value) %>%
#   arrange(p_value) %>%
#   mutate(FDR_BH = p.adjust(p_value, method = "BH"))

# # RMST change (weeks) up to tau, comparator vs reference (BH-adjusted)
# rmst_pretty
###### Following analysis group plastic abundance : "Group" #########
####### Cox survival models group abundance: "Group" #########
# Cluster-robust Cox (jar-level clustering)
fit_cr <- coxph(Surv(week_survival, survival_status) ~ Group + cluster(jar_id),
                data = per_snail, ties = "efron")

# Fixed-effects Cox (for PH tests and baseline comparison)
fit_fix <- coxph(Surv(week_survival, survival_status) ~ Group,
                 data = per_snail, ties = "efron")

# Random-effects Cox (frailty) — optional model comparison
fit_re  <- coxme(Surv(week_survival, survival_status) ~ Group + (1|jar_id),
                 data = per_snail)

(summary(fit_cr))
(summary(fit_fix)) # no jar-level random effect

# boundary issues make LRT conservative; report AIC
(AIC(fit_fix))
(AIC(fit_re))
anova(fit_fix, fit_re)

####### pairwise comparisons HRs group abundance: "Group" #########
# pairwise HRs with robust vcov for coxph (clustered)
Vrob <- sandwich::sandwich(fit_cr)   

emm  <- emmeans(fit_cr, ~ Group, vcov. = Vrob)
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
    perc_reduction = 1 - HR
  ) %>%
  dplyr::select(label, estimate,
                SE, HR, HR.low, HR.high, p.value)

# Pairwise HRs (Tukey-adjusted; robust vcov; B vs A)
(pairs_adj)
write.csv(pairs_adj,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_survival_pairwise_grouped.csv")

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
ggsave(filename = "R2_survival_plot_grouped_toplegend.png", plot = km_plot$plot, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)
(km_plot$table)
ggsave(filename = "R2_survival_plot_grouped_num_at_risk.png", plot = km_plot$table, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)

### Cox predicted curves by treatment from clustered model ----
#(survfit.coxph works with new data to get conditional curves by level)
# a model-based estimate that accounts for clustering within jars
levels_T <- levels(per_snail$Group)
pretty_labels <- levels_T
names(pretty_labels) <- levels_T

build_sf_row <- function(one_level) {
  nd_one <- data.frame(Group = factor(one_level, levels = levels_T))
  sf_one <- survfit(fit_cr, newdata = nd_one)
  ss_one <- survminer::surv_summary(sf_one, data = nd_one)
  ss_one$Group <- pretty_labels[[one_level]]
  ss_one
}
df_mod <- dplyr::bind_rows(lapply(levels_T, build_sf_row))

p_cox <- ggplot(df_mod, aes(x = time, y = surv, color = Group)) +
  geom_step() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Group),
              alpha = 0.20, colour = NA) +
  labs(x = "Weeks", y = "Survival probability",
       color = "Group", fill = "Group",
       title = "Cox-predicted survival (cluster-robust fit)") +
  theme_minimal(base_size = 10)

# check KM vs Cox by treatment
df_km <- survminer::surv_summary(sf_km, data = per_snail) %>%
  mutate(Group = sub("^Group=", "", as.character(strata))) %>%
  transmute(time, surv, lower, upper, Group)

cox_df <- df_mod %>% transmute(time, surv, lower, upper, Group)

km_df <- df_km; km_df$Source  <- "KM"
cx_df <- cox_df; cx_df$Source <- "Cox"
overlay_df <- bind_rows(km_df, cx_df)

p_overlay <- ggplot(overlay_df, aes(time, surv, color = Source, linetype = Source)) +
  geom_step() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Source),
              alpha = 0.12, colour = NA) +
  scale_color_manual(values = c("KM" = "black", "Cox" = "steelblue4")) +
  scale_fill_manual( values = c("KM" = "grey70", "Cox" = "steelblue2")) +
  scale_linetype_manual(values = c("KM" = "solid", "Cox" = "dashed")) +
  facet_wrap(~ Group, ncol = 2) +
  labs(x = "Weeks", y = "Survival probability", color = "", fill = "", linetype = "",
       title = "KM vs Cox (diagnostic only)") +
  theme_minimal(base_size = 10)

# stratified Cox (matches KM so just diagnostic)
sf_strat <- survfit(Surv(week_survival, survival_status) ~ strata(Group), data = per_snail)
df_strat <- survminer::surv_summary(sf_strat, data = per_snail) %>%
  mutate(Group = sub("^strata\\(Group\\)=?", "", strata))

p_strat <- ggplot() +
  geom_step(data = df_km,   aes(time, surv, color = Group), alpha = 0.7) +
  geom_step(data = df_strat, aes(time, surv, color = Group), linetype = "dotted") +
  labs(x = "Weeks", y = "Survival probability", color = "Group",
       title = "Stratified Cox (dotted) over KM (solid) — diagnostic") +
  theme_minimal(base_size = 10)

# (p_cox)
# (p_overlay) # why macro+ micro not starts at week 1? 
# (p_strat) # why macro+ micro not starts at week 1? 

####### Omitted since assumption met: RMST (restricted mean survival time) for "Group" #########
# # RMST since PH not met (primary estimand), τ = 13 weeks
# tau <- 13  # common follow-up window; adjust as needed/wanted
# 
# # pairwise comparisons for RMST
# pair_levels <- combn(levels_T, 2, simplify = FALSE)
# pair_labels <- setNames(levels_T, levels_T)
# 
# robust_extract <- function(fit) {
#   # Extract RMST estimates for each arm and the unadjusted difference row
#   arm0 <- as.data.frame(fit$RMST.arm0$result)
#   arm1 <- as.data.frame(fit$RMST.arm1$result)
#   get_est <- function(df) {
#     idx <- if (length(rownames(df))) which(grepl("RMST", rownames(df), ignore.case = TRUE))[1] else NA_integer_
#     if (is.na(idx)) idx <- 1
#     est_col <- intersect(c("Est.", "Estimate", "est", "RMST"), names(df))
#     if (!length(est_col)) est_col <- names(df)[1]
#     as.numeric(df[idx, est_col[1]])
#   }
#   est0 <- get_est(arm0); est1 <- get_est(arm1)
#   
#   ud <- as.data.frame(fit$unadjusted.result)
#   r_ud <- if (length(rownames(ud))) which(grepl("^\\s*RMST", rownames(ud), ignore.case = TRUE))[1] else NA_integer_
#   if (is.na(r_ud)) r_ud <- 1
#   
#   est_col   <- intersect(c("Est.", "Estimate", "est", "Diff", "Difference"), names(ud)); if (!length(est_col)) est_col <- names(ud)[1]
#   lower_col <- names(ud)[grep("lower|lcl", names(ud), ignore.case = TRUE)][1]; if (is.na(lower_col)) lower_col <- names(ud)[max(2, ncol(ud)-2)]
#   upper_col <- names(ud)[grep("upper|ucl", names(ud), ignore.case = TRUE)][1]; if (is.na(upper_col)) upper_col <- names(ud)[max(3, ncol(ud)-1)]
#   p_col     <- names(ud)[grep("^p($|\\b)|p.value", names(ud), ignore.case = TRUE)][1]; if (is.na(p_col))     p_col     <- names(ud)[ncol(ud)]
#   
#   c(
#     RMST_ref   = est0,
#     RMST_comp  = est1,
#     dRMST      = as.numeric(ud[r_ud, est_col[1]]),
#     dRMST_low  = as.numeric(ud[r_ud, lower_col]),
#     dRMST_high = as.numeric(ud[r_ud, upper_col]),
#     p_value    = as.numeric(ud[r_ud, p_col])
#   )
# }
# 
# rmst_pairwise <- purrr::map_dfr(pair_levels, function(lv) {
#   lv1 <- lv[1]; lv2 <- lv[2]
#   subdat <- per_snail %>%
#     filter(Group %in% c(lv1, lv2)) %>%
#     mutate(grp = factor(Group, levels = c(lv1, lv2)))  # lv1 = 0, lv2 = 1
#   gnum <- as.integer(subdat$grp == lv2)
#   
#   fit <- survRM2::rmst2(time = subdat$week_survival, status = subdat$survival_status, arm = gnum, tau = tau)
#   vals <- robust_extract(fit)
#   
#   tibble::tibble(
#     comparison = sprintf("%s vs %s", pair_labels[[lv2]], pair_labels[[lv1]]),
#     RMST_ref   = vals["RMST_ref"],
#     RMST_comp  = vals["RMST_comp"],
#     dRMST      = vals["dRMST"],
#     dRMST_low  = vals["dRMST_low"],
#     dRMST_high = vals["dRMST_high"],
#     p_value    = vals["p_value"]
#   )
# })
# 
# rmst_pretty <- rmst_pairwise %>%
#   mutate(across(c(RMST_ref, RMST_comp, dRMST, dRMST_low, dRMST_high), ~round(as.numeric(.), 2))) %>%
#   mutate(`95% CI` = sprintf("[%.2f, %.2f]", dRMST_low, dRMST_high)) %>%
#   dplyr::select(comparison, RMST_ref, RMST_comp, dRMST, `95% CI`, p_value) %>%
#   arrange(p_value) %>%
#   mutate(FDR_BH = p.adjust(p_value, method = "BH"))
# 
# # RMST change (weeks) up to tau, comparator vs reference (BH-adjusted)
# rmst_pretty

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
