# ================================
# Load libraries
# ================================
library(glmmTMB)
library(dplyr)
library(ggplot2)
library(emmeans)
library(DHARMa)
library(patchwork)

# ================================
# Data prep (already cleaned upstream)
# ================================
# Source cleaned dataset
source("R1_data_cleaning.R")

all_data_repro <- all_data %>%
  filter(Week_num != 1) %>%
  select(Jar_num, Color, Week_num, Total_egg_mass_num,
         Adult_death_1, Adult_death_2, Adult_death_3, Treatment)

all_data_repro$total_egg_mass_num <- as.numeric(all_data_repro$Total_egg_mass_num)

# Calculate snail number per jar-week
all_data_repro$snails_num <- rowSums(all_data_repro[, c("Adult_death_1", "Adult_death_2", "Adult_death_3")] == 0)
all_data_repro <- all_data_repro[all_data_repro$snails_num > 0, ]  # drop dead jars

# ================================
# Fit model with snail number and random effect of jar
# ================================
fec_RE_snails <- glmmTMB(total_egg_mass_num ~ Treatment + Week_num + snails_num + (1|Jar_num),
                         family = nbinom2, data = all_data_repro)

# ================================
# (A) Raw means ± SE per week-treatment
# ================================
repro_summary <- all_data_repro %>%
  group_by(Treatment, Week_num) %>%
  summarise(
    mean_repro = mean(total_egg_mass_num, na.rm = TRUE),
    sd_repro   = sd(total_egg_mass_num, na.rm = TRUE),
    n          = n(),
    se         = sd_repro / sqrt(n),
    ci_lower   = mean_repro - 1.96 * se,
    ci_upper   = mean_repro + 1.96 * se,
    .groups = "drop"
  )

p_raw <- ggplot(repro_summary, aes(x = Week_num, y = mean_repro, color = Treatment, group = Treatment)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  scale_color_manual(values = c("Control" = "green", "Micro" = "yellow",
                                "Macro" = "orange", "Macro+Micro" = "darkgray")) +
  labs(title = "Fecundity Over Time by Treatment\n(Ao-style: raw means ± SE)",
       x = "Week", y = "Mean Egg Mass Count (per jar)") +
  theme_minimal()

# ================================
# (B) Model-based predictions (holding snail number constant)
# ================================
ref_snails <- round(mean(all_data_repro$snails_num, na.rm = TRUE))

pred_data <- expand.grid(
  Treatment  = levels(all_data_repro$Treatment),
  Week_num   = sort(unique(all_data_repro$Week_num)),
  snails_num = ref_snails,
  Jar_num    = NA
)

preds <- predict(fec_RE_snails, newdata = pred_data,
                 type = "response", se.fit = TRUE, re.form = NA)

pred_data$fit   <- preds$fit
pred_data$lower <- preds$fit - 1.96 * preds$se.fit
pred_data$upper <- preds$fit + 1.96 * preds$se.fit

p_model <- ggplot() +
  geom_jitter(data = all_data_repro,
              aes(x = Week_num, y = total_egg_mass_num, color = Treatment),
              width = 0.2, alpha = 0.2, size = 1) +
  geom_line(data = pred_data,
            aes(x = Week_num, y = fit, color = Treatment, group = Treatment),
            size = 1.2) +
  geom_ribbon(data = pred_data,
              aes(x = Week_num, ymin = lower, ymax = upper, fill = Treatment, group = Treatment),
              alpha = 0.2, linetype = 0) +
  scale_color_manual(values = c("Control" = "green", "Micro" = "yellow",
                                "Macro" = "orange", "Macro+Micro" = "darkgray")) +
  scale_fill_manual(values = c("Control" = "green", "Micro" = "yellow",
                               "Macro" = "orange", "Macro+Micro" = "darkgray")) +
  labs(title = paste0("Fecundity per Week by Treatment\n(Model-based, corrected for ", ref_snails, " snails)"),
       x = "Week", y = "Predicted Egg Mass Count (per jar)") +
  theme_minimal()

p_raw + p_model

# ================================
# Combined Plot: Raw + Model Predictions
# ================================

# Choose reference snail number (mean across jars)
ref_snails <- round(mean(all_data_repro$snails_num, na.rm = TRUE))

# Prediction grid
pred_data <- expand.grid(
  Treatment  = levels(all_data_repro$Treatment),
  Week_num   = sort(unique(all_data_repro$Week_num)),
  snails_num = ref_snails,
  Jar_num    = NA
)

# Predictions (population-level, holding snail number constant)
preds <- predict(fec_RE_snails, newdata = pred_data,
                 type = "response", se.fit = TRUE, re.form = NA)

pred_data$fit   <- preds$fit
pred_data$lower <- preds$fit - 1.96 * preds$se.fit
pred_data$upper <- preds$fit + 1.96 * preds$se.fit

# Combined plot
ggplot() +
  # Raw data distribution by treatment/week
  geom_boxplot(data = all_data_repro,
               aes(x = Week_num, y = total_egg_mass_num, group = interaction(Week_num, Treatment),
                   fill = Treatment),
               alpha = 0.3, position = position_dodge(width = 0.75), outlier.shape = NA, width = 0.6) +
  geom_jitter(data = all_data_repro,
              aes(x = Week_num, y = total_egg_mass_num, color = Treatment),
              width = 0.2, alpha = 0.3, size = 0.8) +
  
  # Model predictions
  geom_line(data = pred_data,
            aes(x = Week_num, y = fit, color = Treatment, group = Treatment),
            size = 1.2) +
  geom_ribbon(data = pred_data,
              aes(x = Week_num, ymin = lower, ymax = upper, fill = Treatment, group = Treatment),
              alpha = 0.15, linetype = 0) +
  
  # Colors
  scale_color_manual(values = c("Control" = "green", "Micro" = "yellow",
                                "Macro" = "orange", "Macro+Micro" = "darkgray")) +
  scale_fill_manual(values = c("Control" = "green", "Micro" = "yellow",
                               "Macro" = "orange", "Macro+Micro" = "darkgray")) +
  
  labs(title = paste0("Fecundity per Week by Treatment (Model vs Raw)\n",
                      "(Model adjusted for ", ref_snails, " snails per jar)"),
       x = "Week", y = "Egg Mass Count (per jar)") +
  theme_minimal()
