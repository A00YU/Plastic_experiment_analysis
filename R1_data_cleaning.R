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
# Load the master data sheet
weekly_data <- read.csv("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/Plastic Tank Experiment Data Sheet Round 1.csv", skip = 2)

# Load the new macro+micro treatment specifics data
macro_micro_treatment <- read.csv("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/Plastic Tank Experiment Data Sheet - Macro+Micro-treatment Specifics.csv", skip = 1)

# Load the macro treatment specifics data
macro_treatment <- read.csv("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/Plastic Tank Experiment Data Sheet - Macro-treatment Specifics.csv", skip = 1)

# Clean and standardize column names for all datasets
colnames(macro_micro_treatment) <- make.names(colnames(macro_micro_treatment))
colnames(macro_treatment) <- make.names(colnames(macro_treatment))
colnames(weekly_data) <- make.names(colnames(weekly_data))

# Clean and prepare treatment data
macro_treatment <- macro_treatment %>%
  rename(jar_num = Jar_num, num_staples = X.) %>%
  mutate(jar_num = as.numeric(jar_num),
         num_staples = as.numeric(num_staples),
         abundance_plastic = ifelse(Abandance == "L", 3, 1))

macro_micro_treatment <- macro_micro_treatment %>%
  rename(jar_num = Jar_num, num_staples = X.) %>%
  mutate(jar_num = as.numeric(jar_num),
         num_staples = as.numeric(num_staples),
         abundance_plastic = ifelse(Abandance == "L", 3, 1))

####### macro df with contamination (week 1-6) ######
data_O <- weekly_data %>%
  filter(Color == "O", Week_num <= 6) %>%
  left_join(macro_treatment, by = c("Jar_num" = "jar_num")) %>%
  rename(Color = Color.x)  %>%
  distinct(Week_num, Jar_num, Color, .keep_all = TRUE)

# table(data_O$Color, data_O$Jar_num)

# Helper fxn to gather the weekly "Adult_death_k" and "Week_of_death_k" into long 
all_long_base_O <- data_O %>%
  select(Color, Jar_num, Week_num, starts_with("Adult_death_"), starts_with("Week_of_death_"))

# per-snail last week (based on Adult_death_k being recorded that week)
last_seen_O <- all_long_base_O %>%
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
first_death_week_O <- all_long_base_O %>%
  pivot_longer(starts_with("Week_of_death_"),
               names_to = "snail_num", values_to = "week_of_death",
               names_pattern = "Week_of_death_(\\d+)") %>%
  group_by(Color, Jar_num, snail_num) %>%
  summarise(
    week_of_death = suppressWarnings(min(as.numeric(week_of_death), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(week_of_death = ifelse(is.finite(week_of_death), week_of_death, NA_real_))

# join + add jar_id, treatment, and define outcome
per_snail_O <- last_seen_O %>%
  full_join(first_death_week_O, by = c("Color","Jar_num","snail_num")) %>%
  left_join(
    data_O %>% distinct(Color, Jar_num, num_staples, abundance_plastic),
    by = c("Color","Jar_num")
  ) %>%
  mutate(
    survival_status = as.integer(!is.na(week_of_death)),
    week_survival   = ifelse(survival_status == 1, week_of_death, last_obs_week)
  ) %>%
  filter(is.finite(week_survival)) %>%
  mutate(
    Color    = factor(Color, levels = c("G","O","W","Y")))

# Snails per jar (should be constant 3)
(per_snail_O %>% count(Jar_num, name = "n_snails") %>% summarise(min=min(n_snails), median=median(n_snails), max=max(n_snails)))

# weeks survival
(per_snail_O %>% count(week_survival) %>% arrange(week_survival) %>% head(10))

####### check if staples have significant impact on O mortality: Yes ###########
# Cox models
# Cluster-robust Cox (jar-level clustering)
fit_cr_O <- coxph(Surv(week_survival, survival_status) ~ num_staples + abundance_plastic + cluster(Jar_num),
                  data = per_snail_O, ties = "efron")

# Fixed-effects Cox (for PH tests and baseline comparison)
fit_fix_O <- coxph(Surv(week_survival, survival_status) ~ num_staples + abundance_plastic,
                   data = per_snail_O, ties = "efron")

# Random-effects Cox (frailty) — optional model comparison
fit_re_O  <- coxme(Surv(week_survival, survival_status) ~ num_staples + abundance_plastic + (1|Jar_num),
                   data = per_snail_O)

(summary(fit_cr_O))
(summary(fit_fix_O))
(summary(fit_re_O)) # contamination directly increases hazard; staple HR = 1.80

# boundary issues make LRT conservative; 
anova(fit_fix_O, fit_re_O) # no diff.; no jar-level random effect 

####### clean micro+macro df with contamination (week 1-6) #########
data_W <- weekly_data %>%
  filter(Color == "W", Week_num <= 6) %>%
  left_join(macro_micro_treatment, by = c("Jar_num" = "jar_num")) %>%
  rename(Color = Color.x)  %>%
  distinct(Week_num, Jar_num, Color, .keep_all = TRUE)

# table(data_W$Color, data_W$Jar_num)

# Helper fxn to gather the weekly "Adult_death_k" and "Week_of_death_k" into long 
all_long_base_W <- data_W %>%
  select(Color, Jar_num, Week_num, starts_with("Adult_death_"), starts_with("Week_of_death_"))

# per-snail last week (based on Adult_death_k being recorded that week)
last_seen_W <- all_long_base_W %>%
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
first_death_week_W <- all_long_base_W %>%
  pivot_longer(starts_with("Week_of_death_"),
               names_to = "snail_num", values_to = "week_of_death",
               names_pattern = "Week_of_death_(\\d+)") %>%
  group_by(Color, Jar_num, snail_num) %>%
  summarise(
    week_of_death = suppressWarnings(min(as.numeric(week_of_death), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(week_of_death = ifelse(is.finite(week_of_death), week_of_death, NA_real_))

# join + add jar_id, treatment, and define outcome
per_snail_W <- last_seen_W %>%
  full_join(first_death_week_W, by = c("Color","Jar_num","snail_num")) %>%
  left_join(
    data_W %>% distinct(Color, Jar_num, num_staples, abundance_plastic),
    by = c("Color","Jar_num")
  ) %>%
  mutate(
    survival_status = as.integer(!is.na(week_of_death)),
    week_survival   = ifelse(survival_status == 1, week_of_death, last_obs_week)
  ) %>%
  filter(is.finite(week_survival)) %>%
  mutate(
    Color    = factor(Color, levels = c("G","O","W","Y")))

# Snails per jar (should be constant 3)
(per_snail_W %>% count(Jar_num, name = "n_snails") %>% summarise(min=min(n_snails), median=median(n_snails), max=max(n_snails)))

# weeks survival
(per_snail_W %>% count(week_survival) %>% arrange(week_survival) %>% head(10))



####### check if staples have significant impact on W mortality: Yes ###########
# Cox models
# Cluster-robust Cox (jar-level clustering)
fit_cr_W <- coxph(Surv(week_survival, survival_status) ~ num_staples + abundance_plastic + cluster(Jar_num),
                  data = per_snail_W, ties = "efron")

# Fixed-effects Cox (for PH tests and baseline comparison)
fit_fix_W <- coxph(Surv(week_survival, survival_status) ~ num_staples + abundance_plastic,
                   data = per_snail_W, ties = "efron")

# Random-effects Cox (frailty) — optional model comparison
fit_re_W  <- coxme(Surv(week_survival, survival_status) ~ num_staples + abundance_plastic + (1|Jar_num),
                   data = per_snail_W)

(summary(fit_cr_W))
(summary(fit_fix_W))
(summary(fit_re_W)) # contamination directly increases hazard, staples HR = 1.97

# boundary issues make LRT conservative; 
anova(fit_fix_W, fit_re_W) # no diff.
####### check if contamination act different in W versus O : no evidence ###########
# staples as numeric: 
per_snail_OW <- bind_rows(per_snail_W, per_snail_O) %>%
  mutate(treatment = factor(Color, levels = c("O","W"))) 

# Cox models
# Fixed-effects Cox (for PH tests and baseline comparison)
fit_fix_OW <- coxph(Surv(week_survival, survival_status) ~ treatment * num_staples + abundance_plastic,
                    data = per_snail_OW, ties = "efron")

(summary(fit_fix_OW)) # no evidence of interactive effect of treatment and staples 

# staples as binary:
per_snail_OW <- bind_rows(per_snail_W, per_snail_O) %>%
  mutate(treatment = factor(Color, levels = c("O","W"))) %>%
  mutate(staples_binary = ifelse(num_staples == 0, 0, 1))

# Cox models
# Fixed-effects Cox (for PH tests and baseline comparison)
fit_fix_OW <- coxph(Surv(week_survival, survival_status) ~ treatment * staples_binary + abundance_plastic,
                    data = per_snail_OW, ties = "efron")

(summary(fit_fix_OW)) # no evidence of interactive effect of treatment and staples 

####### master df without contamination for further analysis ########
# Filter weekly data for relevant groups and jars
filtered_data_GY <- weekly_data %>%
  filter(Color %in% c("G", "Y"), as.numeric(Jar_num) >= 1 & as.numeric(Jar_num) <= 75)

filtered_data_W <- weekly_data %>%
  filter(Color == "W") %>%
  left_join(macro_micro_treatment, by = c("Jar_num" = "jar_num")) %>%
  filter(num_staples == 0) %>%
  rename(Color = Color.x)  %>%
  distinct(Week_num, Jar_num, Color, .keep_all = TRUE)

# table(filtered_data_W$Color, filtered_data_W$Jar_num)

filtered_data_O <- weekly_data %>%
  filter(Color == "O") %>%
  left_join(macro_treatment, by = c("Jar_num" = "jar_num")) %>%
  filter(num_staples == 0) %>%
  rename(Color = Color.x)  %>%
  distinct(Week_num, Jar_num, Color, .keep_all = TRUE)

# table(filtered_data_O$Color, filtered_data_O$Jar_num)

# Combine all filtered data into one dataset
all_data <- bind_rows(filtered_data_GY, filtered_data_O, filtered_data_W) %>%
  mutate(
    Treatment = dplyr::recode(Color,
                       "G" = "Control", "Y" = "Micro", "O" = "Macro", "W" = "Macro+Micro"),
    Treatment = factor(Treatment, levels = c("Control","Micro","Macro","Macro+Micro")),
    # Unique physical jar (color x jar number)
    jar_id = interaction(Color, Jar_num, drop = TRUE)
  )

# table(all_data$Color, all_data$Jar_num)

