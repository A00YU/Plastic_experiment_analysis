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

####### load data#########
source("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/R2_data_cleaning.R") # cleaned master data sheet called all_data_r2
# source("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/R1_growth_analysis.R")
####### measurement error correction due to orientation #########
all_data_growth <- all_data_r2 %>%
  select(Color, Jar_num, Week_num, starts_with("Adult_death_"), starts_with("Snail_length_"), starts_with("snail_mating_"), starts_with("snail_orientation_"), Treatment)

# table(all_data$Color,all_data$Week_num)

# replace NA with 0 for non-mating snails
all_data_growth <- all_data_growth %>%
  mutate(
    snail_mating_1 = ifelse(is.na(snail_mating_1), 0, snail_mating_1),
    snail_mating_2 = ifelse(is.na(snail_mating_2), 0, snail_mating_2),
    snail_mating_3 = ifelse(is.na(snail_mating_3), 0, snail_mating_3),
    snail_orientation_1 = ifelse(is.na(snail_orientation_1), 0, snail_orientation_1),
    snail_orientation_2 = ifelse(is.na(snail_orientation_2), 0, snail_orientation_2),
    snail_orientation_3 = ifelse(is.na(snail_orientation_3), 0, snail_orientation_3)
  )
# if mating, snail length is inaccurate, thus delete
all_data_growth <- all_data_growth %>%
  mutate(
    Snail_length_1 = ifelse(snail_mating_1 == 1, NA, as.numeric(Snail_length_1)),
    Snail_length_2 = ifelse(snail_mating_2 == 1, NA, as.numeric(Snail_length_2)),
    Snail_length_3 = ifelse(snail_mating_3 == 1, NA, as.numeric(Snail_length_3))
  )

# count number of mating snails
num_mating <- sum(all_data_growth$snail_mating_1 == 1, na.rm = TRUE) + sum(all_data_growth$snail_mating_2 == 1, na.rm = TRUE) + sum(all_data_growth$snail_mating_3 == 1, na.rm = TRUE)
num_mating # 130 snails

# if not mating but simply orientation is wrong, correct it with a correction ratio: correct/wrong = 1.030 = 1.02997
# round the product into 3 decimal places
all_data_growth <- all_data_growth %>%
  mutate(
    Snail_length_1 = round(ifelse(snail_orientation_1 == 1 & snail_mating_1 != 1, Snail_length_1 * 1.03, Snail_length_1), 3),
    Snail_length_2 = round(ifelse(snail_orientation_2 == 1 & snail_mating_2 != 1, Snail_length_2 * 1.03, Snail_length_2), 3),
    Snail_length_3 = round(ifelse(snail_orientation_3 == 1 & snail_mating_3 != 1, Snail_length_3 * 1.03, Snail_length_3), 3)
  )

# count number of snails with initial wrong orientation (not including mating as "wrong" orientation)
num_orientation <- sum(all_data_growth$snail_orientation_1 == 1 & all_data_growth$snail_mating_1 != 1, na.rm = TRUE) + 
  sum(all_data_growth$snail_orientation_2 == 1 & all_data_growth$snail_mating_2 != 1, na.rm = TRUE) + 
  sum(all_data_growth$snail_orientation_3 == 1 & all_data_growth$snail_mating_3 != 1, na.rm = TRUE)
num_orientation # 1282 snails

# specifying the 19 jars (20 snails) + 0 jars to eliminate due to missing/duplicate picture or within treatment misplacement
eliminate_jars <- data.frame(
  Week_num = c(1, 1, 1, 1, 1, 4, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 8, 10),  
  Color = c("G", "W-3", "W-3", "W-1", "W-1", "G", "G", "W-3", "G", "O", "O", "O", "W-1", "W-1", "W-1", "Y", "Y", "Y", "0", "0", "0"),  
  Jar_num = c(49, 6, 12, 18, 22, 4, 30, 5, 37, 11, 18, 20, 15, 28, 29, 18, 32, 44, 36, 48, 238)  
)

# Filter the data frame to exclude the specified jars in the specific week, color, and jar number combinations
all_data_growth <- all_data_growth[!with(all_data_growth, paste(Week_num, Color, Jar_num) %in% paste(eliminate_jars$Week_num, eliminate_jars$Color, eliminate_jars$Jar_num)), ]

# count number of snails that have a length measurement
num_length <- sum(!is.na(all_data_growth$Snail_length_1)) + sum(!is.na(all_data_growth$Snail_length_2)) + sum(!is.na(all_data_growth$Snail_length_3))
num_length # 6431 snails have length measurement

# wrong orientation rate
num_orientation / num_length # 0.1993469

# to get the snail_length_avg for each jar at each week (including the ones with snail death)
snail_length_avg_df <- all_data_growth %>%
  mutate(
    Snail_length_1 = ifelse(Adult_death_3 == 1, 0, as.numeric(Snail_length_1)),
    Snail_length_2 = ifelse(Adult_death_2 == 1, 0, as.numeric(Snail_length_2)),
    Snail_length_3 = ifelse(Adult_death_1 == 1, 0, as.numeric(Snail_length_3))) %>% 
  arrange(Jar_num, Color, Week_num) %>%
  mutate(snail_length_sum = Snail_length_1 + Snail_length_2 + Snail_length_3,
         snail_count = rowSums(!is.na(cbind(Adult_death_1, Adult_death_2, Adult_death_3)) & cbind(Adult_death_1, Adult_death_2, Adult_death_3) != 1)) %>%
  arrange(Jar_num, Color, Week_num) %>%
  group_by(Jar_num, Color) %>%
  mutate(prev_snail_length_sum = lag(snail_length_sum),
         prev_snail_count = lag(snail_count),
         weekly_growth = if_else(snail_count == prev_snail_count,
                                 (snail_length_sum - prev_snail_length_sum) / snail_count,
                                 NA_real_),
         prev_snail_length_avg = prev_snail_length_sum/prev_snail_count) %>%
  ungroup() 

####### treatment growth: weekly jar level average length if no snails death that week #########
all_data_growth <- all_data_growth %>%
  mutate(snail_length_sum = Snail_length_1 + Snail_length_2 + Snail_length_3,
         snail_count = rowSums(!is.na(cbind(Adult_death_1, Adult_death_2, Adult_death_3)) & cbind(Adult_death_1, Adult_death_2, Adult_death_3) != 1)) %>%
  arrange(Jar_num, Color, Week_num) %>%
  group_by(Jar_num, Color) %>%
  mutate(prev_snail_length_sum = lag(snail_length_sum),
         prev_snail_count = lag(snail_count),
         weekly_growth = if_else(snail_count == prev_snail_count,
                                 (snail_length_sum - prev_snail_length_sum) / snail_count,
                                 NA_real_)) %>%
  ungroup()

# join the prev_snail_length_avg from snail_length_avg_df to all_data_growth
all_data_growth <- all_data_growth %>%
  left_join(snail_length_avg_df %>% select(Jar_num, Color, Week_num, prev_snail_length_avg), by = c("Jar_num", "Color", "Week_num"))

# Create a new column for combined groups
all_data_growth$Group <- factor(with(all_data_growth, 
                                     ifelse(Color %in% c("Y", "O"), "Seasoned Macro",
                                            ifelse(Color %in% c("W-1", "W-3"), "Virgin Macro","Control"))))

# Verify the new grouping
table(all_data_growth$Group, all_data_growth$Week_num)

# data counts for each senario
na_growth_jar_death <- count(all_data_growth %>%
                               filter(Adult_death_1 == 1 & Adult_death_2 == 1 & Adult_death_3 == 1)) # 131 dead jars (valid NAs)

na_growth <- count(all_data_growth %>%
                     filter(is.na(weekly_growth))) - na_growth_jar_death  # 1192/2704 = 44% NA growth due to snail death in previous/current week, mating, or missing data

negative_growth<- all_data_growth %>%
                          filter(weekly_growth < 0)

negative_growth <-count(all_data_growth %>%
                          filter(weekly_growth < 0)) # 243/2704 = 9% negative growth

count(all_data_growth %>%
        filter(weekly_growth < 0) %>% 
        filter(Color == "O"| Color == "Y")) # 115 in season macro negative growth

count(all_data_growth %>%
        filter(weekly_growth < 0) %>% 
        filter(Color == "W-3"| Color == "W-1")) # 69 in virgin macro negative growth

count(all_data_growth %>%
        filter(weekly_growth < 0) %>% 
        filter(Color == "G")) # 59 in season macro negative growth

# hist(negative_growth$weekly_growth)
table(negative_growth$Color,negative_growth$Week_num)

# save a copy of all_data_growth before filtering NAs and negative growths for later test
all_data_growth_length_covariate_r2 <- all_data_growth

#### 1. Rebuild valid weekly growth data 

# Use corrected measurements saved BEFORE growth filtering.
length_data <- all_data_growth_length_covariate_r2 %>%
  select(
    Color, Jar_num, Week_num, Treatment,
    Snail_length_1, Snail_length_2, Snail_length_3
  )

stopifnot(
  !anyDuplicated(length_data[c("Color", "Jar_num", "Week_num")]),
  !anyDuplicated(all_data_r2[c("Color", "Jar_num", "Week_num")])
)

# Reconstruct numbers alive from the FULL death history.
# Once recorded dead, a snail remains dead in subsequent weeks.
alive_history <- all_data_r2 %>%
  select(Color, Jar_num, Week_num, starts_with("Adult_death_")) %>%
  arrange(Color, Jar_num, Week_num) %>%
  group_by(Color, Jar_num) %>%
  mutate(
    dead1 = cummax(as.integer(coalesce(Adult_death_1 == 1, FALSE))),
    dead2 = cummax(as.integer(coalesce(Adult_death_2 == 1, FALSE))),
    dead3 = cummax(as.integer(coalesce(Adult_death_3 == 1, FALSE))),
    n_alive = 3L - dead1 - dead2 - dead3
  ) %>%
  ungroup() %>%
  select(Color, Jar_num, Week_num, n_alive)

growth_intervals <- length_data %>%
  left_join(
    alive_history,
    by = c("Color", "Jar_num", "Week_num")
  ) %>%
  # Zero is not a valid shell length.
  mutate(
    across(
      starts_with("Snail_length_"),
      ~ ifelse(is.finite(.) & . > 0, ., NA_real_)
    )
  )

# Require a valid length for every living snail.
length_matrix <- as.matrix(
  growth_intervals[c(
    "Snail_length_1", "Snail_length_2", "Snail_length_3"
  )]
)

growth_intervals$n_measured <- rowSums(!is.na(length_matrix))
growth_intervals$length_sum <- rowSums(length_matrix, na.rm = TRUE)

growth_intervals <- growth_intervals %>%
  mutate(
    complete_lengths =
      !is.na(n_alive) & n_alive > 0 & n_measured == n_alive,
    mean_length = if_else(
      complete_lengths,
      length_sum / pmax(n_alive, 1L),
      NA_real_
    )
  ) %>%
  arrange(Color, Jar_num, Week_num) %>%
  group_by(Color, Jar_num) %>%
  mutate(
    previous_week = lag(Week_num),
    prev_snail_length_avg = lag(mean_length),
    n_alive_start = lag(n_alive),
    previous_complete = lag(complete_lengths),
    
    # Consecutive weeks, complete survivor measurements, and no deaths.
    valid_interval = coalesce(
      Week_num - previous_week == 1 &
        complete_lengths &
        previous_complete &
        n_alive == n_alive_start,
      FALSE
    ),
    
    growth_raw = if_else(
      valid_interval,
      mean_length - prev_snail_length_avg,
      NA_real_
    ),
    
    # Preserve your existing negative-growth cleaning decision.
    weekly_growth = pmax(growth_raw, 0)
  ) %>%
  ungroup()

growth_model_data <- growth_intervals %>%
  filter(
    valid_interval,
    !is.na(Treatment),
    is.finite(weekly_growth),
    is.finite(prev_snail_length_avg),
    is.finite(n_alive_start)
  ) %>%
  mutate(
    jar_id = interaction(Color, Jar_num, drop = TRUE),
    Treatment = droplevels(factor(Treatment))
  ) %>%
  arrange(jar_id, Week_num)

stopifnot(
  nrow(growth_model_data) > 0,
  all(growth_model_data$Week_num ==
        floor(growth_model_data$Week_num)),
  !anyDuplicated(growth_model_data[c("jar_id", "Week_num")])
)

# Check the rebuilt dataset.
growth_model_data %>%
  summarise(
    jar_week_intervals = n(),
    jars = n_distinct(jar_id),
    negative_increments_set_to_zero = sum(growth_raw < 0)
  ) %>%
  print()

print(with(growth_model_data, table(Treatment, n_alive_start)))

# Because retained intervals contain no deaths, starting and ending
# live-snail counts are identical.
stopifnot(
  all(growth_model_data$n_alive_start == growth_model_data$n_alive)
)


#### 2. Fit six candidate models on identical observations 

# Without density
formula_base <- weekly_growth ~ Treatment +
  s(prev_snail_length_avg, bs = "cs", k = 5)

# With density: assumes a linear effect per additional living snail
formula_density <- weekly_growth ~ Treatment + n_alive_start +
  s(prev_snail_length_avg, bs = "cs", k = 5)

fit_candidate <- function(add_density, structure) {
  args <- list(
    formula = if (add_density) formula_density else formula_base,
    family = gaussian(link = "identity"),
    method = "ML",
    data = growth_model_data
  )
  
  if (structure %in% c("jar", "both")) {
    args$random <- list(jar_id = ~ 1)
  }
  
  if (structure %in% c("ar1", "both")) {
    args$correlation <- corAR1(form = ~ Week_num | jar_id)
  }
  
  do.call(mgcv::gamm, args)
}

specifications <- data.frame(
  model = paste0("fit", 1:6),
  structure = rep(c("jar", "ar1", "both"), 2),
  density = rep(c(FALSE, TRUE), each = 3)
)

fits <- setNames(vector("list", 6), specifications$model)
fit_status <- rep("OK", 6)
density_varies <- n_distinct(growth_model_data$n_alive_start) > 1

for (i in seq_len(6)) {
  
  if (specifications$density[i] && !density_varies) {
    fit_status[i] <- "Skipped: snail count is constant"
    next
  }
  
  warning_messages <- character()
  
  result <- tryCatch(
    withCallingHandlers(
      fit_candidate(
        add_density = specifications$density[i],
        structure = specifications$structure[i]
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  if (inherits(result, "error")) {
    fit_status[i] <- paste("Error:", conditionMessage(result))
  } else {
    fits[[i]] <- result
    
    if (length(warning_messages) > 0) {
      fit_status[i] <- paste(
        "Review warning:",
        paste(unique(warning_messages), collapse = "; ")
      )
    }
  }
}

# Requested model names:
fit1 <- fits$fit1  # Jar random intercept
fit2 <- fits$fit2  # Within-jar AR(1)
fit3 <- fits$fit3  # Jar random intercept + within-jar AR(1)
fit4 <- fits$fit4  # Jar random intercept + density
fit5 <- fits$fit5  # Within-jar AR(1) + density
fit6 <- fits$fit6  # Jar random intercept + within-jar AR(1) + density


#### 3. Compare AIC using the lme components

aic_comparison <- specifications %>%
  mutate(
    AIC = vapply(
      fits,
      function(x) if (is.null(x)) NA_real_ else AIC(x$lme),
      numeric(1)
    ),
    status = fit_status,
    eligible = status == "OK" & is.finite(AIC)
  )

if (!any(aic_comparison$eligible)) {
  print(aic_comparison)
  stop("No warning-free models available; inspect fitting messages.")
}

best_AIC <- min(aic_comparison$AIC[aic_comparison$eligible])

aic_comparison <- aic_comparison %>%
  mutate(
    delta_AIC = if_else(eligible, AIC - best_AIC, NA_real_)
  ) %>%
  arrange(desc(eligible), AIC)

print(aic_comparison)

# Lowest-AIC model, pending diagnostic checks.
best_model_name <- aic_comparison$model[
  which(aic_comparison$eligible)[1]
]

gamm_model_jar_avgerage <- fits[[best_model_name]]

cat("Lowest-AIC candidate:", best_model_name, "\n")
summary(gamm_model_jar_avgerage$gam) # main results to report here

#export into csv
gam_model_summary <- summary(gamm_model_jar_avgerage$gam)

output_dir <- "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2"

write.csv(
  data.frame(term = rownames(gam_model_summary$p.table),
             gam_model_summary$p.table, row.names = NULL),
  file.path(output_dir, "R2_growth_model_4treatments.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(term = rownames(gam_model_summary$s.table),
             gam_model_summary$s.table, row.names = NULL),
  file.path(output_dir, "R2_growth_model_4treatments_smooth.csv"),
  row.names = FALSE
)
# Random-effect variances and standard deviations
nlme::VarCorr(gamm_model_jar_avgerage$lme)

# Confirm random-effect grouping names
names(gamm_model_jar_avgerage$lme$modelStruct$reStruct) #jar-intercept SD is x cm/week

#### 4. Diagnostics before finalizing the model 

selected_lme <- gamm_model_jar_avgerage$lme
normalized_residuals <- residuals(selected_lme, type = "normalized")

# Homoscedasticity and residual pattern
plot(
  fitted(selected_lme), normalized_residuals,
  xlab = "Fitted growth", ylab = "Normalized residuals"
)
abline(h = 0, lty = 2)

# Residual normality
qqnorm(normalized_residuals)
qqline(normalized_residuals)

# Remaining within-jar autocorrelation, using actual week differences
plot(nlme::ACF(
  selected_lme,
  form = ~ Week_num | jar_id,
  resType = "normalized"
))

### Treatments growth: pairse comparison ----
# 1. Reference shell lengths from the fitted dataset
ref_vals <- quantile(
  growth_model_data$prev_snail_length_avg,
  probs = c(0.25, 0.5, 0.75),
  na.rm = TRUE
)

# 2. Treatment × size grid, holding snail count constant
ref_grid_trt <- expand.grid(
  Treatment = levels(growth_model_data$Treatment),
  prev_snail_length_avg = unname(ref_vals),
  n_alive_start = median(growth_model_data$n_alive_start)
)

Xp_trt <- predict(
  gamm_model_jar_avgerage$gam,
  newdata = ref_grid_trt,
  type = "lpmatrix"
)

## 4. Extract model coefficients and their covariance matrix
##    (these include penalized smooth terms and are required for uncertainty propagation)
beta_trt <- coef(gamm_model_jar_avgerage$gam)
V_trt <- vcov(gamm_model_jar_avgerage$gam)

## 5. Construct an emmeans-compatible object manually
##    (bypasses qdrg() and avoids penalized-smooth incompatibilities)
emm_obj <- emmobj(
  bhat = as.numeric(Xp_trt %*% beta_trt),
  V = Xp_trt %*% V_trt %*% t(Xp_trt),
  levels = list(
    Treatment = levels(all_data_growth$Treatment),
    prev_snail_length_avg = ref_vals
  ),
  linfct = diag(nrow(Xp_trt)),
  df = df.residual(gamm_model_jar_avgerage$gam),
  misc = list(by.vars = "prev_snail_length_avg")
)

# 95% CI (already used for marg_df_1)
ci95 <- as.data.frame(confint(emm_obj, level = 0.95))

# 50% CI
ci50 <- as.data.frame(confint(emm_obj, level = 0.50))

lower95 <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(ci95))[1]
upper95 <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(ci95))[1]

lower50 <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(ci50))[1]
upper50 <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(ci50))[1]
## 6. Perform Tukey-adjusted pairwise comparisons among treatments
##    separately at each prior-week snail length
pairs(
  emm_obj,
  by = "prev_snail_length_avg",
  adjust = "tukey"
)

growth_pairwise <- as.data.frame(
  summary(
    pairs(emm_obj, by = "prev_snail_length_avg", adjust = "tukey"),
    infer = c(TRUE, TRUE)  # Include 95% CIs and p-values
  )
)

write.csv(
  growth_pairwise,
  file = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/R2_growth_pairwise_4treatments.csv",
  row.names = FALSE
)

### Grouped treatments: GAMM with Autoregressive Error Structure ----
# creat group variable
growth_model_data <- growth_model_data %>%
  mutate(
    Group = case_when(
      grepl("^Seasoned Macro\\s*-\\s*(High|Low)$",
            Treatment, ignore.case = TRUE) ~ "Biofouled Macro",
      Treatment == "Control" ~ "Control",
      grepl("^Virgin Macro\\s*-\\s*(High|Low)$",
            Treatment, ignore.case = TRUE) ~ "Virgin Macro",
      TRUE ~ NA_character_
    ),
    Group = factor(
      Group,
      levels = c("Control", "Biofouled Macro", "Virgin Macro")
    )
  )

# Verify the mapping
with(growth_model_data, table(Treatment, Group, useNA = "ifany"))

gamm_model_jar_avgerage_grouped <- gamm(weekly_growth ~ Group + + n_alive_start +
                                  s(prev_snail_length_avg, bs = "cs", k = 5),
                                random = list(jar_id = ~ 1),
                                correlation = nlme::corAR1(form = ~ Week_num | jar_id),
                                family = gaussian(link = "identity"),
                                method = "ML",
                                data = growth_model_data
)

summary(gamm_model_jar_avgerage_grouped$gam)

#export into csv
gam_model_summary <- summary(gamm_model_jar_avgerage_grouped$gam)

output_dir <- "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2"

write.csv(
  data.frame(term = rownames(gam_model_summary$p.table),
             gam_model_summary$p.table, row.names = NULL),
  file.path(output_dir, "R2_growth_model_grouped.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(term = rownames(gam_model_summary$s.table),
             gam_model_summary$s.table, row.names = NULL),
  file.path(output_dir, "R2_growth_model_grouped_smooth.csv"),
  row.names = FALSE
)

# Random-effect variances and standard deviations
nlme::VarCorr(gamm_model_jar_avgerage_grouped$lme)

# Confirm random-effect grouping names
names(gamm_model_jar_avgerage$lme$modelStruct$reStruct) #jar-intercept SD is x cm/week
### Grouped treatments:pairwise comparison ----
# 1. Reference shell lengths from the fitted dataset
ref_grid_grp <- expand.grid(
  Group = levels(growth_model_data$Group),
  prev_snail_length_avg = unname(ref_vals),
  n_alive_start = median(growth_model_data$n_alive_start)
)

Xp_grp <- predict(
  gamm_model_jar_avgerage_grouped$gam,
  newdata = ref_grid_grp,
  type = "lpmatrix"
)

## 4. Extract model coefficients and their covariance matrix
##    (these include penalized smooth terms and are required for uncertainty propagation)
beta_grp <- coef(gamm_model_jar_avgerage_grouped$gam)
V_grp <- vcov(gamm_model_jar_avgerage_grouped$gam)

## 5. Construct an emmeans-compatible object manually
##    (bypasses qdrg() and avoids penalized-smooth incompatibilities)
emm_obj <- emmobj(
  bhat = as.numeric(Xp_grp %*% beta_grp),
  V = Xp_grp %*% V_grp %*% t(Xp_grp),
  levels = list(
    Group = levels(all_data_growth$Group),
    prev_snail_length_avg = ref_vals
  ),
  linfct = diag(nrow(Xp_grp)),
  df = df.residual(gamm_model_jar_avgerage_grouped$gam),
  misc = list(by.vars = "prev_snail_length_avg")
)

# 95% CI (already used for marg_df_1)
ci95 <- as.data.frame(confint(emm_obj, level = 0.95))

# 50% CI
ci50 <- as.data.frame(confint(emm_obj, level = 0.50))

lower95 <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(ci95))[1]
upper95 <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(ci95))[1]

lower50 <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(ci50))[1]
upper50 <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(ci50))[1]
## 6. Perform Tukey-adjusted pairwise comparisons among treatments
##    separately at each prior-week snail length
pairs(
  emm_obj,
  by = "prev_snail_length_avg",
  adjust = "tukey"
)

growth_pairwise <- as.data.frame(
  summary(
    pairs(emm_obj, by = "prev_snail_length_avg", adjust = "tukey"),
    infer = c(TRUE, TRUE)  # Include 95% CIs and p-values
  )
)

write.csv(
  growth_pairwise,
  file = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/R2_growth_pairwise_grouped.csv",
  row.names = FALSE
)

### plot weekly growth ----
R2_growth_plot_grouped <- ggplot(growth_model_data, aes(x = factor(Week_num), y = weekly_growth, color = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width  = 0.15,
    alpha  = 0.4,
    size   = 1
  ) +
  labs(x = "Weeks",
       y = "Mean Shell Length Growth per Snail (cm)",
       color = "Treatment") +
  coord_cartesian(ylim = c(0, 0.1))+
  scale_color_manual(values = c("Control" = "#6baf78",
                                "Biofouled Macro" = "#a5855f",
                                "Virgin Macro" = "#8a8a8a"),
                     labels = c("Control",
                                "Biofouled Macro",
                                "Virgin Macro")) +
  theme_minimal() +
  theme(strip.text.x = element_blank(),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(color = "black", size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(color = "black", size = 11))+
  facet_wrap(~ Group) +
  theme(legend.position = "top",
        legend.direction = "horizontal")

R2_growth_plot_grouped

# Save the plot to a PNG file
ggsave(filename = "R2_growth_plot_grouped_toplegend.png", plot = R2_growth_plot_grouped, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)


R2_growth_plot_4treatment <- ggplot(growth_model_data, aes(x = factor(Week_num), y = weekly_growth, color = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width  = 0.15,
    alpha  = 0.4,
    size   = 1
  ) +
  labs(x = "Weeks",
       y = "Mean Shell Length Growth per Snail (cm)",
       color = "Treatment") +
  coord_cartesian(ylim = c(0, 0.1))+
  scale_color_manual(values = c("Control" = "#6baf78",
                                "Seasoned Macro - Low" = "#d4c2a8",
                                "Seasoned Macro - High" = "#a5855f",
                                "Virgin Macro - Low" = "#bcbcbc",
                                "Virgin Macro - High" = "#8a8a8a"),
                     labels = c("Control",
                                "Biofouled Macro - Low",
                                "Biofouled Macro - High",
                                "Virgin Macro - Low",
                                "Virgin Macro - High")) +
  theme_minimal() +
  theme(strip.text.x = element_blank(),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(color = "black", size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(color = "black", size = 10))+
  facet_wrap(~ Treatment, ncol = 2)+
  theme(legend.position = "top",
        legend.direction = "horizontal")

R2_growth_plot_4treatment
# Save the plot to a PNG file
ggsave(filename = "R2_growth_plot_4treatment_toplegened.png", plot = R2_growth_plot_4treatment, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)

### marginal effect plots ----
# “Treatment effects were visualised using marginal mean predictions from a GAMM, averaging over the empirical distribution of prior-week snail length to account for nonlinear growth dynamics and temporal autocorrelation.”
# (connecting from the pairwise comparison code above)
trt_levels <- levels(growth_model_data$Treatment)
n_trt  <- length(trt_levels)
n_size <- length(ref_vals)
stopifnot(nrow(Xp_trt) == n_trt * n_size)

pred_vec <- as.numeric(Xp_trt %*% beta_trt)

pred_mat <- matrix(pred_vec, nrow = n_trt, ncol = n_size, byrow = FALSE)

size_bins <- cut(
  all_data_growth$prev_snail_length_avg,
  breaks = quantile(
    all_data_growth$prev_snail_length_avg,
    probs = seq(0, 1, length.out = n_size + 1),
    na.rm = TRUE
  ),
  include.lowest = TRUE
)

w <- as.numeric(prop.table(table(size_bins)))
stopifnot(length(w) == n_size)

marginal_means <- as.numeric(pred_mat %*% w)
V_pred <- Xp_trt %*% V_trt %*% t(Xp_trt)
L <- matrix(0, nrow = n_trt, ncol = n_trt * n_size)
for (i in seq_len(n_trt)) {
  cols <- ((i - 1) * n_size + 1):(i * n_size)
  L[i, cols] <- w
}
V_marg <- L %*% V_pred %*% t(L)


emm_marg_r2 <- emmobj(
  bhat   = marginal_means,
  V      = V_marg,
  levels = list(Treatment = trt_levels),
  linfct = diag(n_trt),
  df     = df.residual(gamm_model_jar_avgerage$gam)
)

marg_df_r2 <- as.data.frame(emm_marg_r2)

# Compute 95% CI (already likely present)
alpha95 <- 0.05
z95 <- qnorm(1 - alpha95/2)   # ~1.96
marg_df_r2$lower.CL95 <- marg_df_r2$estimate - z95 * marg_df_r2$SE
marg_df_r2$upper.CL95 <- marg_df_r2$estimate + z95 * marg_df_r2$SE

# Compute 50% CI
alpha50 <- 0.50
z50 <- qnorm(1 - alpha50/2)   # ~0.674
marg_df_r2$lower.CL50 <- marg_df_r2$estimate - z50 * marg_df_r2$SE
marg_df_r2$upper.CL50 <- marg_df_r2$estimate + z50 * marg_df_r2$SE

# Inspect the result
marg_df_r2

plot_r2<- ggplot(marg_df_r2, aes(x = Treatment, y = estimate, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    linewidth = 1, width = 0.15
  ) +
  theme_classic() +
  labs(
    x = "Treatment",
    y = "Marginal mean weekly growth (cm)",
    color = "Treatment") + # Treatment effects marginalised over growth state
  scale_color_manual(values = c("Control" = "#6baf78",
                                "Seasoned Macro - Low" = "#d4c2a8",
                                "Seasoned Macro - High" = "#a5855f",
                                "Virgin Macro - Low" = "#bcbcbc",
                                "Virgin Macro - High" = "#8a8a8a")) +
  theme(strip.text.x = element_blank(),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(color = "black", size = 11),
        legend.position = "none",
        # axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_x_discrete(labels = c("Control",
                              "Biofouled Macro \n- Low",
                              "Biofouled Macro \n- High",
                              "Virgin Macro \n- Low",
                              "Virgin Macro \n- High"))


plot_r2

# would need to read in R file R1_growth_anlysis.R to obtain the df marg_df_r1
source("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/R1_growth_analysis.R") 
marg_df_r1r2 <- bind_rows(
  "Treatment Round 1" = marg_df_r1,
  "Treatment Round 2" = marg_df_r2,
  .id = "source"
)

#label pairwise comparison with compact letter display
letters_df <- data.frame(
  Treatment = c("Control", "Micro", "Macro", "Macro+Micro",
                "Seasoned Macro - Low", "Seasoned Macro - High",
                "Virgin Macro - Low", "Virgin Macro - High"),
  cld = c("a", "ab", "b", "ab", "a", "a*", "a", "a") # Replace with your real letters
)

# Merge letters into your main plotting dataframe
marg_df_r1r2 <- marg_df_r1r2 %>%
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

# Apply this order to your dataframe
marg_df_r1r2$Treatment <- factor(
  marg_df_r1r2$Treatment,
  levels = treatment_order
)

#plot
plot_r1r2<- ggplot(marg_df_r1r2, aes(x = Treatment, y = estimate, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL50, ymax = upper.CL50),
    linewidth = 2.5,      # thicker than 95% CI
    width = 0,        # can be slightly wider
    # alpha = 0.5          # optional: slightly transparent
  ) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    linewidth = 1, width = 0.15
  ) +
  geom_text(aes(label = cld, y = upper.CL),
            vjust = -0.6,           # Push letters above the error bar
            color = "black",        # Make letters black for readability
            fontface = "bold",
            size = 5) +
  theme_classic() +
  coord_cartesian(ylim = c(0.01, 0.03))+
  labs(
    x = NULL,
    y = "Marginal mean weekly growth (cm)",
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

plot_r1r2


# # Save the plot to a PNG file
ggsave(filename = "R1R2_growth_plot_marginalized_mean_seperated.png", plot = plot_r1r2, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)
