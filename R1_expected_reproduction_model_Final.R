rm(list=ls(all=TRUE)) 

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
library(ggsignif)
library(mgcv)
library(boot)
library(ggrepel)

## 0.1. Bootstrap jars for each treatment group----
source("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/R1_data_cleaning.R") # cleaned master data sheet called all_data
set.seed(123)
all_data <- all_data %>% mutate(jar_id = interaction(Color, Jar_num, drop = TRUE))
original_data <- all_data

n_control <- 75 # number of jars to resample for each treatment
n_micro <- 75
n_macro <- 18
n_macro_micro <- 18
B <- 10000 # number of bootstrap samples

# Create bootstrap matrices
control_jars <- unique(as.character(all_data$jar_id[all_data$Treatment == "Control"]))
micro_jars <- unique(as.character(all_data$jar_id[all_data$Treatment == "Micro"]))
macro_jars <- unique(as.character(all_data$jar_id[all_data$Treatment == "Macro"]))
macro_micro_jars <- unique(as.character(all_data$jar_id[all_data$Treatment == "Macro+Micro"]))

n_control <- length(control_jars)
n_micro <- length(micro_jars)
n_macro <- length(macro_jars)
n_macro_micro <- length(macro_micro_jars)

boot_control <- matrix(sample(control_jars, B * n_control, replace = TRUE), nrow = B)
boot_micro <- matrix(sample(micro_jars, B * n_micro, replace = TRUE), nrow = B)
boot_macro <- matrix(sample(macro_jars, B * n_macro, replace = TRUE), nrow = B)
boot_macro_micro <- matrix(sample(macro_micro_jars, B * n_macro_micro, replace = TRUE), nrow = B)

# Check dimensions
print(dim(boot_control)) # should be B rows and n_control columns
print(dim(boot_micro)) 
print(dim(boot_macro))
print(dim(boot_macro_micro))

# Define the treatments and their corresponding bootstrap matrices
treatments <- c("Control", "Micro", "Macro", "Macro+Micro")
bootstrap_matrices <- list(Control = boot_control, Micro = boot_micro, Macro = boot_macro, 'Macro+Micro' = boot_macro_micro)

# Initialize a list to store the bootstrap results for each bootstrap sample
bootstrap_results <- vector("list", length = B)

# Loop through each bootstrap sample
for (b in 1:B) {
  # Initialize an empty dataframe to store the combined results for the current bootstrap sample
  combined_data <- data.frame()
  
  # Loop through each treatment
  for (treatment in treatments) {
    # Get the corresponding bootstrap matrix for the current treatment
    boot_matrix <- bootstrap_matrices[[treatment]]
    
    # Extract jar numbers from the current row of the bootstrap matrix
    jar_numbers <- boot_matrix[b, ]
    
    # Initialize an empty dataframe to store the final results for the current treatment and bootstrap sample
    final_data <- data.frame()
    
    # Loop through each jar number in the current row of the bootstrap matrix
    for (k in seq_along(jar_numbers)) {   jar <- jar_numbers[k]
      # Filter to specific treatment
      treatment_data <- all_data[all_data$Treatment == treatment, ]
      
      # Filter treatment_data to get rows for the current jar number and all weeks (1 to 13)
      jar_data <- treatment_data[treatment_data$jar_id == jar & treatment_data$Week_num %in% 1:13, ]
      jar_data$boot_jar_id <- paste(treatment, k, sep = "_")
      
      # Append the jar_data to the final_data dataframe
      final_data <- rbind(final_data, jar_data)
    }
    
    # Append the final_data dataframe to the combined_data dataframe
    combined_data <- rbind(combined_data, final_data)
  }
  
  # Store the combined_data dataframe in the bootstrap_results list
  bootstrap_results[[b]] <- combined_data
  
}

## 0.2. extract parameters for each bootstrap sample-------
R1_egg_per_mass <- read.csv("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/R1_egg_per_mass_estimates.csv")

# Initialize a list to store the parameters data for each bootstrap sample
parameters_data_list <- vector("list", length = B)

# Loop through each bootstrap sample
for (b in 1:B) {
  # Extract the data for the current bootstrap sample on survival parameter (l_x)
  all_long_base <- bootstrap_results[[b]] %>% select(boot_jar_id, Color, Jar_num, Week_num, starts_with("Adult_death_"), starts_with("Week_of_death_"))
  
  # per-snail last week (based on Adult_death_k being recorded that week)
  last_seen <- all_long_base %>%
    pivot_longer(starts_with("Adult_death_"),
                 names_to = "snail_num", values_to = "death_flag",
                 names_pattern = "Adult_death_(\\d+)") %>%
    group_by(boot_jar_id, Color, Jar_num, snail_num) %>%
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
    group_by(boot_jar_id, Color, Jar_num, snail_num) %>%
    summarise(
      week_of_death = suppressWarnings(min(as.numeric(week_of_death), na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(week_of_death = ifelse(is.finite(week_of_death), week_of_death, NA_real_))
  
  # join + add jar_id, treatment, and define outcome
  per_snail <- last_seen %>%
    full_join(first_death_week, by = c("boot_jar_id", "Color", "Jar_num", "snail_num")) %>%
    left_join(
      original_data %>% distinct(Color, Jar_num, jar_id, Treatment),
      by = c("Color","Jar_num")
    ) %>%
    mutate(
      survival_status = as.integer(!is.na(week_of_death)),
      week_survival   = ifelse(survival_status == 1, week_of_death, last_obs_week)
    ) %>%
    filter(is.finite(week_survival)) %>%
    mutate(
      Color    = factor(Color, levels = c("G","Y","O","W")),
      Treatment = factor(Treatment, levels = c("Control","Micro","Macro","Macro+Micro"))
    )
  
  # Deaths table by week and treatment
  death_table <- per_snail %>%
    mutate(week_of_death = ifelse(is.na(week_of_death), "alive", as.character(week_of_death))) %>% 
    group_by(week_of_death, Treatment) %>%
    summarise(deaths = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Treatment, values_from = deaths, values_fill = list(deaths = 0)) %>%
    mutate(week_of_death = factor(week_of_death, levels = c(as.character(2:13), "alive"))) %>%
    arrange(week_of_death) 
  
  # Convert week_of_deadth to numeric for plotting (excluding "alive")
  death_table <- death_table %>%
    filter(week_of_death != "alive")
  
  # Reshape the data for plotting, retain weeks with and without no deaths
  death_long <- death_table %>% mutate(week_of_death = as.character(week_of_death)) %>%
    pivot_longer(cols = all_of(treatments), names_to = "Treatment", values_to = "deaths") %>%
    complete(week_of_death = as.character(2:13), Treatment = treatments, fill = list(deaths = 0)) %>%
    arrange(Treatment, as.numeric(week_of_death))
  
  Treatment <- c("Control", "Micro", "Macro", "Macro+Micro")
  total <- c(225, 225, 54, 54)
  
  # Create a data frame using these vectors
  total_snails <- data.frame(Treatment, total)
  
  # Merge with total snails to calculate relative deaths
  death_long <- death_long %>%
    left_join(total_snails, by = "Treatment") %>%
    group_by(Treatment) %>%
    mutate(cumulative_death = cumsum(deaths)) %>% 
    filter(week_of_death != 13) # remove week 13 to match the no hatching data for week 13
  
  death_long_r1 <- death_long %>%
    mutate(cumulative_survived = total - cumulative_death,
           l_x = cumulative_survived / total) %>% 
    select(week_of_death, Treatment, l_x) # l_x values is cumulative_relative_survival
  
  parameters_data_r1 <- death_long_r1
  
  # Extract the data for the current bootstrap sample on reproduction parameter (b_x)
  all_data_repro <- bootstrap_results[[b]]
  all_data_repro <- all_data_repro[all_data_repro$Week_num != 1, ]  # Remove week 1 data 
  all_data_repro <- all_data_repro %>% select(Jar_num, Color, Week_num, Total_egg_mass_num, Adult_death_1, Adult_death_2, Adult_death_3, Treatment) %>% 
    filter(Week_num != 13) # remove week 13 to match the no hatching data for week 13
  all_data_repro$total_egg_mass_num <- as.numeric(all_data_repro$Total_egg_mass_num)
  # Add a new column 'snails_num' 
  all_data_repro$snails_num <- rowSums(all_data_repro[, c("Adult_death_1", "Adult_death_2", "Adult_death_3")] == 0)
  # Eliminate rows where snails_num equals 0 (all dead)
  all_data_repro <- all_data_repro[!all_data_repro$snails_num == 0, ]
  
  # simper assumption here assuming same 11 eggs per egg mass for all treatment
  # all_data_repro_r1 <- all_data_repro %>% 
  #   select(Week_num, Treatment, total_egg_mass_num, snails_num) %>%
  #   mutate(eggs_mass_per_snail = total_egg_mass_num / snails_num,
  #          egg_per_snail = round(eggs_mass_per_snail * 11)) # conservative assuming 11 eggs per egg mass 
  
  # instead of egg/mass assumption, we use the average of time-treatment-variant egg/mass to calculate the total viable eggs layed
  all_data_repro_r1 <- all_data_repro %>% 
    select(Week_num, Treatment, total_egg_mass_num, snails_num) %>%
    mutate(eggs_mass_per_snail = total_egg_mass_num / snails_num) %>%
    left_join(
      R1_egg_per_mass %>% 
        select(Week_num, Treatment, egg_per_mass),
      by = c("Week_num", "Treatment")
    ) %>%
    mutate(
      egg_per_snail = round(eggs_mass_per_snail * egg_per_mass)
    )

  all_data_repro_r1 <- all_data_repro_r1 %>%
    group_by(Week_num, Treatment) %>%
    summarise(egg_per_snail = round(mean(egg_per_snail, na.rm = TRUE)), .groups = "drop") %>% 
    mutate(Week_num = as.factor(Week_num))
  
  parameters_data_r1 <- parameters_data_r1 %>%
    left_join(all_data_repro_r1, by = c("week_of_death" = "Week_num", "Treatment" = "Treatment")) 
  
  # Extract the data for the current bootstrap sample on hatch success parameter (h_x)
  all_data_hatch <- bootstrap_results[[b]]
  
  all_data_hatch <- all_data_hatch %>% 
    select(Jar_num, Color, Week_num, Total_egg_mass_num, Enrolled_egg_num, Unhatched_egg_num, Hatch_check_1, Hatch_check_2, Treatment,Adult_death_1, Adult_death_2, Adult_death_3)%>%
    filter(Week_num != 1) %>% # Remove week 1 given no eggs yet 
    filter(Week_num != 13) %>% # & 13 with no hatching check
    mutate(Hatch_check_1 = as.numeric(Hatch_check_1), Hatch_check_2 = as.numeric(Hatch_check_2))
  
  #replace all NA with 0
  all_data_hatch$Enrolled_egg_num[is.na(all_data_hatch$Enrolled_egg_num)] <- 0
  all_data_hatch$Unhatched_egg_num[is.na(all_data_hatch$Unhatched_egg_num)] <- 0
  all_data_hatch$Hatch_check_1[is.na(all_data_hatch$Hatch_check_1)] <- 0
  all_data_hatch$Hatch_check_2[is.na(all_data_hatch$Hatch_check_2)] <- 0
  
  valid_0 <- all_data_hatch %>% 
    filter(Total_egg_mass_num > 0 & Enrolled_egg_num == 0)
  all_data_hatch <- all_data_hatch %>% 
    filter(Enrolled_egg_num >= Unhatched_egg_num & Enrolled_egg_num > 0) 
  all_data_hatch <- rbind(all_data_hatch, valid_0)
  
  # Perform the calculation for hatch success 
  all_data_hatch$hatch <- all_data_hatch$Enrolled_egg_num - all_data_hatch$Unhatched_egg_num + all_data_hatch$Hatch_check_1 + all_data_hatch$Hatch_check_2
  all_data_hatch$unhatch <- pmax(0, all_data_hatch$Unhatched_egg_num - all_data_hatch$Hatch_check_1 - all_data_hatch$Hatch_check_2)
  all_data_hatch$hatch_success <- all_data_hatch$hatch / (all_data_hatch$hatch + all_data_hatch$unhatch) 
  
  # Replace NaN values with 0 in the hatch_success column (due to 0/0 computation but still valid biological 0s)
  all_data_hatch <- all_data_hatch %>%
    mutate(hatch_success = ifelse(is.nan(hatch_success), 0, hatch_success))
  
  all_data_hatch_r1 <- all_data_hatch %>%
    select(Week_num, Treatment, hatch_success) %>%
    group_by(Week_num, Treatment) %>%
    summarise(h_x = mean(hatch_success, na.rm = TRUE), .groups = "drop") %>% 
    filter(Week_num != 13) %>%  # no data collected for week 13 offspring hatching
    mutate(Week_num = as.factor(Week_num))
  
  # week13_rows_r1 <- all_data_hatch_r1 %>%
  #   filter(Week_num == 12) %>%
  #   mutate(Week_num = 13)
  # 
  # all_data_hatch_r1 <- rbind(all_data_hatch_r1, week13_rows_r1) %>% 
  #   mutate(Week_num = as.factor(Week_num))
  # 
  parameters_data_r1 <- parameters_data_r1 %>%
    left_join(all_data_hatch_r1, by = c("week_of_death" = "Week_num", "Treatment" = "Treatment"))
  
  # Store the final combined data for the current bootstrap sample
  parameters_data_list[[b]] <- parameters_data_r1 %>% rename(week = week_of_death,
                                                             # a_x = weeks_to_maturity,
                                                             b_x = egg_per_snail)
}

## 0.3. experimental observation based population growth-----
# Initialize an empty data frame to store the results
end_population_bootstrap_summary <- data.frame()

# Loop through each element in parameters_data_list
for (b in 1:B) {
  # Calculate end_population for the current element
  end_population <- parameters_data_list[[b]] %>%
    group_by(Treatment) %>%
    summarise(end_population = sum(l_x * b_x * h_x), .groups = "drop")
  
  # Add a new variable to keep track of the iteration index
  end_population$iteration <- b
  
  # Combine the results using rbind
  end_population_bootstrap_summary <- rbind(end_population_bootstrap_summary, end_population)
}

summarystats_table<- end_population_bootstrap_summary %>% group_by(Treatment) %>%
  summarise(
    mean_end_population = mean(end_population),
    sd_end_population = sd(end_population)
  )

write.csv(summarystats_table, file = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/R1_end_population_permutation_summarystats_table.csv")

end_population_bootstrap_summary <- end_population_bootstrap_summary %>% select(iteration, Treatment, end_population) %>% pivot_wider(id_cols = iteration, names_from = Treatment, values_from = end_population) %>% arrange(iteration)

write.csv(end_population_bootstrap_summary,  file = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/R1_end_population_bootstrap_summary_giulio.csv") # save time for later use

## 0.4. Prepare original jar histories for randomization ----
boot_saved <- read.csv("G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/R1_end_population_bootstrap_summary_giulio.csv", check.names = FALSE)
names(boot_saved)[names(boot_saved) == "Macro.Micro"] <- "Macro+Micro"

end_population_bootstrap_summary <- boot_saved %>%
  select(iteration, Control, Micro, Macro, `Macro+Micro`) %>%
  pivot_longer(cols = -iteration, names_to = "Treatment", values_to = "end_population")

analysis_weeks <- 2:12
trts <- c("Control", "Micro", "Macro", "Macro+Micro")
block_col <- NULL # Enter the block column name if treatment assignment was blocked.
death_cols <- paste0("Week_of_death_", 1:3)
adult_cols <- paste0("Adult_death_", 1:3)
hatch_cols <- c("Enrolled_egg_num", "Unhatched_egg_num", "Hatch_check_1", "Hatch_check_2")
numeric_cols <- c("Week_num", "Total_egg_mass_num", death_cols, adult_cols, hatch_cols)
required_cols <- c("jar_id", "Treatment", numeric_cols, block_col)
stopifnot(all(required_cols %in% names(original_data)))

raw <- original_data %>% mutate(jar_id = as.character(jar_id), Treatment = as.character(Treatment)) %>%
  mutate(across(all_of(numeric_cols), ~ as.numeric(as.character(.x)))) %>% filter(Week_num %in% 1:13)
raw$.block <- if (is.null(block_col)) "all" else as.character(raw[[block_col]])
stopifnot(!anyNA(raw$jar_id), !anyNA(raw$.block), all(raw$Treatment %in% trts))
jar_table <- raw %>% distinct(jar_id, Treatment, .block) %>% arrange(jar_id)
stopifnot(!anyDuplicated(jar_table$jar_id), !anyDuplicated(raw[c("jar_id", "Week_num")]))
stopifnot(setequal(jar_table$Treatment, trts))
jar_ids <- jar_table$jar_id

# Each matrix row represents one original jar; columns represent analyzed weeks.
# Three initial snails per jar, matching your original denominator.
death_by_jar <- raw %>% group_by(jar_id) %>%
  summarise(across(all_of(death_cols), ~ suppressWarnings(min(.x, na.rm = TRUE))), .groups = "drop")
death_matrix <- as.matrix(death_by_jar[match(jar_ids, death_by_jar$jar_id), death_cols])
stopifnot(all(!is.finite(death_matrix) | death_matrix %in% 2:13))
survival_matrix <- sapply(analysis_weeks, function(w) rowSums(death_matrix > w) / 3)

# Attach original treatment/week egg-conversion estimates before randomization.
egg_rates <- R1_egg_per_mass %>% transmute(Week_num = as.numeric(as.character(Week_num)),
                                           Treatment = as.character(Treatment), egg_per_mass = as.numeric(as.character(egg_per_mass)))
stopifnot(!anyDuplicated(egg_rates[c("Week_num", "Treatment")]))
wk <- raw %>% filter(Week_num %in% analysis_weeks) %>% left_join(egg_rates, by = c("Week_num", "Treatment"))
wk$snails_num <- rowSums(wk[, adult_cols] == 0)
needs_rate <- which(wk$snails_num > 0 & !is.na(wk$Total_egg_mass_num))
stopifnot(all(is.finite(wk$egg_per_mass[needs_rate])))
wk$egg_per_snail <- ifelse(wk$snails_num > 0, round(wk$Total_egg_mass_num / wk$snails_num * wk$egg_per_mass), NA_real_)

# Preserve your existing hatching inclusion and NA-to-zero rules.
wk <- wk %>% mutate(across(all_of(hatch_cols), ~ replace_na(.x, 0)))
valid_hatch <- (wk$Enrolled_egg_num > 0 & wk$Enrolled_egg_num >= wk$Unhatched_egg_num) |
  (!is.na(wk$Total_egg_mass_num) & wk$Total_egg_mass_num > 0 & wk$Enrolled_egg_num == 0)
hatched <- wk$Enrolled_egg_num - wk$Unhatched_egg_num + wk$Hatch_check_1 + wk$Hatch_check_2
unhatched <- pmax(0, wk$Unhatched_egg_num - wk$Hatch_check_1 - wk$Hatch_check_2)
wk$hatch_success <- hatched / (hatched + unhatched)
wk$hatch_success[is.nan(wk$hatch_success)] <- 0
wk$hatch_success[!valid_hatch] <- NA_real_

egg_matrix <- hatch_matrix <- matrix(NA_real_, nrow = length(jar_ids), ncol = length(analysis_weeks))
positions <- cbind(match(wk$jar_id, jar_ids), match(wk$Week_num, analysis_weeks))
egg_matrix[positions] <- wk$egg_per_snail
hatch_matrix[positions] <- wk$hatch_success

## 1. Calculate treatment-level expected reproductive output ----
# ids identifies the original jars assigned to one group.
expected_output <- function(ids) {
  l_x <- colMeans(survival_matrix[ids, , drop = FALSE])
  b_x <- round(colMeans(egg_matrix[ids, , drop = FALSE], na.rm = TRUE))
  h_x <- colMeans(hatch_matrix[ids, , drop = FALSE], na.rm = TRUE)
  if (any(l_x > 0 & !is.finite(b_x))) stop("A group/week has surviving snails but no usable egg-production records.")
  contributes <- l_x > 0 & is.finite(b_x) & b_x > 0
  if (any(contributes & !is.finite(h_x))) stop("A group/week has egg production but no usable hatching records.")
  sum(l_x[contributes] * b_x[contributes] * h_x[contributes])
}

observed_outputs <- data.frame(Treatment = trts, expected_output = vapply(trts,
                                                                          function(t) expected_output(which(jar_table$Treatment == t)), numeric(1)))
print(observed_outputs)

## 2. Randomize assignments among original jars ----
perm_test_jars <- function(groupA, groupB, n_perm = 10000) {
  idx <- which(jar_table$Treatment %in% c(groupA, groupB))
  g <- jar_table$Treatment[idx]
  strata <- jar_table$.block[idx]
  observed <- expected_output(idx[g == groupA]) - expected_output(idx[g == groupB])
  d_perm <- numeric(n_perm)
  
  for (r in seq_len(n_perm)) {
    g_perm <- g
    for (s in unique(strata)) {
      ii <- which(strata == s)
      g_perm[ii] <- sample(g[ii], length(ii), replace = FALSE)
    }
    d_perm[r] <- expected_output(idx[g_perm == groupA]) - expected_output(idx[g_perm == groupB])
  }
  
  p_raw <- (1 + sum(abs(d_perm) >= abs(observed))) / (n_perm + 1)
  list(diff_mean = observed, p_raw = p_raw, perm_dist = d_perm)
}

## 3. Use existing jar-bootstrap estimates for difference intervals ----
# Calculate one difference per bootstrap iteration; do not resample the estimates.
stopifnot(!anyDuplicated(end_population_bootstrap_summary[c("iteration", "Treatment")]))
boot_wide <- end_population_bootstrap_summary %>% select(iteration, Treatment, end_population) %>%
  pivot_wider(names_from = Treatment, values_from = end_population) %>% arrange(iteration)
stopifnot(all(trts %in% names(boot_wide)), all(is.finite(as.matrix(boot_wide[, trts]))))

## 4. Run all six pairwise comparisons ----
set.seed(123)
pair_list <- combn(trts, 2, simplify = FALSE)
pair_results <- lapply(pair_list, function(p) {
  groupA <- p[1]
  groupB <- p[2]
  message("Permuting: ", groupA, " vs ", groupB)
  fit <- perm_test_jars(groupA, groupB, n_perm = 10000)
  boot_difference <- boot_wide[[groupA]] - boot_wide[[groupB]]
  ci <- quantile(boot_difference, c(0.025, 0.975), names = FALSE)
  data.frame(groupA = groupA, groupB = groupB, diff_mean = fit$diff_mean,
             p_raw = fit$p_raw, ci_lower = ci[1], ci_upper = ci[2])
})
res_df <- bind_rows(pair_results)

## 5. Preserve your current multiple-comparison adjustments ----
## 5. Bonferroni correction across all pairwise comparisons ----
res_df$p_Bonferroni <- p.adjust(res_df$p_raw, method = "bonferroni")

## 6. Significance based only on Bonferroni-adjusted p-values ----
alpha <- 0.05
res_df$sig_Bonferroni <- res_df$p_Bonferroni < alpha

## 7. Export results and retain your existing density-plot input ----
print(res_df)
write.csv(res_df, "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2/R1_end_population_bootstrap_summary_final.csv", row.names = FALSE)

jar_data <- end_population_bootstrap_summary %>% transmute(treatment = as.character(Treatment), input_par = end_population)
jar_data$treatment[jar_data$treatment == "Macro+Micro"] <- "Macro.Micro"
jar_data$treatment <- factor(jar_data$treatment, levels = c("Control", "Micro", "Macro", "Macro.Micro"))
## 8. Plot density plot----
mean_values <- jar_data %>% group_by(treatment) %>% summarise(Mean = mean(input_par), .groups = "drop")

density_plot <- ggplot(jar_data, aes(x = input_par, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    x = "Expected Reproductive Output per Snail",
    y = "Density") +
  theme_minimal() +
  theme(strip.text.x = element_blank(),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(color = "black", size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(color = "black", size = 11))+
  scale_fill_manual(name = "Treatment",
                    values = c("Control" = "#6baf78", "Micro" = "#5a9fd6", "Macro" = "#a5855f", "Macro.Micro" = "#9d7ca5"),
                    labels = c("Control", 
                               "Virgin Micro", 
                               "Biofouled Macro", 
                               "Biofouled Macro + Virgin Micro")) +
  scale_color_manual(name = "Treatment",
                     values = c("Control" = "#6baf78", "Micro" = "#5a9fd6", "Macro" = "#a5855f", "Macro.Micro" = "#9d7ca5"),
                     labels = c("Control", 
                                "Virgin Micro", 
                                "Biofouled Macro", 
                                "Biofouled Macro + Virgin Micro")) +
  geom_vline(data = mean_values, aes(xintercept = Mean, color = treatment), linetype = "dashed")+ 
  theme(legend.position = "top",
        legend.direction = "horizontal") 
# +
  # geom_text_repel(data = mean_values, 
  #                 aes(x = Mean, y = 0, label = Mean, color = treatment),  # this labels the mean on the fig
  #                 angle = 0, vjust = 1.5, hjust = 1.5, nudge_y = 0.01, nudge_x = 0,
  #                 size = 4, # Adjust the size of the text
  #                 fontface = "bold", # Make the text bold
  #                 box.padding = 5, # Increase padding around the text
  #                 point.padding = 0.5, # Increase padding around the point
  # )

density_plot

# Save the plot to a PNG file
ggsave(filename = "R1_end_population_permutation_test_plot_biofouled_naming_xaxis.png", plot = density_plot, bg = "white", path = "G:/My Drive/Stanford/Research Projects/Plastic & bulinus Experiment/Snail_Data/FinalOutputs ver2", width = 9, height = 6, dpi = 300)
