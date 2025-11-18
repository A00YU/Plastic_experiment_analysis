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

##################### R1 Final: Bootstrap jars for each treatment group #####################
source("/Users/aoyu/Desktop/Snail_Data/R1_data_cleaning.R") # cleaned master data sheet called all_data
set.seed(123)

n_control <- 75 # number of jars to resample for each treatment
n_micro <- 75
n_macro <- 18
n_macro_micro <- 18
B <- 10000 # number of bootstrap samples

# Create bootstrap matrices
boot_control <- matrix(sample(all_data$Jar_num[all_data$Treatment == "Control"], size = B*n_control, replace = TRUE), nrow = B, ncol = n_control)
boot_micro <- matrix(sample(all_data$Jar_num[all_data$Treatment == "Micro"], size = B*n_micro, replace = TRUE), nrow = B, ncol = n_micro)
boot_macro <- matrix(sample(all_data$Jar_num[all_data$Treatment == "Macro"], size = B*n_macro, replace = TRUE), nrow = B, ncol = n_macro)
boot_macro_micro <- matrix(sample(all_data$Jar_num[all_data$Treatment == "Macro+Micro"], size = B*n_macro_micro, replace = TRUE), nrow = B, ncol = n_macro_micro)

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
    for (jar in jar_numbers) {
      # Filter to specific treatment
      treatment_data <- all_data[all_data$Treatment == treatment, ]
      
      # Filter treatment_data to get rows for the current jar number and all weeks (1 to 13)
      jar_data <- treatment_data[treatment_data$Jar_num == jar & treatment_data$Week_num %in% 1:13, ]
      
      # Append the jar_data to the final_data dataframe
      final_data <- rbind(final_data, jar_data)
    }
    
    # Append the final_data dataframe to the combined_data dataframe
    combined_data <- rbind(combined_data, final_data)
  }
  
  # Store the combined_data dataframe in the bootstrap_results list
  bootstrap_results[[b]] <- combined_data
  
}

############### R1 Final: calculate parameters for each bootstrap sample #########

# Initialize a list to store the parameters data for each bootstrap sample
parameters_data_list <- vector("list", length = B)

# Loop through each bootstrap sample
for (b in 1:B) {
  # Extract the data for the current bootstrap sample on survival parameter (l_x)
  all_long_base <- bootstrap_results[[b]] %>%
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
  
  # join + add jar_id, treatment, and define outcome
  per_snail <- last_seen %>%
    full_join(first_death_week, by = c("Color","Jar_num","snail_num")) %>%
    left_join(
      all_data %>% distinct(Color, Jar_num, jar_id, Treatment),
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
  
  # Reshape the data for plotting
  death_long <- death_table %>%
    pivot_longer(cols = c("Control", "Micro", "Macro", "Macro+Micro"), names_to = "Treatment", values_to = "deaths")
  
  Treatment <- c("Control", "Micro", "Macro", "Macro+Micro")
  total <- c(225, 225, 54, 54)
  
  # Create a data frame using these vectors
  total_snails <- data.frame(Treatment, total)
  
  # Merge with total snails to calculate relative deaths
  death_long <- death_long %>%
    left_join(total_snails, by = "Treatment") %>%
    group_by(Treatment) %>%
    mutate(cumulative_death = cumsum(deaths)) 
  
  death_long_r1 <- death_long %>%
    mutate(cumulative_survived = total - cumulative_death,
           l_x = cumulative_survived / total) %>% 
    select(week_of_death, Treatment, l_x) # l_x values is cumulative_relative_survival
  
  parameters_data_r1 <- death_long_r1
  
  # Extract the data for the current bootstrap sample for growth parameter (week_to_maturity)
  # all_data_growth <- bootstrap_results[[b]] %>%
  #   select(Color, Jar_num, Week_num, 
  #          starts_with("Adult_death_"), starts_with("Snail_length_"), 
  #          starts_with("snail_matching_"), starts_with("Unidentifed_"),
  #          starts_with("snail_mating_"), starts_with("snail_orientation_"), Treatment)
  # 
  # # Replace NA with 0 for non-mating snails
  # all_data_growth <- all_data_growth %>%
  #   mutate(
  #     snail_mating_1 = ifelse(is.na(snail_mating_1), 0, snail_mating_1),
  #     snail_mating_2 = ifelse(is.na(snail_mating_2), 0, snail_mating_2),
  #     snail_mating_3 = ifelse(is.na(snail_mating_3), 0, snail_mating_3),
  #     snail_orientation_1 = ifelse(is.na(snail_orientation_1), 0, snail_orientation_1),
  #     snail_orientation_2 = ifelse(is.na(snail_orientation_2), 0, snail_orientation_2),
  #     snail_orientation_3 = ifelse(is.na(snail_orientation_3), 0, snail_orientation_3)
  #   )
  # 
  # # If mating, snail length is inaccurate, thus delete
  # all_data_growth <- all_data_growth %>%
  #   mutate(
  #     Snail_length_1 = ifelse(snail_mating_1 == 1, NA, as.numeric(Snail_length_1)),
  #     Snail_length_2 = ifelse(snail_mating_2 == 1, NA, as.numeric(Snail_length_2)),
  #     Snail_length_3 = ifelse(snail_mating_3 == 1, NA, as.numeric(Snail_length_3))
  #   )
  # 
  # # If not mating but simply orientation is wrong, correct it with a correction ratio: correct/wrong = 1.030 = 1.02997
  # # Round the product into 3 decimal places
  # all_data_growth <- all_data_growth %>%
  #   mutate(
  #     Snail_length_1 = round(ifelse(snail_orientation_1 == 1 & snail_mating_1 != 1, Snail_length_1 * 1.03, Snail_length_1), 3),
  #     Snail_length_2 = round(ifelse(snail_orientation_2 == 1 & snail_mating_2 != 1, Snail_length_2 * 1.03, Snail_length_2), 3),
  #     Snail_length_3 = round(ifelse(snail_orientation_3 == 1 & snail_mating_3 != 1, Snail_length_3 * 1.03, Snail_length_3), 3)
  #   )
  # 
  # # To get the snail_length_avg for each jar at each week (including the ones with snail death)
  # snail_length_avg_df <- all_data_growth %>%
  #   mutate(
  #     Snail_length_1 = ifelse(Adult_death_3 == 1, 0, as.numeric(Snail_length_1)),
  #     Snail_length_2 = ifelse(Adult_death_2 == 1, 0, as.numeric(Snail_length_2)),
  #     Snail_length_3 = ifelse(Adult_death_1 == 1, 0, as.numeric(Snail_length_3))) %>% 
  #   arrange(Jar_num, Color, Week_num) %>%
  #   mutate(snail_length_sum = Snail_length_1 + Snail_length_2 + Snail_length_3,
  #          snail_count = rowSums(!is.na(cbind(Adult_death_1, Adult_death_2, Adult_death_3)) & cbind(Adult_death_1, Adult_death_2, Adult_death_3) != 1)) %>%
  #   group_by(Jar_num, Color, Week_num) %>%
  #   mutate(row_within_week = row_number()) %>%
  #   ungroup() %>%
  #   arrange(Jar_num, Color, row_within_week, Week_num) %>%
  #   group_by(Jar_num, Color) %>%
  #   mutate(prev_snail_length_sum = if_else(Week_num == 1, NA_real_, lag(snail_length_sum)),
  #          prev_snail_count = if_else(Week_num == 1, NA_real_, lag(snail_count)),
  #          weekly_growth = if_else(snail_count == prev_snail_count,
  #                                  (snail_length_sum - prev_snail_length_sum) / snail_count,
  #                                  NA_real_),
  #          prev_snail_length_avg = prev_snail_length_sum/prev_snail_count) %>%
  #   ungroup() 
  # 
  # all_data_growth <- all_data_growth %>%
  #   mutate(snail_length_sum = Snail_length_1 + Snail_length_2 + Snail_length_3,
  #          snail_count = (!is.na(Snail_length_1)) + (!is.na(Snail_length_2)) + (!is.na(Snail_length_3))) %>%
  #   group_by(Jar_num, Color, Week_num) %>%
  #   mutate(row_within_week = row_number()) %>%
  #   ungroup() %>%
  #   arrange(Jar_num, Color, row_within_week, Week_num) %>%
  #   group_by(Jar_num, Color) %>%
  #   mutate(prev_snail_length_sum = if_else(Week_num == 1, NA_real_, lag(snail_length_sum)),
  #          prev_snail_count = if_else(Week_num == 1, NA_real_, lag(snail_count)),
  #          weekly_growth = if_else(snail_count == prev_snail_count,
  #                                  (snail_length_sum - prev_snail_length_sum) / snail_count,
  #                                  NA_real_)) %>%
  #   ungroup()
  # 
  # # Join the prev_snail_length_avg from snail_length_avg_df to all_data_growth
  # all_data_growth <- all_data_growth %>%
  #   left_join(snail_length_avg_df %>% select(Jar_num, Color, Week_num, prev_snail_length_avg), by = c("Jar_num", "Color", "Week_num"),
  #             relationship = "many-to-many")
  # 
  # all_data_growth <- all_data_growth %>%
  #   mutate(weekly_growth = ifelse(weekly_growth < 0, 0, weekly_growth)) %>% 
  #   filter(!is.na(weekly_growth) & weekly_growth >= 0) 
  # 
  # gamm_model_jar_avgerage <- gamm(weekly_growth ~ Treatment + s(prev_snail_length_avg, bs = "cs", k = 5), correlation = corAR1(), data = all_data_growth)
  # gam_summary <- summary(gamm_model_jar_avgerage$gam)
  # 
  # # Extract the coefficients and standard errors
  # coefficients <- gam_summary$p.coeff
  # 
  # # Create a new data frame with Treatment, estimates, and standard errors
  # model_estimates_df <- data.frame(
  #   Treatment = names(coefficients),
  #   Estimate = coefficients)
  # 
  # # Mutate a new variable with the exponentiated estimates
  # model_estimates_df <- model_estimates_df %>%
  #   mutate(weeks_to_maturity = (0.8*13)/(0.8+Estimate*13)) %>% # 0.8 cm is the threshold size to reach maturity; 13 week is the theoretical growth to maturity
  #   filter(Treatment != "(Intercept)") %>%
  #   mutate(Treatment = dplyr::recode(Treatment, 
  #                                    "TreatmentMicro" = "Micro",
  #                                    "TreatmentMacro" = "Macro",
  #                                    "TreatmentMacro+Micro" = "Macro+Micro")) %>% 
  #   select(Treatment, weeks_to_maturity)
  # 
  # control_reference_estimate_df <- data.frame(
  #   Treatment = "Control",
  #   weeks_to_maturity = 13)
  # 
  # model_estimates_df <- rbind(model_estimates_df, control_reference_estimate_df)
  # 
  # # Combine parameters_data_r1 and model_estimates_df
  # parameters_data_r1 <- parameters_data_r1 %>%
  #   left_join(model_estimates_df, by = "Treatment")
  # 
  # Extract the data for the current bootstrap sample on reproduction parameter (b_x)
  all_data <- bootstrap_results[[b]]
  all_data_repro <- all_data[all_data$Week_num != 1, ]  # Remove week 1 data
  all_data_repro <- all_data_repro %>% select(Jar_num, Color, Week_num, Total_egg_mass_num, Adult_death_1, Adult_death_2, Adult_death_3, Treatment)
  all_data_repro$total_egg_mass_num <- as.numeric(all_data_repro$Total_egg_mass_num)
  # Add a new column 'snails_num' 
  all_data_repro$snails_num <- rowSums(all_data_repro[, c("Adult_death_1", "Adult_death_2", "Adult_death_3")] == 0)
  # Eliminate rows where snails_num equals 0 (all dead)
  all_data_repro <- all_data_repro[!all_data_repro$snails_num == 0, ]
  
  all_data_repro_r1 <- all_data_repro %>% 
    select(Week_num, Treatment, total_egg_mass_num, snails_num) %>%
    mutate(eggs_mass_per_snail = total_egg_mass_num / snails_num,
           egg_per_snail = round(eggs_mass_per_snail * 10)) # conservative assuming 10 eggs per egg mass 
  
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
    filter(Week_num != 13) # no data collected for week 13 offspring hatching
  
  week13_rows_r1 <- all_data_hatch_r1 %>%
    filter(Week_num == 12) %>%
    mutate(Week_num = 13)
  
  all_data_hatch_r1 <- rbind(all_data_hatch_r1, week13_rows_r1) %>% 
    mutate(Week_num = as.factor(Week_num))
  
  parameters_data_r1 <- parameters_data_r1 %>%
    left_join(all_data_hatch_r1, by = c("week_of_death" = "Week_num", "Treatment" = "Treatment"))
  
  # Store the final combined data for the current bootstrap sample
  parameters_data_list[[b]] <- parameters_data_r1 %>% rename(week = week_of_death,
                                                                    # a_x = weeks_to_maturity,
                                                                    b_x = egg_per_snail)
}

################ R1 Final: experimental observation based population growth ##########
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

end_population_bootstrap_summary_giulio <- end_population_bootstrap_summary %>% select(iteration, Treatment, end_population) %>% pivot_wider(id_cols = iteration, names_from = Treatment, values_from = end_population) %>% arrange(iteration)

write.csv(end_population_bootstrap_summary_giulio,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R1_end_population_bootstrap_summary_giulio.csv")

# does the treatment significantly impact the end_population?
# Define a function to perform permutation test
permutation_test <- function(data, treatment1, treatment2, num_permutations = 10000) {
  # Filter data for the two treatments
  data1 <- data %>% filter(Treatment == treatment1) %>% pull(end_population)
  data2 <- data %>% filter(Treatment == treatment2) %>% pull(end_population)
  
  # Calculate observed difference in means
  observed_diff <- mean(data1) - mean(data2)
  
  # Combine the data
  combined_data <- c(data1, data2)
  
  # Initialize a vector to store permutation differences
  perm_diffs <- numeric(num_permutations)
  
  # Perform permutations
  for (i in 1:num_permutations) {
    permuted_data <- sample(combined_data)
    perm_data1 <- permuted_data[1:length(data1)]
    perm_data2 <- permuted_data[(length(data1) + 1):length(combined_data)]
    perm_diffs[i] <- mean(perm_data1) - mean(perm_data2)
  }
  
  # Calculate p-value
  p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
  
  # Return results
  list(observed_diff = observed_diff, p_value = p_value)
}

# Perform pair-wise permutation tests
treatments <- unique(end_population_bootstrap_summary$Treatment)
permutation_results <- list()
num_comparisons <- length(treatments) * (length(treatments) - 1) / 2  # Number of pairwise comparisons

for (i in 1:(length(treatments) - 1)) {
  for (j in (i + 1):length(treatments)) {
    treatment1 <- treatments[i]
    treatment2 <- treatments[j]
    test_result <- permutation_test(end_population_bootstrap_summary, treatment1, treatment2)
    permutation_results[[paste(treatment1, "vs", treatment2)]] <- test_result
  }
}

# Apply Bonferroni correction and create a results table
results_table <- data.frame(
  Comparison = character(),
  Observed_Difference = numeric(),
  P_Value_Raw = numeric(),
  P_Value_Adjusted = numeric(),
  stringsAsFactors = FALSE
)

for (comparison in names(permutation_results)) {
  observed_diff <- permutation_results[[comparison]]$observed_diff
  p_value_raw <- permutation_results[[comparison]]$p_value
  p_value_adjusted <- p_value_raw*num_comparisons
  
  results_table <- rbind(results_table, data.frame(
    Comparison = comparison,
    Observed_Difference = observed_diff,
    P_Value_Raw = p_value_raw,
    P_Value_Adjusted = p_value_adjusted
  ))
}
write.csv(results_table, file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R1_end_population_permutation_test_table.csv")

summarystats_table<- end_population_bootstrap_summary %>% group_by(Treatment) %>%
  summarise(
    mean_end_population = mean(end_population),
    sd_end_population = sd(end_population)
  )

write.csv(summarystats_table, file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R1_end_population_permutation_summarystats_table.csv")

# density plot; plot the distribution of end_population for each treatment
mean_values <- data.frame(
  Treatment = c("Control", "Micro", "Macro", "Macro+Micro"),
  Mean = c(101, 85, 115, 117))

density_plot <- ggplot(end_population_bootstrap_summary, aes(x = end_population, color = Treatment, fill = Treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    x = "Cumulative Viable Offspring Produced per Snail",
    y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "#6baf78", "Micro" = "#5a9fd6", "Macro" = "#a5855f", "Macro+Micro" = "#9d7ca5"),
                    labels = c("Control", 
                               "Virgin Micro", 
                               "Biofouled Macro", 
                               "Biofouled Macro + Virgin Micro")) +
  scale_color_manual(values = c("Control" = "#6baf78", "Micro" = "#5a9fd6", "Macro" = "#a5855f", "Macro+Micro" = "#9d7ca5"),
                     labels = c("Control", 
                                "Virgin Micro", 
                                "Biofouled Macro", 
                                "Biofouled Macro + Virgin Micro")) +
  guides(
    fill = guide_legend(override.aes = list(color = c("#6baf78", "#5a9fd6", "#a5855f", "#9d7ca5"))),
    color = guide_legend(override.aes = list(fill = c("#6baf78", "#5a9fd6", "#a5855f", "#9d7ca5"), alpha = 0.5)))+
  geom_vline(data = mean_values, aes(xintercept = Mean, color = Treatment), linetype = "dashed") +
  geom_text_repel(data = mean_values, aes(x = Mean, y = 0, label = Mean, color = Treatment), angle = 0, vjust = 1.5, hjust = 1.5, nudge_y = 0.01, nudge_x = 0,
                  size = 4, # Adjust the size of the text
                  fontface = "bold", # Make the text bold
                  box.padding = 5, # Increase padding around the text
                  point.padding = 0.5, # Increase padding around the point
                  )

density_plot

# Save the plot to a PNG file
ggsave(filename = "R1_end_population_permutation_test_plot_biofouled_naming_xaxis_mean_labeled.png", plot = density_plot, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)

# histogram plot
ggplot(end_population_bootstrap_summary, aes(x = end_population, fill = Treatment)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 30) +
  labs(title = "Distribution of End Population by Treatment",
       x = "End population",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "#6baf78", "Micro" = "#5a9fd6", "Macro" = "#a5855f", "Macro+Micro" = "#9d7ca5"))


##################### R2 Final: Bootstrap jars for each treatment group #####################
source("/Users/aoyu/Desktop/Snail_Data/R2_data_cleaning.R") # cleaned master data sheet called all_data_r2
set.seed(123)

n_control <- 50 # number of jars to resample for each treatment
n_smacro_low <- 50
n_smacro_high <- 50
n_vmacro_low <- 29
n_vmacro_high <- 29
B <- 10000 # number of bootstrap samples

# Create bootstrap matrices
boot_control <- matrix(sample(all_data_r2$Jar_num[all_data_r2$Treatment == "Control"], size = B*n_control, replace = TRUE), nrow = B, ncol = n_control)
boot_smacro_low <- matrix(sample(all_data_r2$Jar_num[all_data_r2$Treatment == "Seasoned Macro - Low"], size = B*n_smacro_low, replace = TRUE), nrow = B, ncol = n_smacro_low)
boot_smacro_high <- matrix(sample(all_data_r2$Jar_num[all_data_r2$Treatment == "Seasoned Macro - High"], size = B*n_smacro_high, replace = TRUE), nrow = B, ncol = n_smacro_high)
boot_vmacro_low <- matrix(sample(all_data_r2$Jar_num[all_data_r2$Treatment == "Virgin Macro - Low"], size = B*n_vmacro_low, replace = TRUE), nrow = B, ncol = n_vmacro_low)
boot_vmacro_high <- matrix(sample(all_data_r2$Jar_num[all_data_r2$Treatment == "Virgin Macro - High"], size = B*n_vmacro_high, replace = TRUE), nrow = B, ncol = n_vmacro_high)

# Check dimensions
print(dim(boot_control)) # should be B rows and n_control columns
print(dim(boot_smacro_low)) 
print(dim(boot_smacro_high))
print(dim(boot_vmacro_low))
print(dim(boot_vmacro_high))


# Define the treatments and their corresponding bootstrap matrices
treatments <- c("Control", "Seasoned Macro - Low", "Seasoned Macro - High", "Virgin Macro - Low", "Virgin Macro - High")
bootstrap_matrices <- list("Control" = boot_control, "Seasoned Macro - Low" = boot_smacro_low, "Seasoned Macro - High" = boot_smacro_high, "Virgin Macro - Low" = boot_vmacro_low, "Virgin Macro - High" = boot_vmacro_high)

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
    for (jar in jar_numbers) {
      # Filter to specific treatment
      treatment_data <- all_data_r2[all_data_r2$Treatment == treatment, ]
      
      # Filter treatment_data to get rows for the current jar number and all weeks (1 to 13)
      jar_data <- treatment_data[treatment_data$Jar_num == jar & treatment_data$Week_num %in% 1:13, ]
      
      # Append the jar_data to the final_data dataframe
      final_data <- rbind(final_data, jar_data)
    }
    
    # Append the final_data dataframe to the combined_data dataframe
    combined_data <- rbind(combined_data, final_data)
  }
  
  # Store the combined_data dataframe in the bootstrap_results list
  bootstrap_results[[b]] <- combined_data
  
}

############### R2 Final: calculate parameters for each bootstrap sample #########

# Initialize a list to store the parameters data for each bootstrap sample
parameters_data_list <- vector("list", length = B)

# Loop through each bootstrap sample
for (b in 1:B) {
  # Extract the data for the current bootstrap sample on survival parameter (l_x)
  all_long_base <- bootstrap_results[[b]] %>%
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
  
  # Create a new column for combined groups
  per_snail$Group <- factor(with(per_snail, 
                                 ifelse(Color %in% c("Y", "O"), "Seasoned Macro",
                                        ifelse(Color %in% c("W-1", "W-3"), "Virgin Macro","Control"))))
  
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
    filter(week_of_death != "alive") %>%
    mutate(week_of_death = as.numeric(as.character(week_of_death)))
  
  # Reshape the data for plotting
  death_long <- death_table %>%
    pivot_longer(cols = c("Control", "Seasoned Macro - Low", "Seasoned Macro - High",
                          "Virgin Macro - Low", "Virgin Macro - High"), names_to = "Treatment", values_to = "deaths")
  
  Treatment <- c("Control", "Seasoned Macro - Low", "Seasoned Macro - High",
                 "Virgin Macro - Low", "Virgin Macro - High")
  total <- c(150, 150, 150, 87, 87)
  
  # Create a data frame using these vectors
  total_snails <- data.frame(Treatment, total)
  
  # Merge with total snails to calculate relative deaths
  death_long <- death_long %>%
    left_join(total_snails, by = "Treatment") %>%
    group_by(Treatment) %>%
    mutate(cumulative_death = cumsum(deaths)) 
  
  death_long_r1 <- death_long %>%
    mutate(cumulative_survived = total - cumulative_death,
           l_x = cumulative_survived / total) %>% 
    select(week_of_death, Treatment, l_x) # l_x values is cumulative_relative_survival
  
  parameters_data_r1 <- death_long_r1
  
  # Extract the data for the current bootstrap sample for growth parameter (week_to_maturity)
  # all_data_growth <- bootstrap_results[[b]] %>%
  #   select(Color, Jar_num, Week_num, 
  #          starts_with("Adult_death_"), starts_with("Snail_length_"), 
  #          starts_with("snail_matching_"), starts_with("Unidentifed_"),
  #          starts_with("snail_mating_"), starts_with("snail_orientation_"), Treatment)
  # 
  # # Replace NA with 0 for non-mating snails
  # all_data_growth <- all_data_growth %>%
  #   mutate(
  #     snail_mating_1 = ifelse(is.na(snail_mating_1), 0, snail_mating_1),
  #     snail_mating_2 = ifelse(is.na(snail_mating_2), 0, snail_mating_2),
  #     snail_mating_3 = ifelse(is.na(snail_mating_3), 0, snail_mating_3),
  #     snail_orientation_1 = ifelse(is.na(snail_orientation_1), 0, snail_orientation_1),
  #     snail_orientation_2 = ifelse(is.na(snail_orientation_2), 0, snail_orientation_2),
  #     snail_orientation_3 = ifelse(is.na(snail_orientation_3), 0, snail_orientation_3)
  #   )
  # 
  # # If mating, snail length is inaccurate, thus delete
  # all_data_growth <- all_data_growth %>%
  #   mutate(
  #     Snail_length_1 = ifelse(snail_mating_1 == 1, NA, as.numeric(Snail_length_1)),
  #     Snail_length_2 = ifelse(snail_mating_2 == 1, NA, as.numeric(Snail_length_2)),
  #     Snail_length_3 = ifelse(snail_mating_3 == 1, NA, as.numeric(Snail_length_3))
  #   )
  # 
  # # If not mating but simply orientation is wrong, correct it with a correction ratio: correct/wrong = 1.030 = 1.02997
  # # Round the product into 3 decimal places
  # all_data_growth <- all_data_growth %>%
  #   mutate(
  #     Snail_length_1 = round(ifelse(snail_orientation_1 == 1 & snail_mating_1 != 1, Snail_length_1 * 1.03, Snail_length_1), 3),
  #     Snail_length_2 = round(ifelse(snail_orientation_2 == 1 & snail_mating_2 != 1, Snail_length_2 * 1.03, Snail_length_2), 3),
  #     Snail_length_3 = round(ifelse(snail_orientation_3 == 1 & snail_mating_3 != 1, Snail_length_3 * 1.03, Snail_length_3), 3)
  #   )
  # 
  # # To get the snail_length_avg for each jar at each week (including the ones with snail death)
  # snail_length_avg_df <- all_data_growth %>%
  #   mutate(
  #     Snail_length_1 = ifelse(Adult_death_3 == 1, 0, as.numeric(Snail_length_1)),
  #     Snail_length_2 = ifelse(Adult_death_2 == 1, 0, as.numeric(Snail_length_2)),
  #     Snail_length_3 = ifelse(Adult_death_1 == 1, 0, as.numeric(Snail_length_3))) %>% 
  #   arrange(Jar_num, Color, Week_num) %>%
  #   mutate(snail_length_sum = Snail_length_1 + Snail_length_2 + Snail_length_3,
  #          snail_count = rowSums(!is.na(cbind(Adult_death_1, Adult_death_2, Adult_death_3)) & cbind(Adult_death_1, Adult_death_2, Adult_death_3) != 1)) %>%
  #   group_by(Jar_num, Color, Week_num) %>%
  #   mutate(row_within_week = row_number()) %>%
  #   ungroup() %>%
  #   arrange(Jar_num, Color, row_within_week, Week_num) %>%
  #   group_by(Jar_num, Color) %>%
  #   mutate(prev_snail_length_sum = if_else(Week_num == 1, NA_real_, lag(snail_length_sum)),
  #          prev_snail_count = if_else(Week_num == 1, NA_real_, lag(snail_count)),
  #          weekly_growth = if_else(snail_count == prev_snail_count,
  #                                  (snail_length_sum - prev_snail_length_sum) / snail_count,
  #                                  NA_real_),
  #          prev_snail_length_avg = prev_snail_length_sum/prev_snail_count) %>%
  #   ungroup() 
  # 
  # all_data_growth <- all_data_growth %>%
  #   mutate(snail_length_sum = Snail_length_1 + Snail_length_2 + Snail_length_3,
  #          snail_count = (!is.na(Snail_length_1)) + (!is.na(Snail_length_2)) + (!is.na(Snail_length_3))) %>%
  #   group_by(Jar_num, Color, Week_num) %>%
  #   mutate(row_within_week = row_number()) %>%
  #   ungroup() %>%
  #   arrange(Jar_num, Color, row_within_week, Week_num) %>%
  #   group_by(Jar_num, Color) %>%
  #   mutate(prev_snail_length_sum = if_else(Week_num == 1, NA_real_, lag(snail_length_sum)),
  #          prev_snail_count = if_else(Week_num == 1, NA_real_, lag(snail_count)),
  #          weekly_growth = if_else(snail_count == prev_snail_count,
  #                                  (snail_length_sum - prev_snail_length_sum) / snail_count,
  #                                  NA_real_)) %>%
  #   ungroup()
  # 
  # # Join the prev_snail_length_avg from snail_length_avg_df to all_data_growth
  # all_data_growth <- all_data_growth %>%
  #   left_join(snail_length_avg_df %>% select(Jar_num, Color, Week_num, prev_snail_length_avg), by = c("Jar_num", "Color", "Week_num"),
  #             relationship = "many-to-many")
  # 
  # all_data_growth <- all_data_growth %>%
  #   mutate(weekly_growth = ifelse(weekly_growth < 0, 0, weekly_growth)) %>% 
  #   filter(!is.na(weekly_growth) & weekly_growth >= 0) 
  # 
  # gamm_model_jar_avgerage <- gamm(weekly_growth ~ Treatment + s(prev_snail_length_avg, bs = "cs", k = 5), correlation = corAR1(), data = all_data_growth)
  # gam_summary <- summary(gamm_model_jar_avgerage$gam)
  # 
  # # Extract the coefficients and standard errors
  # coefficients <- gam_summary$p.coeff
  # 
  # # Create a new data frame with Treatment, estimates, and standard errors
  # model_estimates_df <- data.frame(
  #   Treatment = names(coefficients),
  #   Estimate = coefficients)
  # 
  # # Mutate a new variable with the exponentiated estimates
  # model_estimates_df <- model_estimates_df %>%
  #   mutate(weeks_to_maturity = (0.8*13)/(0.8+Estimate*13)) %>% # 0.8 cm is the threshold size to reach maturity; 13 week is the theoretical growth to maturity
  #   filter(Treatment != "(Intercept)") %>%
  #   mutate(Treatment = dplyr::recode(Treatment, 
  #                                    "TreatmentMicro" = "Micro",
  #                                    "TreatmentMacro" = "Macro",
  #                                    "TreatmentMacro+Micro" = "Macro+Micro")) %>% 
  #   select(Treatment, weeks_to_maturity)
  # 
  # control_reference_estimate_df <- data.frame(
  #   Treatment = "Control",
  #   weeks_to_maturity = 13)
  # 
  # model_estimates_df <- rbind(model_estimates_df, control_reference_estimate_df)
  # 
  # # Combine parameters_data_r1 and model_estimates_df
  # parameters_data_r1 <- parameters_data_r1 %>%
  #   left_join(model_estimates_df, by = "Treatment")
  # 
  # Extract the data for the current bootstrap sample on reproduction parameter (b_x)
  all_data <- bootstrap_results[[b]]
  all_data_repro <- all_data[all_data$Week_num != 1, ]  # Remove week 1 data
  all_data_repro <- all_data_repro %>% select(Jar_num, Color, Week_num, Total_egg_mass_num, Adult_death_1, Adult_death_2, Adult_death_3, Treatment)
  all_data_repro$total_egg_mass_num <- as.numeric(all_data_repro$Total_egg_mass_num)
  
  # Combine Y&O into SM and W-1&W-3 into VM
  all_data_repro$Group <- factor(with(all_data_repro, 
                                      ifelse(Color %in% c("Y", "O"), "Seasoned Macro",
                                             ifelse(Color %in% c("W-1", "W-3"), "Virgin Macro","Control"))))
  
  # Add a new column 'snails_num' 
  all_data_repro$snails_num <- rowSums(all_data_repro[, c("Adult_death_1", "Adult_death_2", "Adult_death_3")] == 0)
  
  # Eliminate rows where snails_num equals 0 (all dead)
  all_data_repro <- all_data_repro[!all_data_repro$snails_num == 0, ]
  
  all_data_repro_r1 <- all_data_repro %>% 
    select(Week_num, Treatment, total_egg_mass_num, snails_num) %>%
    mutate(eggs_mass_per_snail = total_egg_mass_num / snails_num,
           egg_per_snail = round(eggs_mass_per_snail * 10)) # conservative assuming 10 eggs per egg mass 
  
  all_data_repro_r1 <- all_data_repro_r1 %>%
    group_by(Week_num, Treatment) %>%
    summarise(egg_per_snail = round(mean(egg_per_snail, na.rm = TRUE)), .groups = "drop") 
  
  parameters_data_r1 <- parameters_data_r1 %>%
    left_join(all_data_repro_r1, by = c("week_of_death" = "Week_num", "Treatment" = "Treatment")) 
  
  # Extract the data for the current bootstrap sample on hatch success parameter (h_x)
  all_data_hatch <- bootstrap_results[[1]]
  
  all_data_hatch <- all_data_hatch %>% 
    select(Jar_num, Color, Week_num, Total_egg_mass_num, Enrolled_egg_num, Unhatched_egg_num, Hatch_check_1, Hatch_check_2, Treatment,Adult_death_1, Adult_death_2, Adult_death_3)%>%
    filter(Week_num != 1) %>% # Remove week 1 given no eggs yet 
    mutate(Enrolled_egg_num = as.numeric(Enrolled_egg_num), Unhatched_egg_num = as.numeric(Unhatched_egg_num),
           Hatch_check_1 = as.numeric(Hatch_check_1), Hatch_check_2 = as.numeric(Hatch_check_2))
  
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
  
  # Combine Y&O into SM and W-1&W-3 into VM
  all_data_hatch$Group <- factor(with(all_data_hatch, 
                                      ifelse(Color %in% c("Y", "O"), "Seasoned Macro",
                                             ifelse(Color %in% c("W-1", "W-3"), "Virgin Macro","Control"))))
  
  all_data_hatch_r1 <- all_data_hatch %>%
    select(Week_num, Treatment, hatch_success) %>%
    group_by(Week_num, Treatment) %>%
    summarise(h_x = mean(hatch_success, na.rm = TRUE), .groups = "drop")
  
  parameters_data_r1 <- parameters_data_r1 %>%
    left_join(all_data_hatch_r1, by = c("week_of_death" = "Week_num", "Treatment" = "Treatment"))
  
  # Store the final combined data for the current bootstrap sample
  parameters_data_list[[b]] <- parameters_data_r1 %>% rename(week = week_of_death,
                                                             # a_x = weeks_to_maturity,
                                                             b_x = egg_per_snail)
}

################ R2 Final: experimental observation based population growth ##########
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

# does the treatment significantly impact the end_population?
# Define a function to perform permutation test
permutation_test <- function(data, treatment1, treatment2, num_permutations = 10000) {
  # Filter data for the two treatments
  data1 <- data %>% filter(Treatment == treatment1) %>% pull(end_population)
  data2 <- data %>% filter(Treatment == treatment2) %>% pull(end_population)
  
  # Calculate observed difference in means
  observed_diff <- mean(data1) - mean(data2)
  
  # Combine the data
  combined_data <- c(data1, data2)
  
  # Initialize a vector to store permutation differences
  perm_diffs <- numeric(num_permutations)
  
  # Perform permutations
  for (i in 1:num_permutations) {
    permuted_data <- sample(combined_data)
    perm_data1 <- permuted_data[1:length(data1)]
    perm_data2 <- permuted_data[(length(data1) + 1):length(combined_data)]
    perm_diffs[i] <- mean(perm_data1) - mean(perm_data2)
  }
  
  # Calculate p-value
  p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
  
  # Return results
  list(observed_diff = observed_diff, p_value = p_value)
}

# Perform pair-wise permutation tests
treatments <- unique(end_population_bootstrap_summary$Treatment)
permutation_results <- list()
num_comparisons <- length(treatments) * (length(treatments) - 1) / 2  # Number of pairwise comparisons

for (i in 1:(length(treatments) - 1)) {
  for (j in (i + 1):length(treatments)) {
    treatment1 <- treatments[i]
    treatment2 <- treatments[j]
    test_result <- permutation_test(end_population_bootstrap_summary, treatment1, treatment2)
    permutation_results[[paste(treatment1, "vs", treatment2)]] <- test_result
  }
}

# Apply Bonferroni correction and create a results table
results_table <- data.frame(
  Comparison = character(),
  Observed_Difference = numeric(),
  P_Value_Raw = numeric(),
  P_Value_Adjusted = numeric(),
  stringsAsFactors = FALSE
)

for (comparison in names(permutation_results)) {
  observed_diff <- permutation_results[[comparison]]$observed_diff
  p_value_raw <- permutation_results[[comparison]]$p_value
  p_value_adjusted <- p_value_raw*num_comparisons
  
  results_table <- rbind(results_table, data.frame(
    Comparison = comparison,
    Observed_Difference = observed_diff,
    P_Value_Raw = p_value_raw,
    P_Value_Adjusted = p_value_adjusted
  ))
}
write.csv(results_table, file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_end_population_permutation_test_results.csv")

summarystats_table<- end_population_bootstrap_summary %>% group_by(Treatment) %>%
  summarise(
    mean_end_population = mean(end_population),
    sd_end_population = sd(end_population)
  )

write.csv(summarystats_table, file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R2_end_population_permutation_summarystats_table.csv")


# density plot; plot the distribution of end_population for each treatment
mean_values <- data.frame(
  Treatment = c("Control", "Seasoned Macro - Low", "Seasoned Macro - High", "Virgin Macro - Low", "Virgin Macro - High"),
  Mean = c(164, 271, 324, 198, 198))

density_plot <- ggplot(end_population_bootstrap_summary, aes(x = end_population, color = Treatment, fill = Treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    x = "Cumulative Viable Offspring Produced per Snail",
    y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "#6baf78", "Seasoned Macro - Low" = "#d4c2a8", "Seasoned Macro - High" = "#a5855f", "Virgin Macro - Low" = "#bcbcbc", "Virgin Macro - High" = "#8a8a8a"),
                    labels = c("Control", 
                               "Biofouled Macro - Low", 
                               "Biofouled Macro - High", 
                               "Virgin Macro - Low",
                               "Virgin Macro - High")) +
  scale_color_manual(values = c("Control" = "#6baf78", "Seasoned Macro - Low" = "#d4c2a8", "Seasoned Macro - High" = "#a5855f", "Virgin Macro - Low" = "#bcbcbc", "Virgin Macro - High" = "#8a8a8a"),
                     labels = c("Control", 
                                "Biofouled Macro - Low", 
                                "Biofouled Macro - High", 
                                "Virgin Macro - Low",
                                "Virgin Macro - High"))+
  guides(
                                  fill = guide_legend(override.aes = list(color = c("#6baf78", "#d4c2a8", "#a5855f", "#bcbcbc","#8a8a8a" ))),
                                  color = guide_legend(override.aes = list(fill = c("#6baf78", "#d4c2a8", "#a5855f", "#bcbcbc", "#8a8a8a"), alpha = 0.5)))+
  geom_vline(data = mean_values, aes(xintercept = Mean, color = Treatment), linetype = "dashed") +
  geom_text_repel(data = mean_values, aes(x = Mean, y = 0, label = Mean, color = Treatment), angle = 0, vjust = 1.5, hjust = 1.5, nudge_y = 0.01, nudge_x = 0,
                  size = 4, # Adjust the size of the text
                  fontface = "bold", # Make the text bold
                  box.padding = 2, # Increase padding around the text
                  point.padding = 0.5, # Increase padding around the point
  )


density_plot
# Save the plot to a PNG file
ggsave(filename = "R2_end_population_permutation_test_plot_biofouled_naming_xaxis_mean_labeled.png", plot = density_plot, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)

# histogram plot
ggplot(end_population_bootstrap_summary, aes(x = end_population, fill = Treatment)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 30) +
  labs(title = "Distribution of End Population by Treatment",
       x = "End population",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "green", "Seasoned Macro - Low" = "yellow", "Seasoned Macro - High" = "orange", "Virgin Macro - Low" = "lightgray", "Virgin Macro - High" = "darkgray"))



############### scratch that works: bootstrap that not ideal df structure for next analysis ##############
# set.seed(123)
# 
# n_control <- 75 # number of jars to resample for each treatment
# n_micro <- 75
# n_macro <- 18
# n_macro_micro <- 18
# B <- 2 # number of bootstrap samples
# 
# boot_control <- matrix(sample(all_data$Jar_num[all_data$Treatment == "Control"], size = B*n_control, replace = TRUE), nrow = B, ncol = n_control)
# boot_micro <- matrix(sample(all_data$Jar_num[all_data$Treatment == "Micro"], size = B*n_micro, replace = TRUE), nrow = B, ncol = n_micro)
# boot_macro <- matrix(sample(all_data$Jar_num[all_data$Treatment == "Macro"], size = B*n_macro, replace = TRUE), nrow = B, ncol = n_macro)
# boot_macro_micro <- matrix(sample(all_data$Jar_num[all_data$Treatment == "Macro+Micro"], size = B*n_macro_micro, replace = TRUE), nrow = B, ncol = n_macro_micro)
# 
# # Check dimensions
# print(dim(boot_control)) # should be B rows and n_control columns
# print(dim(boot_micro)) 
# print(dim(boot_macro))
# print(dim(boot_macro_micro))
# 
# # Define the treatments and their corresponding bootstrap matrices
# treatments <- c("Control", "Micro", "Macro", "Macro+Micro")
# bootstrap_matrices <- list(Control = boot_control, Micro = boot_micro, Macro = boot_macro, 'Macro+Micro' = boot_macro_micro)
# 
# # Initialize a list to store the bootstrap results for each treatment
# bootstrap_results <- list()
# 
# # Outer loop to iterate over each treatment
# for (treatment in treatments) {
#   # Get the corresponding bootstrap matrix for the current treatment
#   boot_matrix <- bootstrap_matrices[[treatment]]
#   
#   # Check if boot_matrix is correctly defined and has valid dimensions
#   if (is.null(boot_matrix) || nrow(boot_matrix) == 0) {
#     stop(paste("Invalid bootstrap matrix for treatment:", treatment))
#   }
#   
#   # Initialize a list to store the results for each row of the current bootstrap matrix
#   bootstrap_list <- vector("list", length = nrow(boot_matrix))
#   
#   # Loop through each row of the current bootstrap matrix
#   for (i in 1:nrow(boot_matrix)) {
#     # Extract jar numbers from the current row of the bootstrap matrix
#     jar_numbers <- boot_matrix[i, ]
#     
#     # Initialize an empty dataframe to store the final results for the current row
#     final_data <- data.frame()
#     
#     # Loop through each jar number in the current row of the bootstrap matrix
#     for (jar in jar_numbers) {
#       # Filter to specific treatment
#       treatment_data <- all_data[all_data$Treatment == treatment, ]
#       
#       # Filter treatment_data to get rows for the current jar number and all weeks (1 to 13)
#       jar_data <- treatment_data[treatment_data$Jar_num == jar & treatment_data$Week_num %in% 1:13, ]
#       
#       # Append the jar_data to the final_data dataframe
#       final_data <- rbind(final_data, jar_data)
#     }
#     
#     # Store the final_data dataframe in the bootstrap_list
#     bootstrap_list[[i]] <- final_data
#   }
#   
#   # Store the bootstrap_list in the bootstrap_results list with the treatment name as the key
#   bootstrap_results[[treatment]] <- bootstrap_list
# }

################ scratch that works; growth: week_to maturity ################
# all_data_growth <- bootstrap_results[[1]] %>%
#   select(Color, Jar_num, Week_num, 
#          starts_with("Adult_death_"), starts_with("Snail_length_"), 
#          starts_with("snail_matching_"),starts_with("Unidentifed_"),
#          starts_with("snail_mating_"), starts_with("snail_orientation_"), Treatment)
# 
# # replace NA with 0 for non-mating snails
# all_data_growth <- all_data_growth %>%
#   mutate(
#     snail_mating_1 = ifelse(is.na(snail_mating_1), 0, snail_mating_1),
#     snail_mating_2 = ifelse(is.na(snail_mating_2), 0, snail_mating_2),
#     snail_mating_3 = ifelse(is.na(snail_mating_3), 0, snail_mating_3),
#     snail_orientation_1 = ifelse(is.na(snail_orientation_1), 0, snail_orientation_1),
#     snail_orientation_2 = ifelse(is.na(snail_orientation_2), 0, snail_orientation_2),
#     snail_orientation_3 = ifelse(is.na(snail_orientation_3), 0, snail_orientation_3)
#   )
# # if mating, snail length is inaccurate, thus delete
# all_data_growth <- all_data_growth %>%
#   mutate(
#     Snail_length_1 = ifelse(snail_mating_1 == 1, NA, as.numeric(Snail_length_1)),
#     Snail_length_2 = ifelse(snail_mating_2 == 1, NA, as.numeric(Snail_length_2)),
#     Snail_length_3 = ifelse(snail_mating_3 == 1, NA, as.numeric(Snail_length_3))
#   )
# 
# # if not mating but simply orientation is wrong, correct it with a correction ratio: correct/wrong = 1.030 = 1.02997
# # round the product into 3 decimal places
# all_data_growth <- all_data_growth %>%
#   mutate(
#     Snail_length_1 = round(ifelse(snail_orientation_1 == 1 & snail_mating_1 != 1, Snail_length_1 * 1.03, Snail_length_1), 3),
#     Snail_length_2 = round(ifelse(snail_orientation_2 == 1 & snail_mating_2 != 1, Snail_length_2 * 1.03, Snail_length_2), 3),
#     Snail_length_3 = round(ifelse(snail_orientation_3 == 1 & snail_mating_3 != 1, Snail_length_3 * 1.03, Snail_length_3), 3)
#   )
# 
# # to get the snail_length_avg for each jar at each week (including the ones with snail death)
# snail_length_avg_df <- all_data_growth %>%
#   mutate(
#     Snail_length_1 = ifelse(Adult_death_3 == 1, 0, as.numeric(Snail_length_1)),
#     Snail_length_2 = ifelse(Adult_death_2 == 1, 0, as.numeric(Snail_length_2)),
#     Snail_length_3 = ifelse(Adult_death_1 == 1, 0, as.numeric(Snail_length_3))) %>% 
#   arrange(Jar_num, Color, Week_num) %>%
#   mutate(snail_length_sum = Snail_length_1 + Snail_length_2 + Snail_length_3,
#          snail_count = rowSums(!is.na(cbind(Adult_death_1, Adult_death_2, Adult_death_3)) & cbind(Adult_death_1, Adult_death_2, Adult_death_3) != 1)) %>%
#   group_by(Jar_num, Color, Week_num) %>%
#   mutate(row_within_week = row_number()) %>%
#   ungroup() %>%
#   arrange(Jar_num, Color, row_within_week, Week_num) %>%
#   group_by(Jar_num, Color) %>%
#   mutate(prev_snail_length_sum = if_else(Week_num == 1, NA_real_, lag(snail_length_sum)),
#          prev_snail_count = if_else(Week_num == 1, NA_real_, lag(snail_count)),
#          weekly_growth = if_else(snail_count == prev_snail_count,
#                                  (snail_length_sum - prev_snail_length_sum) / snail_count,
#                                  NA_real_),
#          prev_snail_length_avg = prev_snail_length_sum/prev_snail_count) %>%
#   ungroup() 
# 
# all_data_growth <- all_data_growth %>%
#   mutate(snail_length_sum = Snail_length_1 + Snail_length_2 + Snail_length_3,
#          snail_count = (!is.na(Snail_length_1)) + (!is.na(Snail_length_2)) + (!is.na(Snail_length_3))) %>%
#   group_by(Jar_num, Color, Week_num) %>%
#   mutate(row_within_week = row_number()) %>%
#   ungroup() %>%
#   arrange(Jar_num, Color, row_within_week, Week_num) %>%
#   group_by(Jar_num, Color) %>%
#   mutate(prev_snail_length_sum = if_else(Week_num == 1, NA_real_, lag(snail_length_sum)),
#          prev_snail_count = if_else(Week_num == 1, NA_real_, lag(snail_count)),
#          weekly_growth = if_else(snail_count == prev_snail_count,
#                                  (snail_length_sum - prev_snail_length_sum) / snail_count,
#                                  NA_real_)) %>%
#   ungroup()
# 
# # join the prev_snail_length_avg from snail_length_avg_df to all_data_growth
# all_data_growth <- all_data_growth %>%
#   left_join(snail_length_avg_df %>% select(Jar_num, Color, Week_num, prev_snail_length_avg), by = c("Jar_num", "Color", "Week_num"),
#             relationship = "many-to-many")
# 
# all_data_growth <- all_data_growth %>%
#   mutate(weekly_growth = ifelse(weekly_growth < 0, 0, weekly_growth)) %>% 
#   filter(!is.na(weekly_growth) & weekly_growth >= 0) 
# 
# gamm_model_jar_avgerage <- gamm(weekly_growth ~ Treatment + s(prev_snail_length_avg, bs = "cs", k = 5), correlation = corAR1(), data = all_data_growth)
# gam_summary <- summary(gamm_model_jar_avgerage$gam)
# 
# # Extract the coefficients and standard errors
# coefficients <- gam_summary$p.coeff
# 
# # Create a new data frame with Treatment, estimates, and standard errors
# model_estimates_df <- data.frame(
#   Treatment = names(coefficients),
#   Estimate = coefficients)
# 
# # Mutate a new variable with the exponentiated estimates
# model_estimates_df <- model_estimates_df %>%
#   mutate(weeks_to_maturity = (0.8*13)/(0.8+Estimate*13)) %>% # 0.8 cm is the threshold size to reach maturity; 13 week is the theoretical growth to maturity
#   filter(Treatment != "(Intercept)") %>%
#   mutate(Treatment = dplyr::recode(Treatment, 
#                                    "TreatmentMicro" = "Micro",
#                                    "TreatmentMacro" = "Macro",
#                                    "TreatmentMacro+Micro" = "Macro+Micro")) %>% 
#   select(Treatment, weeks_to_maturity)
# 
# control_reference_estimate_df <- data.frame(
#   Treatment = "Control",
#   weeks_to_maturity = 13)
# 
# model_estimates_df <- rbind(model_estimates_df, control_reference_estimate_df) 
# 
# parameters_data_r1 <- parameters_data_r1 %>%
#   left_join(model_estimates_df, by = "Treatment")

################# scratch that works: fecundity ###############
all_data <- bootstrap_results[[1]]
all_data_repro <- all_data[all_data$Week_num != 1, ]  # Remove week 1 data
all_data_repro <- all_data_repro %>% select(Jar_num, Color, Week_num, Total_egg_mass_num, Adult_death_1, Adult_death_2, Adult_death_3, Treatment)
all_data_repro$total_egg_mass_num <- as.numeric(all_data_repro$Total_egg_mass_num)
# Add a new column 'snails_num' 
all_data_repro$snails_num <- rowSums(all_data_repro[, c("Adult_death_1", "Adult_death_2", "Adult_death_3")] == 0)
# Eliminate rows where snails_num equals 0 (all dead)
all_data_repro <- all_data_repro[!all_data_repro$snails_num == 0, ]

all_data_repro_r1 <- all_data_repro %>% 
  select(Week_num, Treatment, total_egg_mass_num, snails_num) %>%
  mutate(eggs_mass_per_snail = total_egg_mass_num / snails_num,
         egg_per_snail = round(eggs_mass_per_snail * 10)) # conservative assuming 10 eggs per egg mass 

# table(all_data_repro_r1$Week_num, all_data_repro_r1$Treatment) # occational missing jars are due to all 3 snail death

all_data_repro_r1 <- all_data_repro_r1 %>%
  group_by(Week_num, Treatment) %>%
  summarise(egg_per_snail = round(mean(egg_per_snail, na.rm = TRUE))) 

baseline_rows_r1 <- data.frame(
  Week_num = rep(1, 4),
  Treatment = c("Control", "Micro", "Macro", "Macro+Micro"),
  egg_per_snail = 0)

all_data_repro_r1 <- rbind(baseline_rows_r1, all_data_repro_r1)

parameters_data_r1 <- parameters_data_r1 %>%
  left_join(all_data_repro_r1, by = c("week_of_death" = "Week_num", "Treatment" = "Treatment")) %>%
  # rename(week = week_of_death,
  #        l_x = cumulative_relative_survival,
  #        b_x = egg_per_snail) 





################# scratch that works: offspring fitness ############
all_data_hatch <- bootstrap_results[[1]]

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
  summarise(h_x = mean(hatch_success, na.rm = TRUE)) %>% 
  filter(Week_num != 13) # no data collected for week 13 offspring hatching

week13_rows_r1 <- all_data_hatch_r1 %>%
  filter(Week_num == 12) %>%
  mutate(Week_num = 13)

all_data_hatch_r1 <- rbind(all_data_hatch_r1, week13_rows_r1)

parameters_data_r1 <- parameters_data_r1 %>%
  left_join(all_data_hatch_r1, by = c("week_of_death" = "Week_num", "Treatment" = "Treatment")) 
