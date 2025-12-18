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

## -----------------------------------------------------------
## 0.1. Bootstrap jars for each treatment group
## -----------------------------------------------------------
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

## -----------------------------------------------------------
## 0.2. extract parameters for each bootstrap sample
## -----------------------------------------------------------
R1_egg_per_mass <- read.csv("/Users/aoyu/Desktop/Snail_Data/FinalOutputs/R1_egg_per_mass_estimates.csv")

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
    mutate(cumulative_death = cumsum(deaths)) %>% 
    filter(week_of_death != 13) # remove week 13 to match the no hatching data for week 13
  
  death_long_r1 <- death_long %>%
    mutate(cumulative_survived = total - cumulative_death,
           l_x = cumulative_survived / total) %>% 
    select(week_of_death, Treatment, l_x) # l_x values is cumulative_relative_survival
  
  parameters_data_r1 <- death_long_r1
  
  # Extract the data for the current bootstrap sample on reproduction parameter (b_x)
  all_data <- bootstrap_results[[b]]
  all_data_repro <- all_data[all_data$Week_num != 1, ]  # Remove week 1 data 
  all_data_repro <- all_data_repro %>% select(Jar_num, Color, Week_num, Total_egg_mass_num, Adult_death_1, Adult_death_2, Adult_death_3, Treatment) %>% 
    filter(Week_num != 13) # remove week 13 to match the no hatching data for week 13
  all_data_repro$total_egg_mass_num <- as.numeric(all_data_repro$Total_egg_mass_num)
  # Add a new column 'snails_num' 
  all_data_repro$snails_num <- rowSums(all_data_repro[, c("Adult_death_1", "Adult_death_2", "Adult_death_3")] == 0)
  # Eliminate rows where snails_num equals 0 (all dead)
  all_data_repro <- all_data_repro[!all_data_repro$snails_num == 0, ]
  
  # all_data_repro_r1 <- all_data_repro %>% 
  #   select(Week_num, Treatment, total_egg_mass_num, snails_num) %>%
  #   mutate(eggs_mass_per_snail = total_egg_mass_num / snails_num,
  #          egg_per_snail = round(eggs_mass_per_snail * 11)) # conservative assuming 11 eggs per egg mass 
  
  # instead of egg/mass assumption, we use the average of time-treatment-variant egg/mass to calculate the total viable eggs laied
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

## -----------------------------------------------------------
## 0.3. experimental observation based population growth
## -----------------------------------------------------------
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

write.csv(summarystats_table, file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R1_end_population_permutation_summarystats_table.csv")

end_population_bootstrap_summary_giulio <- end_population_bootstrap_summary %>% select(iteration, Treatment, end_population) %>% pivot_wider(id_cols = iteration, names_from = Treatment, values_from = end_population) %>% arrange(iteration)

write.csv(end_population_bootstrap_summary_giulio,  file = "/Users/aoyu/Desktop/Snail_Data/R1_end_population_bootstrap_summary_giulio.csv") # save time for later use

## -----------------------------------------------------------
## 0.4. Load population growth data
## -----------------------------------------------------------

jar_data <- read.csv("R1_end_population_bootstrap_summary_giulio.csv", header = TRUE, stringsAsFactors = FALSE)
jar_data <- jar_data %>% select(Control, Micro, Macro, Macro.Micro) %>% pivot_longer(cols = everything(),
                                                                                     names_to = "treatment",
                                                                                     values_to = "input_par") # 10,000 x4

# --- and inspect them  ---
head(jar_data)
str(jar_data)
summary(jar_data)
jar_data$treatment <- factor(
  jar_data$treatment,
  levels = c("Control", "Micro", "Macro", "Macro.Micro"))
table(jar_data$treatment)

jar_data %>%
  group_by(treatment) %>%
  summarise(mean_input_par = mean(input_par))

## -----------------------------------------------------------
## 1. Helper: difference in summary statistic between two groups
## -----------------------------------------------------------

obs_diff <- function(x, g, groupA, groupB, fun = mean) {
  xa <- x[g == groupA]
  xb <- x[g == groupB]
  fun(xa) - fun(xb)  # see comments below about this "fun()"
}

## -----------------------------------------------------------
## 2. Permutation test for two groups (case #2)
## -----------------------------------------------------------

perm_test_two_groups <- function(x, # numeric data (MOI values)
                                 g, # which treatment each datum belongs to
                                 groupA,  # e.g. "Control", but valid for any pairwise comparison
                                 groupB,  # valid for any pairwise comparison
                                 n_perm = 10000,# number of replicated
                                 fun = mean,  # VERY IMPORTANT, see comment below * 
                                 alternative = c("two.sided", "greater", "less"),
                                 seed = NULL) {
  
  # * above I set: <<fun = mean>> (as arithmetic mean) for simplicity, to be sure
  # that this the same function that is used in the function obs_diff() defined
  # above. This is actually the very complex set of computations that you
  # performed in order to to estimate the expected number of larvae that a snail
  # at the beginning of the experiment during the entire duration of the
  # experiment, i.e., sum by x of: l(x) * f(x) * h(x) - you have already done it,
  # so you know how to properly code it! And you probably have to remove <<fun =
  # mean take>> from the input parameters of this function
  
  
  alternative <- match.arg(alternative) # check that we used the right labels
  if (!is.null(seed)) set.seed(seed)
  
  # Extract the two groups of interest
  # creates a logical vector that is TRUE when g is either groupA or groupB, 
  # and FALSE otherwise, so So idx is a logical “mask” saying which jars belong 
  # to the two treatments we’re comparing.
  idx <- g %in% c(groupA, groupB)
  
  # Then use this mask to subset only the jars in groupA or groupB, 
  # ignoring other treatments
  x_sub <- x[idx] 
  
  # Then, to compelte subsetting, drop unused factor levels from g
  g_sub <- droplevels(g[idx])
  
  # Observed difference in statistic
  D_obs <- obs_diff(x_sub, g_sub, groupA, groupB, fun = fun); D_obs
  
  # Permutation distribution under H0: no treatment effect
  D_perm <- numeric(n_perm)
  for (r in seq_len(n_perm)) {
    g_perm <- sample(g_sub, replace = FALSE) # without replacement
    D_perm[r] <- obs_diff(x_sub, g_perm, groupA, groupB, fun = fun)
  }
  
  # p-value
  if (alternative == "two.sided") {
    p_val <- (sum(abs(D_perm) >= abs(D_obs)) + 1) / (n_perm + 1)
  } else if (alternative == "greater") {
    p_val <- (sum(D_perm >= D_obs) + 1) / (n_perm + 1)
  } else { # "less"
    p_val <- (sum(D_perm <= D_obs) + 1) / (n_perm + 1)
  }
  p_val
  
  list(
    groupA      = groupA,
    groupB      = groupB,
    statistic   = D_obs,
    perm_dist   = D_perm,
    p.value     = p_val,
    alternative = alternative)
}

# --- check whether it works in this mock-up example! ---
res_perm <- perm_test_two_groups(
  x        = jar_data$input_par,
  g        = jar_data$treatment,
  groupA   = "Macro.Micro",
  groupB   = "Macro",
  n_perm   = 10000,
  alternative = "two.sided",
  seed = NULL
)

res_perm$p.value    # permutation p-value
res_perm$statistic  # observed difference in mean Metric of interest (MacroMicro vs - Macroplastic)

## -----------------------------------------------------------
## 3. Bootstrap CI for difference in means between two groups
## -----------------------------------------------------------

boot_diff_ci <- function(x, g,
                         groupA,
                         groupB,
                         B = 10000,
                         fun = mean,
                         conf = 0.95,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  xa <- x[g == groupA]
  xb <- x[g == groupB]
  nA <- length(xa)
  nB <- length(xb)
  
  D_boot <- numeric(B)
  for (b in seq_len(B)) {
    xa_b <- sample(xa, size = nA, replace = TRUE)
    xb_b <- sample(xb, size = nB, replace = TRUE)
    D_boot[b] <- fun(xa_b) - fun(xb_b)
  }
  
  alpha <- 1 - conf
  ci    <- quantile(D_boot, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  
  list(
    statistic = fun(xa) - fun(xb),
    ci_lower  = ci[1],
    ci_upper  = ci[2]
  )
}

## -----------------------------------------------------------
## 4. Run all 6 pairwise comparisons + bootstrap CIs
## -----------------------------------------------------------
trts  <- levels(jar_data$treatment)
pairs <- combn(trts, 2)   # all pairs of treatments

set.seed(123)  # for reproducibility

pair_results <- lapply(seq_len(ncol(pairs)), function(k) {
  A <- pairs[1, k]
  B <- pairs[2, k]
  
  # Permutation test (two-sided)
  perm_res <- perm_test_two_groups(
    x      = jar_data$input_par,
    g      = jar_data$treatment,
    groupA = A,
    groupB = B,
    n_perm = 10000,
    alternative = "two.sided"
  )
  
  # Bootstrap CI for difference
  boot_res <- boot_diff_ci(
    x      = jar_data$input_par,
    g      = jar_data$treatment,
    groupA = A,
    groupB = B,
    B      = 10000,
    conf   = 0.95
  )
  
  data.frame(
    groupA    = perm_res$groupA,
    groupB    = perm_res$groupB,
    diff_mean = perm_res$statistic,   # A - B
    p_raw     = perm_res$p.value,
    ci_lower  = boot_res$ci_lower,
    ci_upper  = boot_res$ci_upper,
    stringsAsFactors = FALSE
  )
})

res_df <- do.call(rbind, pair_results)

## -----------------------------------------------------------
## 5. Apply multiple-comparisons corrections
## -----------------------------------------------------------
res_df$p_BH   <- p.adjust(res_df$p_raw,  method = "BH")    # Benjamini–Hochberg
res_df$p_Holm <- p.adjust(res_df$p_raw,  method = "holm")  # Holm–Bonferroni

## -----------------------------------------------------------
## 6. Mark which comparisons are significant (α = 0.05)
## -----------------------------------------------------------
alpha <- 0.05

# significance according to adjusted p-values
res_df$sig_BH   <- res_df$p_BH   < alpha
res_df$sig_Holm <- res_df$p_Holm < alpha

# optional: significance according to bootstrap CI (CI excludes 0)
res_df$sig_CI <- (res_df$ci_lower > 0) | (res_df$ci_upper < 0)

## -----------------------------------------------------------
## 7. Inspect results
## -----------------------------------------------------------
res_df

write.csv(res_df,  file = "/Users/aoyu/Desktop/Snail_Data/Outputs/R1_end_population_bootstrap_summary_final.csv")

## -----------------------------------------------------------
## 8. Plot density plot
## -----------------------------------------------------------
mean_values <- data.frame(
  treatment = c("Control", "Micro", "Macro", "Macro.Micro"),
  Mean = c(97, 89, 123, 122))

density_plot <- ggplot(jar_data, aes(x = input_par, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    x = "Expected Reproductive Output per Snail",
    y = "Density") +
  theme_minimal() +
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
  geom_vline(data = mean_values, aes(xintercept = Mean, color = treatment), linetype = "dashed") +
  geom_text_repel(data = mean_values, 
                  aes(x = Mean, y = 0, label = Mean, color = treatment), 
                  angle = 0, vjust = 1.5, hjust = 1.5, nudge_y = 0.01, nudge_x = 0,
                  size = 4, # Adjust the size of the text
                  fontface = "bold", # Make the text bold
                  box.padding = 5, # Increase padding around the text
                  point.padding = 0.5, # Increase padding around the point
  )

density_plot

# Save the plot to a PNG file
ggsave(filename = "R1_end_population_permutation_test_plot_biofouled_naming_xaxis_mean_labeled.png", plot = density_plot, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)
