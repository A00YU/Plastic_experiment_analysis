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

####### load data and creat a df for each round called parameters_data_rx #########
# R1 survival rate l_x
source("/Users/aoyu/Desktop/Snail_Data/R1_survival_analysis.R") 

death_long_r1 <- death_long %>%
  mutate(cumulative_survived = total - cumulative_death,
         cumulative_relative_survival = cumulative_survived / total) %>% 
  select(week_of_death, Treatment, cumulative_relative_survival) 

baseline_rows_r1 <- data.frame(
  week_of_death = rep(1, 4),
  Treatment = c("Control", "Micro", "Macro", "Macro+Micro"),
  cumulative_relative_survival = 1)

parameters_data_r1 <- rbind(baseline_rows_r1, death_long_r1) # l_x values in parameters_data$cumulative_relative_survival

# table(parameters_data_r1$week_of_death, parameters_data_r1$Treatment)

# R2 survival rate l_x
source("/Users/aoyu/Desktop/Snail_Data/R2_survival_analysis.R") 

death_long_r2 <- death_long %>%
  mutate(cumulative_survived = total - cumulative_death,
         cumulative_relative_survival = cumulative_survived / total) %>% 
  select(week_of_death, Treatment, cumulative_relative_survival) 

baseline_rows_r2 <- data.frame(
  week_of_death = rep(1, 5),
  Treatment = c("Control", "Seasoned Macro - Low","Seasoned Macro - High", "Virgin Macro - Low", "Virgin Macro - High"),
  cumulative_relative_survival = 1)

parameters_data_r2 <- rbind(baseline_rows_r2, death_long_r2) # l_x values in parameters_data$cumulative_relative_survival

# table(parameters_data_r2$week_of_death, parameters_data_r2$Treatment)

# R1 fecundity b_x
source("/Users/aoyu/Desktop/Snail_Data/R1_fecundity_analysis.R") 

all_data_repro_r1 <- all_data_repro %>% 
  select(Week_num, Treatment, total_egg_mass_num, snails_num) %>%
  mutate(eggs_mass_per_snail = total_egg_mass_num / snails_num,
         egg_per_snail = round(eggs_mass_per_snail * 10)) # conservative assuming 10 eggs per egg mass 

table(all_data_repro_r1$Week_num, all_data_repro_r1$Treatment) # occational missing jars are due to all 3 snail death

all_data_repro_r1 <- all_data_repro_r1 %>%
  group_by(Week_num, Treatment) %>%
  summarise(egg_per_snail = mean(egg_per_snail, na.rm = TRUE)) 

baseline_rows_r1 <- data.frame(
  Week_num = rep(1, 4),
  Treatment = c("Control", "Micro", "Macro", "Macro+Micro"),
  egg_per_snail = 0)

all_data_repro_r1 <- rbind(baseline_rows_r1, all_data_repro_r1)

parameters_data_r1 <- parameters_data_r1 %>%
  left_join(all_data_repro_r1, by = c("week_of_death" = "Week_num", "Treatment" = "Treatment")) %>%
  rename(week = week_of_death,
         l_x = cumulative_relative_survival,
         b_x = egg_per_snail) 

# R2 fecundity b_x
source("/Users/aoyu/Desktop/Snail_Data/R2_fecundity_analysis.R") 

all_data_repro_r2 <- all_data_repro %>% 
  select(Week_num, Treatment, total_egg_mass_num, snails_num) %>%
  mutate(eggs_mass_per_snail = total_egg_mass_num / snails_num,
         egg_per_snail = round(eggs_mass_per_snail * 10)) # conservative assuming 10 eggs per egg mass 

table(all_data_repro_r2$Week_num, all_data_repro_r2$Treatment) # occational missing jars are due to all 3 snail death

all_data_repro_r2 <- all_data_repro_r2 %>%
  group_by(Week_num, Treatment) %>%
  summarise(egg_per_snail = mean(egg_per_snail, na.rm = TRUE)) 

baseline_rows_r2 <- data.frame(
  Week_num = rep(1, 5),
  Treatment = c("Control", "Seasoned Macro - Low","Seasoned Macro - High", "Virgin Macro - Low", "Virgin Macro - High"),
  egg_per_snail = 0)

all_data_repro_r2 <- rbind(baseline_rows_r2, all_data_repro_r2)

parameters_data_r2 <- parameters_data_r2 %>%
  left_join(all_data_repro_r2, by = c("week_of_death" = "Week_num", "Treatment" = "Treatment")) %>%
  rename(week = week_of_death,
         l_x = cumulative_relative_survival,
         b_x = egg_per_snail) 

# R1 offspring fitness h_x
source("/Users/aoyu/Desktop/Snail_Data/R1_offspring_fitness_analysis.R")

all_data_hatch_r1 <- all_data_hatch %>%
  select(Week_num, Treatment, hatch_success) %>%
  group_by(Week_num, Treatment) %>%
  summarise(h_x = mean(hatch_success, na.rm = TRUE)) %>% 
  filter(Week_num != 13) # no data collected for week 13 offspring hatching

baseline_rows_r1 <- data.frame(
  Week_num = rep(1, 4),
  Treatment = c("Control", "Micro", "Macro", "Macro+Micro"),
  h_x = 0)

week13_rows_r1 <- all_data_hatch_r1 %>%
  filter(Week_num == 12) %>%
  mutate(Week_num = 13)

all_data_hatch_r1 <- rbind(baseline_rows_r1, all_data_hatch_r1, week13_rows_r1)

parameters_data_r1 <- parameters_data_r1 %>%
  left_join(all_data_hatch_r1, by = c("week" = "Week_num", "Treatment" = "Treatment")) 

# R2 offspring fitness h_x
source("/Users/aoyu/Desktop/Snail_Data/R2_offspring_fitness_analysis.R")

all_data_hatch_r2 <- all_data_hatch %>%
  select(Week_num, Treatment, hatch_success) %>%
  group_by(Week_num, Treatment) %>%
  summarise(h_x = mean(hatch_success, na.rm = TRUE)) 

baseline_rows_r2 <- data.frame(
  Week_num = rep(1, 5),
  Treatment = c("Control", "Seasoned Macro - Low","Seasoned Macro - High", "Virgin Macro - Low", "Virgin Macro - High"),
  h_x = 0)

all_data_hatch_r2 <- rbind(baseline_rows_r2, all_data_hatch_r2)

parameters_data_r2 <- parameters_data_r2 %>%
  left_join(all_data_hatch_r2, by = c("week" = "Week_num", "Treatment" = "Treatment")) 

# R1 growth to maturity 
source("/Users/aoyu/Desktop/Snail_Data/R1_growth_analysis.R") 

gam_summary <- summary(gamm_model_jar_avgerage$gam)

# Extract the coefficients and standard errors
coefficients <- gam_summary$p.coeff

# Create a new data frame with Treatment, estimates, and standard errors
model_estimates_df <- data.frame(
  Treatment = names(coefficients),
  Estimate = coefficients)

# Mutate a new variable with the exponentiated estimates
model_estimates_df <- model_estimates_df %>%
  mutate(weeks_to_maturity = (0.8*13)/(0.8+Estimate*13)) %>% # 0.8 cm is the threshold size to reach maturity; 13 week is the theoretical growth to maturity
  filter(Treatment != "(Intercept)") %>%
  mutate(Treatment = dplyr::recode(Treatment, 
                            "TreatmentMicro" = "Micro",
                            "TreatmentMacro" = "Macro",
                            "TreatmentMacro+Micro" = "Macro+Micro")) %>% 
  select(Treatment, weeks_to_maturity)

control_reference_estimate_df <- data.frame(
  Treatment = "Control",
  weeks_to_maturity = 13)

model_estimates_df <- rbind(model_estimates_df, control_reference_estimate_df) 

parameters_data_r1 <- parameters_data_r1 %>%
  left_join(model_estimates_df, by = "Treatment") 

# R2 growth to maturity
source("/Users/aoyu/Desktop/Snail_Data/R2_growth_analysis.R") 

gam_summary <- summary(gamm_model_jar_avgerage_abundance$gam)

# Extract the coefficients and standard errors
coefficients <- gam_summary$p.coeff

# Create a new data frame with Treatment, estimates, and standard errors
model_estimates_df <- data.frame(
  Treatment = names(coefficients),
  Estimate = coefficients)

# Mutate a new variable with the exponentiated estimates
model_estimates_df <- model_estimates_df %>%
  mutate(weeks_to_maturity = (0.8*13)/(0.8+Estimate*13)) %>% # 0.8 cm is the threshold size to reach maturity; 13 week is the theoretical growth to maturity
  filter(Treatment != "(Intercept)") %>%
  mutate(Treatment = dplyr::recode(Treatment, 
                                   "TreatmentSeasoned Macro - Low" = "Seasoned Macro - Low",
                                   "TreatmentSeasoned Macro - High" = "Seasoned Macro - High",
                                   "TreatmentVirgin Macro - Low" = "Virgin Macro - Low",
                                   "TreatmentVirgin Macro - High" = "Virgin Macro - High")) %>% 
  select(Treatment, weeks_to_maturity)

control_reference_estimate_df <- data.frame(
  Treatment = "Control",
  weeks_to_maturity = 13)

model_estimates_df <- rbind(model_estimates_df, control_reference_estimate_df) 

parameters_data_r2 <- parameters_data_r2 %>%
  left_join(model_estimates_df, by = "Treatment") 


####### R1:Euler Lotka equation to get treatment dependent lambda#########
# solve for finite growth rate (lambda)
# l_x: the fraction of individuals surviving to age x
# ? Adult growth based on literature increases reproductive output and impact survival, which we have direct measurement of
# b_x:the number of offspring born per individual of age x (assume all viable? if so adjusted offspring needed)
# Adjust birth rate to account for actual hatching success
# b_x_adjust <- b_x * h_x
# h_x <- c()  # hatching success rate
# age classes in terms of experiment weeks; in theory from age 1 to max age of reproduction

# Define the Euler-Lotka equation function
euler_lotka <- function(l_hatched_to_maturity, lambda, l_x, b_x, h_x, ages, weeks_to_maturity) {
  (l_hatched_to_maturity * sum(lambda^(-ages-weeks_to_maturity) * l_x * b_x * h_x)) - 1
}

# Define a wrapper function for uniroot that only takes lambda as input
euler_lotka_root <- function(lambda) {
  euler_lotka(l_hatched_to_maturity,lambda, l_x, b_x, h_x, ages, weeks_to_maturity)
}

# Define a function to calculate lambda and doubling time for a given treatment
calculate_lambda_and_doubling_time <- function(treatment_data) {
  l_x <- treatment_data$l_x
  b_x <- treatment_data$b_x
  h_x <- treatment_data$h_x
  ages <- 1:13
  weeks_to_maturity <- treatment_data$weeks_to_maturity
  l_hatched_to_maturity <- 0.95
  
  euler_lotka_root <- function(lambda) {
    euler_lotka(l_hatched_to_maturity, lambda, l_x, b_x, h_x, ages, weeks_to_maturity)
  }
  
  result <- uniroot(euler_lotka_root, interval = c(0.01, 10))
  lambda <- result$root
  doubling_time <- log(2) / log(lambda)
  
  return(list(lambda = lambda, doubling_time = doubling_time))
}

# Filter data for each treatment
control_r1 <- parameters_data_r1 %>% filter(Treatment == "Control")
micro_r1 <- parameters_data_r1 %>% filter(Treatment == "Micro")
macro_r1 <- parameters_data_r1 %>% filter(Treatment == "Macro")
macro_micro_r1 <- parameters_data_r1 %>% filter(Treatment == "Macro+Micro")

# Calculate lambda and doubling time for each treatment
control_r1_results <- calculate_lambda_and_doubling_time(control_r1)
micro_r1_results <- calculate_lambda_and_doubling_time(micro_r1)
macro_r1_results <- calculate_lambda_and_doubling_time(macro_r1)
macro_micro_r1_results <- calculate_lambda_and_doubling_time(macro_micro_r1)

# Create a results table
results_r1 <- data.frame(
  Treatment = c("Control", "Micro", "Macro", "Macro+Micro"),
  Lambda = c(control_r1_results$lambda, micro_r1_results$lambda, macro_r1_results$lambda, macro_micro_r1_results$lambda),
  DoublingTime = c(control_r1_results$doubling_time, micro_r1_results$doubling_time, macro_r1_results$doubling_time, macro_micro_r1_results$doubling_time),
  stringsAsFactors = FALSE
)

# Print the results table
print(results_r1)

####### R2:Euler Lotka equation to get treatment dependent lambda#########
# Define the Euler-Lotka equation function
euler_lotka <- function(l_hatched_to_maturity, lambda, l_x, b_x, h_x, ages, weeks_to_maturity) {
  (l_hatched_to_maturity * sum(lambda^(-ages-weeks_to_maturity) * l_x * b_x * h_x)) - 1
} # account for 7 weeks to reach maturity # the l_0m probability of snails hatched from eggs surviving to maturity

# Define a function to calculate lambda and doubling time for a given treatment
calculate_lambda_and_doubling_time <- function(treatment_data) {
  l_x <- treatment_data$l_x
  b_x <- treatment_data$b_x
  h_x <- treatment_data$h_x
  ages <- 1:13
  weeks_to_maturity <- treatment_data$weeks_to_maturity
  l_hatched_to_maturity <- 0.95
  
  euler_lotka_root <- function(lambda) {
    euler_lotka(l_hatched_to_maturity, lambda, l_x, b_x, h_x, ages, weeks_to_maturity)
  }
  
  result <- uniroot(euler_lotka_root, interval = c(0.01, 10))
  lambda <- result$root
  doubling_time <- log(2) / log(lambda)
  
  return(list(lambda = lambda, doubling_time = doubling_time))
}

# Filter data for each treatment
control_r2 <- parameters_data_r2 %>% filter(Treatment == "Control")
seasoned_macro_low_r2 <- parameters_data_r2 %>% filter(Treatment == "Seasoned Macro - Low")
seasoned_macro_high_r2 <- parameters_data_r2 %>% filter(Treatment == "Seasoned Macro - High")
virgin_macro_low_r2 <- parameters_data_r2 %>% filter(Treatment == "Virgin Macro - Low")
virgin_macro_high_r2 <- parameters_data_r2 %>% filter(Treatment == "Virgin Macro - High")

# Calculate lambda and doubling time for each treatment
control_r2_results <- calculate_lambda_and_doubling_time(control_r2)
seasoned_macro_low_r2_results <- calculate_lambda_and_doubling_time(seasoned_macro_low_r2)
seasoned_macro_high_r2_results <- calculate_lambda_and_doubling_time(seasoned_macro_high_r2)
virgin_macro_low_r2_results <- calculate_lambda_and_doubling_time(virgin_macro_low_r2)
virgin_macro_high_r2_results <- calculate_lambda_and_doubling_time(virgin_macro_high_r2)

# Create a results table
results_r2_treatment <- data.frame(
  Treatment = c("Control", "Seasoned Macro - Low", "Seasoned Macro - High", "Virgin Macro - Low", "Virgin Macro - High"),
  Lambda = c(control_r2_results$lambda, seasoned_macro_low_r2_results$lambda, seasoned_macro_high_r2_results$lambda, virgin_macro_low_r2_results$lambda, virgin_macro_high_r2_results$lambda),
  DoublingTime = c(control_r2_results$doubling_time, seasoned_macro_low_r2_results$doubling_time, seasoned_macro_high_r2_results$doubling_time, virgin_macro_low_r2_results$doubling_time, virgin_macro_high_r2_results$doubling_time),
  stringsAsFactors = FALSE
)

# Print the results table
print(results_r2_treatment)
