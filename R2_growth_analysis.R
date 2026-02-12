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
source("/Users/aoyu/Desktop/Snail_Data/R2_data_cleaning.R") # cleaned master data sheet called all_data_r2
source("/Users/aoyu/Desktop/Snail_Data/R1_growth_analysis.R")
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

####### calculate growth: weekly jar level average length if no snails death that week #########
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
all_data_growth_length_covariate <- all_data_growth

all_data_growth <- all_data_growth %>%
  mutate(weekly_growth = ifelse(weekly_growth < 0, 0, weekly_growth)) %>% # include negative growth as 0 growth
  filter(!is.na(weekly_growth) & weekly_growth >= 0) # 1381 data points for growth analysis

# Calculate the IQR
Q1 <- quantile(all_data_growth$weekly_growth, 0.25, na.rm = TRUE)
Q3 <- quantile(all_data_growth$weekly_growth, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1

# Define the lower and upper bounds for outliers
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

# Identify outliers
outliers <- all_data_growth %>%
  filter(weekly_growth < lower_bound | weekly_growth > upper_bound)

### Grouped treatments: GAMM with Autoregressive Error Structure ----
gamm_model_jar_avgerage <- gamm(weekly_growth ~ Group + s(prev_snail_length_avg, bs = "cs", k = 5), correlation = corAR1(), data = all_data_growth)
summary(gamm_model_jar_avgerage$gam)

gamm_model_jar_avgerage1 <- gamm(weekly_growth ~ Group + s(Week_num, bs = "cs", k = 5), correlation = corAR1(), data = all_data_growth)
summary(gamm_model_jar_avgerage1$gam)
AIC(gamm_model_jar_avgerage$lme, gamm_model_jar_avgerage1$lme) # model with prev length is better

### Grouped treatments:pairwise comparison ----
## 1. Choose biologically representative values of prior-week snail length
##    (used to condition treatment comparisons on size)
ref_vals <- quantile(
  all_data_growth$prev_snail_length_avg,
  probs = c(0.25, 0.5, 0.75),
  na.rm = TRUE
)

## 2. Create a reference grid with all Treatment × size combinations
##    (this defines where the GAMM will be evaluated)
ref_grid <- expand.grid(
  Group = levels(all_data_growth$Group),
  prev_snail_length_avg = ref_vals
)

## 3. Extract the GAM linear predictor matrix at the reference grid
##    (this includes parametric + smooth basis terms)
Xp <- predict(
  gamm_model_jar_avgerage$gam,
  newdata = ref_grid,
  type = "lpmatrix"
)

## 4. Extract model coefficients and their covariance matrix
##    (these include penalized smooth terms and are required for uncertainty propagation)
beta <- coef(gamm_model_jar_avgerage$gam)
V <- vcov(gamm_model_jar_avgerage$gam)

## 5. Construct an emmeans-compatible object manually
##    (bypasses qdrg() and avoids penalized-smooth incompatibilities)
library(emmeans)

emm_obj <- emmobj(
  ## Predicted marginal means at each Treatment × size combination
  bhat = as.numeric(Xp %*% beta),
  
  ## Full variance–covariance matrix of those predictions
  V = Xp %*% V %*% t(Xp),
  
  ## Factor levels used for inference
  levels = list(
    Group = levels(all_data_growth$Group),
    prev_snail_length_avg = ref_vals
  ),
  
  ## Identity matrix: each row corresponds to one marginal mean
  linfct = diag(nrow(Xp)),
  
  ## Residual degrees of freedom from the GAM
  df = df.residual(gamm_model_jar_avgerage$gam),
  
  ## Specify that contrasts are grouped by size
  misc = list(by.vars = "prev_snail_length_avg")
)

## 6. Perform Tukey-adjusted pairwise comparisons among treatments
##    separately at each prior-week snail length
pairs(
  emm_obj,
  by = "prev_snail_length_avg",
  adjust = "tukey"
)

### Treatments: GAMM with Autoregressive Error Structure ----
gamm_model_jar_avgerage_abundance <- gamm(weekly_growth ~ Treatment + s(prev_snail_length_avg, bs = "cs", k = 5) , correlation = corAR1(), data = all_data_growth)
gamm_model_jar_avgerage_abundance_rm <- gamm(weekly_growth ~ Treatment + s(prev_snail_length_avg, bs = "cs", k = 5), random = list(Jar_num = ~1), correlation = corAR1(), data = all_data_growth)
summary(gamm_model_jar_avgerage_abundance$gam)

gamm_model_jar_avgerage_abundance1 <- gamm(weekly_growth ~ Treatment + s(Week_num, bs = "cs", k = 5), correlation = corAR1(), data = all_data_growth)

AIC(gamm_model_jar_avgerage_abundance$lme, gamm_model_jar_avgerage_abundance1$lme) # model with prev length is better
AIC(gamm_model_jar_avgerage_abundance_rm$lme, gamm_model_jar_avgerage_abundance$lme) # model without jar random effect is better (not that diff)

### Treatments: pairse comparison ----
## 1. Choose biologically representative values of prior-week snail length
##    (used to condition treatment comparisons on size)
ref_vals <- quantile(
  all_data_growth$prev_snail_length_avg,
  probs = c(0.25, 0.5, 0.75),
  na.rm = TRUE
)

## 2. Create a reference grid with all Treatment × size combinations
##    (this defines where the GAMM will be evaluated)
ref_grid <- expand.grid(
  Treatment = levels(all_data_growth$Treatment),
  prev_snail_length_avg = ref_vals
)

## 3. Extract the GAM linear predictor matrix at the reference grid
##    (this includes parametric + smooth basis terms)
Xp <- predict(
  gamm_model_jar_avgerage_abundance$gam,
  newdata = ref_grid,
  type = "lpmatrix"
)

## 4. Extract model coefficients and their covariance matrix
##    (these include penalized smooth terms and are required for uncertainty propagation)
beta <- coef(gamm_model_jar_avgerage_abundance$gam)
V <- vcov(gamm_model_jar_avgerage_abundance$gam)

## 5. Construct an emmeans-compatible object manually
##    (bypasses qdrg() and avoids penalized-smooth incompatibilities)
library(emmeans)

emm_obj <- emmobj(
  ## Predicted marginal means at each Treatment × size combination
  bhat = as.numeric(Xp %*% beta),
  
  ## Full variance–covariance matrix of those predictions
  V = Xp %*% V %*% t(Xp),
  
  ## Factor levels used for inference
  levels = list(
    Treatment = levels(all_data_growth$Treatment),
    prev_snail_length_avg = ref_vals
  ),
  
  ## Identity matrix: each row corresponds to one marginal mean
  linfct = diag(nrow(Xp)),
  
  ## Residual degrees of freedom from the GAM
  df = df.residual(gamm_model_jar_avgerage_abundance$gam),
  
  ## Specify that contrasts are grouped by size
  misc = list(by.vars = "prev_snail_length_avg")
)

## 6. Perform Tukey-adjusted pairwise comparisons among treatments
##    separately at each prior-week snail length
pairs(
  emm_obj,
  by = "prev_snail_length_avg",
  adjust = "tukey"
)

### plot weekly growth ----
R2_growth_plot_grouped <- ggplot(all_data_growth, aes(x = factor(Week_num), y = weekly_growth, color = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width  = 0.15,
    alpha  = 0.4,
    size   = 1
  ) +
  labs(x = "Weeks",
       y = "Mean Shell Length Growth per Snail (cm)",
       color = "Treatment") +
  scale_color_manual(values = c("Control" = "#6baf78",
                                "Seasoned Macro" = "#a5855f",
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
ggsave(filename = "R2_growth_plot_grouped_toplegend.png", plot = R2_growth_plot_grouped, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)


R2_growth_plot_4treatment <- ggplot(all_data_growth, aes(x = factor(Week_num), y = weekly_growth, color = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width  = 0.15,
    alpha  = 0.4,
    size   = 1
  ) +
  labs(x = "Weeks",
       y = "Mean Shell Length Growth per Snail (cm)",
       color = "Treatment") +
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
ggsave(filename = "R2_growth_plot_4treatment_toplegened.png", plot = R2_growth_plot_4treatment, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)

### model assumptions ----
#check model Independence of Errors assumption
gamm_resid <- resid(gamm_model_jar_avgerage$lme, type = "normalized")
acf(gamm_resid, main = "ACF of Residuals") # residuals are approximately independent, sufficient model adjustment

# check Homoscedasticity assumption
plot(gamm_model_jar_avgerage$gam$fitted.values, residuals(gamm_model_jar_avgerage$gam), 
     xlab = "Predicted Values", ylab = "Residuals")
abline(h = 0, col = "red")

#check model normality of residual assumption
qqnorm(gamm_resid, main = "Q-Q Plot of Residuals")

### marginal effect plots ----
# “Treatment effects were visualised using marginal mean predictions from a GAMM, averaging over the empirical distribution of prior-week snail length to account for nonlinear growth dynamics and temporal autocorrelation.”
# (connecting from the pairwise comparison code above)
trt_levels <- levels(all_data_growth$Treatment)
n_trt  <- length(trt_levels)
n_size <- length(ref_vals)
stopifnot(nrow(Xp) == n_trt * n_size)

pred_vec <- as.numeric(Xp %*% beta)  # length = n_trt * n_size

pred_mat <- matrix(
  pred_vec,
  nrow = n_trt,
  ncol = n_size,
  byrow = FALSE
)

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
V_pred <- Xp %*% V %*% t(Xp)   # (n_trt*n_size) × (n_trt*n_size)
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
  cld = c("a", "b", "c", "d", "ab", "b", "c", "c") # Replace with your real letters
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
# ggsave(filename = "R1R2_growth_plot_marginalized_mean_seperated.png", plot = plot_r1r2, bg = "white", path = "/Users/aoyu/Desktop/Snail_Data/Outputs", width = 9, height = 6, dpi = 300)

###check for potential errors or anomalies ----
test<- all_data_growth %>% filter(Group == "Seasoned Macro" | Group == "Virgin Macro") # hist comparison
test<- all_data_growth %>% filter(Group == "Virgin Macro") # single hist
ggplot(all_data_growth, aes(x = weekly_growth, color = Group)) +
  geom_histogram() +
  theme_minimal()  +
  facet_wrap(~ Group) 

ggplot(test, aes(x = weekly_growth, color = Group)) +
  geom_histogram() +
  theme_minimal()  +
  facet_wrap(~ Week_num) 

### initial vs end of the experiment size report #####
source("/Users/aoyu/Desktop/Snail_Data/R2_data_cleaning.R") # cleaned master data sheet called all_data
randomization_check <- all_data_r2 %>% 
  filter(Week_num == 1) %>% 
  select(Color, Jar_num, Snail_length_1, Snail_length_2, Snail_length_3) %>% 
  filter(Color == "G"|Color == "Y" |Color ==  "W-1"|Color ==  "W-3"|Color ==  "O") 

table(randomization_check$Color, randomization_check$Jar_num)

randomization_check_long <- randomization_check %>%
  pivot_longer(cols = starts_with("Snail_length_"),
               names_to = "Snail_ID",
               values_to = "Snail_length") %>% 
  mutate(Snail_length = as.numeric(Snail_length))

# ANOVA to check if initial snail size differs by treatment
anova_result <- aov(Snail_length ~ Color, data = randomization_check_long)
summary(anova_result) # significant difference numerically (could be due to snail picture angle), but acceptable biologically
TukeyHSD(anova_result) # post-hoc test

# Boxplot to visualize the distribution of initial snail sizes by treatment
ggplot(randomization_check_long, aes(x = Color, y = Snail_length, color = Color)) +
  geom_boxplot() +
  labs(title = "Distribution of Initial Snail Sizes by Treatment",
       x = "Treatment Group",
       y = "Snail Length") +
  theme_minimal() 

# final size and whether diff between treatment:
randomization_check <- all_data_r2 %>% 
  filter(Week_num == 13) %>%
  select(Color, Jar_num, Snail_length_1, Snail_length_2, Snail_length_3) %>% 
  filter(Color == "G"|Color == "Y" |Color ==  "W-1"|Color ==  "W-3"|Color ==  "O") 

table(randomization_check$Color, randomization_check$Jar_num)

final_size_long <- randomization_check %>%
  pivot_longer(cols = starts_with("Snail_length_"),
               names_to = "Snail_ID",
               values_to = "Snail_length") %>% 
  mutate(Snail_length = as.numeric(Snail_length))

# ANOVA to check if initial snail size differs by treatment
anova_result <- aov(Snail_length ~ Color, data = final_size_long)
summary(anova_result) # significant difference numerically (could be due to snail picture angle), but acceptable biologically
TukeyHSD(anova_result) # post-hoc test

# Boxplot to visualize the distribution of initial snail sizes by treatment
ggplot(randomization_check_long, aes(x = Color, y = Snail_length, color = Color)) +
  geom_boxplot() +
  labs(title = "Distribution of Initial Snail Sizes by Treatment",
       x = "Treatment Group",
       y = "Snail Length") +
  theme_minimal() 
