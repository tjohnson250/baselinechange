# R Code Conversion from SAS
# eAppendix - Causal Inference Analysis with G-computation

# Load required libraries
library(tidyverse)
library(broom)

# Set seed for reproducibility
set.seed(123)

# Generate data
n <- 1000000
data_a <- data.frame(
  id = 1:n,
  PA = 0.5,
  A = rbinom(n, 1, 0.5),
  PU = 0.3,
  U = rbinom(n, 1, 0.3)
)

# Generate W based on A and U
data_a$PW <- 0.2 + 0.4 * data_a$A + 0.3 * data_a$U
data_a$W <- rbinom(n, 1, data_a$PW)

# Generate S based on W
data_a$PS <- 0.1 + 0.8 * data_a$W
data_a$S <- rbinom(n, 1, data_a$PS)

# Generate Y based on U and A
data_a$PY <- 0.3 + 0.5 * data_a$U - 0.2 * data_a$A
data_a$Y <- rbinom(n, 1, data_a$PY)

# Truth - Full population model
cat("\n========== TRUTH ==========\n")
# Use quasibinomial to avoid warnings with identity link
truth_model <- glm(Y ~ A, data = data_a, family = quasibinomial(link = "identity"))
summary(truth_model)
print(tidy(truth_model))

# Crude - Restricted to S=0
cat("\n========== CRUDE (S=0) ==========\n")
data_s0 <- data_a[data_a$S == 0, ]
crude_model <- glm(Y ~ A, data = data_s0, family = quasibinomial(link = "identity"))
summary(crude_model)
print(tidy(crude_model))

# Stratified on W - Restricted to S=0
cat("\n========== STRATIFIED (S=0) ==========\n")
stratified_model <- glm(Y ~ A + W, data = data_s0, family = quasibinomial(link = "identity"))
summary(stratified_model)
print(tidy(stratified_model))

# G-formula (G-computation)
cat("\n========== G-COMP MODELS ==========\n")

# Model for Y|A,W,S=0
y_model <- glm(Y ~ A + W, data = data_s0, family = quasibinomial(link = "identity"))
yest <- tidy(y_model)
cat("\nG-comp model P(Y|A,W,S=0):\n")
print(yest)

# Model for W|A (using full data)
w_model <- glm(W ~ A, data = data_a, family = quasibinomial(link = "identity"))
west <- tidy(w_model)
cat("\nG-comp model P(W|A):\n")
print(west)

# Extract coefficients for G-computation
# From W model
intw <- west$estimate[west$term == "(Intercept)"]
intw_se <- west$std.error[west$term == "(Intercept)"]
wona <- west$estimate[west$term == "A"]
wona_se <- west$std.error[west$term == "A"]

# From Y model
inty <- yest$estimate[yest$term == "(Intercept)"]
inty_se <- yest$std.error[yest$term == "(Intercept)"]
yona <- yest$estimate[yest$term == "A"]
yona_se <- yest$std.error[yest$term == "A"]
yonw <- yest$estimate[yest$term == "W"]
yonw_se <- yest$std.error[yest$term == "W"]

# G-computation calculation
gcomp0 <- 0
gcomp1 <- 0

for (A in 0:1) {
  for (W in 0:1) {
    # P(W|A)
    pw <- intw + wona * A
    if (W == 0) pw <- 1 - pw
    
    # E[Y|A,W]
    ey <- inty + yona * A + yonw * W
    
    if (A == 0) {
      gcomp0 <- gcomp0 + pw * ey
    } else {
      gcomp1 <- gcomp1 + pw * ey
    }
  }
}

# Calculate the causal effect
gcomp <- gcomp1 - gcomp0

cat("\n========== G-COMPUTATION RESULT ==========\n")
cat("E[Y(1)] - E[Y(0)] =", gcomp, "\n")

# Create summary data frame similar to SAS output
results_df <- data.frame(
  Parameter = c("Intercept_W", "A_on_W", "Intercept_Y", "A_on_Y", "W_on_Y", "G-comp"),
  Estimate = c(intw, wona, inty, yona, yonw, gcomp),
  StdError = c(intw_se, wona_se, inty_se, yona_se, yonw_se, NA)
)

cat("\n========== SUMMARY OF ALL PARAMETERS ==========\n")
print(results_df, row.names = FALSE)

# Alternative approach using predict() for verification
cat("\n========== VERIFICATION USING PREDICT ==========\n")

# Create counterfactual datasets
cf_data <- expand.grid(A = 0:1, W = 0:1)

# Get P(W|A) for each A
pw_a0 <- predict(w_model, newdata = data.frame(A = 0), type = "response")
pw_a1 <- predict(w_model, newdata = data.frame(A = 1), type = "response")

# Calculate E[Y(0)] and E[Y(1)]
ey_a0 <- sum(c(1 - pw_a0, pw_a0) * 
               predict(y_model, newdata = data.frame(A = 0, W = 0:1), type = "response"))
ey_a1 <- sum(c(1 - pw_a1, pw_a1) * 
               predict(y_model, newdata = data.frame(A = 1, W = 0:1), type = "response"))

gcomp_verify <- ey_a1 - ey_a0
cat("G-comp (using predict):", gcomp_verify, "\n")

# ========== REPRODUCING BRESKIN ET AL. PAPER RESULTS ==========
cat("\n\n========== BRESKIN ET AL. (2018) PAPER IMPLEMENTATION ==========\n")

# Following the exact identification strategy from the paper:
# E[Y(a)] = Σ_w E[Y|W=w,S=0,A=a] * P(W=w|A=a)

# This is valid because:
# 1. S(a) ⊥ Y(a) | W(a) (from SWIG)
# 2. W(a) ⊥ A (randomization)
# 3. Y(a) ⊥ A | {W(a), S(a)} (from SWIG)

# Step 1: Estimate P(W|A) from full data (before selection)
cat("\nStep 1: Model P(W|A) from full data\n")
w_full_model <- glm(W ~ A, data = data_a, family = quasibinomial(link = "identity"))
print(tidy(w_full_model))

# Step 2: Estimate E[Y|A,W,S=0] from selected data
cat("\nStep 2: Model E[Y|A,W,S=0] from S=0 data\n")
y_selected_model <- glm(Y ~ A + W, data = data_s0, family = quasibinomial(link = "identity"))
print(tidy(y_selected_model))

# Step 3: Apply the identification formula
cat("\nStep 3: Apply identification formula E[Y(a)] = Σ_w E[Y|W=w,S=0,A=a] * P(W=w|A=a)\n")

# For A=0
p_w0_a0 <- predict(w_full_model, newdata = data.frame(A = 0), type = "response")
p_w1_a0 <- 1 - p_w0_a0
e_y_w0_a0 <- predict(y_selected_model, newdata = data.frame(A = 0, W = 0), type = "response")
e_y_w1_a0 <- predict(y_selected_model, newdata = data.frame(A = 0, W = 1), type = "response")
e_y_a0 <- p_w1_a0 * e_y_w0_a0 + p_w0_a0 * e_y_w1_a0

# For A=1
p_w0_a1 <- predict(w_full_model, newdata = data.frame(A = 1), type = "response")
p_w1_a1 <- 1 - p_w0_a1
e_y_w0_a1 <- predict(y_selected_model, newdata = data.frame(A = 1, W = 0), type = "response")
e_y_w1_a1 <- predict(y_selected_model, newdata = data.frame(A = 1, W = 1), type = "response")
e_y_a1 <- p_w1_a1 * e_y_w0_a1 + p_w0_a1 * e_y_w1_a1

# Causal effect
breskin_effect <- e_y_a1 - e_y_a0
cat("\nBreskin et al. approach - Causal effect E[Y(1)] - E[Y(0)]:", breskin_effect, "\n")

# Verify this matches our G-computation
cat("\nVerification - matches our G-computation result:", gcomp, "\n")
cat("Difference:", abs(breskin_effect - gcomp), "\n")

# ========== INVERSE PROBABILITY WEIGHTING (IPW) ==========
cat("\n\n========== INVERSE PROBABILITY WEIGHTING ==========\n")

# IPW needs to account for BOTH treatment and selection/censoring
# We need to model P(S=0|A,W) and create weights that account for selection

# 1. Model the selection mechanism P(S=0|W) in full data
cat("\nModeling selection mechanism P(S=0|W):\n")
# First create the S=0 indicator (observed)
data_a$observed <- 1 - data_a$S

selection_model <- glm(observed ~ W, data = data_a, family = binomial(link = "logit"))
print(tidy(selection_model))

# Get probability of being observed for everyone
data_a$p_observed <- predict(selection_model, type = "response")

# 2. Model propensity score P(A|W) in full data
ps_model_full <- glm(A ~ W, data = data_a, family = binomial(link = "logit"))
cat("\n\nPropensity Score Model P(A|W) in full data:\n")
print(tidy(ps_model_full))

# Get propensity scores
data_a$ps <- predict(ps_model_full, type = "response")

# 3. Create combined weights for those observed (S=0)
# Weight = 1/[P(A=a|W) * P(S=0|W)]
data_s0 <- data_a[data_a$S == 0, ]
data_s0$treatment_weight <- ifelse(data_s0$A == 1, 
                                   1 / data_s0$ps,
                                   1 / (1 - data_s0$ps))
data_s0$selection_weight <- 1 / data_s0$p_observed
data_s0$combined_ipw <- data_s0$treatment_weight * data_s0$selection_weight

# Check weight distribution
cat("\n\nCombined IPW Weight Summary (accounting for selection):\n")
cat("Treated (A=1):\n")
summary(data_s0$combined_ipw[data_s0$A == 1])
cat("\nControl (A=0):\n")
summary(data_s0$combined_ipw[data_s0$A == 0])

# IPW estimator with selection weights
ipw_selection_model <- glm(Y ~ A, data = data_s0, weights = combined_ipw, 
                           family = quasibinomial(link = "identity"))
cat("\n\nIPW Model Results (with selection adjustment):\n")
print(tidy(ipw_selection_model))
ipw_selection_effect <- coef(ipw_selection_model)["A"]

# 4. Stabilized weights accounting for selection
# Stabilized weight = P(A=a) * P(S=0) / [P(A=a|W) * P(S=0|W)]
p_a <- mean(data_a$A)
p_s0 <- mean(data_a$observed)

data_s0$stabilized_combined <- ifelse(data_s0$A == 1,
                                      (p_a * p_s0) / (data_s0$ps * data_s0$p_observed),
                                      ((1 - p_a) * p_s0) / ((1 - data_s0$ps) * data_s0$p_observed))

cat("\n\nStabilized Combined Weight Summary:\n")
cat("Treated (A=1):\n")
summary(data_s0$stabilized_combined[data_s0$A == 1])
cat("\nControl (A=0):\n")
summary(data_s0$stabilized_combined[data_s0$A == 0])

# Stabilized IPW with selection
sipw_selection_model <- glm(Y ~ A, data = data_s0, weights = stabilized_combined, 
                            family = quasibinomial(link = "identity"))
cat("\n\nStabilized IPW Model Results (with selection adjustment):\n")
print(tidy(sipw_selection_model))
sipw_selection_effect <- coef(sipw_selection_model)["A"]

# 5. Alternative: Model P(S=0|A,W) directly
cat("\n\n========== ALTERNATIVE IPW APPROACH ==========\n")
# This models the selection probability conditional on both A and W
selection_model_aw <- glm(observed ~ A * W, data = data_a, family = binomial(link = "logit"))
cat("\nSelection model P(S=0|A,W):\n")
print(tidy(selection_model_aw))

data_a$p_observed_aw <- predict(selection_model_aw, type = "response")

# Make sure we're working with the correct data_s0 that has all the weights
data_s0$p_observed_aw <- data_a$p_observed_aw[data_a$S == 0]

# Create weights using this model
data_s0$ipw_alt <- ifelse(data_s0$A == 1,
                          1 / (data_s0$ps * data_s0$p_observed_aw),
                          1 / ((1 - data_s0$ps) * data_s0$p_observed_aw))

# Alternative IPW estimator
ipw_alt_model <- glm(Y ~ A, data = data_s0, weights = ipw_alt, 
                     family = quasibinomial(link = "identity"))
cat("\n\nAlternative IPW Model Results:\n")
print(tidy(ipw_alt_model))
ipw_alt_effect <- coef(ipw_alt_model)["A"]

# Note about IPW with selection bias
cat("\n\nNOTE: IPW with selection bias (S=0) requires strong assumptions:\n")
cat("1. No unmeasured confounders of S and Y relationship\n")
cat("2. Positivity: P(S=0|A,W) > 0 for all A,W combinations\n")
cat("3. In this case, we have unmeasured confounder U affecting both W and Y\n")
cat("   This violates the assumptions needed for IPW to recover the true effect\n")

# Original IPW (without selection adjustment) for comparison
# Model propensity score among those with S=0 only
ps_model_s0 <- glm(A ~ W, data = data_s0, family = binomial(link = "logit"))
data_s0$ps_s0 <- predict(ps_model_s0, type = "response")
data_s0$ipw_naive <- ifelse(data_s0$A == 1, 
                            1 / data_s0$ps_s0,
                            1 / (1 - data_s0$ps_s0))

ipw_naive_model <- glm(Y ~ A, data = data_s0, weights = ipw_naive, 
                       family = quasibinomial(link = "identity"))
ipw_naive_effect <- coef(ipw_naive_model)["A"]

# 3. Doubly Robust Estimator (combines IPW and outcome modeling)
cat("\n\n========== DOUBLY ROBUST ESTIMATOR ==========\n")

# Use the outcome model from stratified analysis and stabilized weights
dr_model <- glm(Y ~ A + W, data = data_s0, weights = stabilized_combined,
                family = quasibinomial(link = "identity"))
cat("\nDoubly Robust Model Results:\n")
print(tidy(dr_model))

# For doubly robust, we need to compute the counterfactual predictions
# Predict Y(1) and Y(0) for everyone
data_s0$y1_pred <- predict(dr_model, 
                           newdata = transform(data_s0, A = 1), 
                           type = "response")
data_s0$y0_pred <- predict(dr_model, 
                           newdata = transform(data_s0, A = 0), 
                           type = "response")

# Doubly robust estimator
dr_effect <- mean(data_s0$stabilized_combined * (data_s0$A * data_s0$Y + (1 - data_s0$A) * data_s0$y1_pred)) - 
  mean(data_s0$stabilized_combined * ((1 - data_s0$A) * data_s0$Y + data_s0$A * data_s0$y0_pred))

cat("\nDoubly Robust Effect Estimate:", dr_effect, "\n")
cat("\nDoubly Robust Model Results:\n")
print(tidy(dr_model))

# For doubly robust, we need to compute the counterfactual predictions
# Predict Y(1) and Y(0) for everyone
data_s0$y1_pred <- predict(dr_model, 
                           newdata = transform(data_s0, A = 1), 
                           type = "response")
data_s0$y0_pred <- predict(dr_model, 
                           newdata = transform(data_s0, A = 0), 
                           type = "response")

# Doubly robust estimator
dr_effect <- mean(data_s0$stabilized_combined * (data_s0$A * data_s0$Y + (1 - data_s0$A) * data_s0$y1_pred)) - 
  mean(data_s0$sw * ((1 - data_s0$A) * data_s0$Y + data_s0$A * data_s0$y0_pred))

cat("\nDoubly Robust Effect Estimate:", dr_effect, "\n")

# 4. Summary of all estimates
cat("\n\n========== SUMMARY OF ALL CAUSAL EFFECT ESTIMATES ==========\n")
summary_effects <- data.frame(
  Method = c("Truth (full data)", 
             "Crude (S=0)", 
             "Stratified (S=0)", 
             "G-computation", 
             "IPW (naive - S=0 only)", 
             "IPW (with selection weights)",
             "Stabilized IPW (with selection)", 
             "Alternative IPW (P(S=0|A,W))",
             "Doubly Robust (with selection)"),
  Estimate = c(coef(truth_model)["A"], 
               coef(crude_model)["A"], 
               coef(stratified_model)["A"],
               gcomp, 
               ipw_naive_effect,
               ipw_selection_effect, 
               sipw_selection_effect,
               ipw_alt_effect,
               dr_effect),
  Description = c("True causal effect", 
                  "Naive estimate among S=0",
                  "Adjusted for W among S=0",
                  "Standardization over W",
                  "IPW without selection adjustment",
                  "IPW accounting for P(S=0|W)",
                  "Stabilized IPW with selection",
                  "IPW using P(S=0|A,W)",
                  "Combines IPW and outcome modeling")
)
print(summary_effects, row.names = FALSE)

cat("\n\nKEY INSIGHTS:\n")
cat("1. Naive IPW (among S=0 only) doesn't account for selection bias\n")
cat("2. IPW with selection weights attempts to correct for S=0 selection\n")
cat("3. However, due to unmeasured confounder U, IPW cannot fully recover the true effect\n")
cat("4. G-computation works because of the special independence shown by SWIGs\n")
cat("5. The weights can become extreme when P(S=0|W) is small (strong selection)\n")

# Plot weight distributions
cat("\n\n========== WEIGHT DISTRIBUTION PLOTS ==========\n")
par(mfrow = c(2, 2))

# IPW weights by treatment
boxplot(combined_ipw ~ A, data = data_s0, 
        main = "IPW Weights by Treatment", 
        xlab = "Treatment (A)", ylab = "IPW Weight",
        col = c("lightblue", "lightcoral"))

# Stabilized weights by treatment
boxplot(stabilized_combined ~ A, data = data_s0, 
        main = "Stabilized Weights by Treatment", 
        xlab = "Treatment (A)", ylab = "Stabilized Weight",
        col = c("lightblue", "lightcoral"))

# Propensity score distribution
hist(data_s0$ps[data_s0$A == 1], breaks = 30, col = rgb(1,0,0,0.5),
     main = "Propensity Score Distribution", xlab = "Propensity Score",
     xlim = c(0, 1))
hist(data_s0$ps[data_s0$A == 0], breaks = 30, col = rgb(0,0,1,0.5), add = TRUE)
legend("topright", c("A=1", "A=0"), fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))

# ========== VISUALIZING THE CAUSAL STRUCTURE ==========
cat("\n\n========== VISUALIZING THE CAUSAL STRUCTURE ==========\n")

# Create a simple visualization of the relationships
par(mfrow = c(1, 2))

# Plot 1: Show selection bias impact
selected_effect <- coef(crude_model)["A"]
true_effect <- coef(truth_model)["A"]
barplot(c(true_effect, selected_effect, gcomp), 
        names.arg = c("Truth", "Selected (S=0)", "G-computation"),
        main = "Causal Effect Estimates",
        ylab = "Effect of A on Y",
        col = c("darkgreen", "red", "blue"),
        ylim = c(-0.3, 0))
abline(h = true_effect, lty = 2, col = "darkgreen")

# Plot 2: Show how conditioning affects the data
plot(density(data_a$PY), main = "Distribution of P(Y)", 
     xlab = "P(Y)", xlim = c(0, 1), lwd = 2)
lines(density(data_s0$PY), col = "red", lwd = 2)
legend("topright", c("Full data", "S=0 only"), 
       col = c("black", "red"), lwd = 2)

par(mfrow = c(1, 1))

# Reset plot parameters
par(mfrow = c(1, 1))