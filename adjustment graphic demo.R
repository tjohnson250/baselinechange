set.seed(42)
n <- 3000

# Z (confounder) affects both X and Y
z <- rnorm(n, mean = 50, sd = 10)

# X (exposure) is influenced by Z
x <- 0.7 * z + rnorm(n, mean = 0, sd = 1)

# Y (outcome) is also influenced by Z, and has a true causal effect from X
# We'll make the true causal effect of X on Y relatively small here to highlight confounding
y <- 0.5 * z + 0.1 * x + rnorm(n, mean = 0, sd = 1)

# Combine into a data frame
confounded_data <- data.frame(x, y, z)

# Plot the y vs x without adjusting for z
ggplot(confounded_data, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Relationship between X and Y (Unadjusted for Z)",
       x = "X (without controlling for Z)",
       y = "Y (without controlling for Z)") +
  theme_minimal()

# Calculate residuals after regressing X on Z
residuals_x_on_z <- residuals(lm(x ~ z, data = confounded_data))

# Calculate residuals after regressing Y on Z
residuals_y_on_z <- residuals(lm(y ~ z, data = confounded_data))

# Create a data frame with the residuals
adjusted_data <- data.frame(residuals_x_on_z, residuals_y_on_z)

# Plot the residuals against each other (adjusted for Z)
ggplot(adjusted_data, aes(x = residuals_x_on_z, y = residuals_y_on_z)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Relationship between X and Y (Adjusted for Z)",
       x = "Residuals of X (after controlling for Z)",
       y = "Residuals of Y (after controlling for Z)") +
  theme_minimal()