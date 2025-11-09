# Load required packages
library(haven)
library(ggplot2)
library(stats)
library(dplyr)
library(mfp)
library(lspline)
library(rms)
library(broom)

# Import dataset
data <- read_dta("vitd.dta")

# Summary statistics
summary(data$vitd)

# Plot histogram of Vitamin D levels
ggplot(data, aes(x = vitd)) + 
  geom_histogram(aes(y = after_stat(density)), binwidth = 5, fill = "blue", color = "black") +
  labs(title = "Histogram of Vitamin D", x = "Vitamin D", y = "Density") +
  theme_minimal()

# Generate time variable (t)
data <- data %>%
  mutate(t = as.numeric(date_test - as.Date(paste(year, "01", "01", sep = "-"))) / 365.25)

# Scatter plot of Vitamin D vs time
ggplot(data, aes(x = t, y = vitd)) + 
  geom_point(color = "blue", alpha = 0.6) +
  labs(title = "Scatter plot of Vitamin D over Time", x = "Time of Year (t)", y = "Vitamin D") +
  theme_minimal()

# Filter sample 1
data_1 <- data %>%
  filter(sample == 1)

# LOWESS smoothing with different bandwidths
for (f in c(0.1, 0.5, 0.9)) {
  lowess_fit <- lowess(data_1$t, data_1$vitd, f = f)
  predicted <- data.frame(t = lowess_fit$x, vitd = lowess_fit$y)
  ggplot(data_1, aes(x = t, y = vitd)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_line(data = predicted, aes(x = t, y = vitd), color = "red", linewidth = 1) +
    labs(title = paste("LOWESS (Bandwidth =", f, ")"), x = "Time (t)", y = "Vitamin D") +
    theme_minimal()
}

# Quadratic regression
quad_model <- lm(vitd ~ t + I(t^2), data = data)
data$yhat_quad <- predict(quad_model)
summary(quad_model)

# Plot quadratic regression
ggplot(data, aes(x = t, y = vitd)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_line(aes(y = yhat_quad), color = "red", linewidth = 1) +
  labs(title = "Quadratic Regression", x = "Time (t)", y = "Vitamin D") +
  theme_minimal()

# Cubic regression
cubic_model <- lm(vitd ~ t + I(t^2) + I(t^3), data = data)
data$yhat_cub <- predict(cubic_model)
summary(cubic_model)

# Plot cubic regression
ggplot(data, aes(x = t, y = vitd)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_line(aes(y = yhat_cub), color = "red", linewidth = 1) +
  labs(title = "Cubic Regression", x = "Time (t)", y = "Vitamin D") +
  theme_minimal()

# Fractional polynomial regression
fp_model <- mfp(vitd ~ fp(t), data = data)
summary(fp_model)
data$yhat_fp <- predict(fp_model)

# Plot fractional polynomial fit
ggplot(data, aes(x = t, y = vitd)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_line(aes(y = yhat_fp), color = "red", linewidth = 1) +
  labs(title = "Fractional Polynomial Regression", x = "Time (t)", y = "Vitamin D") +
  theme_minimal()

# Linear spline model
spline_model_linear <- lm(vitd ~ lspline(t, knots = c(0.2, 0.4, 0.6, 0.8)), data = data)
data$yhat_spline_lin <- predict(spline_model_linear)
summary(spline_model_linear)

# Plot linear spline
ggplot(data, aes(x = t, y = vitd)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_line(aes(y = yhat_spline_lin), color = "red", linewidth = 1) +
  labs(title = "Linear Spline Regression", x = "Time (t)", y = "Vitamin D") +
  theme_minimal()

# Cubic spline model (rms)
ddist <- datadist(data)
options(datadist = "ddist")
spline_model_cubic <- ols(vitd ~ rcs(t, 4), data = data)
data$yhat_spline_cub <- predict(spline_model_cubic)

# Plot cubic spline
ggplot(data, aes(x = t, y = vitd)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_line(aes(y = yhat_spline_cub), color = "red", linewidth = 1) +
  labs(title = "Cubic Spline Regression", x = "Time (t)", y = "Vitamin D") +
  theme_minimal()

# Sinusoidal regression
data <- data %>%
  mutate(cos = cos(2 * pi * t),
         sin = sin(2 * pi * t))

sin_cos_model <- lm(vitd ~ cos + sin, data = data)
data$yhat_sin <- predict(sin_cos_model)
summary(sin_cos_model)

# Plot sinusoidal model
ggplot(data, aes(x = t, y = vitd)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_line(aes(y = yhat_sin), color = "red", linewidth = 1) +
  labs(title = "Sinusoidal Regression", x = "Time (t)", y = "Vitamin D") +
  theme_minimal()

# Logistic regression analysis
data <- read_dta("bariatric_surgery.dta")

logit_model <- glm(bariatric ~ bmi + age + factor(education) + sex + CVD + diuretics, 
                   data = data, family = binomial(link = "logit"))
summary(logit_model)

# Compute and display odds ratios
tidy_model <- tidy(logit_model, exponentiate = TRUE, conf.int = TRUE)
print(tidy_model)

# Predicted probabilities
data$p_bariatric <- predict(logit_model, type = "response")
