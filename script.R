# --- Libraries ---
library(BART)
library(bcf)
library(dplyr)
library(tidyr)
library(purrr)
library(progressr)
library(ggplot2)
# library(ggpubr)
library(gridExtra)
library(MASS)  # For robust linear regression

setwd("~/scratch")
options(future.rng.onMisuse = "ignore")
options(progressr.enable = TRUE)
handlers("progress")

# --- Simulation Parameters ---
N <- 250
REPS <- 50
BCF_NBURN <- 200
BCF_NSIM <- 200

# --- Data Generating Function ---
g <- function(x4) ifelse(x4 == 1, 2, ifelse(x4 == 2, -1, -4))  # Define the function g()

generate_data <- function(n, hetero = FALSE, nonlinear = FALSE) {
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  x4 <- rbinom(n, 1, 0.5) + 1; x5 <- sample(1:3, n, TRUE)
  X <- cbind(x1, x2, x3, x4, x5)
  
  g_val <- g(x4)
  mu <- if (nonlinear) -6 + g_val + 6 * abs(x3 - 1) else 1 + g_val + x1 * x3
  
  s <- sd(mu)
  u <- runif(n)
  pi <- 0.8 * pnorm((3 * mu) / s - 0.5 * x1) + 0.05 + u / 10
  z <- rbinom(n, 1, pi)
  tau <- if (hetero) 1 + 2 * x2 * x5 else rep(3, n)
  
  y <- mu + tau * z + rnorm(n)
  list(X = X, y = y, z = z, mu = mu, tau = tau, pi = pi)
}

## 3.  Evaluate ATE for a single fit
## ---------------------------------------------------------
evaluate_metrics <- function(tau_draws, true_tau) {
  ## tau_draws  :  nsim × n   (bcf output)
  ate_samples <- rowMeans(tau_draws)          # one ATE per posterior draw
  ate_hat     <- mean(ate_samples)            # posterior mean ATE
  ci          <- quantile(ate_samples, c(.025, .975))
  list(
    ate   = ate_hat,
    bias  = ate_hat - mean(true_tau),
    rmse  = sqrt((ate_hat - mean(true_tau))^2),
    cover = (ci[1] <= mean(true_tau) && ci[2] >= mean(true_tau)),
    len   = diff(ci),
    p_val = mean(ate_samples > 0)             # Pr(ATE > 0)
  )
}

# --- Core Simulation per DGP ---
run_simulation <- function(n, hetero, nonlinear) {
  data <- generate_data(n, hetero, nonlinear)
  X <- data$X; y <- data$y; z <- data$z; pi <- data$pi; true_tau <- data$tau
  
  results <- list()
  
  # --- BCF ---
  ps_model <- glm(z ~ ., family = binomial(), data = as.data.frame(X))
  ps_hat <- predict(ps_model, type = "response")
  
  bcf_fit <- bcf(
    y = y, z = z,
    x_control = X, x_moderate = X,
    pihat = ps_hat,
    nburn = 200, nsim = 500,
    ntree_control = 200, ntree_moderate = 50,
    base_control = 0.95, power_control = 2,
    base_moderate = 0.25, power_moderate = 3,
    use_muscale = TRUE, sd_control = 2 * sd(y),
    use_tauscale = TRUE, sd_moderate = sd(y),
    no_output = TRUE
  )
  
  # Extract ATE and CATE
  bcf_ate_samples <- colMeans(bcf_fit$tau)
  bcf_cate <- rowMeans(bcf_fit$tau)
  
  # Evaluate performance metrics
  bcf_eval <- evaluate_metrics(bcf_ate_samples, bcf_cate, true_tau)
  bcf_eval$method <- "BCF"
  results[[length(results) + 1]] <- bcf_eval
  
  # --- BART ---
  bart_fit <- wbart(x.train = cbind(X, z), y.train = y, ndpost = BCF_NSIM, nskip = BCF_NBURN)
  y1 <- predict(bart_fit, newdata = cbind(X, z = 1))
  y0 <- predict(bart_fit, newdata = cbind(X, z = 0))
  bart_ate_samples <- rowMeans(y1 - y0)
  bart_cate <- colMeans(y1 - y0)
  bart_eval <- evaluate_metrics(bart_ate_samples, bart_cate, true_tau)
  bart_eval$method <- "BART"
  results[[length(results) + 1]] <- bart_eval
  
  # --- PCA ---
  pcs <- prcomp(X, scale. = TRUE)$x[, 1:3]
  pca_model <- lm(y ~ ., data = data.frame(z = z, pcs))
  pca_ate <- coef(pca_model)["z"]
  pca_se <- summary(pca_model)$coeff["z", "Std. Error"]
  pca_samples <- rnorm(BCF_NSIM, pca_ate, pca_se)
  pca_eval <- evaluate_metrics(pca_samples, rep(pca_ate, n), true_tau)
  pca_eval$method <- "PCA"
  results[[length(results) + 1]] <- pca_eval
  
  # --- GLM ---
  glm_model <- glm(y ~ ., data = data.frame(y = y, z = z, as.data.frame(X)))
  glm_ate <- coef(glm_model)["z"]
  glm_se <- summary(glm_model)$coeff["z", "Std. Error"]
  glm_samples <- rnorm(BCF_NSIM, glm_ate, glm_se)
  glm_eval <- evaluate_metrics(glm_samples, rep(glm_ate, n), true_tau)
  glm_eval$method <- "GLM"
  results[[length(results) + 1]] <- glm_eval
  
  # --- Robust Linear Regression ---
  rlm_model <- rlm(y ~ ., data = data.frame(z = z, X))
  rlm_ate <- coef(rlm_model)["z"]
  rlm_se <- summary(rlm_model)$coeff["z", "Std. Error"]
  rlm_samples <- rnorm(BCF_NSIM, rlm_ate, rlm_se)
  rlm_eval <- evaluate_metrics(rlm_samples, rep(rlm_ate, n), true_tau)
  rlm_eval$method <- "Robust Linear Regression"
  results[[length(results) + 1]] <- rlm_eval
  
  bind_rows(results)
}

# --- Full Simulation Loop ---
run_all <- function(n, reps = 50) {
  configs <- expand.grid(hetero = c(FALSE, TRUE), nonlinear = c(FALSE, TRUE))
  sim_list <- rep(split(configs, seq(nrow(configs))), each = reps)
  
  with_progress({
    p <- progressor(steps = length(sim_list))
    results <- purrr::map_dfr(sim_list, function(cfg) {
      out <- run_simulation(n, cfg$hetero, cfg$nonlinear)
      out <- mutate(out, hetero = cfg$hetero, nonlinear = cfg$nonlinear, n = n)
      p(); gc(); out
    })
  })
  
  results %>%
    group_by(method, hetero, nonlinear, n) %>%
    summarise(across(c(ate, rmse, bias, p_value, cover, len), mean), .groups = "drop")
}

# --- Run Main ---
set.seed(42)
results_250 <- run_all(250, reps = REPS)
write.csv(results_250, "results_combined.csv", row.names = FALSE)
print(results_250)

# --- Visualizations ---
results <- read.csv("results_combined.csv")
results <- results %>%
  mutate(
    Method = factor(method, levels = c("BCF", "BART", "GLM", "PCA", "Robust Linear Regression")),
    Heterogeneity = ifelse(hetero, "Heterogeneous", "Homogeneous"),
    Model = ifelse(nonlinear, "Nonlinear", "Linear")
  )

# ATE Estimates Boxplot
ggplot(results, aes(x = Method, y = ate, fill = Method)) +
  geom_boxplot() +
  facet_grid(Heterogeneity ~ Model) +
  labs(title = "Distribution of ATE Estimates", y = "ATE", x = "") +
  theme_bw() +
  theme(legend.position = "none")

# CI Coverage Bar Plot
ggplot(results, aes(x = Method, fill = as.factor(cover))) +
  geom_bar(position = "fill") +
  facet_grid(Heterogeneity ~ Model) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "95% CI Coverage", y = "Coverage Rate", fill = "Covers True ATE") +
  theme_minimal()

# RMSE Boxplot
ggplot(results, aes(x = Method, y = rmse, fill = Method)) +
  geom_boxplot() +
  facet_grid(Heterogeneity ~ Model) +
  labs(title = "Root Mean Squared Error (RMSE)", y = "RMSE") +
  theme_light() +
  theme(legend.position = "none")

# Confidence Interval Length Boxplot
ggplot(results, aes(x = Method, y = len, fill = Method)) +
  geom_boxplot() +
  facet_grid(Heterogeneity ~ Model) +
  labs(title = "Length of 95% Confidence Intervals", y = "Interval Length") +
  theme_minimal()

# Performance Trade-offs (Bias vs RMSE vs CI Length)
summary_data <- results %>%
  group_by(Method, Heterogeneity, Model) %>%
  summarise(
    RMSE = mean(rmse),
    Bias = mean(bias),
    Len = mean(len),
    .groups = 'drop'
  )

ggplot(summary_data, aes(x = Bias, y = RMSE, size = Len, color = Method)) +
  geom_point(alpha = 0.7) +
  facet_grid(Heterogeneity ~ Model) +
  labs(title = "Performance Trade-offs: Bias vs RMSE vs CI Length",
       x = "Bias", y = "RMSE", size = "CI Length") +
  theme_classic()

# Save the plots
ggsave("plot_ate_boxplot.png", width = 8, height = 6)
ggsave("plot_ci_coverage.png", width = 8, height = 6)
ggsave("plot_rmse.png", width = 8, height = 6)
ggsave("plot_ci_length.png", width = 8, height = 6)

# CATE RMSE (PEHE) Comparison Boxplot
ggplot(results, aes(x = method, y = cate_rmse, fill = method)) +
  geom_boxplot() +
  facet_grid