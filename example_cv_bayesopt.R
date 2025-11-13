# Example: Using desurv_cv_bayesopt() to optimize hyperparameters
# This script demonstrates how to use Bayesian Optimization to efficiently
# tune hyperparameters for desurv_cv()

library(DeSurv)

# Load your data
data <- readRDS("/home/naimrashid/Downloads/data.rds")
X <- data$ex
y <- data$sampInfo$time
d <- data$sampInfo$event
dataset <- data$sampInfo$dataset
samp_keeps <- data$samp_keeps

# Define the hyperparameter search space
# You can customize which parameters to include
bo_bounds <- desurv_cv_bo_default_bounds(
  include_factor_penalties = TRUE,   # Include lambdaW_grid and lambdaH_grid
  include_n_starts = TRUE,            # Include n_starts
  include_ngene = TRUE                # Include ngene (for preprocessing)
)

# Alternatively, define custom bounds:
# bo_bounds <- list(
#   k_grid = list(lower = 2, upper = 12, type = "integer"),
#   alpha_grid = list(lower = 0, upper = 0.95, type = "continuous"),
#   lambda_grid = list(lower = 1e-5, upper = 1e5, scale = "log10"),
#   nu_grid = list(lower = 0, upper = 1, type = "continuous"),
#   lambdaW_grid = list(lower = 1e-5, upper = 1e5, scale = "log10"),
#   lambdaH_grid = list(lower = 1e-5, upper = 1e5, scale = "log10"),
#   n_starts = list(lower = 10, upper = 100, type = "integer"),
#   ngene = list(lower = 1000, upper = 5000, type = "integer")
# )

# Run Bayesian Optimization
# Adjust n_init and n_iter based on your computational budget
bo_result <- desurv_cv_bayesopt(
  X = X,
  y = y,
  d = d,
  dataset = dataset,
  samp_keeps = samp_keeps,
  preprocess = TRUE,
  method_trans_train = "rank",
  engine = "cold",              # Use "warmstart" for faster alpha path search
  # Fixed CV/optimization parameters
  nfolds = 5,
  tol = 1e-4,
  maxit = 100,
  # BO parameters
  bo_bounds = bo_bounds,
  n_init = 15,                  # Initial LHS samples (default: max(5, 3*d))
  n_iter = 20,                  # Sequential BO iterations
  exploration_weight = 0.01,    # Exploration vs exploitation tradeoff
  seed = 2025,
  verbose = TRUE,
  cv_verbose = FALSE            # Set to TRUE to see CV progress
)

# View results
print(bo_result)

# Access best parameters
best_params <- bo_result$best$params
print("Best hyperparameters:")
print(best_params)

# Access full history
history <- bo_result$history
print("First few evaluations:")
print(head(history))

# Plot optimization progress (if you have ggplot2)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  history_ok <- history[history$status == "ok", ]

  p <- ggplot(history_ok, aes(x = eval_id, y = mean_cindex)) +
    geom_point(aes(color = stage)) +
    geom_line(alpha = 0.3) +
    geom_hline(yintercept = max(history_ok$mean_cindex),
               linetype = "dashed", color = "red") +
    labs(
      title = "Bayesian Optimization Progress",
      x = "Evaluation",
      y = "Mean C-index",
      color = "Stage"
    ) +
    theme_minimal()

  print(p)
  ggsave("bo_progress.pdf", p, width = 8, height = 5)
}

# Refit final model with best parameters
# Option 1: Use desurv_cv() with best parameters to get CV results + final model
final_cv <- desurv_cv(
  X = X,
  y = y,
  d = d,
  dataset = dataset,
  samp_keeps = samp_keeps,
  preprocess = TRUE,
  method_trans_train = "rank",
  k_grid = best_params["k_grid"],
  alpha_grid = best_params["alpha_grid"],
  lambda_grid = best_params["lambda_grid"],
  nu_grid = best_params["nu_grid"],
  lambdaW_grid = best_params["lambdaW_grid"],
  lambdaH_grid = best_params["lambdaH_grid"],
  n_starts = best_params["n_starts"],
  ngene = best_params["ngene"],
  nfolds = 5,
  tol = 1e-4,
  maxit = 100,
  engine = "cold",
  cv_only = FALSE,  # Fit final model
  verbose = TRUE
)

# Option 2: Just fit final model directly without CV
# (if you want to train on full dataset with best params)
# processed_data <- preprocess_data(
#   X = X,
#   dataset = dataset,
#   samp_keeps = samp_keeps,
#   ngene = best_params["ngene"],
#   method_trans_train = "rank"
# )
#
# final_fit <- desurv_fit(
#   X = processed_data$X,
#   y = y,
#   d = d,
#   k = best_params["k_grid"],
#   alpha = best_params["alpha_grid"],
#   lambda = best_params["lambda_grid"],
#   nu = best_params["nu_grid"],
#   lambdaW = best_params["lambdaW_grid"],
#   lambdaH = best_params["lambdaH_grid"],
#   n_starts = best_params["n_starts"],
#   tol = 1e-4,
#   maxit = 100,
#   verbose = TRUE
# )

# Save results
saveRDS(bo_result, "bo_optimization_results.rds")
saveRDS(final_cv, "final_cv_model.rds")

cat("\nBayesian Optimization complete!\n")
cat("Best C-index:", bo_result$best$mean_cindex, "\n")
cat("Total evaluations:", nrow(bo_result$history), "\n")
cat("Results saved to bo_optimization_results.rds and final_cv_model.rds\n")
