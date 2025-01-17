test_that("generate_rules works as expected!", {
  # Generate sample data
  set.seed(3784)
  dataset_cont <- generate_cre_dataset(n = 300, rho = 0, n_rules = 2, p = 10,
                                       effect_size = 0.5,
                                       binary_outcome = FALSE)
  y <- dataset_cont[["y"]]
  z <- dataset_cont[["z"]]
  X <- dataset_cont[["X"]]


  ite_method <- "bart"
  include_ps <- "TRUE"
  ps_method <- "SL.xgboost"
  or_method <- NA
  ntrees_rf <- 100
  ntrees_gbm <- 50
  node_size <- 20
  max_nodes <- 5
  max_depth <- 15
  replace <- FALSE
  intervention_vars <- NULL

  # Step 1: Split data
  X <- as.matrix(X)
  y <- as.matrix(y)
  z <- as.matrix(z)

  # Step 2: Estimate ITE
  ite <- estimate_ite(y, z, X, ite_method,
                           include_ps = include_ps,
                           ps_method = ps_method,
                           or_method = or_method)

  # Correct outputs
  rules <- generate_rules(X, ite, intervention_vars, ntrees_rf, ntrees_gbm,
                          node_size, max_nodes, max_depth, replace)

  expect_true(class(rules) == "character")



  ite_method <- "bart"
  include_ps <- "TRUE"
  ps_method <- "SL.xgboost"
  or_method <- NA
  ntrees_rf <- 100
  ntrees_gbm <- 0
  node_size <- 20
  max_nodes <- 5
  max_depth <- 15
  replace <- FALSE
  intervention_vars <- c()

  # Check for binary outcome
  binary <- ifelse(length(unique(y)) == 2, TRUE, FALSE)

  # Step 1: Split data
  X <- as.matrix(X)
  y <- as.matrix(y)
  z <- as.matrix(z)

  # Step 2: Estimate ITE
  ite <- estimate_ite(y, z, X, ite_method,
                           include_ps = include_ps,
                           ps_method = ps_method,
                           or_method = or_method)

  # Correct outputs
  rules <- generate_rules(X, ite, intervention_vars, ntrees_rf, ntrees_gbm,
                          node_size, max_nodes, max_depth, replace)

  expect_true(class(rules) == "character")

})
