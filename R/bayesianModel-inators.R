#' @title TauInator
#'
#' @description Calculates population level variation for single-feature models under an empirical Bayesian framework, followed by fitting a full random slopes model for different features in a Bayesian model using the `brms` package.
#'
#' @param data A data frame containing the data to be analyzed.
#' @param feature_col A string specifying the name of the column in `data` that contains the feature identifiers (e.g., gene names).
#' @param response_col A string specifying the name of the column in `data` that contains the response variable (e.g., expression levels).
#' @param group_col A string specifying the name of the column in `data` that contains the grouping variable for random intercepts (e.g., sample IDs).
#' @param covariate_cols A vector of strings specifying the names of the columns in `data` that contain covariates to be included in the model.
#' @param prior_list A list of prior distributions to be used in the Bayesian model.
#' @param iter An integer specifying the number of iterations for the Bayesian model fitting (default: 2000).
#' @param chains An integer specifying the number of chains for the Bayesian model fitting (default: 4).
#' @param cores An integer specifying the number of CPU cores to use for parallel processing (default: 4).
#' @param seed An integer specifying the random seed for reproducibility (default: 123).
#' @param ... Additional arguments to be passed to the `brm` function from the `brms` package.
#' @param family A string specifying the family of the response variable (default: "gaussian").
#' @param offset_term A string specifying the name of the column in `data` to be used as an offset term in the model (default: NULL).
#' @param covariate_whitelist A vector of strings specifying which covariates to include in the tau calculation. If NULL, all covariates will be used (default: NULL).
#' @param make_models_directory Logical, if TRUE (default) creates a directory called "brms_models" in the working directory to store model files.
#' @return A list containing tau values for each covariate specified in `covariate_whitelist`.
#' @import brms
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @importFrom stats as.formula
#' @export
#' @examples
#'  # Example usage:
#'  # create fake `df`, a data frame with columns: "gene", "expression", "sample", "covariate1", "covariate2"
#'  df <- data.frame(
#'    gene = rep(paste0("gene", 1:10), each = 20),
#'    expression = unlist(map(.f = function(x) rnorm(mean = x, n = 20), sort(1:10))),
#'    sample = rep(paste0("sample", 1:20), times = 10),
#'    covariate1 = unlist(map(.f = function(x) rnorm(mean = x, n = 20), sort(1:10))),
#'    covariate2 = unlist(map(.f = function(x) rnorm(mean = x, n = 20), sort(1:10)))
#'  )
#' prior_vector = c(set_prior("normal(0, 2", class = "b"), set_prior("normal(0, 2", class = "Intercept"))
#' results <- TauInator(data = df, feature_col = "gene", response_col = "expression",
#'                      group_col = "sample", covariate_cols = c("covariate1", "covariate2"),
#'                      prior_vector = prior_list, iter = 2000, chains = 4, cores = 4, seed = 123,
#'                      covariate_whitelist = c("covariate1","covariate2"))
#'
#'

TauInator <- function(data, feature_col, response_col, group_col, covariate_cols = NULL,
                      family = "gaussian",
                      prior_vector = c(set_prior("normal(0, 2", class = "b"),
                                       set_prior("normal(0, 2", class = "Intercept")),
                      iter = 2000, chains = 4, cores = 4, seed = 123,
                      offset_term = NULL,
                      covariate_whitelist = NULL,
                      make_models_directory = TRUE,
                      ...) {
  # Ensure required packages are loaded
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("The 'brms' package is required but not installed. Please install it.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("The 'dplyr' package is required but not installed. Please install it.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("The 'tidyr' package is required but not installed. Please install it.")
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("The 'purrr' package is required but not installed. Please install it.")
  }

  #make a directory for the models if it doesn't already exist
  if (make_models_directory) {
    if (!dir.exists(here::here("brms_models"))) {
      dir.create(here::here("brms_models"))
    }
  }
  # Convert feature_col, response_col, group_col, and covariate_cols to symbols
  feature_sym <- rlang::sym(feature_col)
  response_sym <- rlang::sym(response_col)
  group_sym <- rlang::sym(group_col)
  covariate_syms <- rlang::syms(covariate_cols)
  if (!is.null(offset_term)){
    offset_sym <- rlang::sym(offset_term)
  }
  # Nest data by feature
  nested_data <- data %>%
    dplyr::group_by(!!feature_sym) %>%
    tidyr::nest()


  #create formula using some common modeling approaches.
  #generally, this can probably be simplistic, as a more complicated model will be fit downstream
  if (!is.null(covariate_cols) & !is.null(offset_term)){
    formula_string <- paste(response_col, "~", paste(covariate_cols, collapse = " + "), "+ (1 |", group_col, ") + offset(", offset_term, ")")
  } else if (!is.null(covariate_cols) & is.null(offset_term)){
    formula_string <- paste(response_col, "~", paste(covariate_cols, collapse = " + "), "+ (1 |", group_col, ")")
  } else if (is.null(covariate_cols) & !is.null(offset_term)){
    formula_string <- paste(response_col, "~ 1 + (1 |", group_col, ") + offset(", offset_term, ")")
  } else {
    formula_string <- paste(response_col, "~ 1 + (1 |", group_col, ")")
  }

  #loop to fit a bayesian model for each feature using brms
  model_results <- nested_data %>%
    dplyr::group_by(!!feature_sym) %>%
    dplyr::do({
      cur_group <- as.character(.data[1,1][[1]])
      file_name <- ifelse(make_models_directory,
                     here::here("brms_models", paste0(as.character(cur_group), "_brms_model.rds")),
                     NULL)
      fit <- brms::brm(formula_string ,
                       data = .data[1,2][[1]],
                       family = family,
                       prior = c(set_prior("normal(0, 2", class = "b"),
                                 set_prior("normal(0, 2", class = "Intercept")
                       ),
                       chains = 4,
                       cores = 4,
                       backend = "cmdstanr",
                       file = file_name,
                       file_refit = 'on_change')
      # extract fixed effects estimates and standard errors for a subset of beta coefficients
      if (!is.null(covariate_whitelist)){
        print("Using covariate whitelist for tau")
        #check covariates
        for (i in 1:length(covariate_whitelist)){
          if (!covariate_whitelist[i] %in% rownames(fixef(fit))){
            stop(paste0("Covariate ", covariate_whitelist[i], " not found in model fixed effects. Check spelling?"))
          }
        }
        betas <- data.frame()
        #loop over covariates
        for (covariate in covariate_whitelist) {
          betas_tmp <- fixef(fit)[covariate, ]
          betas_tmp <- data.frame(t(betas_tmp))
          betas_tmp$feature <- cur_group
          betas_tmp$covariate <- covariate
          betas <- rbind(betas, betas_tmp)
        }
        #otherwise, get all betas
      } else {
        print("Pooling across all covariates for tau. (Do you really want this?)")
        betas <- data.frame(fixef(fit))
        betas$feature <- cur_group
        betas$covariate <- rownames(betas)
      }
      #return betas
      betas
    })
  #use ashr to estimate tau (beta coefficient variation across features) from the fixed effects estimates and standard errors
  #multi-beta case
  if (!is.null (covariate_whitelist) & length(covariate_whitelist) > 1){
    model_results <- model_results %>%
      dplyr::filter(covariate %in% covariate_whitelist)
    #initialize a tau_df to append to
    tau_df <- data.frame(covariate = character(),
                         tau = numeric())

    for (covariate_name in covariate_whitelist) {
      model_results_subset <- model_results %>%
        dplyr::filter(covariate == covariate_name)
      tau_df_subset <- .fit_ashr_model(model_results_subset, covariate)
      tau_df <- rbind(tau_df, tau_df_subset)
    }
    #case when there's just one covariate to care about (probably the typical case)
  } else if (!is.null (covariate_whitelist) & length(covariate_whitelist) == 1){
    tau_df <- .fit_ashr_model(model_results, covariate_whitelist)
    #case when you want to pool across all covariates (probably not a good idea)
  } else {
    tau_df <- .fit_ashr_model(model_results, "all")
  }
  return(tau_df)
}


.fit_ashr_model <- function(model_results, covariate_name) {
  ash_fit <- ashr::ash(betahat = model_results$Estimate,
                       sebetahat = model_results$Est.Error,
                       method = "fdr",
                       mixcompdist = 'normal')
  g <- ashr::get_fitted_g(ash_fit)
  means <- g$mean
  sds <- g$sd
  weights <- g$pi
  tau2_hat <- sum(weights * (sds^2 + means^2)) - (sum(weights * means))^2
  tau_hat <- sqrt(tau2_hat)
  return(data.frame(tau = tau_hat, covariate = covariate_name))
}


