% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesianModel-inators.R
\name{TauInator}
\alias{TauInator}
\title{TauInator}
\usage{
TauInator(
  data,
  feature_col,
  response_col,
  group_col,
  covariate_cols = NULL,
  family = "gaussian",
  prior_vector = c(set_prior("normal(0, 2", class = "b"), set_prior("normal(0, 2", class
    = "Intercept")),
  iter = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  offset_term = NULL,
  covariate_whitelist = NULL,
  make_models_directory = TRUE,
  ...
)
}
\arguments{
\item{data}{A data frame containing the data to be analyzed.}

\item{feature_col}{A string specifying the name of the column in \code{data} that contains the feature identifiers (e.g., gene names).}

\item{response_col}{A string specifying the name of the column in \code{data} that contains the response variable (e.g., expression levels).}

\item{group_col}{A string specifying the name of the column in \code{data} that contains the grouping variable for random intercepts (e.g., sample IDs).}

\item{covariate_cols}{A vector of strings specifying the names of the columns in \code{data} that contain covariates to be included in the model.}

\item{family}{A string specifying the family of the response variable (default: "gaussian").}

\item{iter}{An integer specifying the number of iterations for the Bayesian model fitting (default: 2000).}

\item{chains}{An integer specifying the number of chains for the Bayesian model fitting (default: 4).}

\item{cores}{An integer specifying the number of CPU cores to use for parallel processing (default: 4).}

\item{seed}{An integer specifying the random seed for reproducibility (default: 123).}

\item{offset_term}{A string specifying the name of the column in \code{data} to be used as an offset term in the model (default: NULL).}

\item{covariate_whitelist}{A vector of strings specifying which covariates to include in the tau calculation. If NULL, all covariates will be used (default: NULL).}

\item{make_models_directory}{Logical, if TRUE (default) creates a directory called "brms_models" in the working directory to store model files.}

\item{...}{Additional arguments to be passed to the \code{brm} function from the \code{brms} package.}

\item{prior_list}{A list of prior distributions to be used in the Bayesian model.}
}
\value{
A list containing tau values for each covariate specified in \code{covariate_whitelist}.
}
\description{
Calculates population level variation for single-feature models under an empirical Bayesian framework, followed by fitting a full random slopes model for different features in a Bayesian model using the \code{brms} package.
}
\examples{
 # Example usage:
 # create fake `df`, a data frame with columns: "gene", "expression", "sample", "covariate1", "covariate2"
 df <- data.frame(
   gene = rep(paste0("gene", 1:10), each = 20),
   expression = unlist(map(.f = function(x) rnorm(mean = x, n = 20), sort(1:10))),
   sample = rep(paste0("sample", 1:20), times = 10),
   covariate1 = unlist(map(.f = function(x) rnorm(mean = x, n = 20), sort(1:10))),
   covariate2 = unlist(map(.f = function(x) rnorm(mean = x, n = 20), sort(1:10)))
 )
prior_vector = c(set_prior("normal(0, 2", class = "b"), set_prior("normal(0, 2", class = "Intercept"))
results <- TauInator(data = df, feature_col = "gene", response_col = "expression",
                     group_col = "sample", covariate_cols = c("covariate1", "covariate2"),
                     prior_vector = prior_list, iter = 2000, chains = 4, cores = 4, seed = 123,
                     covariate_whitelist = c("covariate1","covariate2"))


}
