#' Truncated Infinite Factor Analysis Imputation Model
#'
#' @description
#' Applies the tIFA imputation method to the input dataset.
#'
#' @param input_data The dataset with missing value requiring imputation in matrix format.
#' @param coding The coding used to indicate a data entry is missing. Defaults to `NA`.
#' @param n.iters The number of MCMC iterations desired. Defaults to 10000.
#' @param k.star Parameter \eqn{k^{*}} in the tIFA model. Defaults to 5.
#' @param verbose Is a readout desired? Defaults to `TRUE`.
#' @param burn The number of MCMC iterations to be discarded in an initial burn. Defaults to 5000.
#' @param thin The level of thinning to take place in the MCMC chain. Defaults to 5, meaning only every fifth draw will be kept.
#' @param mu_varphi Parameter \eqn{\varphi} in the tIFA model. Defaults to 0.1.
#' @param kappa_1 Parameter \eqn{\kappa_1} in the tIFA model. Defaults to 3.
#' @param kappa_2 Parameter \eqn{\kappa_2} in the tIFA model. Defaults to 2.
#' @param a_sigma Parameter \eqn{a_{\sigma}} in the tIFA model. Defaults to 1.
#' @param b_sigma Parameter \eqn{b_{\sigma}} in the tIFA model. Defaults to 0.3.
#' @param a_1 Parameter \eqn{a_1} in the tIFA model. Defaults to 2.1.
#' @param a_2 Parameter \eqn{a_2} in the tIFA model. Defaults to 3.1.
#' @param return_chain Should resulting draws from the MCMC chain be returned? Defaults to `FALSE`.
#'
#' @return Returns a list containing entries:
#'
#' * imputed_dataset: The imputed dataset in matrix format.
#'
#' * imputation_info: A dataframe containing information on imputed values and associated uncertainty for each missing entry in the original dataset.
#'
#' * chain: Returned only if return_chain = `TRUE`. Provides all draws from the MCMC chain and further information.
#'
#'     * store_data: MCMC draws of the imputed data entries.
#'
#'     * store_alpha: MCMC draws of the \eqn{\alpha} parameter.
#'
#'     * store_sigma_inv: MCMC draws of the \eqn{\sigma^{-2}} entries.
#'
#'     * store_lambda: MCMC draws of the loadings matrix.
#'
#'     * store_eta: MCMC draws of the latent factor scores.
#'
#'     * store_mu: MCMC draws of the mean vector.
#'
#'     * store_k_t: MCMC draws of the effective factor number (unchanging).
#'
#'     * store_phi: MCMC draws of the local shrinkage parameters.
#'
#'     * store_delta: MCMC draws of the multiplicative gamma process parameters.
#'
#'     * store_b_sigma: MCMC draws of the \eqn{b_{\sigma}} parameter (unchanging).
#'
#'     * store_Z: MCMC draws of the missingness designation indicator for each imputed entry.
#'
#'     * input_data: The dataset input into the tIFA model.
#'
#'     * input_missingness: The missingness pattern input into the tIFA model.
#'
#' @importFrom Rfast transpose mat.mult colsums
#' @importFrom TruncatedNormal rtmvnorm
#' @importFrom truncdist rtrunc
#' @importFrom stats rgamma rbeta rbinom runif prcomp cov
#' @importFrom truncnorm ptruncnorm rtruncnorm
#' @importFrom Rcpp cppFunction
#' @importFrom softImpute softImpute
#'
#' @export
#'
#' @examples
#' # generate example data
#' example_data <- matrix(abs(rnorm(100)), nrow = 5)
#'
#' # add missingness to example data coded as 0.001
#' example_data[4, 2] <- 0.001
#' example_data[2, 18] <- 0.001
#'
#' # run tIFA model
#' # short chain for example
#' res <- tIFA_model(input_data = example_data, coding = 0.001, n.iters = 100, k.star = 3,
#'                   burn = 50, thin = 5)
#'
#' # checkout imputed dataset
#' res$imputed_dataset
#'
#' # checkout further information on imputed entries
#' res$imputation_info
#'
tIFA_model <- function(input_data, coding = NA, n.iters = 10000, k.star = 5,
                       verbose = TRUE, burn = 5000, thin = 5,
                       mu_varphi = 0.1, kappa_1 = 3L, kappa_2 = 2L,
                       a_sigma = 1L, b_sigma = 0.3, a_1 = 2.1, a_2 = 3.1,
                       return_chain = FALSE, verbose_frequency = 100) {

  Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }", depends = "RcppArmadillo")

  # code data as expected
  # extract missingness information
  prepped_data <- tifa_prep(input_data, coding = coding)

  missingness <- prepped_data$missingness_pattern

  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~
  # initialisation
  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~

  # data dimensions
  n <- nrow(prepped_data$data)
  p <- ncol(prepped_data$data)

  # seeding function for rcpp functions
  get_seed <- function() {

    sample.int(.Machine$integer.max, 1)

  }

  # truncation point
  trunc.point <- min(prepped_data$data, na.rm = TRUE)

  # indices of missing data
  # col 1 contains row indices in data
  # col 2 contains col indices in data
  # so miss_ind[1, ] are the indices of the first missing point in the data
  miss_ind <- which(missingness == 0, arr.ind = TRUE)
  miss_vec <- which(missingness == 0)

  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~
  # def starting values
  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~

  # initialise imputation with svd imputation
  start_svd <- softImpute(prepped_data$data, rank.max = k.star, type = "svd")

  start_imp <- start_svd$u %*% diag(start_svd$d) %*% t(start_svd$v)

  initial_dataset <- prepped_data$data
  initial_dataset[is.na(initial_dataset)] <- start_imp[is.na(initial_dataset)]

  # data for tIFA
  data_working <- initial_dataset

  # initialise remaining parameters
  mu <- matrix(apply(prepped_data$data, 2, mean, na.rm = TRUE))
  mu_tilde <- mu
  alpha <- runif(1, 0, 1)

  # use pca to initialise factors and loadings
  pca_res <- prcomp(data_working, scale. = FALSE, center = FALSE)

  lambda <- pca_res$rotation[ , 1:k.star]  # p x k.star
  rownames(lambda) <- NULL
  colnames(lambda) <- NULL
  lambda[which(lambda < 0)] <- abs(lambda[which(lambda < 0)])

  eta <- matrix(NA, nrow = n, ncol = k.star)

  for (i in 1:n) {

    eta[i , ] <- TruncatedNormal::rtmvnorm(n = 1,
                                           mu = rep(0, k.star),
                                           sigma = diag(k.star),
                                           lb = rep(0, k.star))

  }  # n x k.star

  eps <- matrix(NA, nrow = n, ncol = p)

  for (i in 1:n) {

    eps[i, ] <- data_working[i, ] - mu - (lambda %*% matrix(eta[i ,]))

  }

  sigma_inv <- matrix(1 / diag(cov(eps, use = "complete")))

  Sigma_inv <- diag(sigma_inv[ , 1])

  phi <- matrix(rgamma(n = p * k.star, shape = kappa_1, rate = kappa_2),  # p x k.star
                nrow = p)

  delta <- matrix(rep(0, k.star), nrow = k.star)  # k.star x 1
  delta[1, 1] <- rgamma(n = 1, shape = a_1, rate = 1)
  delta[-1, 1] <- truncdist::rtrunc(spec = "gamma",
                                    n = k.star - 1,
                                    shape = a_2,
                                    rate = 1,
                                    a = 1)

  tau <- matrix(cumprod(delta))  # k.star x 1

  rm(prepped_data)

  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~
  # set up storage
  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~

  # imputed dataset
  store_data <- list()

  # alpha - prob of MAR
  store_alpha <- list()

  # sigma
  store_sigma_inv <- list()

  # lambda - factor loadings
  store_lambda <- list()

  # eta
  store_eta <- list()

  # mu - mean
  store_mu <- list()

  # effective factors
  store_k_t <- list()

  # phi
  store_phi <- list()

  # delta
  store_delta <- list()

  # Z
  store_Z <- list()

  for (m in seq_len(n.iters)) {

    # create a readout
    if (verbose) {
      if (m %% verbose_frequency == 0) {
        statement <- paste("tIFA process running. Now on iteration ", m, " of ", n.iters)
        print(statement)
      }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # lambda update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    eta_t_eta <- mat.mult(Rfast::transpose(eta) , eta)  # k.star x k.star

    for (j in 1:p) {

      if (k.star == 1) {

        lambda_var <- armaInv(sigma_inv[j] * eta_t_eta + diag(c(phi[j, ] * tau), nrow = 1))  # k.star x k.star

      } else {

        lambda_var <- armaInv(sigma_inv[j] * eta_t_eta + diag(c(phi[j, ] * tau)))  # k.star x k.star

      }

      lambda_mean <- mat.mult(lambda_var, mat.mult(
        Rfast::transpose(eta), as.matrix(data_working[ , j] - mu[j])) * sigma_inv[j])  # k.star x 1


      lambda[j, ] <- TruncatedNormal::rtmvnorm(n = 1, mu = as.vector(lambda_mean),
                                               sigma = lambda_var, lb = rep(0, k.star))

    }

    # set NaN to small number
    # as approach zero from above
    lambda[which(lambda == "NaN" | lambda == Inf)] <- 1e-8

    # lambda is a p x k.star matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # eta update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    lambda_t_Sigma <- mat.mult(Rfast::transpose(lambda), Sigma_inv)

    eta_var <- armaInv(diag(k.star) + mat.mult(lambda_t_Sigma, lambda))  # k.star x k.star

    for (i in 1:n) {

      eta_mean <- mat.mult(eta_var, mat.mult(lambda_t_Sigma, as.matrix(data_working[i, ] - mu)))  # k.star x 1

      eta[i , ] <-  TruncatedNormal::rtmvnorm(n = 1, mu = as.vector(eta_mean),
                                              sigma = eta_var, lb = rep(0, k.star))

    }

    # in case of computation difficulty on first iter
    eta[which(eta == "NaN" | eta == Inf)] <- 0

    # eta is a n x k.star matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # sigma_inv update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    lambda_eta_t <- mat.mult(lambda, Rfast::transpose(eta))  # p x n

    sigma_inv <- matrix(rgamma(p, shape = a_sigma + n/2,
                               rate = b_sigma + 0.5 * colsums((sweep(data_working, 2, mu, FUN = "-") -
                                                                 Rfast::transpose(lambda_eta_t))^2)))

    # rate correct here as per relationship between Ga with shape and rate and IG with shape and scale

    Sigma_inv <- diag(sigma_inv[ , 1])

    # sigma_inv is a 100 x 1 matrix
    # Sigma_inv is a 100 x 100 diag matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # mu update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    mu_temp1 <- mu_varphi * diag(p)

    mu_var <- diag(1 / (mu_varphi + n * sigma_inv[ , 1]))

    # note the colsum here gives a large negative value
    mu_temp2 <- mat.mult(Sigma_inv , matrix(colsums(data_working - mat.mult(eta, t(lambda)))))

    mu_mean <- mat.mult(mu_var, mu_temp2 + mat.mult(mu_temp1, mu_tilde))

    mu <- TruncatedNormal::rtmvnorm(n = 1, mu = as.vector(mu_mean),
                                    sigma = mu_var, lb = rep(0, p))


    mu <- as.matrix(mu) # p x 1

    rm(mu_temp1)
    rm(mu_temp2)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # phi update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for (h in 1:k.star) {

      phi[ , h] <- rgamma(p, shape = rep(0.5 + kappa_1, p), rate = kappa_2 + 0.5 * tau[h] * lambda[ , h]^2)

    }

    # phi is a p x k.star matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # delta and tau update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    delta[1, 1] <- rgamma(1,
                          shape = a_1 + 0.5 * k.star * p,
                          rate = 1 + 0.5 * sum((tau / delta[1, 1]) * colsums(lambda^2 * phi)))

    tau <- matrix(cumprod(delta))

    if (k.star > 1) {

      for (h in 2:k.star) {

        delta_shape <- a_2 + 0.5 * p * (k.star - h + 1)
        delta_rate <- 1 + 0.5 * sum(((tau / delta[h, 1]) * colsums(lambda^2 * phi))[h:k.star])

        delta[h, 1] <- truncdist::rtrunc(spec = "gamma",
                                         n = 1,
                                         shape = delta_shape,
                                         rate = delta_rate,
                                         a = 1)

        tau <- matrix(cumprod(delta))

      }

    }


    # delta is a k.star x 1 matrix
    # tau is a k.star x 1 matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # alpha update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    alpha_shape1 <- sum((1 - missingness) * (data_working >= trunc.point)) + 1

    alpha_shape2 <- sum((missingness) * (data_working >= trunc.point)) + 1

    alpha <- rbeta(1, shape1 = alpha_shape1, shape2 = alpha_shape2)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # imputation
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Z <- rep(NA, length(miss_vec))

    if (dim(miss_ind)[1] != 0) {

      impute_mean <- sweep(mat.mult(lambda, Rfast::transpose(eta)), 1, mu, "+")  # p x n
      impute_var <- 1 / sigma_inv  # p x p

      for (point_ind in 1:dim(miss_ind)[1]) {

        lower_prob <- truncnorm::ptruncnorm(trunc.point,
                                            a = 0,
                                            mean = impute_mean[miss_ind[point_ind, ][2], miss_ind[point_ind, ][1]],
                                            sd = sqrt(impute_var[miss_ind[point_ind, ][2], 1]))

        higher_prob <- 1 - lower_prob

        zij_prob <- (lower_prob) / (lower_prob + alpha * higher_prob)

        zij <- rbinom(1, 1, zij_prob)

        Z[point_ind] <- zij

        if (zij == 0) {

          # impute as MAR
          data_working[miss_ind[point_ind, ][1], miss_ind[point_ind, ][2]] <-
            truncnorm::rtruncnorm(1, a = trunc.point,
                                  b = Inf,
                                  mean = impute_mean[miss_ind[point_ind, ][2], miss_ind[point_ind, ][1]],
                                  sd = sqrt(impute_var[miss_ind[point_ind, ][2], 1]))

        } else {

          # impute as MNAR
          data_working[miss_ind[point_ind, ][1], miss_ind[point_ind, ][2]] <-
            truncnorm::rtruncnorm(1, a = 0, b = trunc.point,
                                  mean = impute_mean[miss_ind[point_ind, ][2], miss_ind[point_ind, ][1]],
                                  sd = sqrt(impute_var[miss_ind[point_ind, ][2], 1]))

        }

      }

    }

    # ~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~
    # update storage
    # ~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~

    if (m %% thin == 0) {

      # imputed dataset
      store_data[[m / thin]] <- data_working[miss_vec]

      # alpha - prob of MAR
      store_alpha[[m / thin]] <- alpha

      # sigma
      store_sigma_inv[[m / thin]] <- sigma_inv

      # lambda - factor loadings
      store_lambda[[m / thin]] <- lambda

      # eta
      store_eta[[m / thin]] <- eta

      # mu - mean
      store_mu[[m / thin]] <- mu

      # effective factors
      store_k_t[[m / thin]] <- k.star

      # phi
      store_phi[[m / thin]] <- phi

      # delta
      store_delta[[m / thin]] <- delta

      # Z
      store_Z[[m / thin]] <- Z

    }

  }

  tifa_res <- list("store_data" = store_data,
                   "store_alpha" = store_alpha,
                   "store_sigma_inv" = store_sigma_inv,
                   "store_lambda" = store_lambda,
                   "store_eta" = store_eta,
                   "store_mu" = store_mu,
                   "store_k_t" = store_k_t,
                   "store_phi" = store_phi,
                   "store_delta" = store_delta,
                   "store_b_sigma" = b_sigma,
                   "store_Z" = store_Z,
                   "input_data" = input_data,
                   "input_missingness" = missingness
  )

  formatted_res <- tifa_res_format(tifa_res, burn = burn, thin = thin)

  output <- list("imputed_dataset" = formatted_res$imputed_dataset,
                 "imputation_info" = formatted_res$imputation_info)

  if (return_chain) {

    output$chain <- tifa_res

  }

  return(output)
}
