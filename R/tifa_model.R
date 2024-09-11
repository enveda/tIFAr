#' Truncated Infinite Factor Analysis Imputation Model
#'
#' @description
#' Applies the tIFA imputation method to the input dataset.
#'
#'
#' @param data The data on which to apply imptutation. Should be processed using [FA_process()] function.
#' @param missingness The missingness pattern of the input data. This is returned from the [FA_process()] function applied to the input data.
#' @param M The number of MCMC iterations desired. A multiple of 5 is expected.
#' @param initial_dataset Starting imputation for the input data. E.g., the imputed dataset resulting from mean imputation.
#' @param k The number of "effective-factors" desired for use in the tIFA model. Defaults to 5.
#'
#' @return Returns a list containing 11 entries.
#'
#' * store_data: MCMC draws of the imputed data entries.
#'
#' * store_alpha: MCMC draws of the alpha parameter.
#'
#' * store_sigma_inv: MCMC draws of the sigma^-2 entries.
#'
#' * store_lambda: MCMC draws of the loadings matrix.
#'
#' * store_eta: MCMC draws of the latent factor scores.
#'
#' * store_mu: MCMC draws of the mean vector.
#'
#' * store_k_t: MCMC draws of the effective factor number (unchanging).
#'
#' * store_phi: MCMC draws of the local shrinkage parameters.
#'
#' * store_delta: MCMC draws of the multiplicative gamma process parameters.
#'
#' * store_b_sigma: MCMC draws of the b_sigma parameter (unchanging).
#'
#' * store_Z: MCMC draws of the missingness designation indicator for each imputed entry.
#'
#' * input_data: The dataset input into the tIFA model.
#'
#' * input_missingness: The missingness pattern input into the tIFA model.
#'
#' @importFrom Rfast transpose mat.mult colsums
#' @importFrom TruncatedNormal rtmvnorm
#' @importFrom truncdist rtrunc
#' @importFrom stats rgamma rbeta rbinom runif prcomp cov
#' @importFrom truncnorm ptruncnorm rtruncnorm
#' @importFrom Rcpp cppFunction
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
#' # apply data pre-processing function
#' input_data <- FA_process(example_data, coding = 0.001)
#'
#' # perform zero imputation to use as initial imputed data value in tIFA
#' zero_imp <- input_data$data
#' zero_imp[is.na(zero_imp)] <- 0
#'
#' # run tIFA model
#' tIFA_model(data = input_data$data, missingness = input_data$missingness, M = 50,
#'            initial_dataset = zero_imp, k = 5)
#'
tIFA_model <- function(data, missingness, M, initial_dataset, k = 5) {

  Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }", depends = "RcppArmadillo")

  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~
  # initialisation
  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~

  # data dimensions
  n <- nrow(data)
  p <- ncol(data)

  # seeding function for rcpp functions
  get_seed <- function() {

    sample.int(.Machine$integer.max, 1)

  }

  # truncation point
  trunc.point <- min(data, na.rm = TRUE)

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

  data_working <- initial_dataset
  mu <- matrix(apply(data, 2, mean, na.rm = TRUE))
  mu_tilde <- mu
  mu_varphi <- 0.1
  kappa_1 <- 3L
  kappa_2 <- 2L
  a_sigma <- 1L
  b_sigma <- 0.3
  a_1 <- 2.1
  a_2 <- 3.1
  alpha <- runif(1, 0, 1)

  # use pca to initialise factors and loadings
  pca_res <- prcomp(data_working, scale. = FALSE, center = FALSE)

  lambda <- pca_res$rotation[ , 1:k]  # p x k
  rownames(lambda) <- NULL
  colnames(lambda) <- NULL
  lambda[which(lambda < 0)] <- abs(lambda[which(lambda < 0)])

  eta <- matrix(NA, nrow = n, ncol = k)

  for (i in 1:n) {

    eta[i , ] <- TruncatedNormal::rtmvnorm(n = 1,
                                           mu = rep(0, k),
                                           sigma = diag(k),
                                           lb = rep(0, k))

  }  # n x k

  eps <- matrix(NA, nrow = n, ncol = p)

  for (i in 1:n) {

    eps[i, ] <- data_working[i, ] - mu - (lambda %*% matrix(eta[i ,]))

  }

  sigma_inv <- matrix(1 / diag(cov(eps, use = "complete")))

  Sigma_inv <- diag(sigma_inv[ , 1])

  phi <- matrix(rgamma(n = p * k, shape = kappa_1, rate = kappa_2),  # p x k
                nrow = p)

  delta <- matrix(rep(0, k), nrow = k)  # k x 1
  delta[1, 1] <- rgamma(n = 1, shape = a_1, rate = 1)
  delta[-1, 1] <- truncdist::rtrunc(spec = "gamma",
                                    n = k - 1,
                                    shape = a_2,
                                    rate = 1,
                                    a = 1)

  tau <- matrix(cumprod(delta))  # k x 1

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

  for (m in seq_len(M)) {

    # create a readout
    statement <- paste(Sys.time(), "iteration:", m, ", k_t:", k)
    print(statement)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # lambda update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    eta_t_eta <- mat.mult(Rfast::transpose(eta) , eta)  # k x k

    for (j in 1:p) {

      if (k == 1) {

        lambda_var <- armaInv(sigma_inv[j] * eta_t_eta + diag(c(phi[j, ] * tau), nrow = 1))  # k x k

      } else {

        lambda_var <- armaInv(sigma_inv[j] * eta_t_eta + diag(c(phi[j, ] * tau)))  # k x k

      }

      lambda_mean <- mat.mult(lambda_var, mat.mult(
        Rfast::transpose(eta), as.matrix(data_working[ , j] - mu[j])) * sigma_inv[j])  # k x 1


      lambda[j, ] <- TruncatedNormal::rtmvnorm(n = 1, mu = as.vector(lambda_mean),
                                               sigma = lambda_var, lb = rep(0, k))

    }

    # set NaN to small number
    # as approach zero from above
    lambda[which(lambda == "NaN" | lambda == Inf)] <- 1e-8

    # lambda is a p x k matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # eta update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    lambda_t_Sigma <- mat.mult(Rfast::transpose(lambda), Sigma_inv)

    eta_var <- armaInv(diag(k) + mat.mult(lambda_t_Sigma, lambda))  # k x k

    for (i in 1:n) {

      eta_mean <- mat.mult(eta_var, mat.mult(lambda_t_Sigma, as.matrix(data_working[i, ] - mu)))  # k x 1

      eta[i , ] <-  TruncatedNormal::rtmvnorm(n = 1, mu = as.vector(eta_mean),
                                              sigma = eta_var, lb = rep(0, k))

    }

    eta[which(eta == "NaN" | eta == Inf)] <- 0

    # eta is a n x k matrix

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

    for (h in 1:k) {

      phi[ , h] <- rgamma(p, shape = rep(0.5 + kappa_1, p), rate = kappa_2 + 0.5 * tau[h] * lambda[ , h]^2)

    }

    # phi is a p x k matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # delta and tau update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    delta[1, 1] <- rgamma(1,
                          shape = a_1 + 0.5 * k * p,
                          rate = 1 + 0.5 * sum((tau / delta[1, 1]) * colsums(lambda^2 * phi)))

    tau <- matrix(cumprod(delta))

    if (k > 1) {

      for (h in 2:k) {

        delta_shape <- a_2 + 0.5 * p * (k - h + 1)
        delta_rate <- 1 + 0.5 * sum(((tau / delta[h, 1]) * colsums(lambda^2 * phi))[h:k])

        delta[h, 1] <- truncdist::rtrunc(spec = "gamma",
                                         n = 1,
                                         shape = delta_shape,
                                         rate = delta_rate,
                                         a = 1)

        tau <- matrix(cumprod(delta))

      }

    }


    # delta is a k x 1 matrix
    # tau is a k x 1 matrix

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

    if (m %% 5 == 0) {

      # imputed dataset
      store_data[[m / 5]] <- data_working[miss_vec]

      # alpha - prob of MAR
      store_alpha[[m / 5]] <- alpha

      # sigma
      store_sigma_inv[[m / 5]] <- sigma_inv

      # lambda - factor loadings
      store_lambda[[m / 5]] <- lambda

      # eta
      store_eta[[m / 5]] <- eta

      # mu - mean
      store_mu[[m / 5]] <- mu

      # effective factors
      store_k_t[[m / 5]] <- k

      # phi
      store_phi[[m / 5]] <- phi

      # delta
      store_delta[[m / 5]] <- delta

      # Z
      store_Z[[m / 5]] <- Z

    }

  }

  output <- list("store_data" = store_data,
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
                 "input_data" = data,
                 "input_missingness" = missingness
  )

  return(output)
}
