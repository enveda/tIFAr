#' Extract information from input data in preparation for imputation.
#'
#' @description
#' Extracts information from a dataset for use in a tIFA model.
#'
#'
#' @param data The data to be processed. If not a matrix, will be transformed in to a matrix.
#' @param scaled `TRUE`/`FALSE.` Should the input data be scaled? Defaults to `FALSE.`
#' @param coding The coding used to indicate a data entry is missing. Defaults to `NA.`
#' @param pareto_scale `TRUE`/`FALSE.` Should the input data be Pareto scaled? Defaults to `FALSE.`
#' @param log_trans `TRUE`/`FALSE.` Should the input data be log transformed? Defaults to `FALSE.`
#'
#' @return A list containing three entries.
#'
#' * data: Contains the processed data in matrix format with missing entries coded as `NA`.
#'
#' * missingness_pattern: Contains the missingness pattern of the data for use in the tIFA model.
#'
#' * divisors: Contains a vector of the divisors used to scale the input data, if applicable. If unscaled, the value will be `NULL`.
#'
#' @importFrom stats sd
#'
#' @noRd
tifa_prep <- function(data, scaled = FALSE, coding = NA, pareto_scale = FALSE,
                       log_trans = FALSE) {

  # convert to matrix as required
  if (!(is.matrix(data))){
    data <- as.matrix(data)
  }

  # how missing data is encoded
  data[data == coding] <- NA

  # centre or scale as required
  if (scaled){
    data <- scale(data, center = FALSE, scale = scaled)
  }

  divisors <- NULL

  if (pareto_scale) {

    divisors <- sqrt(apply(data, 2, sd, na.rm = TRUE))

    data <- sweep(data, MARGIN = 2, STATS = divisors, FUN = "/")
  }

  if (log_trans) {
    data <- log(data)
  }

  # generate missingness pattern
  missingness_pattern <- 1 - is.na(data)

  return(list("data" = data, "missingness_pattern" = missingness_pattern, "divisors" = divisors))
}

#' Mode
#'
#' @param x A vector.
#'
#' @return The mode of x. If x is multi-modal, the first mode (by index in x) is returned.
#'
#' @noRd
Mode <- function(x) {

  options <- unique(x)

  return(options[which.max(tabulate(match(x, options)))])

}

#' Results formatting for tIFA results
#'
#' @param tifa_res Draws from the MCMC chain of the tIFA model
#' @param burn The number of MCMC iterations to be discarded in an initial burn.
#' @param thin The level of thinning to take place in the MCMC chain. E.g., `thin = 5` will retain every fifth draw, post-burn.
#'
#' @return Returns a list containing two entries.
#'
#' * imputed_dataset: The returned dataset in matrix format with the imputed values inserted.
#'
#' * imputation_info: A dataframe containing information on imputed values and associated uncertainty for each missing entry in the original dataset.
#'
#' @importFrom stats quantile
#'
#' @noRd
tifa_res_format <- function(tifa_res, burn, thin) {

  # extract input data from results object
  original_data <- tifa_res$input_data

  n <- dim(original_data)[1]
  p <- dim(original_data)[2]

  M <- length((tifa_res$store_data[-(1:(burn/thin))]))

  # extract missingness information from results object
  missingness <- tifa_res$input_missingness

  miss_vec <- which(missingness == 0)
  miss_ind <- which(missingness == 0, arr.ind = TRUE)

  # Z designation
  store_Z <- sapply(tifa_res$store_Z[-(1:(burn/thin))], FUN = matrix)

  # imputed data
  data_mat <- array(NA, dim = c(length((tifa_res$store_data[-(1:(burn/thin))])[[1]]),
                                M))

  print(dim(data_mat))

  for (i in 1:M) {

    data_mat[ , i] <- (tifa_res$store_data[-(1:(burn/thin))])[[i]]

  }

  # imputed values
  just_imputed <- data_mat

  # modal designations
  mode_Z <- apply(store_Z, 1, Mode)

  MNAR_designations <- which(mode_Z == 1)
  MAR_designations <- which(mode_Z == 0)

  designations_vec <- rep("MNAR", dim(store_Z)[1])
  designations_vec[MAR_designations] <- "MAR"

  MNAR_post_mean <- c()
  MNAR_95_cred_upper <- c()
  MNAR_95_cred_lower <- c()

  for (mnar_ind in MNAR_designations) {

    relevant_draws <- just_imputed[mnar_ind, which(store_Z[mnar_ind, ] == mode_Z[mnar_ind])]

    post_mean <- mean(relevant_draws)
    cred_int <- quantile(relevant_draws, probs = c(0.025, 0.975))
    MNAR_post_mean <- c(MNAR_post_mean, post_mean)
    MNAR_95_cred_lower <- c(MNAR_95_cred_lower, cred_int[1])
    MNAR_95_cred_upper <- c(MNAR_95_cred_upper, cred_int[2])

  }

  MAR_post_mean <- c()
  MAR_95_cred_upper <- c()
  MAR_95_cred_lower <- c()

  for (mar_ind in MAR_designations) {

    relevant_draws <- just_imputed[mar_ind, which(store_Z[mar_ind, ] == mode_Z[mar_ind])]

    post_mean <- mean(relevant_draws)
    cred_int <- quantile(relevant_draws, probs = c(0.025, 0.975))
    MAR_post_mean <- c(MAR_post_mean, post_mean)
    MAR_95_cred_lower <- c(MAR_95_cred_lower, cred_int[1])
    MAR_95_cred_upper <- c(MAR_95_cred_upper, cred_int[2])

  }


  tifa_imp <- original_data
  tifa_imp[miss_vec][MNAR_designations] <- MNAR_post_mean
  tifa_imp[miss_vec][MAR_designations] <- MAR_post_mean

  cred_int_upper <- matrix(NA, nrow = n, ncol = p)
  cred_int_upper[miss_vec][MNAR_designations] <- MNAR_95_cred_upper
  cred_int_upper[miss_vec][MAR_designations] <- MAR_95_cred_upper

  cred_int_lower <- matrix(NA, nrow = n, ncol = p)
  cred_int_lower[miss_vec][MNAR_designations] <- MNAR_95_cred_lower
  cred_int_lower[miss_vec][MAR_designations] <- MAR_95_cred_lower

  # get prop of each point's Zs that are MNAR
  z_props <- apply(store_Z, 1, sum) / dim(store_Z)[2]

  # flip designated MAR to 1 - MNAR certainty
  z_props[MAR_designations] <- 1 - z_props[MAR_designations]

  # for ease of indexing
  z_uncertainty_tifa <- matrix(NA, nrow = n, ncol = p)
  z_uncertainty_tifa[miss_vec][MNAR_designations] <- (1 - z_props[MNAR_designations])
  z_uncertainty_tifa[miss_vec][MAR_designations] <- (1 - z_props[MAR_designations])

  # create dataframe with uncertainty information
  imputation_info <- data.frame(entry_row = miss_ind[ , 1],
                                entry_col = miss_ind[ , 2],
                                imputed_val = tifa_imp[miss_vec],
                                cred_int_upper = cred_int_upper[miss_vec],
                                cred_int_lower = cred_int_lower[miss_vec],
                                miss_mech = designations_vec,
                                miss_mech_unc = z_uncertainty_tifa[miss_vec]
  )


  return(list("imputed_dataset" = tifa_imp,
              "imputation_info" = imputation_info))

}
