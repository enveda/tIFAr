#' Post-processing for tIFA results
#'
#' @param tifa_res The output object of the [tIFA_model()] function.
#' @param burn The desired burn-in. Should be a multiple of five.
#'
#' @return Returns a list containing two entries.
#'
#' * imputed_dataset: The returned dataset in matrix format with the imputed values inserted.
#'
#' * imputation_info: A dataframe containing information on imputed values and associated uncertainty for each missing entry in the original dataset.
#'
#' @importFrom stats quantile
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
#' # run tIFA model (short chain for example)
#' tifa_res <- tIFA_model(data = input_data$data, missingness = input_data$missingness, M = 1000,
#'                        initial_dataset = zero_imp, k = 5)
#'
#' # process tIFA results (short burn for example)
#' processed_res <- tifa_post(tifa_res, burn = 10)
tifa_post <- function(tifa_res, burn) {

  # extract input data from results object
  original_data <- tifa_res$input_data

  n <- dim(original_data)[1]
  p <- dim(original_data)[2]

  M <- length((tifa_res$store_data[-(1:(burn/5))]))

  # extract missingness information from results object
  missingness <- tifa_res$input_missingness

  miss_vec <- which(missingness == 0)
  miss_ind <- which(missingness == 0, arr.ind = TRUE)

  # Z designation
  store_Z <- sapply(tifa_res$store_Z[-(1:(burn/5))], FUN = matrix)

  # imputed data
  data_mat <- array(NA, dim = c(length((tifa_res$store_data[-(1:(burn/5))])[[1]]),
                                M))

  print(dim(data_mat))

  for (i in 1:M) {

    data_mat[ , i] <- (tifa_res$store_data[-(1:(burn/5))])[[i]]

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
