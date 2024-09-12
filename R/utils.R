#' Data preprocessing
#'
#' @description
#' Processes a dataset for use in a tIFA model.
#'
#'
#' @param data The data to be processed. If not a matrix, will be transformed in to a matrix.
#' @param scaled TRUE/FALSE. Should the input data be scaled? Defaults to FALSE.
#' @param coding The coding used to indicate a data entry is missing. Defaults to NA.
#' @param pareto_scale TRUE/FALSE. Should the input data be Pareto scaled? Defaults to FALSE.
#' @param log_trans TRUE/FALSE. Should the input data be log transformed? Defaults to FALSE.
#'
#' @return A list containing three entries.
#'
#' * "data" contains the processed data in matrix format with missing entries coded as NA.
#'
#' * "missingness_pattern" contains the missingness pattern of the data for use in the tIFA model.
#'
#' * "divisors" contains a vector of the divisors used to scale the input data, if applicable. If unscaled, the value will be NULL.
#' @export
#'
#' @importFrom stats sd
#'
#' @examples
#' example_data <- matrix(rnorm(20), nrow = 5)
#' example_data[4, 2] <- 0.001  # add missingness to example data coded as 0.001
#' FA_process(example_data, coding = 0.001)
#'
#' @noRd
FA_process <- function(data, scaled = FALSE, coding = NA, pareto_scale = FALSE,
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
#' @export
#'
#' @examples
#' x <- c(1, 2, 3, 4, 4, 5)
#' Mode(x)
#'
#' @noRd
Mode <- function(x) {

  options <- unique(x)

  return(options[which.max(tabulate(match(x, options)))])

}
