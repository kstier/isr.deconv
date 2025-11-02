#' @title Serum Volume Calculator
#'
#' @description
#' Calculates the estimated serum volume based on patient sex, weight, and
#' height for use in [isr.deconv()].
#'
#' @param subject.sex String for patient sex, `"m"` or `"f"`.
#' @param subject.weight Numeric for subject weight in kilograms.
#' @param subject.height Numeric for subject height in centimeters.
#'
#' @returns
#' Numeric value.
#'
#' @export isr.volume
#'
#' @examples
#' isr.volume("m", 86.2, 181.5)
#'
#' @seealso [isr.deconv()]
#'
isr.volume <- function(subject.sex = c("m", "f"), subject.weight, subject.height) {
  subject.bsa <- subject.weight^0.425 * subject.height^0.725 * 71.84
  if (subject.sex == "m") {
    1.92 * subject.bsa + 0.64
  } else if (subject.sex == "f") {
    1.11 * subject.bsa + 2.04
  } else {
    stop("Error: subject.sex must be one of 'm' or 'f'.")
  }
}
