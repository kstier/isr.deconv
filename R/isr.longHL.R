#' @title C-Peptide Long Half Life Calculator
#'
#' @description
#' Outputs the long half life based on patient age for use in [isr.deconv] per
#' the Van Cauter method of estimating insulin secretion rate.
#'
#' @param subject.age Numeric for patient age in years.
#'
#' @returns
#' Numeric value.
#'
#' @export isr.longHL
#'
#' @examples
#' isr.longHL(18.08)
#'
#' @seealso [isr.deconv()]
#'
isr.longHL <- function(subject.age) 0.14 * subject.age + 29.2
