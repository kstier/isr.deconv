#' @title C-Peptide Fraction Attributable to Short Half Life Lookup
#'
#' @description
#' Outputs the relevant fraction attributable to short half life for use in
#' [isr.deconv] based on patient type (normal, obese, or non-insulin dependent
#' diabetes mellitus).
#'
#' @param subject.type String for patient type.
#'
#' @returns
#' Numeric value.
#'
#' @export isr.fraction
#'
#' @examples
#' isr.fraction("normal")
#'
#' isr.fraction("obese")
#'
#' isr.fraction("niddm")
#'
#' @seealso [isr.deconv()]
#'
isr.fraction <- function(subject.type = c("normal", "obese", "niddm")) {
  if (is.character(subject.type)) {
    switch(subject.type,
      "normal" = 0.76,
      "obese" = 0.78,
      "niddm" = 0.78,
      stop("Error: subject.type must be one of 'normal', 'obese', or 'niddm'.")
    )
  } else {
    stop("Error: subject.type must be one of 'normal', 'obese', or 'niddm'.")
  }
}
