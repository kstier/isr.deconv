#' @title C-Peptide Short Half Life Calculator
#'
#' @description
#' Outputs the relevant long half life of c-peptide for use in [isr.deconv()]
#' based on patient type (normal, obese, or non-insulin dependent diabetes
#' mellitus).
#'
#' @param subject.type String for patient type.
#'
#' @returns
#' Numeric value.
#'
#' @export isr.shortHL
#'
#' @examples
#' isr.shortHL("normal")
#'
#' isr.shortHL("obese")
#'
#' isr.shortHL("niddm")
#'
#' @seealso [isr.deconv()]
#'
isr.shortHL <- function(subject.type = c("normal", "obese", "niddm")) {
  if (is.character(subject.type)) {
    switch(subject.type,
      "normal" = 4.95,
      "obese" = 4.55,
      "niddm" = 4.52,
      stop("Error: subject.type must be one of 'normal', 'obese', or 'niddm'.")
    )
  } else {
    stop("Error: subject.type must be one of 'normal', 'obese', or 'niddm'.")
  }
}
