#' @title Insulin Secretion Rate Deconvolution
#'
#' @description
#' Estimates insulin secretion rate based on a time series of c-peptide values
#' based on the Van Cauter method. C-peptide values are interpolated using
#' either linear approximation or a cubic spline; the linear method (default)
#' will output a time series of insulin secretory rates at time points between
#' each input time, whereas the spline method will output a function which can
#' be called to return insulin secretory rates for specific time points.
#'
#' Typically, c-peptide values are provided in pmol/mL and time is in minutes,
#' resulting in insulin secretion rate outputs in pmol/min.
#'
#' @returns If `method = "linear"` is selected (default), returns a list of time
#' points and insulin secretory rates. If `method = "spline"` is selected, returns
#' a function which can be called to return insulin secretory rates at specified
#' time points.
#'
#' @param timeseries Vector of numeric time values corresponding to the c-peptide values in `cpepseries`
#' @param cpepseries Vector of numeric c-peptide values at each time point in `timeseries`
#' @param vol Volume of distribution in the main compartment (i.e., serum volume), which can be calculated with the [isr.volume] function
#' @param shl Short half-life, which can be calculated with the [isr.shortHL] function
#' @param lhl Long half-life, which can be calculated with the [isr.longHL] function
#' @param frc Fraction attributable to the short half-life, which can be calculated with the [isr.fraction] function
#' @param method Determines whether C-peptide values are interpolated with linear connections between the points or a cubic spline
#' @param plotspline Dictates whether to print a plot of the C-peptide interpolation
#' @param plotisr Dictates whether to print a plot of the ISR output values
#' @param isr.validated.vals Used for debugging; validated ISR values to compare to past analyses
#' @param isr.validated.time Used for debugging; time values corresponding to validated ISR values to compare to past analyses
#'
#' @importFrom npreg ss
#' @importFrom graphics points
#' @importFrom stats approx integrate predict
#'
#' @seealso [isr.volume()], [isr.shortHL()], [isr.longHL()], [isr.fraction()]
#'
#' @examples
#' isr.deconv(
#'   timeseries = c(-30, 0, 30, 60, 90, 120),
#'   cpepseries = c(1.72, 1.72, 5.40, 5.23, 2.71, 1.79),
#'   vol = 6104,
#'   shl = 4.55,
#'   lhl = 31.05,
#'   frc = 0.78,
#'   method = "linear",
#'   plotspline = TRUE,
#'   plotisr = TRUE,
#' )
#'
#' @export isr.deconv

isr.deconv <- function(timeseries, cpepseries, vol, shl, lhl, frc, method = c("linear", "spline"), plotspline = FALSE, plotisr = FALSE, isr.validated.vals, isr.validated.time) {
  ## Validate inputs
  method <- match.arg(method)
  if(length(timeseries) != length(cpepseries)) {
    stop("Error: timeseries and cpepseries must have the same number of values.")
  }

  ## Sort values
  ord <- order(timeseries)
  timeseries <- timeseries[ord]
  cpepseries <- cpepseries[ord]

  ## Define coefficients
  A <- frc / vol
  a <- log(2) / shl
  B <- (1 - frc) / vol
  b <- log(2) / lhl
  k2 <- (A * b + a * B) / (A + B)
  k3 <- a * b / k2
  k1 <- a + b - (k2 + k3)
  C <- cpepseries * vol
  t0 <- min(timeseries)
  tmax <- max(timeseries)

  if (method == "spline") {
    ## Define spline
    f_spline <- ss(x = timeseries, y = C, spar = 0.3)

    ## Value of the spline at a point t
    spline_value <- function(t) predict(f_spline, t)$y

    ## Derivative of the spline at a point t; requires scaling adjustment
    spline_derivative <- function(t) predict(f_spline, t, deriv = 1)$y / (tmax - t0)

    ## Integral of the spline from t0 to t
    spline_integral <- function(t0, t) integrate(function(s) predict(f_spline, s)$y * exp(k2 * s), lower = t0, upper = t)$value

    ## ISR function
    isr <- function(t0, t) -exp(-k2 * t) * (k1 * spline_value(t0) * exp(k2 * t0) + k1 * k2 * spline_integral(t0, t)) + spline_derivative(t) + (k1 + k3) * spline_value(t)

    ## Spline and ISR value plot
    if (plotspline | plotisr) {
      t_grid <- seq(t0, tmax, length.out = 100)
      if (plotspline) {
        plot_vals <- vapply(t_grid, spline_value, numeric(1))
        plot(t_grid, plot_vals,
          type = "l", lwd = 2,
          xlab = "time", ylab = "total c-peptide",
          main = "Spline Fit to c-peptide*volume Time Series"
        )
        points(timeseries, C, pch = 19, cex = 0.6)
      }
      if (plotisr) {
        isr_vals <- vapply(t_grid, function(s) isr(t0, s), numeric(1))
        plot(t_grid, isr_vals,
          type = "l", lwd = 2,
          xlab = "time", ylab = "ISR",
          main = "Estimated Insulin Secretory Rate"
        )
        if (!(missing(isr.validated.vals) | missing(isr.validated.time))) points(isr.validated.time, isr.validated.vals, pch = 4, col = "green", lwd = 2)
      }
    }

    ## Output function
    return(function(t) vapply(t, function(s) isr(t0, s), numeric(1)))
  } else if (method == "linear") {
    timeseries.right <- timeseries[-1]
    timeseries.left <- timeseries[-length(timeseries)]
    C.right <- C[-1]
    C.left <- C[-length(C)]

    midpoint.times <- (timeseries.right + timeseries.left) / 2
    midpoint.times.t0 <- append(t0, midpoint.times)

    C.slope <- (C.right - C.left) / (timeseries.right - timeseries.left)

    f_linear <- function(t) approx(x = timeseries, y = C, xout = t, method = "linear", rule = 1)$y

    f_linear.vals <- f_linear(midpoint.times)

    f_linear_integrand_fun1 <- function(s) f_linear(s) * exp(k2 * s)

    f_linear_integrand_vals <- vapply(midpoint.times.t0, f_linear_integrand_fun1, numeric(1))

    f_linear_integrand_fun2 <- function(t) approx(x = midpoint.times.t0, y = f_linear_integrand_vals, xout = t, method = "linear", rule = 1)$y

    linear.integral <- vapply(midpoint.times, function(t) integrate(f_linear_integrand_fun2, lower = t0, upper = t)$value, numeric(1))

    linear.isr.vals <- -exp(-k2 * midpoint.times) * (k1 * f_linear(t0) * exp(k2 * t0) + k1 * k2 * linear.integral) + C.slope + (k1 + k3) * f_linear.vals

    ## Spline and ISR value plot
    if (plotspline | plotisr) {
      if (plotspline) {
        plot(timeseries, C,
          type = "o", lwd = 2, pch = 19,
          xlab = "time", ylab = "total c-peptide",
          main = "Linear Fit to c-peptide*volume Time Series"
        )
      }
      if (plotisr) {
        plot(midpoint.times, linear.isr.vals,
          type = "o", lwd = 2, pch = 19, xlim = c(t0, tmax),
          xlab = "time", ylab = "ISR",
          main = "Estimated Insulin Secretory Rate"
        )
        if (!(missing(isr.validated.vals) | missing(isr.validated.time))) points(isr.validated.time, isr.validated.vals, pch = 4, col = "green", lwd = 2)
      }
    }

    ## Output values
    return(list(x = midpoint.times, y = linear.isr.vals))
  }
}
