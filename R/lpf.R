#' Apply a low-pass filter to some time-series data.
#'
#' @param y time series values.
#' @param t time.
#' @param f frequency,
#' @return Filtered values that represent a rolling average, more or less...
#' @examples
#' y <- sin(seq(1,30,0.1))
#' y <- sin(seq(1,30,0.1))
#' y <- jitter(y,2,2)
#' plot(y)
#' points(lpf(y,1,1),col="red",pch=3)
lpf <- function( y, t, f ) {
  rc <- 1 / ( 2 * pi * f )
  a  <- t / ( t + rc )
  n  <- length( y )
  yf <- y
  for( i in 2:length(y) ) {
    yf[i] <- a * y[i] + (1-a) * yf[i-1]
  }
  return( yf )
}
