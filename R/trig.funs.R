#' convert radians to degrees
#'
#' @param x Numeric; value in radians
#' @return A single value
#' @export
#' @seealso \code{\link{rad}}
deg <- function(x){180/pi*x}


#' convert degrees to radians
#'
#' @param x Numeric; value in degrees
#' @return A single value
#' @export
#' @seealso \code{\link{deg}}
rad <- function(x){pi/180*x}

#' compute the heading between to cratesian points (radians counter clockwise from north or vertical)
#'
#' @param x1 Numeric; point A x value
#' @param y1 numeric; point A y value
#' @param x2 numeric; point B x value
#' @param y2 numeric; point B y value
#' @return A single value in radians
#' @export
#'
bearing.xy <- function(x1,x2,y1,y2){
  theta <- atan((x2-x1)/(y1-y2))
  return(theta)
}

#example

#' @examples
#' A <- c(0,0)
#' B <- c(-5,5)
#' thet <- bearing.xy(A[1],B[1],A[2],B[2])
#' deg(thet)
