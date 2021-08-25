#' @title Computes new position of a point rotating about an orgin.
#' @description Computes the coordinate position of a point after rotation about an origin.
#' @param p numeric vector of length 2, the x and y coordinates of the point that will rotate
#' @param o numeric vector of length 2, the x and y coordinates of the origin
#' @param theta numeric, the angle of rotation
#' @return A vector of length 2 containing the resulting x,y coordinates
#' @export
#' @examples
#'p <- c(2,2)
#'o <- c(0,0)
#'
#'plot(c(-4,4),col=NULL)
#'points(p[1],p[2],col="red")
#'points(o[1],o[2])
#'
#'new.p <- point.ang.orig(p,o,pi/2*-1)
#'points(new.p[1],new.p[2],col="blue")
#'
point.ang.orig<- function(p,o,theta){
  xrot<-cos(theta)*(p[1]-o[1])-sin(theta)*(p[2]-o[2])+o[1]
  yrot<-sin(theta)*(p[1]-o[1])+cos(theta)*(p[2]-o[2])+o[2]
  return(c(xrot,yrot))
}

#' @title Computes orthogonal distance between a point and a line.
#' @description Computes 2D orthogonal distance between a point and a line given the points coordinates and the line's slope and intercept or model formula
#' @param x numeric, the x coordinate of the point
#' @param y numeric, the y coordinate of the point
#' @param slope numeric, the slope of the line
#' @param intercept numeric, the slope of the line
#' @param form  formula of type \code{lm} describing the line, ignored if \code{slope} and \code{intercept} are specified
#' @return The distance between the point and line
#' @details if \code{slope} and \code{intercept} are missing, a model formula of form \code{lm(y~x)} can be passed to \code{form}.
#' @export
#' @examples
#'x <- runif(1:10)
#'y <- x*0.5
#'
#'xy.lm <- lm(y~x)
#'
#'x.p <- y.p <- -3
#'
#'#one way . . .
#'dist.2d.line(x.p,y.p,form=xy.lm)
#'
#'#or another
#'dist.2d.line(x.p,y.p,slope=coef(xy.lm)[2],intercept=coef(xy.lm)[1])
#'
dist.2d.line <- function(x=NULL, y=NULL, slope=NULL, intercept=NULL,form=NULL) {
  
  if(all(sapply(list(slope,intercept,form), is.null))) stop("must specify slope and intercept of model formula")
  
  if(all(sapply(list(slope,intercept), is.null))&!is.null(form)) {
    intercept <- coef(form)[1]
    slope <- coef(form)[2]
  }
  
  b = c(1, intercept + slope)
  c = c(-intercept / slope, 0)
  a = c(x, y)
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1, v2)
  return(abs(det(m)) / sqrt(sum(v1 * v1)))
}


#' @title Computes distance between two points in Cartesian space.
#'
#' @description Computes distance between two points in Cartesian space using simple trigonometry functions
#'
#' @param x1 Numeric; x position of coordinate 1
#' @param y1 numeric; y position of coordinate 1
#' @param x2 numeric; x position of coordinate 2
#' @param y2 numeric; y position of coordinate 2
#' @return A single value of the distance between p[x1,y1] and p[x2,y2]
#' @export
#' @examples
#' #Find the lengths of the sides of a tringle and print to plot
#' x <- c(0,3,2)
#' y <- c(0,3,0)
#' plot(x,y)
#' lines(x,y)
#' lines(x[c(1,3)],y[c(1,3)])
#' hyp <- dist.2d(x[1],x[2],y[1],y[2])
#' s1 <- dist.2d(x[1],x[3],y[1],y[3])
#' s2 <- dist.2d(x[2],x[3],y[2],y[3])
#' text(mean(x[1:2],mean(y[2:3])),labels=round(hyp,1))
#' text(mean(x[c(1,3)]),y[1]+0.25,labels=round(s1,1))
#' text(mean(x[c(2:3)]),mean(y[2:3]),labels=round(s2,1))

dist.2d <- function(x1, x2, y1, y2) {
  sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2)
}

#' @title Computes angle between two segments sharing a point.
#' @description Computes angle between two segments of a triangle using law of cosines.
#' @param l Numeric; length of segment to the left
#' @param r Numeric; length of segment to the right
#' @param o Numeric; length of opposite segment
#' @return A single value of the angle in radians
#' @export
#' @examples
#'#a right triangle
#'L=3
#'R=3
#'O=sqrt(L^2+R^2)
#'cosine.ang(L,R,O)

cosine.ang <- function(l,r,o) {
  if(o==(l+r)){return(pi)}else{
    a <- acos((o^2-l^2-r^2)/(-2*l*r))
    return(a)}
}

#' @title Converts radians to degrees
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

#' @title Computes the heading between to cartesian points
#' @description Computes the heading (radians counter clockwise from north or vertical)
#'
#' @param x1 Numeric; point A x value
#' @param y1 numeric; point A y value
#' @param x2 numeric; point B x value
#' @param y2 numeric; point B y value
#' @return A single value in radians
#' @export
#' @return a single value in radians
#'
#' @examples
#' #example
#' A <- c(0,0)
#' B <- c(-5,5)
#' thet <- bearing.xy(A[1],B[1],A[2],B[2])
#' deg(thet)
#'
#'
bearing.xy <- function(x1,x2,y1,y2){
  theta <- atan((x2-x1)/(y1-y2))
  return(theta)
}

