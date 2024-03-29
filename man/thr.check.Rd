% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image.R
\name{thr.check}
\alias{thr.check}
\title{Evaluate a range of threshold values for creating binary images}
\usage{
thr.check(img, min = NULL, max = NULL, otsu = TRUE, plot.grid = TRUE)
}
\arguments{
\item{img}{character, a path to an image file.}

\item{min}{numeric, the minimum threshold value (0-1).}

\item{max}{numeric, the maximum threshold value (0-1).}

\item{otsu}{logical, should the automatically determined threshold value be printed. See Details.}

\item{plot.grid}{logical, should the grid of images be plotted.}
}
\value{
If 'otsu=TRUE', a single automated threshold value is returned using Otsu's method. See \code{\link{otsu}}.
}
\description{
Converts an image to grayscale, applies a user defined threshold range to segment the binary the image, and then plots segmented images that reflect 9 discrete values in the range.
}
\details{
Displays a grid of 9 images as an R raster graphic, each the result of 9 discrete values within the range describe by \code{min} and \code{max}. If both \code{min} and \code{max} are \code{NULL}, then the threshold range is defined by the default \code{seq(0.1,0.9,.1)}. This function should help users refine threshold values for detecting ROIs with \code{kin.search} and \code{kin.simple} if automatic thresholding using Otsu's method is not satisfactory.
}
\examples{


#access image in system
 y <-system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")

# default values
thr.check(y)

#custom values
thr.check(y,min=0.2,max=1.5)


}
\seealso{
\code{\link{kin.simple}}, \code{\link{kin.search}}, \code{\link{otsu}}
}
