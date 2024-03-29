% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kin.R
\name{fin.kin}
\alias{fin.kin}
\title{Tracking of fin-like extensions of body contours}
\usage{
fin.kin(
  kin,
  out = "cont",
  frames = NULL,
  fin.pos = NULL,
  smooth.n = 5,
  ml.meth = "hull",
  ml.smooth = 0.9,
  red = NULL
)
}
\arguments{
\item{kin, }{a named list from one of \code{trackter}s \code{kin} functions.}

\item{out}{character, specifying the contour data in 'kin'}

\item{frames}{integer, the frames within 'kin' to analyze. If NULL then all frames are included.}

\item{fin.pos}{numeric, a vector of length 2 indicating the start and end of the contour region that contains the fins of interest as a proportion of the body length.}

\item{smooth.n}{numeric, the number of smoothing operations undertaken by \link{coo_smooth} on the contour described by 'x'. See Details.}

\item{ml.meth}{character, the midline detection method. One of 'ang' for bisection using \code{\link{free.ml.ang}} or 'hull' for bisection using \code{\link{free.ml.hull}}. Delaunay triangulation using \code{\link{free.ml.del}} is not supported. See Details.}

\item{ml.smooth}{numeric (0-1), the smoothing value for the midline. See Details.}

\item{red}{numeric, between 0-1 the proportion of contour coordinates to sample. Will speed up fin position and midline estimations. If 'NULL', the full contour in \code{out} will be used See Details.}
}
\value{
A list with the following components:

\code{cont} a data table consisting of x,y coordinates of the body contour

\code{fin} a data table describing the contour of the fins consisting of the following:

\itemize{
\item 'side': fin side, 'a' or 'b'
\item 'n': the position of fin coordinates where min(n) is closest to the head
\item x,y coordinates within the range of \code{fin.pos}

}

\code{fin.pts} a data table describing fin position consisting of the following:

\itemize{
\item 'side': fin side, 'a' or 'b'
\item 'n': The index matching that of the body contour coordinates 
\item  x,y coordinates of the fin tips, start, and end within the range of \code{fin.pos}.
\item 'n': The index matching that of the body contour coordinates 
\item 'pos': description  of the coordinates' positions, 'start', 'end' or 'tip' or 'tip2'.
}

\code{comp} a data table describing the composite contour of the body minus the fins.
\itemize{
\item 'n': The index matching that of the body contour coordinates 
\item  x,y coordinates of the body except the range of x values within \code{fin.pos}. These values take on a straight line described by the prediction of \code{lm()} based on the start and end of the fin. See Details.
}

\code{midline} a data table describing the estimated midline, 'x', 'y', the smoothed x and y positions, respectively

\code{bl} the body length in pixels

\code{amp} a data table describing the estimated fin amplitudes based on each method of finding the fin tip. 'amp1', the amplitude of the fin positions based on the maximum euclidean distance from the fin base and 'amp2', the distance of the fin points (based on distance inflections) from the fin base.
}
\description{
Estimates the amplitudes of regions along a body contour that are protruding. Useful in computing paired-fin amplitudes from contour data produced from  \link{kin.simple}, \link{kin.search}, or \link{kin.free}. Also computes a smoothed midline based on the body outline with the fin region removed.
}
\details{
If \code{red} is specified, the contour \code{out} is sampled with \code{\link{coo_interpolate}} from the \code{Momocs} package. The number of points sampled (n) equals \code{red} times the number of points in \code{out}

To establish the contour positions that are within \code{fin.pos}, a midline is estimated using one of two methods specified by \code{ml.meth}, "ang" for \code{\link{free.ml.ang}} or "hull" for  \code{\link{free.ml.hull}}. Midline points are indexed by position along the body length by calculating the cumulative distance between midline coordinates in pixels. This midline distance is then used to estimate the position of the fin about the contour using the parameter \code{fin.pos}.

The positions of the tip of the fin appendages is estimated in two ways. The first is simply the point in each appendage that is farthest from the base of the fin. The base is estimated as a straight line between the contour coordinates that match \code{fin.pos}.  The second is a little more complicated and starts with calculation the distance between each fin contours coordinates and the midpoint of the fin base. \code{\link{features}} from the \code{features} package is then use to calculate an inflection in this distance and the point of this inflection is used to estimate the fin position. Amplitudes of each method are calculated based on the orthogonal euclidean distance from the fin bases.

In addition to fin amplitude and contour extraction, this function also produces a composite contour of the body minus the fin area described by \code{fin.pos}. Fin contours are replaced by a simple linear prediction constructed from the coordinates of the first and last values covered by \code{fin.pos}, that is, the fin bases. The result is a straight line between the start and end of each fin. 

From this composite body contour, a midline prediction is made based on the mean position of contour coordinates from each side of the contour sharing a the same position along body using the \code{free.ml} functions. 

#' For midline smoothing, \code{\link{smooth_spline}} from the \code{smoothr} package is used to interpolate points between a reduced number of vertices using piecewise cubic polynomials. The number of vertices is calculated based on the number of midline coordinates times the value of \code{ml.smooth}.
}
\examples{
###plot pectoral-fin amplitudes of a swimming sunfish
\dontrun{
require(ggplot2)

#download example avi video
f <- "https://github.com/ckenaley/exampledata/blob/master/sunfish_pect.avi?raw=true"
download.file(f,"sunfish.avi")

#create directories
 ti <-paste0(tempdir(),"/images")
 dir.create(ti)

#extract images with ffmpeg operations and reduce them to 600 px wide with a filter
filt.red <- " -vf scale=600:-1 " #filter
vid.to.images2(vid.path="sunfish.avi",filt = filt.red,out.dir=ti) #extract

#number of frames
fr <- list.files(ti,full.names=TRUE)
thr.check(fr[3])

#extract contours and other data
kin <- kin.search(image.dir = ti,thr=0.9,ant.per = 0.25,save=FALSE)

#fin data by frame
fin.pos <- c(0.25,.55)
fin.dat <- fin.kin(kin=kin,fin.pos = fin.pos,smooth.n=1,ml.smooth=0.3)

p <- ggplot(dat=fin.dat$amp,aes(x=frame,y=amp2,col=side))+geom_line()+theme_classic(15)
print(p)


## plot body and fin contours of a frame

#plot body contour and fins
fr <- 8
p <- qplot(data=fin.dat$cont[frame==fr],x=x,y=y)
p <- p+geom_point(data=fin.dat$fin[frame==fr],aes(x,y),col="red",size=3)
p <- p+geom_point(data=fin.dat$fin.pts[frame==fr],aes(x,y,shape=pos))
p+xlim(c(0,kin$dim[1]))+ylim(c(0,kin$dim[2]))


#plot body contour minus fins and the body midline
p <- qplot(data=fin.dat$comp[frame==fr],x=x,y=y)
p <- p++geom_point(data=fin.dat$midline[frame==fr],aes(x,y),col="red",size=2)
p+xlim(c(0,kin$dim[1]))+ylim(c(0,kin$dim[2]))

# now the whole sequence with gg.overlay()
 gg.overlay(dat=fin.dat,
 under="cont",
 over="fin.pts",
 size=5,
 animate=TRUE,
 fps=5,
 alph=0.05)
 
unlink(ti,recursive=TRUE)
}


}
\seealso{
\code{\link{kin.simple}}, \code{\link{kin.search}}, \code{\link{kin.free}}, \code{\link{coo_sample}}, \code{\link{coo_smooth}}, \code{\link{smooth_spline}}
}
