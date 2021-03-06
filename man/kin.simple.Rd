% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kin.R
\name{kin.simple}
\alias{kin.simple}
\title{Simplified midline tracking over image sequences}
\usage{
kin.simple(
  image.dir = NULL,
  frames = NULL,
  thr = 0.7,
  size.min = 0.05,
  ant.per = 0.2,
  tips = 0.02,
  smoothing = "loess",
  smooth = 0.25,
  smooth.points = 200,
  save = TRUE,
  out.qual = 1,
  out.dir = NULL,
  plot.pml = TRUE,
  image.type = "orig",
  flip = TRUE,
  show.prog = FALSE
)
}
\arguments{
\item{image.dir}{character, directory containing images to analyze.}

\item{frames}{numeric, vector indicating which images to process.}

\item{thr}{numeric or character ('otsu') threshold to determine binary image. See Details.}

\item{size.min}{numeric, indicating the minimum size of ROIs as a proportion of the pixel field to be considered in analysis. May be useful if smaller unimportant ROIs appear in the frame. Default is 0.05.}

\item{ant.per}{numeric; left-most proportion of ROI that establishes the horizontal reference for the midline displacement.}

\item{tips, }{numeric, the proportion the the midline data to use in calculation of the head and tail position.}

\item{smoothing}{character, the midline smoothing method, either 'loess' or 'spline'.}

\item{smooth}{numeric; if \code{smoothing} is set to 'loess', passed to 'span' parameter of \code{\link{loess}}. If \code{smoothing} is set to 'spline', passed to 'spar' parameter of \code{\link{smooth.spline}}}

\item{smooth.points}{numeric, number of equally spaced points along the ROI midline on which the smoothed midline is computed.}

\item{save}{logical, value indicating if images should be outputted with midline and predicted midline based on the \code{lm()} predictions from \code{ant.per}overlaying original or binary images.}

\item{out.qual, }{numeric, a value between 0-1 representing the quality of outputted images. Ignored if \code{save=FALSE}.}

\item{out.dir}{character, the directory to which outputted images should be saved.}

\item{plot.pml}{logical, value indicating if outputted images should include the predicted midline (in blue) and the points according to \code{ant.per} used to construct the predicted midline (in green).}

\item{image.type}{character; the type of image to be outputted, either 'orig' or 'bin' representing the original or binary images, respectively. Ignored if 'save=FALSE'.}

\item{flip}{logical, indicating if binary image should be flipped.}

\item{show.prog}{logical, indicating if outputted image should be displayed during analysis.}
}
\value{
A list with the following components:

\code{kin.dat} a data table consisting of frame-by-frame position parameters for the ROI determined by LDA analysis.
\itemize{
\item the frame number
\item 'x' and 'y': the position of the tail (rightmost or posteriormost)
\item 'head.x' and 'head.y': the x and y position of the head (leftmost or anteriormost)
\item 'amp': the amplitude (\code{amp}) of the tail relative to the theoretical midline determined by the \code{lm()} predictions from \code{ant.per}
\item 'roi': a character indicating the ROI ranked by size ('a' being the largest)
\item 'head.pval': p values of the \code{lm()} fit that describes the position of the head as determined by \code{ant.per} (green points in the outputted images/video)}

\code{midline} A data table containing, for each frame described by \code{frames}, the following: \itemize{
\item 'x' and 'y.m': x and y positions of the midline of the ROI
#' \item 'y.min' and 'y.max': min and max y positions ROI's contour used in y.m calculation
\item 'mid.pred': the predicted linear midline based on the points/pixels defined by \code{ant.per} (green points in the outputted images/video if 'plot.pml=TRUE')
\item 'y.pred': midline points fit to a smooth spline or loess model with spar or span equal to \code{smooth} (red curve in the outputted images/video)
\item 'wave.y': midline points 'y.pred' relative to 'mid.pred'
\item 'roi': a character indicating ROI size ('a' being the largest)
}

\code{cont} A data table containing x and y positions of the contours used to calculate the data in 'kin.dat'. Contains the following: 
\itemize{
\item 'frame': the frame
\item 'x' and 'y': the x and y positions of the contours
}

\code{all.classes} A data table containing the following for all ROIs detected:  
\itemize{
\item 'frame': the frame
\item 'roi': the name of each ROI found in a frame.
\item 'size': the size of each ROI
}

\code{dim} the x and y dimensions of the images analyzed
}
\description{
Automatically retrieves the midline of a detected ROI based on size. Assumes the ROI of interest is the largest detected and not intersecting the edges of the image frame, conditions often met in kinematic studies. For each ROI of interest, finds the y-value midpoint along the x-value array of the ROI and fits a midline according to a chosen smoothing method (loess or spline). Also outputs the midline amplitude relative to a reference line determined by an anterior section of the ROI and outputs contours ROIs in each frame for subsequent analysis. Supported image formats are jpeg, png, and tiff.
}
\details{
The algorithm assumes a left-right orientation, i.e., the head of the ROI is positioned left, the tail right. ffmpeg operations or even imageJ can rotate images not in this orientation. The \code{ant.per} value therefore establishes the reference line (theoretical straight midline) based on that portion of the head. The midline is calculated as the midpoints between the y extrema for each x position.  

If 'save=TRUE', images are saved as binary or the original with a body midline  overlay and, if chosen, with the theoretical midline (based on \code{ant.per}). 

Thresholding operations can be performed with an arbitrary (user defined) numeric value or with Otsu's method ('thr="otsu"'). The latter chooses a threshold value by minimizing the combined intra-class variance. See \code{\link{otsu}}.
}
\examples{
#### plot caudal amplitude and produce a classic midline waveform plot of swimming fish	
##A very long example.
\dontrun{	

#download example images and place in 'example' subdirectory	
f <- "https://github.com/ckenaley/exampledata/blob/master/example.zip?raw=true"	

download.file(f, paste0(tempdir(),"/temp.zip"))
unzip(paste0(tempdir(),"/temp.zip"), exdir=tempdir())	
unlink(paste0(tempdir(),"/temp.zip"))	

dir.create(paste0(tempdir(),"/processed_images"))
kin <- kin.simple(image.dir =paste0(tempdir(),"/example"),	
      smoothing = "loess",frames=1:50,
      out.dir=paste0(tempdir(),"/processed_images"),
      show.prog = FALSE,thr = "otsu",	
      image.type="bin",smooth=0.4)	

#plot instantaneous amplitude of tail (last/rightmost point) over frames 	
p <- ggplot(dat=kin$kin.dat,aes(x=frame,y=amp))+geom_line()+geom_point()+theme_classic(15)	
print(p)	
	
# midline plot	
ml <- kin$midline	
#leftmost x starts at 0	
ml <- ml[,x2:=x-x[1],by=frame]	

ml <- merge(ml,kin$kin.dat[,list(frame,amp)],by="frame") #merge these	

pal <- wes_palette("Zissou1", 100, type = "continuous") #"Zissou" color palette

p <- ggplot(dat=ml,aes(x=x2,y=wave.y))+theme_classic(15)+scale_color_gradientn(colours = pal)	
p <- p+geom_line(aes(group=frame,color=amp),stat="smooth",method = "loess", size = 1.5)	
print(p)	
	
#Make a video of processed frames	

images.to.video2(image.dir=paste0(tempdir(),"/processed_images"),	
vid.name="trout_test",out.dir=tempdir(),frame.rate=5,qual=100,raw=FALSE)	
file.exists(paste0(tempdir(),"/trout_test_red.mp4"))

}

## A very short example.

#retrieve image 
i <- EBImage::readImage(system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter"))
#create directory and write image to it
t <-tempdir()
dir.create(paste0(t,"/images"))
EBImage::writeImage(i,paste0(t,"/images/sunfish001.jpg"),type = "jpeg")

#run kin.search and save output image to directory
kin.i<- kin.simple(image.dir = paste0(t,"/images"),save = TRUE,out.dir = t)

#plot midline
with(kin.i$midline,plot(x,wave.y))
i2 <- EBImage::readImage(paste0(t,"/sunfish001_000.jpg"))
EBImage::display(i2,method="raster")
#clean up
unlink(paste0(t,"/images"),recursive=TRUE)
}
\seealso{
\code{\link{kin.search}}
}
