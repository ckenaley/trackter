#' @title Evaluate a range of threshold values for creating binary images

#' @description Converts an image to grayscale, applies a user defined threshold range to segment the binary the image, and then plots segmented images that reflect 9 discrete values in the range.
#'
#' @param img character, a path to an image file.
#' @param min numeric, the minimum threshold value (0-1).
#' @param max numeric, the maximum threshold value (0-1).
#' @param otsu logical, should the automatically determined threshold value be printed. See Details.
#' @param plot.grid logical, should the grid of images be plotted.
#'
#' @export
#'
#' @importFrom graphics text
#'
#' @details
#' Displays a grid of 9 images as an R raster graphic, each the result of 9 discrete values within the range describe by \code{min} and \code{max}. If both \code{min} and \code{max} are \code{NULL}, then the threshold range is defined by the default \code{seq(0.1,0.9,.1)}. This function should help users refine threshold values for detecting ROIs with \code{kin.search} and \code{kin.simple} if automatic thresholding using Otsu's method is not satisfactory.
#'
#' @return If 'otsu=TRUE', a single automated threshold value is returned using Otsu's method. See \code{\link{otsu}}.
#'
#' @seealso \code{\link{kin.simple}}, \code{\link{kin.search}}, \code{\link{otsu}}
#'
#'
#' @examples
#'
#'
#' #access image in system
#'  y <-system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")
#'
#' # default values
#' thr.check(y)
#'
#' #custom values
#' thr.check(y,min=0.2,max=1.5)
#'
#'

thr.check <-
  function(img,
           min = NULL,
           max = NULL,
           otsu = TRUE,
           plot.grid = TRUE) {
    if (!grepl(".jpeg|.jpg|.png|.tiff", img, ignore.case = TRUE))
      stop("file in file path doesn't appear to be an image.")
    
    x <- EBImage::readImage(img, all = FALSE)
    EBImage::colorMode(x) = EBImage::Grayscale
    
    ot <- EBImage::otsu(x)[1]
    
    if (plot.grid) {
      if (all(sapply(list(min, max), is.null)))
        thr <- seq(0.1, 0.9, .1) #sequence of threshold values
      
      if (length(which(sapply(list(min, max), is.null) == F)) == 1)
        stop("both 'min' and 'max' must have value=NULL or numeric 0-1")
      
      if (all(!sapply(list(min, max), is.null)) &&
          min > max)
        stop("'min' must be < 'max'")
      
      if (all(!sapply(list(min, max), is.null)) &&
          min == max)
        stop("'min' must be value different from 'max'")
      
      if (all(!sapply(list(min, max), is.null)))
        thr <- seq(min, max, length.out = 9)
      
      img.l <- list()
      for (i in thr) {
        y <- x > i #threshold the image using value i
        img.l[[paste(i)]] <- EBImage::getFrame(y, 1)
      }
      
      img.t <-   EBImage::combine(img.l)
      
      nx <- ny <-  3
      width <- dim(x)[1]
      height <-  dim(x)[2]
      x_offset <-  y_offset <-  20
      n <- EBImage::numberOfFrames(img.t, 'render')
      
      EBImage::display(img.t,
                       method = "raster",
                       all = TRUE,
                       nx = nx)
      text(
        x = rep(seq(
          from = 0,
          by = width,
          length.out = nx
        ), ny) + x_offset,
        y = rep(seq(
          from = 0,
          by = height,
          length.out = ny
        ), each = nx) + y_offset,
        label = thr,
        adj = c(0, 1),
        col = "red",
        cex = 1.5
      )
    }
    if (all(!sapply(list(min, max), is.null))) {
      if (ot > max |
          ot < min)
        message("Otsu value is outside defined threshold range")
    }
    
    if (otsu)
      return(ot)
    
  }

#' @title Crop an image

#' @description Crops an image with a rectangle and saves it
#'
#' @param img character, a path to an image file.
#' @param ul numeric, the upper-left xy position of the rectangle in pixels.
#' @param br numeric, the bottom-right xy position in pixel.
#' @param out.dir character, the directory to which the image will be saved.
#' @param type character, image type (optional). Must be: "jpeg", "png", or "tiff (case ignored)." If missing, file format is automatically determined by file name extension.
#' @param locate logical, if TRUE, the original image is printed to the graphics device and the user picks \code{ul} and \code{br} (in that order) with the cursor and two clicks of the mouse. If TRUE, the coordinates of the points picked by the user are returned.
#' 
#' @export
#'
#' @examples
#'
#'
#' #retrieve image in system
#'  y <-system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")
#' od <- paste0(tempdir(),"/cropimg")
#' dir.create(od)
#' 
#' #display original
#' EBImage::display(EBImage::readImage(y),method="raster")
#' 
#' 
#' # crop and save as orignal format
#' crop.img(img=y,ul=c(5,30),br=c(100,200),out.dir=od)
#' 
#' # crop and save as pnng
#' crop.img(img=y,ul=c(5,30),br=c(200,200),out.dir=od,type="png")
#' 
#'#display cropped images
#'EBImage::display(EBImage::readImage(paste0(od,"/sunfish_BCF.jpg")),method="raster")
#'EBImage::display(EBImage::readImage(paste0(od,"/sunfish_BCF.png")),method="raster")
#'
#' #clean up
#' unlink(od,recursive=TRUE)
#'
#'# use graphics device to choose crop margins
#'\dontrun{
#'crop.img(img=y,out.dir=od,locate=TRUE)
#'}


crop.img <-
  function(img = NULL,
           ul = NULL,
           br = NULL,
           out.dir = NULL,
           type = NULL,
           locate = FALSE) {
    if (!grepl(".jpeg|.jpg|.png|.tiff", img, ignore.case = TRUE))
      stop("file in file path doesn't appear to be an image.")
    
    if (is.null(out.dir))
      stop("`out.dir` is NULL")
    
    if (!dir.exists(out.dir))
      stop("`out.dir` doesn't exhist")
    
    if(any(sapply(list(ul,br), is.null)) && !locate)
      stop("both 'ul' and 'br' must be specified or 'locate=TRUE'")
    
    
    x <- EBImage::readImage(img, all = FALSE)
    
    
    if(!locate) xcrop <- x[ul[1]:br[1], ul[2]:br[2], ]
    
    if (locate) {
      message("select two points, upper left then lower right")
      x2 <- x
      x2[1:dim(x)[1], c(1, dim(x)[2]), ] <-  0
      x2[c(1, dim(x)[1]), 1:dim(x)[2], ]  <- 0
      
      EBImage::display(x2, method = "raster")
      pts <- graphics::locator(2,type = "p")
      message("two points selected. . . well done")
      
      pts <- lapply(pts, function(x)
        round(x, 0))
      xcrop <- x[pts$x[1]:pts$x[2], pts$y[1]:pts$y[2], ]
      
        return(pts)
    }

      if (is.null(type))
        EBImage::writeImage(xcrop, files = paste0(out.dir, "/", basename(img)))
      if (!is.null(type)) {
        n <- gsub('\\..[^\\.]*$', '', basename(img))
        if (!grepl("jpeg|jpg|png|tiff", type, ignore.case = TRUE))
          stop("invalid 'type'")
        
        EBImage::writeImage(xcrop,
                            type = type,
                            files = paste0(out.dir, "/", n, ".", type))
      }
    
  }


#' @title 2D Manual digitization of image data

#' @description  Uses the active graphics device for interactive digitization of shape data. 
#' 
#' @param img character, a one-dimension vector of image file paths.
#' @param n numeric, the max number of points to digitize for each frame (i.e., image file) in \code{img}. The number of points added to each frame may be less than this value. See details.
#' @param review logical. If 'TRUE', the user inputs 'return/enter' after points have been added for each frame. If 'FALSE', the function advances to the next frame without review after the last point has been added. 
#' @param smoothing character, the smoothing method: either 'loess' or 'spline' for line smoothing or 'outline' for contours. See details.
#' @param smooth numeric; if \code{smoothing} is set to 'loess', this values is passed to 'span' parameter of \code{\link{loess}}. If \code{smoothing} is set to 'spline', this value is ignored and the smooth points are determined by the default values of \code{\link{spline}}. See details.
#' @export
#' @importFrom graphics lines locator points
#' @importFrom stats complete.cases fitted lm loess  predict smooth.spline
#'
#' @details
#' Users interact with each frame using the \code{\link{locator}} function. Although the maximum number of points for each frame is determined by \code{n}, the user may advance to the next frame by double-clicking on the last point. If \code{review=TRUE}, the user is asked to review the points and, if satisfied, a "return"/"enter" key stroke will advance the graphics device to the next frame. After double clicking in the last frame and a "return"/"enter" key stroke if \code{review=TRUE}, the results are returned. If \code{smooth=NULL}, only the clicked points are returned. If \code{smooth="spline"} or \code{smooth="loess"}, interpolated smoothed points are returned.
#' 
#'Smoothed points for outlines are produced with \code{\link{coo_smooth}}. Lines are smoothed with \code{\link{smooth.spline}} with \code{smooth} passed to the \code{spar} parameter or with \code{\link{smooth.spline}} with \code{smooth} passed to the \code{spar}
#'
#' @return A list of one or two data tables, alwyas \code{pts} with the picked points and, optionally, \code{smooth.pts} containing the the smoothed points if \code{smoothing} is chosen. Each contains the following columns:
#'
#' \itemize{
#' \item 'frame', the frame file name
#' \item 'frame.n', the frame number
#' \item 'pt', the point number
#' \item x,y coordinates of the selected points
#' }
#'
#' @seealso \code{\link{locator}}, \code{\link{loess}}, \code{\link{smooth.spline}}
#' @export
#'
#' @importFrom graphics lines
#' @importFrom stats complete.cases fitted lm loess  predict smooth.spline pt spline
#' @importFrom utils head  tail flush.console
#'
#' @examples
#' 
#' \dontrun{
#'#download example avi video
#' f <- "https://github.com/ckenaley/exampledata/blob/master/sunfish_pect.avi?raw=true"
#' download.file(f,"sunfish.avi")
#' 
#'#extract images with ffmpeg operations and reduce them to 600 px wide with a filter
#'filt.red <- " -vf scale=600:-1 " #filter
#'dir.create(paste0(tempdir(),"/out"))
#'vid.to.images2(vid.path="sunfish.avi",filt = filt.red,out.dir = paste0(tempdir(),"/out")) 
#'y <- list.files(paste0(tempdir(),"/out"),full.names = TRUE)[1:3]
#'pick.pts(y,review=F)
#'}

pick.pts <- function(img,n=1000,review=TRUE,smoothing=NULL,smooth=NULL){
  y <- NULL
  
xy.l <- list()

xysm.l <- list()

sm <- TRUE

if(is.null(smoothing)) sm <- FALSE

if(sm)if(!smoothing %in% c("loess","spline","outline")) stop("'smoothing' must be set to either 'loess', 'spline', or 'outline'")
i <- 1
while(i<=length(img)){

    x <- EBImage::readImage(img[i], all = FALSE)
    EBImage::display(x, all = FALSE,method="raster")
    
    d <- 10
    xy.lx <- list()
    for(ii in 1:n){
      
      
      xy.i <- unlist(locator(1))
      points(xy.i[1],xy.i[2],col="red",pch=16)
      
      
      xy.lx[[ii]]<- xy.i
      
      if(ii>1){
        pts.x <- data.frame(do.call(rbind,xy.lx))
        d <- dist.2d(pts.x$x[ii],pts.x$x[ii-1],pts.x$y[ii],pts.x$y[ii-1])
        
        with(pts.x,lines(x,y,col="red"))
        if(d<0.5) break
      }
      
   
      }  
  
    fr <- basename(img[i])  
 pt.i <- data.table(frame=fr,frame.n=i,pt=1:ii,do.call(rbind,xy.lx))
 

 if(sm) pt.sm <- pt.i[pt<ii,{
   n.o <- length(min(x):max(x))
   s <- spline(x,y,n=n.o);
   list(
     frame=fr,
     frame.n=i,
     pt=1:length(s$x),
     x=s$x,
     y=s$y
     
   )
 }]
 

 
 if (sm) pt.sm <- pt.i[pt<ii,{
    x.s <- seq(min(x),max(x),length.out = max(x)-min(x));
     s <- predict(loess(y~x,span=smooth),newdata=data.frame(x=x.s));
     list(
       frame=fr,
       frame.n=i,
       pt=1:length(x.s),
       x=x.s,
       y=s
     )
   }]
 


 if (sm){  
   x.s <- with(pt.i,seq(min(x),max(x),length.out = max(x)-min(x)))
   
   pt.out <-Momocs::Out(list(as.matrix(pt.i[pt<ii, list(x, y)])))[[1]][[1]]
   pt.cc <- Momocs::coo_perimcum(pt.out)
   s <- spline(pt.cc,pt.out[,1],n=length(x.s))
   
   pt.spline <- cbind(spline(pt.cc, pt.out[, 1], length(x.s)/2)$y,
                      spline(pt.cc, pt.out[,2], length(x.s)/2)$y)
   
   
   pt.sm <- data.table(
     frame=fr,
     frame.n=i,
     pt=1:length(pt.spline[,1]),
     x=pt.spline[,1],
     y=pt.spline[,2]
   )
 }
 
 pt.i<- pt.i[pt<ii]
 
 if(sm){with(pt.sm,lines(x,y,col="yellow",lwd=1))
   xysm.l[[i]] <- pt.sm
 }else{ xysm.l[[i]] <- NULL}
 

 xy.l[[i]] <- pt.i
  

if(review){  
cat("Are you happy with these points? 'return/enter' to advance, any other key to go back.") 
  ln <- readline()
  if(ln=="") {i <- i+1; flush.console()} else{i <- i}
}else{i <- i+1}
 
}

pts <- data.table(do.call(rbind,xy.l))
pts.sm <- data.table(do.call(rbind,xysm.l))
return(list(pts=pts,smooth.pts=pts.sm) )
}
