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

#' @description Crops an image or images with a rectangle
#'
#' @param img character, a path to an image file or directory containing images.
#' @param ul numeric, the upper-left xy position of the rectangle in pixels.
#' @param br numeric, the bottom-right xy position in pixel.
#' @param out.dir character, the directory to which the image will be saved.
#' @param type character, image type (optional). Must be: "jpeg", "png", or "tiff (case ignored)." If missing, file format is automatically determined by file name extension(s).
#' @param locate logical, if TRUE and the 'img' is a single file path, the original image is printed to the graphics device and the user picks \code{ul} and \code{br} (in that order) with the cursor and two clicks of the mouse. If TRUE, the coordinates of the points picked by the user are returned. Ignored is 'img' is a directory containing images.
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
#' # crop and save as orignal format
#' crop.img(img=y,ul=c(5,30),br=c(100,200),out.dir=od)
#' 
#' # crop and save as png
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
    
    one <- length(list.files(img))==0
    many <- length(list.files(img))>0
    
    if (is.null(out.dir))
      stop("`out.dir` is NULL")
    
    if (!dir.exists(out.dir))
      stop("`out.dir` doesn't exist")
    
    if(one){
    if (!grepl(".jpeg|.jpg|.png|.tiff", img, ignore.case = TRUE))
      stop("file in file path doesn't appear to be an image.")
  
    
    if(any(sapply(list(ul,br), is.null)) && !locate)
      stop("both 'ul' and 'br' must be specified or 'locate=TRUE'")
  
    x <- EBImage::readImage(img, all = FALSE)
    
    
    if(!locate){
      
      dim.x <- dim(x)
      
      if(!all(ul[1]<=dim.x[1],ul[2]<=dim.x[2],br[1]<=dim.x[1],br[2]<=dim.x[2]))
        stop("cropping region set by 'ul' or 'br' beyond image bounds")
      
      xcrop <- x[ul[1]:br[1], ul[2]:br[2], ]
    }
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
    
    if(many){
      if (any(!grepl(".jpeg|.jpg|.png|.tiff", list.files(img), ignore.case = TRUE)))
        stop("not all files in the directory appear to be an image.")
     
      
      if(any(sapply(list(ul,br), is.null)))
        stop("both 'ul' and 'br' must be specified")
     
       imgs <- list.files(img,full.names = TRUE)
      
      for(i in imgs){
      x <- EBImage::readImage(i, all = FALSE)
      dim.x <- dim(x)
      if(!all(ul[1]<=dim.x[1],ul[2]<=dim.x[2],br[1]<=dim.x[1],br[2]<=dim.x[2]))
        stop("cropping region set by 'ul' or 'br' beyond image bounds")
      
      xcrop <- x[ul[1]:br[1], ul[2]:br[2], ]
      
      
      if (is.null(type))
        EBImage::writeImage(xcrop, files = paste0(out.dir, "/", basename(i)))
      if (!is.null(type)) {
        n <- gsub('\\..[^\\.]*$', '', basename(i))
        if (!grepl("jpeg|jpg|png|tiff", type, ignore.case = TRUE))
          stop("invalid 'type'")
        
        EBImage::writeImage(xcrop,
                            type = type,
                            files = paste0(out.dir, "/", n, ".", type))
      }
      
    }
    }
}

#' @title Change image contrast

#' @description Changes the contrast of an image or images
#'
#' @param img character, a path to an image file or directory containing images.
#' @param c numeric, the factor by which to change the contrast. 'c'>1 enhances, 'c'<1 reduces.
#' @param out.dir character, the directory to which the image will be saved.
#' @param type character, image type (optional). Must be: "jpeg", "png", or "tiff (case ignored)." If missing, file format is automatically determined by file name extension(s).
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
#' # adjust contrast and save as original format
#' contrast.img(img=y,c=0.5,out.dir=od)
#' 
#'#display modified image
#'EBImage::display(EBImage::readImage(paste0(od,"/sunfish_BCF.jpg")),method="raster")
#'
#' #clean up
#' unlink(od,recursive=TRUE)


contrast.img <-
  function(img = NULL,
           c = 2,
           out.dir = NULL,
           type = NULL) {
    
    one <- length(list.files(img))==0
    many <- length(list.files(img))>0
    
    if (is.null(out.dir))
      stop("`out.dir` is NULL")
    
    if (!dir.exists(out.dir))
      stop("`out.dir` doesn't exist")
    
    if(one){
      if (!grepl(".jpeg|.jpg|.png|.tiff", img, ignore.case = TRUE))
        stop("file in file path doesn't appear to be an image.")

      
      x <- EBImage::readImage(img, all = FALSE)
      
      x2 <- x*c
      
  
      if (is.null(type))
        EBImage::writeImage(x2, files = paste0(out.dir, "/", basename(img)))
      if (!is.null(type)) {
        n <- gsub('\\..[^\\.]*$', '', basename(img))
        if (!grepl("jpeg|jpg|png|tiff", type, ignore.case = TRUE))
          stop("invalid 'type'")
        
        EBImage::writeImage(x2,
                            type = type,
                            files = paste0(out.dir, "/", n, ".", type))
      }
      
    }
    
    if(many){
      if (any(!grepl(".jpeg|.jpg|.png|.tiff", list.files(img), ignore.case = TRUE)))
        stop("not all files in the directory appear to be an image.")
      
      
      imgs <- list.files(img,full.names = TRUE)
      
      for(i in imgs){
        x <- EBImage::readImage(i, all = FALSE)
        
        x2 <- x*c
        
        if (is.null(type))
          EBImage::writeImage(x2, files = paste0(out.dir, "/", basename(i)))
        if (!is.null(type)) {
          n <- gsub('\\..[^\\.]*$', '', basename(i))
          if (!grepl("jpeg|jpg|png|tiff", type, ignore.case = TRUE))
            stop("invalid 'type'")
          
          EBImage::writeImage(x2,
                              type = type,
                              files = paste0(out.dir, "/", n, ".", type))
        }
        
      }
    }
  }


#' @title 2D Manual digitization of image data

#' @description  Uses the active graphics device for interactive digitization of shape and position data. 
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
#' @return A list of one or two data tables, always \code{pts} with the picked points and, optionally, \code{smooth.pts} containing the the smoothed points if \code{smoothing} is chosen. Each contains the following columns:
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



#' @title Plot data over an image

#' @description  Simple wrapper for \code{\link{points}} to plot data derived from an image over that image.
#' 
#' @param img character, an image file path.
#' @param over a data frame, data table, or matrix with two columns, the first representing x and the second y coordinates. 
#' @param ... other arguments to be passed to \code{\link{points}}
#' 
#' @details 
#' Simply plots 2D dimensional data over the image specified in \code{img}. May be useful in plotting output from \code{trackter}'s kin functions.
#' 
#' @return The image is plotted in the graphics device with points given in \code{x}
#'
#' @seealso \code{\link{kin.free}}, \code{\link{kin.search}}, \code{\link{kin.simple}}

#' @export
#'
#' @importFrom graphics points
#'
#' @examples
#' 
#' f <-system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")
#' 
#' d <- dim(EBImage::readImage(f))
#' 
#' x <- runif(10,1,d[1])
#' y <- runif(10,1,d[2])
#' pts <- cbind(x,y)
#' 
#' data.overlay(img=f,over=pts,col="red",type="p")

data.overlay <- function(img,over,...){
  img<-  EBImage::readImage(img)
  suppressMessages(EBImage::display(img, method = "raster"))
  points(unlist(over[,1]),unlist(over[,2]),...)
}

#' @title Plot output from \code{kin} functions

#' @description  Plots data output from one of \code{trackter}'s \code{kin} functions using \code{ggplot2}. Useful for quickly visualizing results from \code{kin} functions. Data and geometries include body-contour polygons with midline points, head or fin points, both midline and head points, or none or all of these. Optionally builds and saves an animation of data with reference to frame from \code{kin} data. 
#' 
#' @param kin a list of data returned from \code{kin} functions. See Details.
#' @param frames integer, a vector for frames within the data (a column common to all \code{kin} functions' output). Will subset the data or, if NULL, will not. If more than one frame is specified, plot will be faceted by frame with \code{facet_wrap} or animated across frames.
#' @param under character, the data to plot under \code{over}, the contour(s) from \code{kin} functions, i.e, 'cont' or 'cont.sm'. See Details.
#' @param over character, the data to plot over the contour(s). Must be either 'kin.dat', 'midline', 'both', 'fin', 'fin.pt', 'all', or 'none'. See Details.
#' @param zoom logical, should the plotted area zoom to the extents of the contours. See Detail.
#' @param animate logical, should an animation be plotted to the graphics device with state changes reflected by \code{frame}.
#' @param shadow logical,should \code{\link{shadow_mark}} be implemented to draw geometries in previous frames. Ignored if \code{animate=FALSE}.
#' @param alpha numeric, the opacity of geometries in previous frames during animation. Ignored if \code{animate=FALSE}.
#' @param fps integer, the play back speed of animation in frames per second. Ignored if \code{animate=FALSE}.
#' @param save logical, should the animation be saved as a GIF. Ignored if \code{animate=FALSE}. See Details.
#' @param filename character, the name given to the GIF file. Ignored if \code{animate=FALSE} or \code{save=FALSE}.
#' @param out.dir character, the file path of the directory to which the GIF file is saved. Ignored if \code{animate=FALSE} or \code{save=FALSE}.
#' @param ... other arguments to be passed to midline and head point geometries, e.g. 'size', 'color'.
#' 
#' @details 
#' Simply plots 2D dimensional data over contours retrieved from \code{kin} functions. May be useful in quickly assesses their results. The \code{under} layer can be one of the named contour data tables in lists returned by the \code{kin} functions \code{\link{kin.search}},\code{\link{kin.simple}}, \code{\link{kin.free}}, or \code{\link{fin.kin}}. The overlayed data layer specified by \code{over} are non-contour data from the \code{kin} functions.
#' 
#' If the list specified by 'kin' is from \code{\link{kin.search}},\code{\link{kin.simple}}, or  \code{\link{kin.free}}, 'over' must be 'cont' or 'cont.sm' and under one of 'midline','kin', 'none' or 'both'. 'midline' will produce the smoothed midline coordinates.
#' 
#' If the list specified by 'kin' is from \code{\link{fin.kin}}, 'over' must be 'cont' or 'comp' and under one of 'midline','fin', 'fin.pts', 'none' or 'all'.
#' 
#' 
#' Animations are produced with \code{gganimate} and optionally saved as a GIF with \code{anim_save} using the default \code{gifski_renderer()}.
#' 
#' 
#' 
#' @return A \code{ggplot} printed to the graphics device, a \code{gganimate} object plotted to graphics device, or a GIF file saved to a local directory specified in \code{out.dir}.
#'
#' @seealso \code{\link{geom_point}}, \code{\link{geom_polygon}},\code{\link{gganimate}}, \code{\link{kin.free}}, \code{\link{kin.search}}, \code{\link{kin.simple}},\code{\link{fin.kin}}

#' @export
#' 
#'
#' @import ggplot2 gganimate animation
#'
#' @examples
#' 
#' \dontrun{
#' #animate a ropefish swimming with its midline
#' 
#' #download example avi and place in subdirectory
#' f <- "https://github.com/ckenaley/exampledata/blob/master/ropefish.avi?raw=true"
#'
#' download.file(f, paste0(tempdir(),"/ropefish.avi"))
#'
#' dir.create(paste0(tempdir(),"/images"))
#' 
#' vid.to.images(paste0(tempdir(),"/ropefish.avi"), out.dir = paste0(tempdir(),"/images"))
#' 
#' kin <- kin.free(image.dir =paste0(tempdir(),"/images"),
#' par=TRUE,
#'       ml.smooth=list("spline",0.9),
#'       thr = "otsu",
#'       size.min=0.01,
#'       red=0.5
#'       )
#'
#' gg.overlay(kin=kin,
#' frames=seq(20,320,20),
#' under="cont.sm",
#' over="midline",
#' size=1,
#' animate=TRUE,
#' col="red",
#' fps=10)
#' 
#' #clean up 
#' unlink(paste0(tempdir(),"/images"),recursive=TRUE)
#'}
#'
gg.overlay <- function(kin=NULL,frames=NULL,under="cont",over="midline",zoom=FALSE,animate=FALSE,shadow=TRUE,alpha=0.1,fps=10,save=FALSE,filename=NULL,out.dir=NULL,...){
  
  
  frame <- x <- y <- head.x <- head.y <- x.sm <- y.sm <-pos <- side <- NULL
  
  if(is.null(frames)) frames <-unique(kin[[1]]$frame)
  
if(!"fin" %in% names(kin))  stopifnot(class(kin)=="list",c("midline","kin.dat","cont.sm") %in% names(kin),frames %in% kin$cont.sm$frame,!is.null(kin))

  if("fin" %in% names(kin))  stopifnot(class(kin)=="list",c("cont","fin","comp","fin.pts") %in% names(kin),frames %in% kin$cont$frame,!is.null(kin))
  an <- animate
  
  if(save&is.null(animate)&is.null(out.dir)) stop("'save=TRUE' but 'out.dir' is not specified")
  if(!is.null(out.dir))if(!dir.exists(out.dir)) stop("'out.dir' does not exists")
  
  z <- kin$cont[frame%in%frames,list(x=range(x)*c(0.9,1.1),y=range(y)*c(0.9,1.1))]
  
  if(is.null(kin$dim)) kin$dim <-c(max(z$x)*1.25,max(z$y)*1.25)
  cart <- list(xlim(c(0,kin$dim[1])),ylim(c(kin$dim[2],0)))
  
  if(zoom & !animate) cart <- coord_cartesian(xlim=z$x,ylim=z$y[2:1]) 
  c<- kin$cont[frame%in%frames]
  
  #midline data
  m <- kin$midline[frame%in%frames]
  if(!"x.sm" %in% colnames(m)) m[,x.sm:=x]
  
  if(!"fin" %in% names(kin)) {
  c <- kin$cont[frame%in%frames]
  if(under=="cont.sm") c <- kin$cont.sm[frame%in%frames]
  k <- kin$kin.dat[frame%in%frames]
  g.k <- geom_point(data=k,aes(head.x,head.y),...)
  g.m <- geom_point(data=m,aes(x.sm,y.sm),...)

  if(over=="kin.dat") g <- g.k
  if(over=="midline") g <- g.m
 
  if(over=="both") { g.m2 <- geom_point(data=m,aes(x.sm,y.sm),col="blue",...)
  g.k2 <- geom_point(data=k,aes(head.x,head.y),col="red",size=g.m2$aes_params$size+1)
  g <- list(g.m2,g.k2)
  }
  
  }
  
  if("fin" %in% names(kin) ) {
    c <- kin$cont[frame%in%frames]
    if(under=="comp") c <- kin$comp[frame%in%frames]
    
    fp <- kin$fin.pts[frame%in%frames]
    f <- kin$fin[frame%in%frames]
    g.f <- geom_point(data=f,aes(x,y,col=side),...)
    g.fp <- geom_point(data=fp,aes(x,y,col=side,shape=pos),...)
    g.m <- geom_point(data=m,aes(x,y),...)
    
    if(over=="fin") g <- g.f
    if(over=="fin.pts") g <- g.fp
    if(over=="midline") g <- g.m
    
    if(over=="all") { g.m2 <- geom_point(data=m,aes(x,y),col="blue",...)
    g.f2 <- geom_point(data=f,aes(head.x,head.y))
    g.fp2 <- geom_point(data=fp,aes(head.x,head.y))
    g <- list(g.m2,g.f2,g.fp2)
    }
    
  }
  
  if(length(frames)>1 | is.null(frames)) facet <- facet_wrap(frame~.) else facet <- NULL

  
  if(over=="none") g <- NULL
  p<- ggplot()+geom_polygon(data=c,aes(x,y))+g+cart+theme_void()+facet
  
  ta <- alpha
  
  sh <- NULL
  if(shadow) sh<- gganimate::shadow_mark(alpha = ta)
  p.anim <-  ggplot()+geom_polygon(data=c,aes(x,y))+g+theme_void()+cart+sh+transition_time(as.integer(frame)) +labs(title = "frame: {frame_time}")
  
  FPS <- fps
  fn <- filename
  if(animate) {p <- gganimate::animate(p.anim, fps=FPS,height = kin$dim[2], width = kin$dim[1]) 
  
  if(save) gganimate::anim_save(filename = fn,animation = p,path=out.dir)
  }
  print(p)
  
}


