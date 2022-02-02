#' @title  Midline and outline tracking over image sequences

#' @description  Wrapper functions for \code{\link{find.roi}} that automatically retrieve the contour and midline coordinates of a detected ROI in each image of a sequence through thresholding and segmentation. Also outputs the midline amplitude relative to a reference line determined by an anterior section of the ROI. Supported image formats are jpeg, png, and tiff.
#' 
#' 
#'\code{kin.search} and \code{kin.simple} find the y-value midpoint along the x-value array of the ROI and fits a midline according to a chosen smoothing method (loess or spline). Thus, these assume a horizontal position (see Details).
#' 
#' \code{kin.search} and \code{kin.free} include arguments for flexible ROI selection.
#' 
#' \code{kin.simple} is itself a wrapper for \code{kin.search}, finding the largest ROI in field using Otsu thresholding for segmentation.
#' 
#' \code{kin.free} does not assume any particular orientation and is intended for finding ROIs freely moving within the image field. This function estimates midlines by various methods and supports parallel processing of frames (see Details). 
#' 
#' @param image.dir character, directory containing images to analyze.
#' @param frames numeric, vector indicating which images to process. NULL, the default, will result in all images in \code{image.dir} processed.
#' @param par logical, should the frames be processed in parallel using \code{cores.n}.
#' @param cores.n numeric, the number of CPU cores to use if \code{par=TRUE}. If \code{cores.n=NULL} (the default), the total number of cores minus 1 are used.
#' @param ant.per numeric; left-most percentage of ROI that establishes the horizontal reference for the midline displacement.
#' @param ant.pos character, one of NULL, "l","r","u",or "d" to specify the position of the anterior of the ROI. If not NULL, the default algorithm to find the anterior is overridden. See Details.
#' @param tips, numeric, the proportion the the midline data to use in calculation of the head and tail position.
#' @param tips, numeric, the proportion the the midline data to use in calculation of the head and tail position.
#' @param smooth.n, numeric, the number of contour smoothing iterations. See Details.
#' @param red numeric, between 0-1 the proportion of contour coordinates to sample for midline estimates. Ignored if \code{ml.meth} is not 'del'. Will speed up midline estimations with Delaunay triangulation. If 'NULL', the full contour retrieved from the ROI will be passed to \code{\link{free.ml.del}}. See Details.
#' @param red numeric, between 0-1 the proportion of contour coordinates to sample for midline estimates. Ignored if \code{ml.meth} is not 'del'. Will speed up midline estimations with Delaunay triangulation. If 'NULL', the full contour retrieved from the ROI will be passed to \code{\link{free.ml.del}}. See Detail.
#' @param ml.meth character, the midline detection method. One of 'ang' for bisection using \code{\link{free.ml.ang}}, 'hull' for bisection using \code{\link{free.ml.hull}}, or 'del' for Delaunay triangulation using \code{\link{free.ml.del}}. See Details.
#' @param smooth.n, numeric, the number of contour smoothing iterations. See Details.
#' @param ml.meth character, the midline detection method. One of 'ang' for bisection using \code{\link{free.ml.ang}}, 'hull' for bisection using \code{\link{free.ml.hull}}, or 'del' for Delaunay triangulation using \code{\link{free.ml.del}}. See Details
#' @param ml.smooth a list of length two with unnamed components including a character string specifying the midline smoothing method, either 'loess' or "spline", and a numeric value specifying the amount of smoothing. See Details.
#' @param save logical, value indicating if images should be outputted with midline and predicted midline based on the \code{ant.per} \code{lm()} overlaying original or binary images.
#' @param plot.pml logical, value indicating if outputted images should include an overlay of the theoretical midline based on \code{ant.per}.
#' @param out.qual, numeric, a value between 0-1 representing the quality of outputted images. Ignored if \code{save=FALSE}.
#' @param out.dir character, the directory to which outputted images should be saved.
#' @param ..., other parameters passed to \code{link{find.roi}}.
#'
#' @export
#'
#' @details
#' 
#'The algorithms in \code{kin.simple} and \code{kin.search} assume a left-right horizontal orientation, i.e., the head of the ROI is positioned left, the tail right. If this is not the case, consider using \code{kin.free} or rotating images before processing. The \code{ant.per} value therefor establishes the reference line (theoretical straight midline) based on that portion of the head. The midline is calculated as the midpoints between the y extrema for each x position. 
#'
#'\code{kin.search} and \code{kin.free} choose ROIs based on relative ROI size or position according to \code{\link{find.roi}}. Parameters  for this function are passed through additional arguments with \code{...}. Thresholding operations can be performed with an arbitrary (user defined) numeric value or with Otsu's method ('thr="otsu"'). The latter chooses a threshold value by minimizing the combined intra-class variance. See \code{\link{otsu}}. 
#'
#' \code{kin.simple} is more streamlined than \code{kin.search}. It attempts to find the largest ROI using Otsu thresholding and invokes other default values of \code{\link{find.roi}}.
#'
#'
#' With \code{kin.free}, the position of the anterior of the ROI (that which is moving forward in the field) is determined by the displacement of the ROI between the first two frames. Thus, \code{frames} must be >1. For analyses of relatively static ROIs in the field (e.g., steadily swimming animals in flumes, etc.), automatically determining the anterior of the ROI may be spurious. In this case, the default automatic determination of the anterior should be overridden by specifying 'l', 'r', 'u', 'd' with \code{ant.pos}. These values specify that the anterior region of the ROI is leftmost, rightmost, upmost, or downmost in the field, respectively, and assumes that the origin of the field (0,0) is the upper left corner of each frame. 
#' 
#' Midline estimation in \code{kin.free} is pursued by one of three algorithms: bisection of contours across the long axis defined by the tips using \code{\link{free.ml.ang}} or \code{\link{free.ml.hull}}  or by Delaunay triangulation using \code{\link{free.ml.del}}. The default is 'hull' This choice is not arbitrary. The use of \code{free.ml.ang} and \code{free.ml.hull} can be faster, but perform poorly for tips that snake back on themselves (i.e., a high degree curvature). The use of \code{free.ml.del} can be slower for high resolution outlines, but  produces better results when contour regions overlap (i.e, those that snake back on themselves), but produces less precise midlines for complicated contours. Using Delaunay triangulation can be hastened (but possibly with a trade off in precision) by reducing the the complexity of the contour with the 'red' argument. For example, a contour of 1000 coordinates would be reduced to one of 500 with 'red=0.5'.
#' 
#' For midline smoothing, if \code{ml.smooth} contains 'spline', \code{\link{smooth_spline}} from the \code{smoothr} package is used to interpolate points between a reduced number of vertices using piecewise cubic polynomials. The number of vertices is calculated based on the number of midline coordinates times numeric value of the list in \code{ml.smooth}. If \code{ml.smooth} contains 'loess', \code{loess} is used to fit a polynomial surface. For contours that have a complicated midline with non-unique x values, say an orginisms swimming vertically in the file, loess smoothing can produce poor results. Thus, spline smoothing is usually the advisable option. 
#' 
#' For contour smoothing before midline estimate in \code{kin.free}, \code{smooth.n} is passed to the \code{smooth.n} parameter of \code{\link{free.ml.ang}}, \code{\link{free.ml.hull}}, or \code{\link{free.ml.del}}, which smooths coordinates using a simple moving average. Contours are similarly smoothed in \code{kin.search} and \code{kin.simple} by invoking \code{\link{coo_smooth}} from the \code{\link{Momocs}} package.  Users should be wary of oversmoothing by smoothing both the contour (from which the midline is calculated) and the midline.
#'
#'
#' @return A list with the following components:
#'
#' \code{kin.dat} a data table consisting of frame-by-frame position parameters for the ROI determined by \code{search.for}.
#' \itemize{
#' \item the frame number
#' \item 'roi': a character indicating the ROI ranked by size ('a' being the largest)
#' \item 'x' and ''y': the position of the tail (rightmost or posteriormost)
#' \item 'head.x' and 'head.y': the x and y position of the head (leftmost or anteriormost)
#' \item 'amp': the amplitude (\code{amp}) of the tail relative to thr theoretical midline determined by the \code{lm()} predictions from \code{ant.per}
#' \item 'head.pval': p values of the \code{lm()} fit that describes the position of the head as determined by \code{ant.per} (green points in the outputted images/video)
#' }
#'
#' \code{midline} A data table containing, for each frame described by \code{frames}, the following:
#' \itemize{
#' \item the frame number
#' \item 'roi': a character indicating ROI size ('a' being the largest)
#' \item 'n': the index of the points where n=1 is headmost
#' \item 'x' and 'y': unsmoothed x and y positions of the midline of the ROI
#' \item 'x.sm' and 'y.sm': midline points fit to a smooth spline or loess model with spar or span equal to \code{smooth} (red curve in the outputted images/video)
#' \item 'wave.y': midline points 'y.pred' relative to 'mid.pred'
#' \item 'per.bl': the percentage of 'x.sm' along the body length calculated as the cumulative sum of distances between points
#' }
#'
#' \code{cont} A data table containing x and y positions of the contours used to calculate the data in 'kin.dat'. Contains the following:
#' \itemize{
#' \item 'frame': the frame
#' \item 'x' and 'y': the x and y positions of the contours
#' }
#'
#' \code{cont.sm} A data table containing the smoothed x and y positions of the contours used to calculate the data in 'kin.dat'. Contains the following:
#' \itemize{
#' \item 'frame': the frame
#' \item 'n': the position of the coordinate. n=1 and n=max(n) are adjacent at the head
#' \item 'x' and 'y': the x and y positions of the contours
#' }
#'
#' \code{all.classes} A data table containing the following for all ROIs detected:
#' \itemize{
#' \item 'frame': the frame
#' \item 'roi': the name of each ROI found in a frame.
#' \item 'edge': indicating whether ROI was on the edge of the image field
#' \item 'size': size of the ROI in pixels^2
#' \item 'offset.x': ROI distance from horizontal center
#' \item 'offset.y': ROI distance from vertical center
#' \item 'offset': linear distance of ROI's centroid to image center
#' }
#'
#' \code{mid.pred} the theoretical midline based on a linear model established by the anterior section of the smoothed midline established by \code{ant.per}. Used to calculate \code{midline$wave.y} as the orthogonal distance between the line defined by 'x' and 'mid.pred' and each coordinate defined by '\code{midline$x.sm} and \code{midline$y.sm}. A data table that contains the following:
#' 
#' \itemize{
#' \item 'frame': the frame
#' \item 'x': x position of the predicted midline
#' \item 'mid.pred': the y position of the predicted midline
#' }
#' 
#' \code{dim} the x and y dimensions of the images analyzed
#'
#' @export
#'
#' @import data.table
#'
#' @importFrom graphics lines frame
#' @importFrom stats complete.cases fitted lm loess predict smooth.spline
#' @importFrom utils head tail
#' @importFrom grDevices dev.off jpeg
#' @importFrom stats complete.cases fitted lm loess predict smooth.spline coef dist 
#' @importFrom utils head tail flush.console
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel 
#' @importFrom foreach foreach %do% %dopar%
#'
#' @seealso \code{\link{kin.free}}, \code{\link{find.roi}}
#' 
#' @examples
#'
#' #### plot caudal amplitude and produce a classic midline waveform plot of a swimming rainbow trout
#' library(data.table)
#' ##A very long example using kin.search()
#' \dontrun{
#' 
#' #download example images and place in 'example' subdirectory
#' f <- "https://github.com/ckenaley/exampledata/blob/master/example.zip?raw=true"
#'
#' download.file(f, paste0(tempdir(),"/temp.zip"))
#' unzip(paste0(tempdir(),"/temp.zip"), exdir=tempdir())
#' unlink(paste0(tempdir(),"/temp.zip"))
#'
#' dir.create(paste0(tempdir(),"/processed_images"))
#' kin <- kin.search(image.dir =paste0(tempdir(),"/example"),
#'       frames=1:50,
#'       out.dir=paste0(tempdir(),"/processed_images"))
#'       
#' 
#' #plot instantaneous amplitude of tail (last/rightmost point) over frames
#' library(ggplot2)
#' p <- ggplot(dat=kin$kin.dat,aes(x=frame,y=amp))+geom_line()+geom_point()+theme_classic(15)
#' print(p)
#'
#' # midline plot
#' ml <- kin$midline
#' 
#' #leftmost x starts at 0
#' ml <- ml[,x2:=x-x[1],by=frame]
#'
#' ml <- merge(ml,kin$kin.dat[,list(frame,amp)],by="frame") #merge these
#'
#' p <- ggplot(dat=ml,aes(x=x2,y=wave.y))+theme_classic(15)
#' p <- p+geom_line(aes(group=frame,color=amp),stat="smooth",method = "loess", size = 1.5)
#' print(p)
#'
#'}
#'
#' ## A very quick example using kin.simple() and kin.search().
#'
#' #retrieve image with arguments passed to find.roi()
#' i <- EBImage::readImage(system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter"))
#' 
#' #create directory and write image to it
#' t <- tempdir()
#' dir.create(paste0(t,"/images"))
#' 
#' EBImage::writeImage(i,paste0(t,"/images/sunfish001.jpg"),type = "jpeg")
#'
#' fi <- list.files(paste0(t,"/images"),full.names=TRUE)
#' #run kin.search and save output image to directory
#' 
#' kin.srch<- kin.search(image.dir = paste0(t,"/images"),
#' save = TRUE,out.dir = t,search.for="largest",size.min=0.01)
#'
#' kin.simp<- kin.simple(image.dir = paste0(t,"/images"),
#' save = TRUE,out.dir = t)
#' 
#' #plot similar results
#' library(ggplot2)
#' 
#' kin.both <- rbind(data.table(kin.srch$midline,fun="search"),
#' data.table(kin.simp$midline,fun="simple"))
#' 
#' qplot(data=kin.both,x=x,y=y.sm,col=fun)
#' 
#' #' #plot midline over original image from kin.simple()
#' i2 <- EBImage::readImage(paste0(t,"/sunfish001.jpg"))
#' EBImage::display(i2,method="raster")
#'
#'
#' ##A somewhat long example using kin.free()
#' #### plot midline waveform on images of swimming ropefish
#' \dontrun{
#' library(data.table)
#' #download example video and establish directories
#' f <- "https://github.com/ckenaley/exampledata/blob/master/ropefish.avi?raw=true"
#' download.file(f, paste0(tempdir(),"/ropefish.avi"))
#'
#' dir.create(paste0(tempdir(),"/images"))
#' dir.create(paste0(tempdir(),"/out"))
#' #extract images
#' vid.to.images(paste0(tempdir(),"/ropefish.avi"), out.dir = paste0(tempdir(),"/images"))
#' 
#' #run kin.free() 
#' 
#' kin <- kin.free(image.dir =paste0(tempdir(),"/images"),
#'       par=TRUE,
#'       save=TRUE,
#'       out.dir=paste0(tempdir(),"/out"),
#'       ml.smooth=list("spline",0.5),
#'       thr = "otsu",
#'       ml.meth="ang",
#'       ant.pos="l",
#'       red=0.5,
#'       size.min=0.01
#'       )
#'
#'
#' #see results
#' #images
#' fi <- list.files(paste0(tempdir(),"/out"),full.names=TRUE)
#' EBImage::display(EBImage::readImage(fi[320]),"raster")
#' 
#' #with gg.overlay() on first 300 frames
#' 
#' gg.overlay(kin=kin,
#' under="cont.sm",
#' over="midline",
#' frames=0:299,
#' size=.2,
#' animate=TRUE,
#' zoom=FALSE,
#' alpha=0.01,
#' col="red",
#' fps=10)
#' 
#' 
#' 
#'}
#'
kin.search <-function(image.dir = NULL,frames = NULL, ant.per = 0.10, tips = 0.02, smooth.n=0,ml.meth="hull",ml.smooth = list("loess",0.25),save = FALSE,plot.pml = TRUE,out.qual = 1,out.dir = NULL, ...) {
    
    size <-
      x <-
      y.pred <-
      wave.y <-
      mid.pred <- roi <- y <-n <- y.sm <- bl <- per.bl <- amp <- head.pval <- NULL # to avoid NSE errors on R CMD check
    

    if (!file.exists(image.dir))
      stop("Directory specified by 'image.dir' (",
           paste0(image.dir),
           ") does not exist")
    
    if (!data.table::between(out.qual, 0, 1))
      stop("'out.qual' must be >=0 and <=1")


    if (save) {
      if (is.null(out.dir))
        stop("'out.dir' not specified")
      if (!file.exists(out.dir))
        stop("Directory specified by 'out.dir' (",
             paste0(out.dir),
             ") does not exist")
      
    }

    
    images <- list.files(image.dir, full.names = TRUE)
    
    if (!length(images) > 0)
      stop("no images in image.dir")
    
    if (any(frames > length(images)))
      stop("variable 'frames' out of range of image sequence")
    if (!is.null(frames))
      images <- images[frames]
    
    if(!is.null(ml.smooth)){
      if(!any(sapply(ml.smooth,class)=="character") | !any(sapply(ml.smooth,class)=="numeric") | !is.list(ml.smooth) ) stop("'ml.smooth' must be a list of length 2 consisting of a character and numeric value")}
    
    smoothing <- unlist(ml.smooth[which(sapply(ml.smooth,class)=="character")])
    smooth <- unlist(ml.smooth[which(sapply(ml.smooth,class)=="numeric")])
    
    trial <- gsub("\\.[^.]*$", "", basename(images[1]))
    
    kin.l <- list()
    midline.l <- list()
    classes.l <- list()
    lms <- list()
    conts <- list()
    conts.sm <- list()
  
    
    roi.outs <- list() #store the rois for each image
    
    for (im in images) {
      frame <- which(im == images) - 1
      
      #find roi
      roi <- find.roi(img=im,...)[[1]]
      #roi <- find.roi(img=im,thr=0.1,search.for="largest",size.min=0.01)[[1]]
      
      best.cont <- best.cont.orig<- roi$best
      best.n <- roi$best.class
      best.class <- roi$classes[roi==best.n]
      
      if(!is.null(smooth.n)){
        if(smooth.n>1|smooth.n<0) stop("'smooth.n' must be 0-1")
        if(smooth.n>0) best.cont <- data.table(Momocs::coo_smooth(as.matrix(best.cont),smooth.n))
      }
      
      colnames(best.cont) <- c("x","y")
     
      
      conts[[paste0(frame)]] <- data.table(frame = frame, best.cont.orig[,n:=1:.N])
      conts.sm[[paste0(frame)]] <- data.table(frame = frame, best.cont[,n:=1:.N])
      classes.l[[paste0(frame)]] <- data.table(frame = frame,roi$classes)
      
      y.df <-
        best.cont[, list(y.min = min(y),
                         y.max = max(y),
                         y = mean(y)), by = list(x)]
      setkey(y.df, "x")
      
      ends <- ceiling(nrow(y.df) * tips)
      tip.y <-
        mean(tail(y.df$y[!is.na(y.df$y)], ends))#tip is mean  of last 30 pixels
      tip.x <-
        mean(tail(y.df$x[!is.na(y.df$y)], ends))#tip is mean y of last 30 pixels
      
      head.y <-
        mean(head(y.df$y[!is.na(y.df$y)], ends))#tip is mean y of first 30 pixels
      head.x <-
        mean(head(y.df$x[!is.na(y.df$y)], ends))#tip is mean y of first 30 pixels
  
      midline <-y.df #two hundred points on midline
      
      midline <- midline[complete.cases(midline)]
      midline <- data.table(frame, midline)
      
      
      ####which type of lines to be fitted, spline or loess
      if (!any(c("spline", "loess") == smoothing))
        stop("'ml.smooth' must contain 'loess' or 'spline'")
      
      if (smoothing == "loess")
        ml.pred <-
        fitted(loess(
          midline$y ~ midline$x,
          span = smooth,
          degree = 1
        ))
      if (smoothing == "spline")
        ml.pred <-
        smooth.spline(x = midline$x,
                      y = midline$y,
                      spar = smooth)$y
      
      midline[, y.sm := ml.pred]#add smoothed predictions
 
      #head section
      head.dat <- midline[1:(ant.per * nrow(midline)), ]
      head.lm <- lm(y.sm ~ x, head.dat)
      
      head.p <- summary(head.lm)$r.squared #how well does head lm fit
      
      midline$mid.pred <-
        predict(head.lm, newdata = midline)#add lm prediction to midline df
      
      
      midline <- midline[complete.cases(midline), ]
      
      midline[, wave.y := y.sm - mid.pred] #wave y based on midline y and straight head.lm pred points
      
      midline[, roi := best.n]
  
      #add bl
      midline[, bl := dist.2d(x, dplyr::lag(x), y.sm, dplyr::lag(y.sm))]
      midline[, per.bl := c(0, cumsum(bl[!is.na(bl)]) / sum(bl, na.rm = TRUE))][, bl := NULL]
      midline[,n:=1:.N]
      
      kin.l[[paste(frame)]] <-
        data.table(
          frame,
          x = tip.x,
          y = tip.y,
          head.x,
          head.y,
          amp = last(midline$wave.y),
          head.pval = head.p,
          best.class
        )
      midline.l[[paste(frame)]] <- midline
      
      
      if (save) {
        jpeg(filename = file.path(out.dir,basename(im)),
          quality = out.qual * 100,
          width = roi$dim[1],
          height = roi$dim[2]
        )
          img <- EBImage::readImage(im)
          suppressMessages(EBImage::display(img, method = "raster"))
        
        
        if (plot.pml)
          lines(predict(lm(mid.pred ~ x, midline)),
                x = midline$x,
                col = "blue",
                lwd = 4)
        with(midline, lines(y.sm ~ x, col = "red", lwd = 4))
        if (plot.pml)
          with(midline[1:ceiling(ant.per * nrow(midline)), ],
               points(
                 x,
                 mid.pred,
                 col = "green",
                 pch = 16,
                 cex = 0.75
               ))
        
        dev.off()
      }
      
    }
    
    classes.dat <- do.call(rbind, classes.l)
    kin.dat <- do.call(rbind, kin.l)
    midline.dat <- data.table(do.call(rbind, midline.l))
    cont.dat <- do.call(rbind, conts)
    cont.sm.dat <- do.call(rbind, conts.sm)
    
    return(
      list(
        kin.dat = kin.dat[,list(frame,roi,x,y,head.x,head.y,amp,head.pval)],
        midline = midline.dat[,list(frame,roi,n,x,y,y.sm,wave.y,per.bl)],
        cont = cont.dat[,list(frame,x,y)],
        cont.sm=cont.sm.dat[,list(frame,n,x,y)],
        mid.pred = midline.dat[,list(frame,x,mid.pred)],
        all.classes = classes.dat,
        dim = roi$dim
      )
    )
  }

NULL

#' @rdname kin.search
#' @export

kin.simple <-function(image.dir = NULL,frames = NULL,ant.per = 0.20,tips = 0.02,smooth.n=0,ml.meth = "hull", ml.smooth = list("loess",0.25), save = FALSE,plot.pml = TRUE,out.qual = 1, out.dir = NULL) {
   
      r <- kin.search(image.dir=image.dir,frames=frames,ant.per=ant.per,tips=tips,ml.meth = ml.meth,ml.smooth=ml.smooth,save=save,out.qual=1,out.dir=out.dir,plot.pml=plot.pml)
  
      return(r)

  }

NULL

#' @rdname kin.search
#' @export
#' 
kin.free <-
  function(image.dir = NULL,frames=NULL,par=FALSE,cores.n=NULL,ant.per = 0.10,ant.pos=NULL,tips = 0.02,smooth.n=1,red=NULL,ml.meth="hull",ml.smooth = list("spline",0.25),save = FALSE,out.qual = 1,out.dir = NULL,plot.pml = TRUE,...) {
    
    #to prevent NSE warnings and NSB notes
    size <-x <-y.pred <-wave.y <-mid.pred <-roi <- prev.dist <-prev.distF <-nxt.dist <-nxt.distF <- x.c <-y.c <- x.sm <-min.x <-range.x <-max.x <-m <-b <- y.sm <-max.n <- n <- keep <-  l <-  head.pval <-  bl <-  per.bl <-  amp <- y <- im <- or.dist <- head.dist <- closer <- n.shift <- x<- dist.next <- dist.current <- close.n <- edge <-flip <- NULL
  
    if (!is.null(frames) &
        length(frames) <= 1)
      stop("number of frames must be >1")
    
    if (!file.exists(image.dir))
      stop("Directory specified by 'image.dir' (",
           paste0(image.dir),
           ") does not exist")
    
    if (!data.table::between(out.qual, 0, 1))
      stop("'out.qual' must be >=0 and <=1")
    
    if (!save &
        !is.null(out.dir))
      stop("'out.dir' specified but 'save=FALSE'. To save processed images, 'save' must be 'TRUE'")
    
    if (save) {
      if (is.null(out.dir))
        stop("'out.dir' not specified")
      if (!file.exists(out.dir))
        stop("Directory specified by 'out.dir' (",
             paste0(out.dir),
             ") does not exist")
      
    }
    
    proc.dir <- out.dir
    
    images <- list.files(image.dir, full.names = TRUE)
    
    if (!length(images) > 0)
      stop("no images in image.dir")
    
    if (any(frames > length(images)))
      stop("variable 'frames' out of range of image sequence")
    if (!is.null(frames))
      images <- images[frames]
    
    
    if (!ml.meth %in% c("ang", "hull", "del"))
      stop("'ml.meth' must be set to 'ang', 'hull', or 'del')")
    
    
    if(!is.null(ml.smooth)){
      if(!any(sapply(ml.smooth,class)=="character") | !any(sapply(ml.smooth,class)=="numeric") | !is.list(ml.smooth) ) stop("'ml.smooth' must be a list of length 2 consisting of a character and numeric value")}
    
    smoothing <- unlist(ml.smooth[which(sapply(ml.smooth,class)=="character")])
    smooth <- unlist(ml.smooth[which(sapply(ml.smooth,class)=="numeric")])
    
    
    trial <- gsub("\\.[^.]*$", "", basename(images[1]))
    
    
    `%fun%` <- `%do%`
    if (par == TRUE) {
      if (!is.null(cores.n))
        cores <- cores.n
      else
        cores <- detectCores()
      cl <- makeCluster(cores[1] - 1) #not to overload your computer
      registerDoParallel(cl)
      `%fun%` <- `%dopar%`
    }
    
    #im=images[1]
    
    kin.res <-
      foreach(
        im = images,
        .combine = append,
        .packages = c("trackter","data.table")
      ) %fun% {
        frame <- which(im == images) - 1
        
        #find roi
        roi <- find.roi(img=im,...)[[1]]
        #roi <- find.roi(img=im,results=FALSE)
        
        best.cont <- roi$best
        best.n <- roi$best.class
        best.class <- roi$classes[roi==best.n]
        
        cont.im <- data.table(frame = frame, best.cont)
        
        ##ml methods
       if(ml.meth=="hull") fml <- free.ml.hull(out = as.matrix(best.cont), smooth.n = smooth.n,red=red)
        if(ml.meth=="ang") fml <- free.ml.ang(out = as.matrix(best.cont), smooth.n = smooth.n,red=red)
      if(ml.meth=="del") fml <- free.ml.del(out = as.matrix(best.cont), smooth.n = smooth.n,red=red)
        
        cont.sm <- data.table(frame = frame, fml$cont.sm)
        
        ml<- fml$ml
        
        #EBImage::display(z,"raster")
        #with(ml,points(x,y,cex=0.1,col="red"))
        
        ends <- ceiling(nrow(ml) * tips)
        tip.y <-
          mean(tail(ml$y[!is.na(ml$y)], ends))
        tip.x <-
          mean(tail(ml$x[!is.na(ml$y)], ends))
        
        head.y <-
          mean(head(ml$y[!is.na(ml$y)], ends))
        head.x <-
          mean(head(ml$x[!is.na(ml$y)], ends))
        
        ml <- data.table(frame, ml)
        
        #with(data.table(head.x,head.y),points(head.x,head.y,col="red"))
        #with(data.table(tip.x,tip.y),points(tip.x,tip.y,col="red"))
        
        ml[, roi := roi$classes$roi]
        
        kin.im <-
          data.table(
            frame,
            x = tip.x,
            y = tip.y,
            head.x,
            head.y,
            amp = NA,
            head.pval = NA,
            best.class
          )
        
        par.res <- list(
          kin = kin.im,
          mid = ml,
          cont = cont.im,
          cont.sm = cont.sm,
          class = data.frame(frame=frame,roi$classes),
          dim = roi$dim
        )
        
        par.res <- list(par.res)
        #names(par.res=im)
        
        return(par.res)
      }
    
    if (par)
      stopCluster(cl)
    
    kin.dat <- kin.dat1 <- do.call(rbind, lapply(kin.res, function(x)
      x$kin))
    midline.dat <- do.call(rbind, lapply(kin.res, function(x)
      x$mid))
    cont.dat <- do.call(rbind, lapply(kin.res, function(x)
      x$cont))
    #cont.dat[,n=1:.N,by=frame] #add index
    cont.sm.dat <- do.call(rbind, lapply(kin.res, function(x)
      x$cont.sm))
    class.dat <- do.call(rbind, lapply(kin.res, function(x)
      x$class))
    dim <-  do.call(rbind, lapply(kin.res, function(x)
      x$dim))[1, ]
    
    #so begins making head and tails
    #find ends closest in next frame
    
    ends <- copy(kin.dat[, list(frame, x, y, head.x, head.y)])
    
    
    ends[, c("prev.dist", "prev.distF") := list(
      dist.2d(head.x,  dplyr::lag(head.x), head.y, dplyr::lag(head.y)),
      dist.2d(head.x,  dplyr::lag(x), head.y,  dplyr::lag(y))
    )][, flip := ifelse(prev.dist > prev.distF, "flip", "unflip")]
    
    ends[1:2, c("nxt.dist", "nxt.distF") := list(
      dist.2d(head.x,  dplyr::lead(head.x), head.y, dplyr::lead(head.y)),
      dist.2d(head.x,  dplyr::lead(x), head.y,  dplyr::lead(y))
    )][1, flip := ifelse(nxt.dist > nxt.distF, "flip", "unflip")]
    
    
    ends[1, c("prev.dist", "prev.distF") := list(nxt.dist, nxt.distF)]
    
    
    ends2 <- copy(ends)
    
    for (i in 1:nrow(ends)) {
      if (with(ends2[i,], prev.dist > prev.distF)) {
        ends2[i, c("x", "y", "head.x", "head.y") := list(head.x, head.y, x, y)]
        
      }
      
      ends2[, c("prev.dist", "prev.distF") := list(
        dist.2d(head.x,  dplyr::lag(head.x), head.y, dplyr::lag(head.y)),
        dist.2d(head.x,  dplyr::lag(x), head.y,  dplyr::lag(y))
      )][, flip := ifelse(prev.dist > prev.distF, "flip", "unflip")]
      
      if (i == 1 &
          with(ends2[1,], nxt.dist > nxt.distF))
        ends2[1, flip := "flip"]
      else
        ends2[1, flip := "unflip"]
      
    }
    
    
    flip.f <- which(ends2[, list(x)] != ends[, list(x)])
    
    midline.dat2 <- copy(midline.dat)[, n := 1:.N, by = list(frame)]
    if (length(flip.f) != 0){
      midline.dat2[frame %in% (flip.f - 1), n := max(n):1, by = list(frame)]
    }
    
    kin.dat <- ends2[,list(frame,x,y,head.x,head.y)]
    #now flip if head is the wrong direction
    
    setkeyv(midline.dat2, c("frame", "n"))
    
    #qplot(d=cont.dat[frame==0],x,y)+geom_point(d=midline.dat2[frame==38],aes(x,y,col=n))
    #centroid in next frame
    
    keep.n <- midline.dat2[,list(max.n=round(max(n)*ant.per)),by=frame]
    midline.dat2 <- midline.dat2[keep.n, on="frame"] 
    
    head.ml <- copy(midline.dat2[n <= max.n])
    head.ml[, c("x.c","y.c"):=list(mean(x), y.c = mean(y)), by = list(frame)]
    
    ends3 <- copy(head.ml)[,keep:=n==min(n)|n==max(n),by=frame][keep==TRUE]
    
    distMax.cent <- ends3[n==max.n][,c("dist.next","dist.current"):=list(dist.2d(x,dplyr::lead(x.c),y,dplyr::lead(y.c)),dist.2d(x,x.c,y,y.c))][frame==0][,diff:=dist.next-dist.current]
    dist0.cent <-ends3[n!=max.n][,c("dist.next","dist.current"):=list(dist.2d(x,dplyr::lead(x.c),y,dplyr::lead(y.c)),dist.2d(x,x.c,y,y.c))][frame==0][,diff:=dist.next-dist.current]
    
    ## if head is moving toward max(n) position, it's tailward
    direction <-ifelse(distMax.cent$diff<dist0.cent$diff, "tailward", "headward")
    
    #add aligned x,y to kin.dat
    kin.dat[,c("x","y","head.x","head.y"):=ends2[,list(x,y,head.x,head.y)]]
    
    if (direction == "tailward"){
      midline.dat2[, n := max(n):1, by = list(frame)]
      kin.dat[,c("x","y","head.x","head.y"):=ends2[,list(head.x,head.y,x,y)]]
    }
    # qplot(data=midline.dat2[frame>30&frame<40],x=x,y=y,col=n)+facet_wrap(frame~.)
    
   
    #override head position
    if(!is.null(ant.pos)){
    
      if(!ant.pos%in% c("l","r","u","d")) stop("'ant.pos' must be one of 'l','r','u', or 'd'")
      
      ends4 <- copy(midline.dat2)[,keep:=n%in%c(max(n),min(n)),by=frame][keep==TRUE]
      
      #head is "l"
      if(ant.pos=="l"){
        head.check <- ends4[,list(head=n[which.min(x)]==1),by=frame]
        midline.dat2 <- midline.dat2[frame%in%head.check[head==FALSE]$frame, n := rev(n), by = list(frame)]
      }
      
        #head is "r"
        if(ant.pos=="r"){
          head.check <- ends4[,list(head=n[which.max(x)]==1),by=frame]
          midline.dat2 <- midline.dat2[frame%in%head.check[head==FALSE]$frame, n := rev(n), by = list(frame)] 
        }
          
          #head is "u"
          if(ant.pos=="u"){
            head.check <- ends4[,list(head=n[which.min(y)]==1),by=frame]
            midline.dat2 <- midline.dat2[frame%in%head.check[head==FALSE]$frame, n := rev(n), by = list(frame)] 
          }
            #head is "d"
            if(ant.pos=="d"){
              head.check <- ends4[,list(head=n[which.max(y)]==1),by=frame]
              midline.dat2 <- midline.dat2[frame%in%head.check[head==FALSE]$frame, n := rev(n), by = list(frame)] 
            
      }
    }
    #rerun kin, head, and ml calculations with new data
    
    ####which type of lines to be fitted, spline or loess
    if (!any(c("spline", "loess") == smoothing))
      stop("'ml.smooth' must contain 'loess' or 'spline'")
    
    if (smoothing == "loess")
      midline.dat2[, y.pred :=
                     fitted(loess(y ~ x,
                                  span = smooth,
                                  degree = 0)), by = list(frame)]
    
    
    if (smoothing == "spline") {
      if (smooth >= 1)
        stop("'smooth' must <1")
      
      sm.spline <- function(x, sm = smooth) {
        s <-  round(seq(1, nrow(x), length.out = nrow(x) * (1 - smooth)), 0)
        x <- as.matrix(x)
        d <- data.table(smoothr::smooth_spline(x[s, ], n = nrow(x)))
        colnames(d) <- c("x.sm", "y.sm")
        return(d)
      }
      
      
      
      midline.dat2[, c("x.sm", "y.sm") := sm.spline(x = data.frame(x = x, y =
                                                                     y), sm = smooth), by = frame]
    }
    
    if (is.null(ml.smooth))
      midline.dat2[, c("x.sm", "y.sm") := list(x, y)]
    
    
    setkeyv(midline.dat2, c("frame", "n"))
    
    head.dat <-
      copy(midline.dat2[n <=max.n, ])
    
    head.l <- list()
    for (h in unique(head.dat$frame))
      head.l[[paste(h)]] <- lm(y.sm ~ x.sm, head.dat[frame == h])
    
    head.p <- lapply(head.l, function(x)
      summary(x)$r.squared)
    
    head.x.fit <- lapply(head.l, fitted)
    
    
    head.dir   <-
      head.dat[, list(dir = ifelse(last(x.sm) < first(x.sm), TRUE, FALSE)), by =
                 frame]
    head.x.m <- lapply(head.l, function(x)
      coef(x)[2])
    head.x.b <- lapply(head.l, function(x)
      coef(x)[1])
    
    
    fit.ab <- function(m, x, b) {
      m * x + b
    }
    
    setkeyv(midline.dat2, c("frame", "n"))
    
    x.range <-
      midline.dat2[n <=max.n, ][, list(x.sm, frame)][, list(
        min.x = dplyr::first(x.sm),
        max.x = dplyr::last(x.sm),
        range.x = diff(range(x.sm))
      ), by = frame][, c("m", "b", "dir") := list(unlist(head.x.m), unlist(head.x.b), head.dir$dir)][, c("min.x", "max.x") :=
                                                                                                       list(ifelse(dir, min.x - range.x * 5, min.x),
                                                                                                            ifelse(dir, min.x, max.x * 5))]
    
    in.circ <- function(x, y, pt.x, pt.y, r) {
      (x - pt.x) ^ 2 + (y - pt.y) ^ 2 <= r ^ 2
    }
    
    
    mid.pred2 <-
      x.range[, list(x = seq(min.x, max.x, 0.5),
                     mid.pred = fit.ab(m, seq(min.x, max.x, 0.5), b)), by = frame]
    
    bls <-
      midline.dat2[, list(l = dist.2d(x.sm[1], dplyr::last(x.sm), y.sm[1], dplyr::last(y.sm))), by =
                     frame]
    
    mid.pred2 <- mid.pred2[bls, on = "frame"]
    mid.pred2 <-
      mid.pred2[head.dat[n == 1, list(frame, x.sm, y.sm)], on = "frame"]
    mid.pred2[, keep := in.circ(x, mid.pred, x.sm, y.sm, l)]
    mid.pred2 <- mid.pred2[keep == TRUE, ]
    
    
    kin.dat[, head.pval := head.p]
    
    
    midline.dat2[, wave.y := dist.2d.line(x.sm, y.sm, unlist(head.x.m[paste(frame)]), unlist(head.x.b[paste(frame)])), by =
                   list(frame, n)]
    midline.dat2[, bl := dist.2d(x.sm, dplyr::lag(x.sm), y.sm, dplyr::lag(y.sm)), by =
                   frame]
    midline.dat2[, per.bl := c(0, cumsum(bl[!is.na(bl)]) / sum(bl, na.rm = TRUE)), by =
                   frame][, bl := NULL]
    
    
   kin.dat[, amp := midline.dat2[, list(wave.y = last(wave.y)), by = frame]$wave.y]
    
    #add roi data back in
    
    kin.dat <- kin.dat[kin.dat1[,list(frame,roi)],on="frame"]
    #qplot(d=cont.sm.dat[frame==9],x,y,col=n)
    
    #re index cont.sm data
    
    #kin.dat <- kin$kin.dat
    #cont.sm.dat <- kin$cont.sm
    cont.sm.dat[,n:=1:.N,by=frame]
    cont.sm.dat2 <-copy(cont.sm.dat[kin.dat[,list(head.x,head.y,frame)],on="frame"])
    
    cont.flip <- cont.sm.dat2[,head.dist:=dist.2d(x,head.x,y,head.y)][,list(close.n=n[which.min(head.dist)]),by=frame]
    
    
    cont.sm.dat2 <-  cont.sm.dat2[cont.flip, on="frame"]
    #reset n=1 to head
    
    cont.sm.dat2 <- cont.sm.dat2[,{c <- Momocs::coo_slide(as.matrix(data.frame(x=x,y=y)),close.n[1]);list(x=c[,1],y=c[,2])},by=frame]
    
    cont.sm.dat2[,n:=1:.N,by=frame]
    
    setkeyv(cont.sm.dat,c("frame","n"))
    
    ##save
    
    if (save) {
      for (xx in images) {
        frame2 <- which(xx == images)
        img2 <-  EBImage::readImage(xx)
        jpeg(file.path(proc.dir,basename(xx)),
             quality = out.qual * 100)
        suppressMessages(EBImage::display(img2, method = "raster"))
        
        
        if (plot.pml) {
          with(mid.pred2[frame == (frame2 - 1)], lines(
            x = x,
            mid.pred,
            col = "blue",
            lwd = 4
          ))
          
          with(midline.dat2[frame == (frame2 - 1)], lines(y.sm ~ x.sm, col = "red", lwd =
                                                            4))
          with(midline.dat2[frame == (frame2 - 1)][n<=max.n,],
               points(
                 x,
                 y.sm,
                 col = "green",
                 pch = 16,
                 cex = 0.75
               ))
        }
        dev.off()
      }
    }
    
    return(
      list(
        kin.dat = kin.dat[,list(frame,roi,x,y,head.x,head.y,amp,head.pval)],
        midline = midline.dat2[,list(frame,roi,n,x,y,x.sm,y.sm,wave.y,per.bl)],
        cont = cont.dat[,list(frame,n,x,y)],
        cont.sm = cont.sm.dat2[,list(frame,n,x,y)],
        all.classes = class.dat,
        mid.pred = mid.pred2[, list(frame, x, mid.pred)],
        dim = dim
      )
    )
  }
NULL

#' @title Tracking of fin-like extensions of body contours

#' @description  Estimates the amplitudes of regions along a body contour that are protruding. Useful in computing paired-fin amplitudes from contour data produced from  \link{kin.simple}, \link{kin.search}, or \link{kin.free}. Also computes a smoothed midline based on the body outline with the fin region removed.
#'
#' @param kin, a named list from one of \code{trackter}s \code{kin} functions. 
#' @param out character, specifying the contour data in 'kin'
#' @param frames integer, the frames within 'kin' to analyze. If NULL then all frames are included.
#' @param fin.pos numeric, a vector of length 2 indicating the start and end of the contour region that contains the fins of interest as a proportion of the body length.
#' @param smooth.n numeric, the number of smoothing operations undertaken by \link{coo_smooth} on the contour described by 'x'. See Details.
#' @param ml.smooth numeric (0-1), the smoothing value for the midline. See Details. 
#' @param red numeric, between 0-1 the proportion of contour coordinates to sample. Will speed up fin position and midline estimations. If 'NULL', the full contour in \code{out} will be used See Details.
#' @param ml.meth character, the midline detection method. One of 'ang' for bisection using \code{\link{free.ml.ang}} or 'hull' for bisection using \code{\link{free.ml.hull}}. Delaunay triangulation using \code{\link{free.ml.del}} is not supported. See Details.
#' 
#' @export
#' @importFrom graphics lines
#' @importFrom stats complete.cases fitted lm loess  predict smooth.spline
#' @importFrom utils head setTxtProgressBar tail txtProgressBar
#'
#' @details
#' If \code{red} is specified, the contour \code{out} is sampled with \code{\link{coo_interpolate}} from the \code{Momocs} package. The number of points sampled (n) equals \code{red} times the number of points in \code{out}
#'
#' To establish the contour positions that are within \code{fin.pos}, a midline is estimated using one of two methods specified by \code{ml.meth}, "ang" for \code{\link{free.ml.ang}} or "hull" for  \code{\link{free.ml.hull}}. Midline points are indexed by position along the body length by calculating the cumulative distance between midline coordinates in pixels. This midline distance is then used to estimate the position of the fin about the contour using the parameter \code{fin.pos}.
#' 
#' The positions of the tip of the fin appendages is estimated in two ways. The first is simply the point in each appendage that is farthest from the base of the fin. The base is estimated as a straight line between the contour coordinates that match \code{fin.pos}.  The second is a little more complicated and starts with calculation the distance between each fin contours coordinates and the midpoint of the fin base. \code{\link{features}} from the \code{features} package is then use to calculate an inflection in this distance and the point of this inflection is used to estimate the fin position. Amplitudes of each method are calculated based on the orthogonal euclidean distance from the fin bases.
#'
#'In addition to fin amplitude and contour extraction, this function also produces a composite contour of the body minus the fin area described by \code{fin.pos}. Fin contours are replaced by a simple linear prediction constructed from the coordinates of the first and last values covered by \code{fin.pos}, that is, the fin bases. The result is a straight line between the start and end of each fin. 
#'
#'From this composite body contour, a midline prediction is made based on the mean position of contour coordinates from each side of the contour sharing a the same position along body using the \code{free.ml} functions. 
#'
#'#' For midline smoothing, \code{\link{smooth_spline}} from the \code{smoothr} package is used to interpolate points between a reduced number of vertices using piecewise cubic polynomials. The number of vertices is calculated based on the number of midline coordinates times the value of \code{ml.smooth}. 
#'
#'
# 'smooth.n' is passed to the 'n' parameter of \code{\link{coo_smooth}}, which smooths coordinates using a simple moving average. Users should be careful not to oversmooth. If the input contour has few points (say just a 100 or so extracted from \code{kin} functions run on low resolution images), much detail will be lost. In general, \code{smooth.n} should be <5.
#'
#'
#' @return A list with the following components:
#'
#' \code{cont} a data table consisting of x,y coordinates of the body contour
#'
#' \code{fin} a data table describing the contour of the fins consisting of the following:
#'
#' \itemize{
#' \item 'side': fin side, 'a' or 'b'
#' \item 'n': the position of fin coordinates where min(n) is closest to the head
#' \item x,y coordinates within the range of \code{fin.pos}
#'
#' }
#'
#' \code{fin.pts} a data table describing fin position consisting of the following:
#'
#' \itemize{
#' \item 'side': fin side, 'a' or 'b'
#' \item 'n': The index matching that of the body contour coordinates 
#' \item  x,y coordinates of the fin tips, start, and end within the range of \code{fin.pos}.
#' \item 'n': The index matching that of the body contour coordinates 
#' \item 'pos': description  of the coordinates' positions, 'start', 'end' or 'tip' or 'tip2'.
#' }
#'
#' \code{comp} a data table describing the composite contour of the body minus the fins.
#' \itemize{
#' \item 'n': The index matching that of the body contour coordinates 
#' \item  x,y coordinates of the body except the range of x values within \code{fin.pos}. These values take on a straight line described by the prediction of \code{lm()} based on the start and end of the fin. See Details.
#' }
#'
#' \code{midline} a data table describing the estimated midline, 'x', 'y', the smoothed x and y positions, respectively
#' 
#' \code{bl} the body length in pixels
#' 
#' \code{amp} a data table describing the estimated fin amplitudes based on each method of finding the fin tip. 'amp1', the amplitude of the fin positions based on the maximum euclidean distance from the fin base and 'amp2', the distance of the fin points (based on distance inflections) from the fin base.
#'
#' @seealso \code{\link{kin.simple}}, \code{\link{kin.search}}, \code{\link{kin.free}}, \code{\link{coo_sample}}, \code{\link{coo_smooth}}, \code{\link{smooth_spline}}
#' @export
#'
#' @importFrom graphics lines
#' @importFrom stats complete.cases fitted lm loess  predict smooth.spline
#' @importFrom utils head tail 
#'
#' @examples
#' ###plot pectoral-fin amplitudes of a swimming sunfish
#' \dontrun{
#' require(ggplot2)
#'
#' #download example avi video
#' f <- "https://github.com/ckenaley/exampledata/blob/master/sunfish_pect.avi?raw=true"
#' download.file(f,"sunfish.avi")
#' 
#' #create directories
#'  ti <-paste0(tempdir(),"/images")
#'  dir.create(ti)
#'
#' #extract images with ffmpeg operations and reduce them to 600 px wide with a filter
#' filt.red <- " -vf scale=600:-1 " #filter
#' vid.to.images2(vid.path="sunfish.avi",filt = filt.red,out.dir=ti) #extract
#'
#' #number of frames
#' fr <- list.files(ti,full.names=TRUE)
#' thr.check(fr[3])
#' 
#' #extract contours and other data
#' kin <- kin.search(image.dir = ti,thr=0.9,ant.per = 0.25,save=FALSE)
#' 
#' #fin data by frame
#' fin.pos <- c(0.25,.55)
#' fin.dat <- fin.kin(kin=kin,fin.pos = fin.pos,smooth.n=1,ml.smooth=0.3)
#' 
#' p <- ggplot(dat=fin.dat$amp,aes(x=frame,y=amp2,col=side))+geom_line()+theme_classic(15)
#'print(p)
#'
#'
#' ## plot body and fin contours of a frame
#'
#' #plot body contour and fins
#' fr <- 8
#' p <- qplot(data=fin.dat$cont[frame==fr],x=x,y=y)
#' p <- p+geom_point(data=fin.dat$fin[frame==fr],aes(x,y),col="red",size=3)
#' p <- p+geom_point(data=fin.dat$fin.pts[frame==fr],aes(x,y,shape=pos))
#' p+xlim(c(0,kin$dim[1]))+ylim(c(0,kin$dim[2]))
#' 
#' 
#' #plot body contour minus fins and the body midline
#' p <- qplot(data=fin.dat$comp[frame==fr],x=x,y=y)
#' p <- p++geom_point(data=fin.dat$midline[frame==fr],aes(x,y),col="red",size=2)
#' p+xlim(c(0,kin$dim[1]))+ylim(c(0,kin$dim[2]))
#'
#'# now the whole sequence with gg.overlay()
#'  gg.overlay(dat=fin.dat,
#'  under="cont",
#'  over="fin.pts",
#'  size=5,
#'  animate=TRUE,
#'  fps=5,
#'  alph=0.05)
#'  
#'unlink(ti,recursive=TRUE)
#' }
#' 
#'
fin.kin <-
  function(kin,out="cont",frames=NULL,fin.pos = NULL,smooth.n = 5, ml.meth="hull",ml.smooth = 0.9,red=NULL) {
    
    y <- x <- x.sm <- y.sm <- n <- m <- b <- dist2 <- x.c <- y.c <- tip1 <- tip2 <- method <- amp2 <-pos <- y.pred <-side <- ends <- head.dist <-frame <- gap <- bl <-size <-   NULL # due to NSE notes in R CMD checks
    # 
    if(is.null(frames)) frames <-unique(kin[[1]]$frame)
    if(!is.null(frames) & any(!frames %in% kin[[1]]$frame)) stop("not all frames are in 'kin' list")
    
     if (is.null(fin.pos))               
      stop("'fin.pos' not defined")
    
    if (length(fin.pos) != 2)
      stop("length of 'fin.pos' argument must be 2")
    
    if(!out %in% names(kin)) stop("outline data not in kin list. Did you mean to enter 'cont'")
    
    
    if (!ml.meth %in% c("ang", "hull"))
      stop("'ml.meth' must be set to 'ang' or 'hull")
    
    sm.spline <- function(x, sm = 0.5) {
      s <-  round(seq(1, nrow(x), length.out = nrow(x) * (1 - sm)), 0)
      x <- as.matrix(x)
      d <- data.table(smoothr::smooth_spline(x[s, ], n = nrow(x)))
      colnames(d) <- c("x.sm", "y.sm")
      return(d)
    }
    
    dist.2d.line <- function(x, y, slope, intercept) {
      b = c(1, intercept + slope)
      c = c(-intercept / slope, 0)
      a = c(x, y)
      v1 <- b - c
      v2 <- a - b
      m <- cbind(v1, v2)
      return(abs(det(m)) / sqrt(sum(v1 * v1)))
    }
    

    k <- kin$kin.dat
    c <- kin$cont
    if(out=="cont.sm") c <- kin$cont.sm
    
 #f=0
    r <- lapply(frames, function(f,...){
      o <- c[frame==f,list(x,y)]
      
    if (!is.matrix(o))
     o<- as.matrix(o)
    if (is.null(colnames(o)))
      colnames(o) <- c("x", "y")
    
    
    if(ml.meth=="ang") fml <- free.ml.ang(out=as.matrix(o),smooth.n =smooth.n,red=red)
    if(ml.meth=="hull") fml <- free.ml.hull(out=as.matrix(o),smooth.n =smooth.n,red=red)
    
    #with(fml$cont.sides,plot(x,y))
    
    ## shift indext to head position
    
    #the head
    h <- k[frame==f,]
  
    cont.sm <- fml$cont.sm[,n:=1:.N]
    #qplot(d=cont.sm,x,y,col=n)
    
    cont.flip <- cont.sm[,head.dist:=dist.2d(x,h$head.x,y,h$head.y)][,list(close.n=n[which.min(head.dist)])]
    

    #reset n=1 to head
    cont.sm <- data.table(Momocs::coo_slide(as.matrix(cont.sm[,list(x,y)]),cont.flip$close.n))
    cont.sm[,n:=1:.N]
    #qplot(d=cont.sm,x,y,col=n)
    
    
    setkeyv(cont.sm,c("n"))
    
    ml2 <- sm.spline(x=as.matrix(fml$ml[,list(x,y)]),sm=0.99)[,n:=1:.N]
    
    #qplot(d=ml2,x.sm,y.sm,col=n)
    min.dist <- with(ml2,dist.2d(x.sm,h$head.x,y.sm,h$head.y))
    
    #    #flip cont.sides and ml2 if needed
    if(which.min(min.dist)>which.max(min.dist)){ ml2[,n:=max(n):min(n)]
      fml$cont.sides[,n:=max(n):min(n),by=side]
    }
    
    setkeyv(fml$cont.sides,c("side","n"))
    setkeyv(ml2,c("n"))
    
    #qplot(d=ml2,x.sm,y.sm,col=n)
    
    #qplot(d=fml$cont.sides,x,y,col=n)
    rot <- ml2[,{pa <- point.ang.orig(c(x.sm,y.sm),c(ml2$x.sm[1],ml2$y.sm[1]),pi/2);list(x=pa[1],y=pa[2])},by=n][c(range(n)),]
    rot.lm <- lm(y~x,rot)
    #rot$pred <- predict(rot.lm)
    
   # qplot(d=fml$cont.sides,x,y,col=dist)+geom_point(d=rot,aes(x,),col="red")

    
    fml$cont.sides[,dist:=dist.2d.line(x,y,coef(rot.lm)[2],coef(rot.lm)[1]),by=list(n,side)][,bl:=dist/max(dist),by=side]
    
    fins <- fml$cont.sides[data.table::between(bl,fin.pos[1],fin.pos[2]),]
    
    #fml$cont.sides$bl
    #qplot(data= fins,x=x,y=y,col=side)
    
    lm.coefs <- fins[,ends:=n%in%c(min(n),max(n)),by=side][ends==TRUE,{l <- lm(y~x);list(b=coef(l)[1],m=coef(l)[2])},by=side]
    
    cents <- fins[ends==TRUE,list(x.c=mean(x),y.c=mean(y)),by=side]
    fins <- fins[lm.coefs,on="side"][cents,on="side"][,dist:=dist.2d.line(x,y,m,b),by=list(side,n)][,dist2:=dist.2d(x,x.c,y,y.c),by=side]
    
    
    fins[,tip1:=dist==dist[which.max(dist)],by=side]

    
    fins[tip1==TRUE,tip1:=dist2==max(dist2)&n==max(n),by=side] #&n==max(n) for duplicated coords
    
    #qplot(d=fins,x,y,col=side)+geom_point(d=fins[tip1==TRUE],aes(x,y),col="black")
    
    cp.a <- with(fins[side=='a'],features(n,dist2,smoother = "smooth.spline",spar=0.3))
    cp.an <- round(cp.a$cpts,0)
    cp.b <- with(fins[side=='b'],features(n,dist2,smoother = "smooth.spline",spar=0.3))
    cp.bn <- round(cp.b$cpts,0)
    
    
    cp.bn <- fins[side=='b'][n%in%cp.bn,][which.max(dist2),]$n
    cp.an <- fins[side=='a'][n%in%cp.an,][which.max(dist2),]$n
    
    if(length(cp.an)==0) cp.an <- fins[side=='a'][which.min(dist2)]$n
    if(length(cp.bn)==0) cp.an <- fins[side=='b'][which.min(dist2)]$n
    cp <- data.table(side=c("a","b"),n=c(cp.an,cp.bn),tip2=TRUE)
    
    
  
    #qplot(d=fml$cont.sides,x,y,col=n)+facet_grid(side~.)
    
    fins <- cp[fins,on=c("side","n")]

    
    lms <- list()
    for(i in c("a","b")){
      lms[[i]] <- lm(y~x,fins[ends==TRUE & side==i])
    }
    
    fins[, y.pred := predict(lms[[side]], newdata = data.frame(x = x)), by = list(side)]
    
  
    #number of coordinates not in fins

    comp <- merge(fml$cont.sides, fins[,list(y.pred,n,side)], by = c("n","side"), all.x = TRUE)
    comp[,y:=as.numeric(y)]
    comp[!is.na(y.pred), y := y.pred]
    
    comp[,dist:=dist.2d(x,cont.sm[n==1]$x,y,cont.sm[n==1]$y),by=side]
    setkeyv(comp,c("side","dist"))
  
    #qplot(data= cont.sm,x=x,y=y,col=n)
    #qplot(data= comp,x=x,y=y,col=n)


    comp[,n:=1:.N,by=side]
    
    setkeyv(comp, c("n"))
    
    
   #with(comp,plot(x,y))
    
    fins2 <- copy(fins)
    fins2[, `:=`(pos, "pt"),by=side]
    
    fins2[, `:=`(pos, ifelse(n==min(n),"start",pos)),by=side]
    fins2[, `:=`(pos, ifelse(n==max(n),"end",pos)),by=side]
    fins2[tip1==TRUE,pos:="tip1"]
    fins2[tip2==TRUE,pos:="tip2"]
    #add tips are both 1 and 2
    dup.tip <- fins2[tip1==TRUE & tip2==TRUE & pos=="tip2"][,pos:="tip1"]
    fins2 <- rbind(fins2,dup.tip)
    
    
    finPts <- fins2[pos!="pt",list(side,n,x,y,pos)]
    finPts <- finPts[!duplicated(finPts),]
    setkeyv(finPts,c("side","pos"))
    
    # qplot(data= fml$cont.sides,x=x,y=y)+geom_point(data=fins[tip2==TRUE],aes(x,y),size=3,col="red")
    
    amp <- fins2[!is.na(pos)&tip1|tip2,][,method:=pos][,list(amp1=dist[tip1],amp2=dist[tip2]),by=side][!is.na(amp2)]
  
    #with(comp,plot(x,y))
    
 
    ### get points that are inside the composite 
    
    cont.mat <- as.matrix(cont.sm[,list(x,y)])
    comp.mat <- as.matrix(comp[,list(x,y)])
    sr <- SpatialPolygons(list(Polygons(list(Polygon(comp.mat)),ID=1)))
    ins <-  over(SpatialPoints(cont.mat),sr)
    
    ### select the points not on contour and redefine composite
    cont.pts <-  data.table(cont.mat[!is.na(ins),])
    comp <-cont.sm[cont.pts, on=c("x","y")]
    setkeyv(comp,"n")
    
  
    
    #get rid of points that sneak in
    sr2 <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(comp[,list(x,y)]))),ID=1)))
    ins2 <-  over(SpatialPoints(as.matrix(comp[,list(x,y)])),SpatialPoints(comp.mat))
    
    comp2 <-  comp[!is.na(ins2),]
    comp2[,gap:=dplyr::lead(n)-n>5|abs(dplyr::lag(n)-n)>5]
   #qplot(d=comp2,x,y)
    
    # plot(cont.pts)
    # with(comp,points(x,y,col="red"))
    # points(comp.mat,col="red")
    
    if(nrow(comp2[gap==TRUE])>4) stop(paste0("'fin.pos' appears not to capture all of the fin in frame ",f," . Try expanding the upper limit."))
  
    comp2[gap==TRUE,gaps:=c("a","a","b","b")]
    
    #qplot(d=comp2,x,y,col=gap)
    #predict gap points
    
    gaps <- comp2[gap==TRUE, list(
      n=c(n[1]:n[2])[-1],
      x = seq(x[1],x[2],length.out = n[2]-n[1]),
      y = predict(lm(y~n), newdata = data.frame(n = seq(n[1],n[2],length.out = n[2]-n[1])))), by = list(gaps)]
    
    #redefine comp again, include predicted gap points
    comp <- rbind(comp[,list(x,y,n)],gaps[,list(x,y,n)])
setkeyv(comp,"n")

    #qplot(data= comp,x=x,y=y,col=n)

    if(ml.meth=="ang") fml.comp <- free.ml.ang(as.matrix(comp[,list(x,y)]))
    if(ml.meth=="hull") fml.comp <- free.ml.hull(as.matrix(comp[,list(x,y)]))
    
  ml.comp <- fml.comp$ml 

    #qplot(data= ml.comp,x=x,y=y,col=n)+geom_point(d=fml.comp$cont.sides,aes(x,y))
  
    
    ml.comp.s <- ml.comp[,{s <- predict(smooth.spline(data.frame(x,y),spar = ml.smooth));list(x=s$x,y=s$y)}]
    
   # qplot(data= comp,x=x,y=y,col=n)+geom_point(dat=ml.comp.s,aes(x,y,col=n))
    
  bl <- sum(copy(ml.comp.s[,dist:=dist.2d(x,dplyr::lead(x),y,dplyr::lead(y))])$dist,na.rm = TRUE)
     
     # p <- qplot(data= fml$cont.sides,x=x,y=y)+geom_point(dat=fins,aes(x,y,col=side))+geom_point(data=fins[tip2==TRUE],aes(x,y),size=3,col="black")+geom_point(d=cents,aes(x.c,y.c),col="blue")+geom_point(d= ml.comp,aes(x,y))
  
  #qplot(data= fml$cont.sides,x=x,y=y)+geom_point(d=finPts,aes(x,y,col=pos))
  
     # print(p)
    # 
    return(list(
      cont = data.table(frame=f,cont.sm), #smoothed contour
      fin = data.table(frame=f,fins[,list(side,n,x,y)]),
      fin.pts = data.table(frame=f,finPts),
      comp = data.table(frame=f,comp[,list(n,side,x,y)]),
      midline = data.table(frame=f,ml.comp.s[,list(x,y)]),
      bl = data.table(frame=f,bl) ,
      amp = data.table(frame=f,amp)
    ))
    }
    )
    
    
    cont.dat <- do.call(rbind, lapply(r, function(x)
      x$cont))
    fin.dat <- do.call(rbind, lapply(r, function(x)
      x$fin))
    fin.pts <- do.call(rbind, lapply(r, function(x)
      x$fin.pts))
    comp.dat <- do.call(rbind, lapply(r, function(x)
      x$comp))
    mid.dat <- do.call(rbind, lapply(r, function(x)
      x$midline))
    bl.dat <- do.call(rbind, lapply(r, function(x)
      x$bl))
    amp.dat <- do.call(rbind, lapply(r, function(x)
      x$amp))
    
    r2 <- list(cont=cont.dat,
         fin=fin.dat,
         fin.pts=fin.pts,
         comp=comp.dat,
         midline=mid.dat,
         amp=amp.dat,
         bl=bl.dat
    )
    

    #qplot(data= r2$comp[frame==1],x=x,y=y,col=side)+geom_point(dat=ml.comp.s,aes(x,y,col=n))
    
    return(r2)
}
