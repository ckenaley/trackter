#' @title  Midline and outline tracking over image sequences

#' @description  Wrapper functions for \code{\link{find.roi}} that automatically retrieve the contour and midline coordinates of a detected ROI in each image of a sequence through thresholding and segmentation.  Functions find the y-value midpoint along the x-value array of the ROI and fits a midline according to a chosen smoothing method (loess or spline). Also outputs the midline amplitude relative to a reference line determined by an anterior section of the ROI. Supported image formats are jpeg, png, and tiff.
#' 
#' \code{kin.search} includes arguments for flexible ROI selection.
#' 
#' \code{kin.simple} is itself a wrapper for \code{kin.search}, finding the largest ROI in field using Otsu thresholding for segmentation.
#' 
#' @param image.dir character, directory containing images to analyze.
#' @param frames numeric, vector indicating which images to process. NULL, the default, will result in all images in \code{image.dir} processed.
#' @param ant.per numeric; left-most percentage of ROI that establishes the horizontal reference for the midline displacement.
#' @param tips, numeric, the proportion the the midline data to use in calculation of the head and tail position.
#' @param plot.pml logical, value indicating if outputted images should include an overlay of the theoretical midline based on \code{ant.per}.
#' @param smoothing character, the midline smoothing method, either 'loess' or "spline".
#' @param smooth numeric; if \code{smoothing} if set to 'loess', smoothing parameter value for plotted midline.
#' @param save logical, value indicating if images should be outputted with midline and predicted midline based on the \code{ant.per} \code{lm()} overlaying original or binary images.
#' @param out.qual, numeric, a value between 0-1 representing the quality of outputted images. Ignored if \code{save=FALSE}.
#' @param out.dir character, the directory to which outputted images should be saved.
#' @param ... other parameters passed to \code{link{find.roi}} (\code{kin.search} only).
#'
#' @export
#'
#'
#' @details
#'The algorithm assumes a left-right orientation, i.e., the head of the ROI is positioned left, the tail right. If this is not the case, consider using \code{\link{kin.free}} or rotating images before processing. The \code{ant.per} value therefor establishes the reference line (theoretical straight midline) based on that portion of the head. The midline is calculated as the midpoints between the y extrema for each x position. 
#'
#'\code{kin.search} chooses ROIs based on relative ROI size or position according too to \code{\link{find.roi}}. Thresholding operations can be performed with an arbitrary (user defined) numeric value or with Otsu's method ('thr="otsu"'). The latter chooses a threshold value by minimizing the combined intra-class variance. See \code{\link{otsu}}. Other search arguments can be adjusted by passing arguments to \code{link{find.roi}}.
#'
#' \code{kin.search} is more streamline. It attempts to find the largest ROI using Outsu thresholding and invokes other default values of \code{\link{find.roi]}.
#'
#'
#' @return A list with the following components:
#'
#' \code{kin.dat} a data table consisting of frame-by-frame position parameters for the ROI determined by \code{search.for}.
#' \itemize{
#' \item the frame number
#' \item 'x' and ''y': the position of the tail (rightmost or posteriormost)
#' \item 'head.x' and 'head.y': the x and y position of the head (leftmost or anteriormost)
#' \item 'amp': the amplitude (\code{amp}) of the tail relative to thr theoretical midline determined by the \code{lm()} predictions from \code{ant.per}
#' \item 'head.pval': p values of the \code{lm()} fit that describes the position of the head as determined by \code{ant.per} (green points in the outputted images/video)
#' \item 'roi': a character indicating the ROI ranked by size ('a' being the largest)
#' \item 'edge': indicating whether ROI was on the edge of the image field
#' \item 'size': size of the ROI in pixels^2
#' \item 'offset.x': ROI distance from horizontal center
#' \item 'offset.y': ROI distance from vertical center
#' \item 'offset': linear distance of ROI's centroid to image center
#' }
#'
#' \code{midline} A data table containing, for each frame described by \code{frames}, the following:
#' \itemize{
#' \item 'x' and 'y.m': x and y positions of the midline of the ROI
#' #' \item 'y.min' and 'y.max': min and max y positions ROI's contour used in y.m calculation
#' \item 'mid.pred': the predicted linear midline based on the points/pixels defined by \code{head.per} (green points in the outputted images/video)
#' \item 'y.pred': midline points fit to a smooth spline or loess model with spar or span equal to \code{smooth} (red curve in the outputted images/video)
#' \item 'wave.y': midline points 'y.pred' relative to 'mid.pred'
#' \item 'roi': a character indicating ROI size ('a' being the largest)
#' }
#'
#' \code{cont} A data table containing x and y positions of the contours used to calculate the data in 'kin.dat'. Contains the following:
#' \itemize{
#' \item 'frame': the frame
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
#' \code{dim} the x and y dimensions of the images analyzed
#'
#' @return A list with the following components:
#'
#' \code{kin.dat} a data frame consisting of frame-by-frame position parameters for the ROI
#' \itemize{
#' \item the frame number
#'
#' \item 'head.x' and 'head.y': the x and y position of the head (leftmost or anteriormost)
#' \item 'x' and 'y': the position of the tail (rightmost or posteriormost)
#' \item 'amp': the amplitude (\code{amp}) of the tail
#' \item 'cent.x' and 'cent.y': centroid coordinate of ROI
#' \item 'roi': a character indicating ROI size ('a' being the largest)
#' \item 'head.pval': p values of the \code{lm()} fit that describes the position of the head as determined by \code{ant.per} (green points in the outputted images/video)}
#'
#' \code{midline} A data frame containing, for each frame described by \code{frames}, the following: \itemize{
#' \item 'x' and 'y.m': x and y positions of the midline of the ROI
#' \item 'roi': a character indicating ROI size ('a' being the largest)
#' \item 'mid.pred': the predicted linear midline based on the points/pixels defined by \code{head.per} (green points in the outputted images/video)
#' \item 'y.pred': midline points fit to a smooth spline or loess model with spar or span equal to \code{smooth} (red curve in the outputted images/video)
#' \item 'wave.y': midline points 'y.pred' normalized to 'mid.pred'
#' }
#'
#' \code{dim} the x and y dimensions of the images analyzed
#'
#' @export
#'
#' @importFrom graphics lines
#' @importFrom stats complete.cases fitted lm loess  predict smooth.spline
#' @importFrom utils head tail
#'
#'
#' @seealso \code{\link{kin.simple}}, \code{\link{kin.free}}, \code{\link{find.roi}}
#' @examples
#'
#' #### plot caudal amplitude and produce a classic midline waveform plot of a swimming rainbow trout
#' ##A very long example.
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
#'       smoothing = "loess",frames=1:50,
#'       out.dir=paste0(tempdir(),"/processed_images"),
#'       smooth=0.4)
#'       
#' 
#' #plot instantaneous amplitude of tail (last/rightmost point) over frames
#' library(ggplot2)
#' p <- ggplot(dat=kin$kin.dat,aes(x=frame,y=amp))+geom_line()+geom_point()+theme_classic(15)
#' print(p)
#'
#' # midline plot
#' ml <- kin$midline
#' #leftmost x starts at 0
#' ml <- ml[,x2:=x-x[1],by=frame]
#'
#' ml <- merge(ml,kin$kin.dat[,list(frame,amp)],by="frame") #merge these
#'
#' p <- ggplot(dat=ml,aes(x=x2,y=wave.y))+theme_classic(15)
#' p <- p+geom_line(aes(group=frame,color=amp),stat="smooth",method = "loess", size = 1.5)
#' print(p)
#'
#' #Make a video of processed frames
#'
#' images.to.video2(image.dir=paste0(tempdir(),"/processed_images"),
#' vid.name="trout_test",out.dir=tempdir(),frame.rate=5,qual=100,raw=FALSE)
#' file.exists(paste0(tempdir(),"/trout_test_red.mp4"))
#'
#'}
#'
#' ## A very quick example using kin.simple() and kin.search().
#'
#' #retrieve image with arguments passed to find.roi()
#' i <- EBImage::readImage(system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter"))
#' #create directory and write image to it
#' t <- tempdir()
#'
#'
#' dir.create(paste0(t,"/images"))
#' 
#' EBImage::writeImage(i,paste0(t,"/images/sunfish001.jpg"),type = "jpeg")
#'
#' fi <- list.files(paste0(t,"/images"),full.names=TRUE)
#' #run kin.search and save output image to directory
#' 
#' kin.srch<- kin.search(image.dir = paste0(t,"/images"),smooth=0.2,
#' save = TRUE,out.dir = t,thr=0.8,search.for="largest",size.min=0.01)
#'
#' kin.simp<- kin.simple(image.dir = paste0(t,"/images"),smooth=0.2,
#' save = TRUE,out.dir = t)
#' 
#' #plot similar results
#' library(ggplot2)
#' 
#' kin.both <- rbind(data.table(kin.srch$midline,fun="search"),
#' data.table(kin.simp$midline,fun="simple"))
#' 
#' qplot(data=kin.both,x=x,y=y.pred,col=fun)
#' 
#' #' #plot midline over original image from kin.simple()
#' i2 <- EBImage::readImage(paste0(t,"/sunfish001_000.jpg"))
#' EBImage::display(i2,method="raster")
#'
#' #clean up
#' unlink(paste0(t,"/images"),recursive=TRUE)
#'

kin.search <-function(image.dir = NULL,frames = NULL,plot.pml = TRUE, ant.per = 0.10, tips = 0.02, smoothing = "loess", smooth = 0.25,save = TRUE,out.qual = 1,out.dir = NULL, ...) {
    
    size <-
      x <-
      y.pred <-
      wave.y <-
      mid.pred <- roi <- NULL # to avoid NSE errors on R CMD check
    
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
    
    trial <- gsub("\\.[^.]*$", "", basename(images[1]))
    
    kin.l <- list()
    midline.l <- list()
    classes.l <- list()
    lms <- list()
    conts <- list()
   
    
    roi.outs <- list() #store the rois for each image
    
    for (im in images) {
      frame <- which(im == images) - 1
      
      #find roi
      roi <- find.roi(img=im,...)[[1]]
      #roi <- find.roi(img=im,thr=0.1,search.for="largest",size.min=0.01)[[1]]
      
      best.cont <- roi$best
      best.n <- roi$best.class
      best.class <- roi$classes[roi==best.n]
     
      
      conts[[paste0(frame)]] <- data.table(frame = frame, best.cont)
      classes.l[[paste0(frame)]] <- data.table(frame = frame,roi$classes)
      
      y.df <-
        best.cont[, list(y.min = min(y),
                         y.max = max(y),
                         y.m = mean(y)), by = list(x)]
      setkey(y.df, "x")
      
      ends <- ceiling(nrow(y.df) * tips)
      tip.y <-
        mean(tail(y.df$y.m[!is.na(y.df$y.m)], ends))#tip is mean y.m of last 30 pixels
      tip.x <-
        mean(tail(y.df$x[!is.na(y.df$y.m)], ends))#tip is mean y.m of last 30 pixels
      
      head.y <-
        mean(head(y.df$y.m[!is.na(y.df$y.m)], ends))#tip is mean y.m of first 30 pixels
      head.x <-
        mean(head(y.df$x[!is.na(y.df$y.m)], ends))#tip is mean y.m of first 30 pixels
  
      midline <-y.df #two hundred points on midline
      
      midline <- midline[complete.cases(midline)]
      midline <- data.table(frame, midline)
      
      
      ####which type of lines to be fitted, spline or loess
      if (!any(c("spline", "loess") == smoothing))
        stop("'smoothing' must = 'loess' or 'spline'")
      
      if (smoothing == "loess")
        ml.pred <-
        fitted(loess(
          midline$y.m ~ midline$x,
          span = smooth,
          degree = 1
        ))
      if (smoothing == "spline")
        ml.pred <-
        smooth.spline(x = midline$x,
                      y = midline$y.m,
                      spar = smooth)$y
      
      midline[, y.pred := ml.pred]#add smoothed predictions
      
      #head section
      head.dat <- midline[1:(ant.per * nrow(midline)), ]
      head.lm <- lm(y.pred ~ x, head.dat)
      
      head.p <- summary(head.lm)$r.squared #how well does head lm fit
      
      midline$mid.pred <-
        predict(head.lm, newdata = midline)#add lm prediction to midline df
      
      midline <- midline[complete.cases(midline), ]
      
      midline[, wave.y := y.pred - mid.pred] #wave y based on midline y and straight head.lm pred points
      
      midline[, roi := best.n]
  
      
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
        jpeg(
          paste0(proc.dir, "/", trial, "_", sprintf("%03d", frame), ".jpg"),
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
        with(midline, lines(y.pred ~ x, col = "red", lwd = 4))
        if (plot.pml)
          with(midline[1:ceiling(ant.per * nrow(midline)), ],
               points(
                 x,
                 y.pred,
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
    
    return(
      list(
        kin.dat = kin.dat,
        midline = midline.dat,
        cont = cont.dat,
        all.classes = classes.dat,
        dim = roi$dim
      )
    )
  }

#' @rdname kin.search
#' @export

kin.simple <-function(image.dir = NULL,frames = NULL,ant.per = 0.20,tips = 0.02, smoothing = "loess", smooth = 0.25, save = TRUE,out.qual = 1, out.dir = NULL,plot.pml = TRUE) {
   
      r <- kin.search(image.dir=image.dir,frames=frames,ant.per=ant.per,tips=tips,smoothing=smoothing,save=save,out.qual=1,out.dir=out.dir,plot.pml=plot.pml)
  
      return(r)

  }

NULL

#' @title  Contour and midline tracking of free-moving ROIs over image sequences

#' @description  A wrapper function for \code{find.roi} that automatically ROIs that are free to move in the spatial field of an image sequence. Does so through thresholding and segmentation. Estimates midlines by various methods. Also outputs the midline amplitude relative to a reference line determined by an anterior section of the ROI. Supported image formats are jpeg, png, and tiff. Supports parallel processing of frames.
#' 
#' @param image.dir character, directory containing images to analyze.
#' @param frames numeric, vector indicating which images to process. Must be >1. See Details.
#' @param par logical, should the frames be processed in parallel using \code{cores.n}.
#' @param cores.n numeric, the number of CPU cores to use if \code{par=TRUE}. If \code{cores.n=NULL} (the default), the total number of cores minus 1 are used.
#' @param thr numeric or character ('otsu'); the threshold to determine binary image. See Details.
#' @param ant.per numeric; anterior percentage of ROI that establishes the reference for the midline displacement.
#' @param ant.pos character, one of , NULL, "l","r","u",or "d" to specify the position of the anterior of the ROI. If not NULL, the default algorithm to find the anterior is overridden. See Details.
#' @param tips, numeric, the proportion the the midline data to use in calculation of the head and tail position.
#' @param smooth.n, numeric, the number of contour smoothing iterations. See Details.
#' @param red numeric, between 0-1 the proportion of contour coordinates to sample for midline estimates. Ignored if \code{ml.meth} is not 'del'. Will speed up midline estimations with Delauny triangulation. If 'NULL', the full contour retrieved from the ROI will be passed to \code{\link{free.ml.del}}. See Detail.
#' @param ml.meth character, the midline detection method. One of 'ang' for bisection using \code{\link{free.ml.ang}}, 'hull' for bisection using \code{\link{free.ml.hull}}, or 'del' for Delaunay triangulation using \code{\link{free.ml.del}}. See Details.
#' @param ml.smooth a list of length two with unnamed components including a character string specifying the midline smoothing method, either 'loess' or "spline", and a numeric value specifying the amount of smoothing. See Details.
#' @param save logical, value indicating if images should be outputted with midline and predicted midline based on the \code{ant.per} \code{lm} overlaying original or binary images.
#' @param out.qual, numeric, a value between 0-1 representing the quality of outputted images. Ignored if \code{save=FALSE}.
#' @param out.dir character, the directory to which outputted images should be saved.
#' @param plot.pml logical, value indicating if outputted images should include an overlay of the midline, head region and theoretical midline based on \code{ant.per}.
#' @param . . . , other arguments passed to \code{\link{find.roi}}
#'
#' @export
#'
#' @details
#' By default, the position of the anterior of the ROI (that which is moving forward in the field) is determined by the displacement of the ROI between the first two frames. Thus, \code{frames} must be >1. For analyses of relatively static ROIs in the field (e.g., steadily swimming animals in flumes, etc.), automatically determining the anterior of the ROI may be spurious. In this case, the automatic determination of the anterior should be overridden by specifying 'l', 'r', 'u', 'd' with \code{ant.pos}. These values specify that the anterior region of the ROI is leftmost, rightmost, upmost, or downmost in the field, respectively, and assumes that the origin of the field (0,0) is the upper left corner of each frame. 
#' 
#'Thresholding operations are preformed with \code{link{find.roi}}. Parameters for this function are passed through addition arguments with ...
#'
#'
#' Midline estimation is pursued by one of three algorithms: bisection of contours across the long axis defined by the tips using \code{\link{free.ml.ang}} or \code{\link{free.ml.hull}}  or by Delaunay triangulation using \code{\link{free.ml.del}}. The default is 'hull' This choice is not arbitrary. The use of \code{free.ml.ang} and \code{free.ml.hull} can be faster, but perfom poorly for tips that snake back on themselves (i.e., a high degree curvature). The use of \code{free.ml.del} can be slowe for high resolution outlines, but  produces better results when contour regions overlap (i.e, those that snake back on themselves), but produces less precise midlines for complicated contours.
#' 
#' For midline smoothing, if \code{ml.smooth} contains 'spline' (the default), \code{\link{smooth_spline}} from the \code{smoothr} package is used to interpolate points between a reduced number of vertices using piecewise cubic polynomials. The number of vertices is calculated based on the number of midline coordinates times numeric value of the list in \code{ml.smooth}. If \code{ml.smooth} contains 'loess', \code{loess} is used to fit a polynomial surface. For contours that have a complicated midline with non-unique x values, loess smoothing can produce poor results. Thus, spline smoothing is usually the advisable option. 
#' 
#' For contour smoothing \code{smooth.n} is passed to the \code{n} parameter of \code{\link{free.ml.ang}}, \code{\link{free.ml.hull}}, or \code{\link{free.ml.del}}, which smooths coordinates using a simple moving average. Users should be wary of oversmoothing by smoothing both the contour (from which the midline is calculated) and the midline.
#'
#' @return A list with the following components:
#'
#' \code{kin.dat} a data table consisting of frame-by-frame position parameters for the ROI
#' \itemize{
#' \item the frame number
#' \item 'x' and 'y': the position of the tail (rightmost or posteriormost)
#' \item 'head.x' and 'head.y': the x and y position of the head (leftmost or anteriormost)
#' \item 'amp': the amplitude of the tail (taken from 'wave.y' of \code{midline})
#' #' \item 'head.pval': p values of the \code{lm()} fit that describes the position of the head as determined by \code{ant.per} (green points in the outputted images)
#' \item 'roi': a character indicating ROI size ('a' being the largest)
#' \item 'edge': was the roi on the edge of the frame}
#'
#' \code{midline} A data table containing, for each frame described by \code{frames}, the following: \itemize{
#' \item 'x' and 'y': x and y positions of the midline of the ROI
#' \item 'x.sm' and 'y.sm': the smoothed midline positions predicted by \code{ml.smooth}.
#' \item 'wave.y': orthogona distance of midline points from predicted midline (see below)
#' \item 'per.bl': the percentage of 'x.sm' along the body length calculated as the cumulative sum of distances between points
#' }
#' 
#' \code{cont} A data table containing x and y positions of the contours used to calculate the data in 'kin.dat'. Contains the following:
#' \itemize{
#' \item 'frame': the frame
#' \item 'x' and 'y': the x and y positions of the contours
#' }
#' 
#' \code{cont.sm} A data table containing x and y positions of the smooth contours. Contains the following:
#' \itemize{
#' \item 'frame': the frame
#' \item 'n': the position of the coordinate. n=1 and max(n) are adjacent at the head
#' \item 'x' and 'y': the x and y positions of the contours
#' }
#'
#'
#'#' \code{all.classes} A data table containing the following for all ROIs detected:
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
#' \code{mid.pred} the theoretical midline based on a linear model established by the anterior section of of the smoothed midline established by \code{ant.per}. Used to calculate \code{midline$wave.y} as the orthogonal distance between the line defined by 'x' and 'mid.pred' and each coordinate defined by '\code{midline$x.sm} and \code{midline$y.sm}. A data table that contains the following:
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
#' @importFrom graphics lines
#' @importFrom stats complete.cases fitted lm loess predict smooth.spline coef dist 
#' @importFrom utils head tail flush.console
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel 
#' @importFrom foreach foreach %do% %dopar%
#' 
#'
#' @seealso \code{\link{kin.search}}, \code{\link{kin.simple}},\code{\link{free.ml.bis}}, \code{\link{free.ml.del}}
#' @examples
#'
#' ##A somewhat long example
#' #### plot midline waveform on images of swimming ropefish
#' \dontrun{
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
#' gg.overlay(kin=kin,data="midline",frames=0:300,size=.2,animate=TRUE,zoom=FALSE,alpha=0.01,col="red",fps=10)
#' 
#' #clean up
#' unlink(paste0(tempdir(),"/images"),recursive=TRUE)
#' unlink(paste0(tempdir(),"/out"),recursive=TRUE)
#' 
#' 
#'}
#'
#'
#' 
kin.free <-
  function(image.dir = NULL,frames=NULL,par=FALSE,cores.n=NULL,ant.per = 0.10,ant.pos=NULL,tips = 0.02,smooth.n=1,red=NULL,ml.meth="hull",ml.smooth = list("spline",0.25),save = FALSE,out.qual = 1,out.dir = NULL,plot.pml = TRUE,flip = TRUE,...) {
    
    #to prevent NSE warnings and NSB notes
    size <-x <-y.pred <-wave.y <-mid.pred <-roi <- prev.dist <-prev.distF <-nxt.dist <-nxt.distF <- x.c <-y.c <- x.sm <-min.x <-range.x <-max.x <-m <-b <- y.sm <-max.n <- n <- keep <-  l <-  head.pval <-  bl <-  per.bl <-  amp <- y <- im <- or.dist <- head.dist <- closer <- n.shift <- x<- dist.next <- dist.current <- close.n <- NULL
  
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
        
        
        ##using free.ml
        
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
            best.class[, list(roi, edge)]
          )
        
        par.res <- list(
          kin = kin.im,
          mid = ml,
          cont = cont.im,
          cont.sm = cont.sm,
          class = roi$classes,
          dim = roi$dim
        )
        
        par.res <- list(par.res)
        #names(par.res=im)
        
        return(par.res)
      }
    
    if (par)
      stopCluster(cl)
    
    kin.dat <- do.call(rbind, lapply(kin.res, function(x)
      x$kin))
    midline.dat <- do.call(rbind, lapply(kin.res, function(x)
      x$mid))
    cont.dat <- do.call(rbind, lapply(kin.res, function(x)
      x$cont))
    cont.sm.dat <- do.call(rbind, lapply(kin.res, function(x)
      x$cont.sm))
    class.dat <- do.call(rbind, lapply(kin.res, function(x)
      x$class))
    dim <-  do.call(rbind, lapply(kin.res, function(x)
      x$dim))[1, ]
    
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
    
    ## if head is moving toward in max(n) position, it's tailward
    direction <-ifelse(distMax.cent$diff<dist0.cent$diff, "tailward", "headward")
    
    #add aligned x,y to kin.dat
    kin.dat[,c("x","y","head.x","head.y"):=ends2[,list(x,y,head.x,head.y)]]
    
    if (direction == "tailward"){
      midline.dat2[, n := max(n):1, by = list(frame)]
      kin.dat[,c("x","y","head.x","head.y"):=ends2[,list(head.x,head.y,x,y)]]
    }
    # qplot(data=midline.dat2[frame>30&frame<40],x=x,y=y,col=n)+facet_wrap(frame~.)
    
   
    #overide head position
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
    
    
    #qplot(d=cont.dat[frame==36],x,y)+geom_point(d=midline.dat2[frame==36],aes(x,y,col=n))+geom_point(d=mid.pred2[frame==36],aes(x,mid.pred))
    #add head.p to kin.dat
    
    kin.dat[, head.pval := head.p]
    
    #dist btw point and line
    
    dist.2d.line <- function(x, y, slope, intercept) {
      b = c(1, intercept + slope)
      c = c(-intercept / slope, 0)
      a = c(x, y)
      v1 <- b - c
      v2 <- a - b
      m <- cbind(v1, v2)
      return(abs(det(m)) / sqrt(sum(v1 * v1)))
    }
    
    
    
    midline.dat2[, wave.y := dist.2d.line(x.sm, y.sm, unlist(head.x.m[paste(frame)]), unlist(head.x.b[paste(frame)])), by =
                   list(frame, n)]
    midline.dat2[, bl := dist.2d(x.sm, dplyr::lag(x.sm), y.sm, dplyr::lag(y.sm)), by =
                   frame]
    midline.dat2[, per.bl := c(0, cumsum(bl[!is.na(bl)]) / sum(bl, na.rm = TRUE)), by =
                   frame][, bl := NULL]
    
    
    kin.dat[, amp := midline.dat2[, list(wave.y = last(wave.y)), by = frame]$wave.y]
    
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
    #qplot(d=cont.sm.dat2,x,y,col=n)+facet_wrap(frame~.)
    
    #qplot(d=cont.sm.dat2[frame==7],x,y,col=n)+geom_point(d=cont.sm.dat2[frame==7&flip==T],aes(x,y),col="re")
    
    setkeyv(cont.sm.dat,c("frame","n"))
    
    ##save
    
    if (save) {
      for (xx in images) {
        frame2 <- which(xx == images)
        trial2 <- gsub("(.*)_\\d+", "\\1", trial)
        img2 <-  EBImage::readImage(xx)
        jpeg(paste0(proc.dir, "/", trial2, "_", sprintf("%03d", frame2), ".jpg"),
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
        kin.dat = kin.dat,
        midline = midline.dat2,
        cont = cont.dat,
        cont.sm = cont.sm.dat2[,list(frame,n,x,y)],
        all.classes = class.dat,
        mid.pred = mid.pred2[, list(frame, x, mid.pred)],
        dim = dim
      )
    )
  }


#' @title Tracking of fin-like extensions of body contours

#' @description  Estimates the amplitudes of regions along a body contour that are protruding. Useful in computing paired-fin amplitudes from contour data produced from  \link{kin.simple} and \link{kin.search}. Also computes a smoothed midline based on the body outline with the fin region removed.
#'
#' @param out a data frame or matrix with 'x' and 'y' data as columns. The first row should be a point near the head, the last near the tail
#' @param fin.pos numeric, a vector of length 2 indicating the start and end of the contour region that contains the fins of interest as a proportion of the body length.
#' @param smooth.n numeric, the number of smoothing operations undertaken by \link{coo_smooth} on the contour described by 'x'. See Details.
#' @param ml.smooth numeric (0-1), the smoothing value for the midline. See Details. 
#' @param red numeric, between 0-1 the proportion of contour coordinates to sample. Will speed up fin position and midline estimations. If 'NULL', the full conour in \code{out} will be used See Details.
#' @param ml.meth character, the midline detection method. One of 'ang' for bisection using \code{\link{free.ml.ang}}, 'hull' for bisection using \code{\link{free.ml.hull}}, or 'del' for Delaunay triangulation using \code{\link{free.ml.del}}. See Details.
#' 
#' @export
#' @importFrom graphics lines
#' @importFrom stats complete.cases fitted lm loess  predict smooth.spline
#' @importFrom utils head setTxtProgressBar tail txtProgressBar
#'
#' @details
#' If \code{red} is specified, the contour \code{out} is sampled with \code{\link{coo_interpolate}} from the \code{Momocs} package. The number of points sampled (n) equals \code{red} times the number of points in \code{out}
#'
#' To establish the contour positions that are within \code{fin.pos}, a midline is estimated using one of two methods specified by \code{ml.meth}, "ang" for \code{\link{free.ml.ang}} or "hull" for  \code{\link{free.ml.hull}}. Midline points are indexed by position along the body length by calculating the cumulative distance between midline coordinates in pixels. This midline distance is then used to estimate the position of the fin about the contour using the paramter \code{fin.pos}.
#' 
#' The positions of the tip of the fin appendages is estimated in two ways. The first is simply the point in each appendage that is farthest from the base of the fin. The base is estimated as a straight line between the contour coordinates that match \code{fin.pos}.  The second is a little more complicated and starts with calculation the distance between each fin contours coordinates and the midpoint of the fin base. \code{\link{features}} from the \code{features} package is then use to calculate an inflection in this distance and the point of this inflection is used to estimate the fin position. Amplitudes of each method are caculated based on the orthogonal euclidean distance from the fin bases.
#'
#'In addition to fin amplitude and contour extraction, this function also produces a composite contour of the body minus the fin area described by \code{fin.pos}. Fin contours are replaced by a simple linear prediction constructed from the coordinates of the first and last values covered by \code{fin.pos}, that is, the fin bases. The result is a straight line between the start and end of each fin. 
#'
#'From this composite body contour, a midline prediction is made based on the  \code{ml.smooth}. The midline is calculated as the midpoints defined between all pairs of coordinates with the same index value
#'
# 'smooth.n' is passed to the 'n' parameter of \code{\link{coo_smooth}}, which smooths coordinates using a simple moving average. Users should be careful not to oversmooth. If the input contour has few points (say just a 100 or so extracted from \code{kin} functions run on low resolution images), much detail will be lost. In general, \code{smooth.n} should be <5.
#'
#' For midline smoothing, \code{\link{smooth_spline}} from the \code{smoothr} package is used to interpolate points between a reduced number of vertices using piecewise cubic polynomials. The number of vertices is calculated based on the number of midline coordinates times \code{1-smooth}. 
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
#' \code{midline} a data table describing the estimated midline, 'x.sm', 'y.sm', the smooth x and y positions, respectively
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
#' @importFrom utils head setTxtProgressBar tail txtProgressBar
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
#' thr.check(fr[1])
#' 
#' #extract contours and other data
#' kin <- kin.free(image.dir = ti,thr=0.9,ant.per = 0.25,red=0.5,smooth.n=2)
#' 
#' #fin data by frame
#' fin.pos <- c(0.25,.5)
#' fin.dat <- fin.kin2(kin=kin,fin.pos = fin.pos,smooth.n=0,red=0.9)
#' 
#' p <- ggplot(dat=fin.dat$amp,aes(x=frame,y=amp2,col=side))+geom_line()+theme_classic(15)
#'print(p)
#'
#'
#' ## plot body and fin contours of frame 8
#' cont <- kin$cont.sm[frame==8,list(x,y)]
#'
#' #plot body contour and fins
#' p <- qplot(data=cont,x=x,y=y)+geom_point(data=fin.dat$fin[frame==8],aes(x,y),col="red",size=3)
#' p+geom_point(data=fin.dat$fin.pts[frame==8],aes(x,y,shape=pos))+xlim(c(0,kin$dim[1]))+ylim(c(0,kin$dim[2]))
#' 
#' #plot body contour minus fins and the body midline
#' p <- qplot(data=fin.dat$comp[frame==8],x=x,y=y)+geom_point(data=fin.dat$midline[frame==8],aes(x,y),col="red",size=2)
#' p+xlim(c(0,kin$dim[1]))+ylim(c(0,kin$dim[2]))
#'
#'unlink(ti,recursive=TRUE)
#' }
#' 
#' out=as.matrix(cont[,list(x,y)])
#'
fin.kin <-
  function(kin,out="cont",frames=NULL,fin.pos = NULL,smooth.n = 5, ml.meth="hull",ml.smooth = 0.9,red=NULL) {
    
    y <- x <- x.sm <- y.sm <- n <- m <- b <- dist2 <- x.c <- y.c <- tip1 <- tip2 <- method <- amp2 <-pos <- y.pred <-side <- ends <- head.dist <- NULL # due to NSE notes in R CMD checks
   
    if(is.null(frames)) frames <-unique(kin[[1]]$frame)
    if(!is.null(frames) & any(!frames %in% c$frame)) stop("not all frames are in 'kin' list")
    
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
    
    #return position of a point rotated theta about another point
    point.ang.orig<- function(p,o,theta){
      xrot<-cos(theta)*(p[1]-o[1])-sin(theta)*(p[2]-o[2])+o[1]
      yrot<-sin(theta)*(p[1]-o[1])+cos(theta)*(p[2]-o[2])+o[2]
      return(c(xrot,yrot))
    }
    
    k <- kin$kin.dat
    c <- kin$cont
    if(out=="cont.sm") c <- kin$cont.sm
    
 
    r <- lapply(frames, function(f,...){
      o <- c[frame==f,list(x,y)]
      
    if (!is.matrix(o))
     out<- as.matrix(o)
    if (is.null(colnames(o)))
      colnames(o) <- c("x", "y")
    
    
    if(ml.meth=="ang") fml <- free.ml.ang(out=as.matrix(o),smooth.n =smooth.n,red=red)
    if(ml.meth=="hull") fml <- free.ml.hull(out=as.matrix(o),smooth.n =smooth.n,red=red)
    
    #with(fml$cont.sides,plot(x,y))
    
    
    ## shift points to match order in input
    
    out2 <- data.table(o,n=1:nrow(o))
    
    #the head
    h <- k[frame==f,]
  
    cont.sm <- fml$cont.sm[,n:=1:.N]
    #qplot(d=cont.sm,x,y,col=n)
    
    cont.flip <- cont.sm[,head.dist:=dist.2d(x,h$head.x,y,h$head.y)][,list(close.n=n[which.min(head.dist)])]
    

    #reset n=1 to head
    cont.sm2 <- data.table(Momocs::coo_slide(as.matrix(cont.sm[,list(x,y)]),cont.flip$close.n))
    cont.sm2[,n:=1:.N]
    #qplot(d=cont.sm2,x,y,col=n)
    
    
    setkeyv(cont.sm2,c("n"))
    
    ml2 <- sm.spline(x=as.matrix(fml$ml[,list(x,y)]),sm=.8)[,n:=1:.N]
    
    #qplot(d=ml2,x.sm,y.sm,col=n)
    min.dist <- with(ml2,dist.2d(x.sm,head$x,y.sm,head$y))
    
    #    #flip cont.sides  and ml2 if needed
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
    #p <- qplot(data= comp,x=x,y=y,col=n)+geom_point(dat=fins,aes(x,y.pred,col=side))
    #print(p)

    comp[,n:=1:.N,by=side]
    
    setkeyv(comp, c("n"))
    
    
   #with(comp,plot(x,y))
    
    fins2 <- copy(fins)
    fins2[, `:=`(pos, "pt"),by=side]
    
    fins2[, `:=`(pos, ifelse(n==min(n),"start",pos)),by=side]
    fins2[, `:=`(pos, ifelse(n==max(n),"end",pos)),by=side]
    fins2[tip1==TRUE,pos:="tip1"]
    fins2[tip2==TRUE,pos:="tip2"]
    
    finPts <- fins2[pos!="pt",list(side,n,x,y,pos)]
    
    amp <- fins2[!is.na(pos)&tip1|tip2,][,method:=pos][,list(amp1=dist[tip1],amp2=dist[tip2]),by=side][!is.na(amp2)]
  
    #with(comp,plot(x,y))
    
    ml.comp <- comp[, list(x = sum(x) / 2, y = sum(y) / 2), by = list(n)]

    #ml.comp.s <- sm.spline(x=as.matrix(ml.comp[,list(x,y)]),sm=ml.smooth)
    ml.comp.s <- ml.comp[,{s <- predict(smooth.spline(data.frame(x,y),spar = ml.smooth));list(x=s$x,y=s$y)}]
    
   # qplot(data= comp,x=x,y=y,col=n)+geom_point(dat=ml.comp.s,aes(x,y,col=n))
    
  bl <- sum(copy(ml.comp.s[,dist:=dist.2d(x,dplyr::lead(x),y,dplyr::lead(y))])$dist,na.rm = TRUE)
    
    # p <- qplot(data= fml$cont.sides,x=x,y=y)+geom_point(dat=fins,aes(x,y,col=side))+geom_point(data=fins[tip1==TRUE],aes(x,y),size=3,col="black")+geom_point(d=cents,aes(x.c,y.c),col="blue")+geom_point(d= ml2,aes(x.sm,y.sm))
    # print(p)
    # 
    return(list(
      cont = data.table(frame=f,cont.sm2), #smoothed contour
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
    
    
    return(list(cont=cont.dat,
                fin=fin.dat,
                fin.pts=fin.pts,
                comp=comp.dat,
                midline=mid.dat,
                amp=amp.dat,
                bl=bl.dat
                ))
}
