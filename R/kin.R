#' @title  Midline tracking over image sequences

#' @description  Automatically retrieves the midline of a detected ROI in each image of a sequence through thresholding and segmentation; finds the y-value midpoint along the x-value array of the ROI and fits a midline according to a chosen smoothing method (loess or spline). Also outputs the midline amplitude relative to a reference line determined by an anterior section of the ROI. Supported image formats are jpeg, png, and tiff.
#'
#' @param image.dir character, directory containing images to analyze.
#' @param frames numeric, vector indicating which images to process.
#' @param thr numeric or character ('otsu') threshold to determine binary image. See Details.
#' @param ant.per numeric; left-most percentage of ROI that establishes the horizontal reference for the midline displacement.
#' @param tips, numeric, the proportion the the midline data to use in calculation of the head and tail position.
#' @param plot.pml logical, value indicating if outputted images should include an overlay of the theoretical midline based on \code{ant.per}.
#' @param smoothing character, the midline smoothing method, either 'loess' or "spline".
#' @param smooth numeric; if \code{smoothing} is set to 'loess', smoothing parameter value for plotted midline.
#' @param smooth.points numeric, number of equally spaced points along the ROI midline on which the smoothed midline is computed.
#' @param flip logical, indicating if binary should be flipped.
#' @param show.prog logical value indicating if outputted image should be displayed during analysis.
#' @param size.min numeric, indicating the minimum size of ROIs as a proportion of the pixel field to be considered in analysis. May be useful if smaller unimportant ROIs appear in the frame. Default is 0.02.
#' @param save logical, value indicating if images should be outputted with midline and predicted midline based on the \code{ant.per} \code{lm()} overlaying original or binary images.
#' @param out.qual, numeric, a value between 0-1 representing the quality of outputted images. Ignored if \code{save=FALSE}.
#' @param out.dir character, the directory to which outputted images should be saved.
#' @param image.type character; the type of image to be outputted, either 'orig' or 'bin' representing the original or binary images, respectively. Ignored if 'save=FALSE'.
#' @param search.for character, the search parameter. See Details.
#' @param edges logical, should ROIs on image edges be evaluated. See Details.
#' @param border if \code{edges=TRUE}, size of border to add in pixels. Dee details.
#'
#' @export
#'
#'
#'
#' @details
#'The algorithm assumes a left-right orientation, i.e., the head of the ROI is positioned left, the tail right. The \code{ant.per} value therefor establishes the reference line (theoretical straight midline) based on that portion of the head. The midline is calculated as the midpoints between the y extrema for each x position. Chooses ROIs based on relative ROI size or position.
#'
#'Thresholding operations can be performed with an arbitrary (user defined) numeric value or with Otsu's method ('thr="otsu"'). The latter chooses a threshold value by minimizing the combined intra-class variance. See \code{\link{otsu}}.
#'
#'If 'edges=TRUE', it is best to add an artificial border so that any part of the ROI in contact with the edge can be distinguished from it.
#'
#' \code{search.for} determines how ROIs are chosen:
#' \itemize{
#' \item "offset", the ROI with a centroid that is the shortest linear distance to the center of the field
#' \item "offset.x", the ROI with a centroid x position that is closest to the x position of the center of the field
#' \item "offset.y", the ROI with a centroid y position that is closest to the y position of the center of the field
#' \item "largest", the largest ROI.
#' }
#'
#' These choices will be made on ROI sets that are not on the edge of the field if 'edges=FALSE'.
#'
#' \code{edges} Set by default to 'FALSE'. It is not advisable to include shapes that are on the edge of any frame and are therefore incomplete.	Yet, if set to 'TRUE', the \code{border} adds a black border to the image so that the intended ROI may be distinguished from the edge.
#'
#'\code{image.type} Can be set as "orig" or "bin". "orig" plots midline and reference lines over the original video frames, "bin" over binary images.
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
#' \code{kin.dat} a data frame consisting of frame-by-frame position parameters for the ROI indicated by \code{n.blob}:
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
#' @importFrom utils head setTxtProgressBar tail txtProgressBar
#'
#'
#' @seealso \code{\link{kin.simple}}
#' @examples
#'
#' #### plot lot caudal amplitude and produce a classic midline waveform plot of swimming fish
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
#'        search.for = "largest",
#'       smoothing = "loess",frames=1:50,
#'       out.dir=paste0(tempdir(),"/processed_images"),
#'       show.prog = FALSE,thr = "otsu",
#'       image.type="bin",smooth=0.4)
#'
#' #plot instantaneous amplitude of tail (last/rightmost point) over frames
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
#' pal <- wes_palette("Zissou1", 100, type = "continuous") #"Zissou" color palette
#'
#' p <- ggplot(dat=ml,aes(x=x2,y=wave.y))+theme_classic(15)+scale_color_gradientn(colours = pal)
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
#' ## A very short example.
#'
#' #retrieve image
#' i <- EBImage::readImage(system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter"))
#' #create directory and write image to it
#' t <- tempdir()
#'
#'
#' dir.create(paste0(t,"/images"))
#' EBImage::writeImage(i,paste0(t,"/images/sunfish001.jpg"),type = "jpeg")
#'
#' list.files(paste0(t,"/images"))
#' #run kin.search and save output image to directory
#' kin.i<- kin.search(image.dir = paste0(t,"/images"),smooth=0.7,save = TRUE,out.dir = t)
#'
#' #plot midline over original image
#' with(kin.i$midline,plot(x,wave.y))
#'
#' i2 <- EBImage::readImage(paste0(t,"/sunfish001_000.jpg"))
#' EBImage::display(i2,method="raster")
#'
#' #clean up
#' unlink(paste0(t,"/images"),recursive=TRUE)
#'

kin.search <-
  function(image.dir = NULL,
           frames = NULL,
           thr = "otsu",
           plot.pml = TRUE,
           show.prog = FALSE,
           ant.per = 0.10,
           tips = 0.02,
           smoothing = "loess",
           smooth = 0.25,
           smooth.points = 200,
           image.type = "orig",
           save = TRUE,
           out.qual = 1,
           out.dir = NULL,
           flip = TRUE,
           size.min = 0.02,
           search.for = "largest",
           edges = FALSE,
           border = 5) {
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
    
    if (!search.for %in% c("offset.x", "offset.y", "largest", "offset"))
      stop("'search.for' must be set to 'offset', 'offset.x','offset.y', or 'largest')")
    
    trial <- gsub("\\.[^.]*$", "", basename(images[1]))
    
    kin.l <- list()
    midline.l <- list()
    classes.l <- list()
    lms <- list()
    conts <- list()
    
    
    
    roi.outs <- list() #store the rois for each image
    
    for (im in images) {
      frame <- which(im == images) - 1
      
      img <- EBImage::readImage(im, all = FALSE)
      
      img.dim <- dim(img)[1:2]
      
      
      # computes binary mask
      if (thr != "otsu" &
          !is.numeric(thr))
        stop("'thr' must be set to 'otsu' or a numeric value=0-1")
      if (thr == "otsu") {
        EBImage::colorMode(img) = EBImage::Grayscale
        thr <- EBImage::otsu(img)[1]
      }
      
      y = img > thr #contrast threshold
      
      if (flip) {
        #flip binary
        y[y == 1] <- 5
        y[y == 0] <- 1
        y[y == 5] <- 0
      }
      z = EBImage::bwlabel(y)
      
      if (edges) {
        bord <- border * 2
        
        z1 <-
          EBImage::resize(
            z,
            w = dim(z)[1],
            h = dim(z)[2],
            output.dim = c(dim(z)[1] + bord, dim(z)[2] + bord),
            bg.col = "white"
          )
        
        EBImage::translate(z1, v = c(border, border), bg.col = "black")
        
      }
      
      rois <- tabulate(z)
      
      pix <- dim(z[, , 1])[1] * dim(z[, , 1])[2]
      w <- dim(z[, , 1])[1] #width of image
      h <- dim(z[, , 1])[2] #height of image
      per <- rois / (w * h) #how big are rois compared to pixel field
      
      c.roi <-
        which(per >= size.min) #candidate rois, filtered by size of of pixel field
      
      names(c.roi) <-
        as.factor(letters[order(rois[c.roi], decreasing = TRUE)])
      
      z.l <- list()
      out.l <- list()
      
      for (r in c.roi) {
        r.name <- as.character(names(c.roi)[c.roi == r])
        z.r <- z
        z.r[z != r] <- 0
        z.r[z == r] <- 1
        z.m <- z.r[, , 1]
        z.m[1, 1] <- 0 #this gets a 1 when
        z.l[[r.name]] <- z.m
        
        z.c <- EBImage::ocontour(z.m)
        
        wall <-
          any(z.c[[1]][, 1] > dim(z)[1] - 2 |
                z.c[[1]][, 1] < 2  | z.c[[1]][, 2] > dim(z)[2] - 2 |
                z.c[[1]][, 2] < 2)
        
        r.out <- Out(EBImage::ocontour(z.m))
        if (wall)
          edge <- TRUE
        if (!wall)
          edge <- FALSE
        r.out$fac <-
          data.frame(
            shape = paste0("roi-", r.name),
            type = paste0("roi"),
            edge = edge
          )
        out.l[[r.name]] <- r.out
        rois[c.roi[r.name]]
        
      }
      
      #don't combine if only one ROI
      if (length(out.l) == 1) {
        roi.out2 <- out.l[[1]]
      } else{
        roi.out2 <- Momocs::combine(out.l)
      }
      
      #centroids
      cent <- coo_centpos(roi.out2) #centroids
      offset.x <- abs((img.dim[1] / 2) - cent[, "x"])
      offset.y <- abs((img.dim[2] / 2) - cent[, "y"])
      
      offset <-
        apply(cent, 1, function(x)
          dist.2d(x["x"], img.dim[1] / 2, x["y"], img.dim[2] / 2))
      
      classes <-
        data.table(
          roi = gsub("roi-", "", roi.out2$fac$shape),
          edge = roi.out2$fac$edge,
          size = rois[c.roi],
          offset.x = offset.x,
          offset.y = offset.y,
          offset = offset
        )
      
      classes.l[[paste0(frame)]] <- data.table(frame = frame, classes)
      
      if (edges == FALSE) {
        classes <- classes[edge == FALSE]
        if (nrow(classes) == 0)
          stop("no ROI found that is not on edge")
      }
      
      if (search.for == "offset") {
        z.best <- z.l[[classes[which.min(offset), ]$roi]] #best segment
        r.name <- classes[which.min(offset), ]$roi
        best.class <- classes[which.min(offset), ]
      }
      
      if (search.for == "offset.x") {
        z.best <- z.l[[classes[which.min(offset.x), ]$roi]] #best segment
        r.name <- classes[which.min(offset.x), ]$roi
        best.class <- classes[which.min(offset.x), ]
      }
      
      if (search.for == "offset.y") {
        z.best <- z.l[[classes[which.min(offset.y), ]$roi]] #best segment
        r.name <- classes[which.min(offset.y), ]$roi
        best.class <- classes[which.min(offset.y), ]
      }
      
      if (search.for == "largest") {
        z.best <- z.l[[classes[which.max(size), ]$roi]] #best segment
        r.name <- classes[which.max(size), ]$roi
        best.class <- classes[which.max(size), ]
      }
      
      
      if (show.prog) {
        suppressMessages(EBImage::display(z.best, method = "raster"))
      }
      
      best.cont <- data.table(EBImage::ocontour(z.best)[[1]])
      colnames(best.cont) <- c("x", "y")
      
      conts[[paste0(frame)]] <- data.table(frame = frame, best.cont)
      
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
      
      #n midline points
      if (is.null(smooth.points))
        smooth.points <- nrow(y.df)
      midline <-
        y.df[seq(1, nrow(y.df), length.out = smooth.points), ] #two hundred points on midline
      
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
      head.dat <- midline[1:(ant.per * smooth.points), ]
      head.lm <- lm(y.pred ~ x, head.dat)
      
      head.p <- summary(head.lm)$r.squared #how well does head lm fit
      
      midline$mid.pred <-
        predict(head.lm, newdata = midline)#add lm prediction to midline df
      
      midline <- midline[complete.cases(midline), ]
      
      midline[, wave.y := y.pred - mid.pred] #wave y based on midline y and straight head.lm pred points
      
      midline[, roi := r.name]
  
      
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
          width = w,
          height = h
        )
        if (image.type == "bin")
          suppressMessages(EBImage::display(z, method = "raster"))
        if (image.type == "orig")
          suppressMessages(EBImage::display(img, method = "raster"))
        
        
        if (plot.pml)
          lines(predict(lm(mid.pred ~ x, midline)),
                x = midline$x,
                col = "blue",
                lwd = 4)
        with(midline, lines(y.pred ~ x, col = "red", lwd = 4))
        if (plot.pml)
          with(midline[1:ceiling(ant.per * smooth.points), ],
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
        dim = img.dim
      )
    )
  }

#' @title  Simplified midline tracking over image sequences

#' @description  Automatically retrieves the midline of a detected ROI based on size. Assumes the ROI of interest is the largest detected and not intersecting the edges of the image frame, conditions often met in kinematic studies. For each ROI of interest, finds the y-value midpoint along the x-value array of the ROI and fits a midline according to a chosen smoothing method (loess or spline). Also outputs the midline amplitude relative to a reference line determined by an anterior section of the ROI and outputs contours ROIs in each frame for subsequent analysis. Supported image formats are jpeg, png, and tiff.
#'
#'
#' @param image.dir character, directory containing images to analyze.
#' @param frames numeric, vector indicating which images to process.
#' @param thr numeric or character ('otsu') threshold to determine binary image. See Details.
#' @param size.min numeric, indicating the minimum size of ROIs as a proportion of the pixel field to be considered in analysis. May be useful if smaller unimportant ROIs appear in the frame. Default is 0.05.
#' @param ant.per numeric; left-most proportion of ROI that establishes the horizontal reference for the midline displacement.
#' @param tips, numeric, the proportion the the midline data to use in calculation of the head and tail position.
#' @param smoothing character, the midline smoothing method, either 'loess' or 'spline'.
#' @param smooth numeric; if \code{smoothing} is set to 'loess', passed to 'span' parameter of \code{\link{loess}}. If \code{smoothing} is set to 'spline', passed to 'spar' parameter of \code{\link{smooth.spline}}
#' @param smooth.points numeric, number of equally spaced points along the ROI midline on which the smoothed midline is computed.
#' @param save logical, value indicating if images should be outputted with midline and predicted midline based on the \code{lm()} predictions from \code{ant.per}overlaying original or binary images.
#' @param out.qual, numeric, a value between 0-1 representing the quality of outputted images. Ignored if \code{save=FALSE}.
#' @param out.dir character, the directory to which outputted images should be saved.
#' @param plot.pml logical, value indicating if outputted images should include the predicted midline (in blue) and the points according to \code{ant.per} used to construct the predicted midline (in green).
#' @param image.type character; the type of image to be outputted, either 'orig' or 'bin' representing the original or binary images, respectively. Ignored if 'save=FALSE'.
#' @param flip logical, indicating if binary image should be flipped.
#' @param show.prog logical, indicating if outputted image should be displayed during analysis.
#'
#' @export
#'
#' @details
#'The algorithm assumes a left-right orientation, i.e., the head of the ROI is positioned left, the tail right. ffmpeg operations or even imageJ can rotate images not in this orientation. The \code{ant.per} value therefore establishes the reference line (theoretical straight midline) based on that portion of the head. The midline is calculated as the midpoints between the y extrema for each x position.
#'
#'If 'save=TRUE', images are saved as binary or the original with a body midline  overlay and, if chosen, with the theoretical midline (based on \code{ant.per}).
#'
#'Thresholding operations can be performed with an arbitrary (user defined) numeric value or with Otsu's method ('thr="otsu"'). The latter chooses a threshold value by minimizing the combined intra-class variance. See \code{\link{otsu}}.
#'
#' @return A list with the following components:
#'
#' \code{kin.dat} a data table consisting of frame-by-frame position parameters for the ROI determined by LDA analysis.
#' \itemize{
#' \item the frame number
#' \item 'x' and 'y': the position of the tail (rightmost or posteriormost)
#' \item 'head.x' and 'head.y': the x and y position of the head (leftmost or anteriormost)
#' \item 'amp': the amplitude (\code{amp}) of the tail relative to the theoretical midline determined by the \code{lm()} predictions from \code{ant.per}
#' \item 'roi': a character indicating the ROI ranked by size ('a' being the largest)
#' \item 'head.pval': p values of the \code{lm()} fit that describes the position of the head as determined by \code{ant.per} (green points in the outputted images/video)}
#'
#' \code{midline} A data table containing, for each frame described by \code{frames}, the following: \itemize{
#' \item 'x' and 'y.m': x and y positions of the midline of the ROI
#' #' \item 'y.min' and 'y.max': min and max y positions ROI's contour used in y.m calculation
#' \item 'mid.pred': the predicted linear midline based on the points/pixels defined by \code{ant.per} (green points in the outputted images/video if 'plot.pml=TRUE')
#' \item 'y.pred': midline points fit to a smooth spline or loess model with spar or span equal to \code{smooth} (red curve in the outputted images/video)
#' \item 'wave.y': midline points 'y.pred' relative to 'mid.pred'
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
#' \item 'size': the size of each ROI
#' }
#'
#' \code{dim} the x and y dimensions of the images analyzed
#' @seealso \code{\link{kin.search}}
#' @export
#'
#' @importFrom graphics lines
#' @importFrom stats complete.cases fitted lm loess predict smooth.spline
#' @importFrom utils head setTxtProgressBar tail txtProgressBar
#' @importFrom grDevices dev.off jpeg
#' @importFrom EBImage bwlabel otsu
#'
#' @examples
#' #### plot caudal amplitude and produce a classic midline waveform plot of swimming fish
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
#' kin <- kin.simple(image.dir =paste0(tempdir(),"/example"),
#'       smoothing = "loess",frames=1:50,
#'       out.dir=paste0(tempdir(),"/processed_images"),
#'       show.prog = FALSE,thr = "otsu",
#'       image.type="bin",smooth=0.4)
#'
#' #plot instantaneous amplitude of tail (last/rightmost point) over frames
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
#' pal <- wes_palette("Zissou1", 100, type = "continuous") #"Zissou" color palette
#'
#' p <- ggplot(dat=ml,aes(x=x2,y=wave.y))+theme_classic(15)+scale_color_gradientn(colours = pal)
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
#' ## A very short example.
#'
#' #retrieve image
#' i <- EBImage::readImage(system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter"))
#' #create directory and write image to it
#' t <-tempdir()
#' dir.create(paste0(t,"/images"))
#' EBImage::writeImage(i,paste0(t,"/images/sunfish001.jpg"),type = "jpeg")
#'
#' #run kin.search and save output image to directory
#' kin.i<- kin.simple(image.dir = paste0(t,"/images"),save = TRUE,out.dir = t)
#'
#' #plot midline
#' with(kin.i$midline,plot(x,wave.y))
#' i2 <- EBImage::readImage(paste0(t,"/sunfish001_000.jpg"))
#' EBImage::display(i2,method="raster")
#' #clean up
#' unlink(paste0(t,"/images"),recursive=TRUE)


kin.simple <-
  function(image.dir = NULL,
           frames = NULL,
           thr = 0.7,
           size.min = 0.05,
           ant.per = 0.20,
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
           show.prog = FALSE) {
    size <-
      x <-
      y.pred <-
      wave.y <-
      mid.pred <- roi <- NULL # to avoid NSE erros or R CD check
    
    if (!file.exists(image.dir))
      stop("Directory specified by 'image.dir' (",
           paste0(image.dir),
           ") does not exist")
    
    if (!save &
        !is.null(out.dir))
      stop("'out.dir' specified but 'save=FALSE'. To save processed images, 'save' must be 'TRUE'")
    
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
      
      img <-
        EBImage::readImage(im, all = FALSE) #if don't add package, others use "display"
      
      img.dim <- dim(img)[1:2]
      
      # computes binary mask
      if (thr != "otsu" &
          !is.numeric(thr))
        stop("'thr' must be set to 'otsu' or a numeric value=0-1")
      if (thr == "otsu") {
        EBImage::colorMode(img) = EBImage::Grayscale
        thr <- EBImage::otsu(img)[1]
      }
      
      y = img > thr #contrast threshold
      
      if (flip) {
        #flip binary
        y[y == 1] <- 5
        y[y == 0] <- 1
        y[y == 5] <- 0
      }
      z = EBImage::bwlabel(y)
      
      rois <- tabulate(z)
      
      pix <- dim(z[, , 1])[1] * dim(z[, , 1])[2]
      w <- dim(z[, , 1])[1] #width of image
      h <- dim(z[, , 1])[2] #height of image
      per <- rois / (w * h) #how big are rois compared to pixel field
      
      c.roi <-
        which(per >= size.min) #candidate rois, filtered by size of of pixel field
      
      names(c.roi) <-
        as.factor(letters[order(rois[c.roi], decreasing = TRUE)])
      
      z.l <- list()
      out.l <- list()
      
      for (r in c.roi) {
        r.name <- as.character(names(c.roi)[c.roi == r])
        z.r <- z
        z.r[z != r] <- 0
        z.r[z == r] <- 1
        z.m <- z.r[, , 1]
        z.m[1, 1] <- 0 #this gets a 1 when
        z.l[[r.name]] <- z.m
        
        z.c <- EBImage::ocontour(z.m)
        
        
        wall <-
          any(z.c[[1]][, 1] > dim(z)[1] - 2 |
                z.c[[1]][, 1] < 2  | z.c[[1]][, 2] > dim(z)[2] - 2 |
                z.c[[1]][, 2] < 2)
        
        
        r.out <- Out(list(EBImage::ocontour(z.m)))
        if (wall)
          edge <- TRUE
        if (!wall)
          edge <- FALSE
        r.out$fac <-
          data.frame(
            shape = paste0("roi-", r.name),
            type = paste0("roi"),
            edge = edge
          )
        out.l[[r.name]] <- r.out
        rois[c.roi[r.name]]
        
      }
      #don't combine if only one ROI
      if (length(out.l) == 1) {
        roi.out2 <- out.l[[1]]
      } else{
        roi.out2 <- Momocs::combine(out.l)
      }
      
      
      
      classes <-
        data.table(
          roi = gsub("roi-", "", roi.out2$fac$shape),
          edge = roi.out2$fac$edge,
          size = rois[c.roi]
        )
      
      classes.l[[paste0(frame)]] <- data.table(frame = frame, classes)
      z.best <- z.l[[classes[edge == FALSE, ][which.max(size)]$roi]]
      r.name <- classes[edge == FALSE, ][which.max(size)]$roi
      best.class <- classes[edge == FALSE, ][which.max(size)]
      
      
      if (show.prog) {
        EBImage::display(z.best, method = "raster")
      }
      
      best.cont <- data.table(EBImage::ocontour(z.best)[[1]])
      colnames(best.cont) <- c("x", "y")
      
      conts[[paste0(frame)]] <- data.table(frame = frame, best.cont)
      
      y.df <-
        best.cont[, list(y.min = min(y),
                         y.max = max(y),
                         y.m = mean(y)), by = list(x)]
      setkey(y.df, "x")
      
      ends <- ceiling(nrow(y.df) * tips)
      tip.y <-
        mean(tail(y.df$y.m[!is.na(y.df$y.m)], ends))
      tip.x <-
        mean(tail(y.df$x[!is.na(y.df$y.m)], ends))
      
      head.y <-
        mean(head(y.df$y.m[!is.na(y.df$y.m)], ends))
      
      head.x <-
        mean(head(y.df$x[!is.na(y.df$y.m)], ends))
      
      #n midline points
      if (is.null(smooth.points))
        smooth.points <- nrow(y.df)
      midline <-
        y.df[seq(1, nrow(y.df), length.out = smooth.points), ] #two hundred points on midline
      
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
        ml.pred <- smooth.spline(x = midline$x,
                                 y = midline$y.m,
                                 spar = smooth)$y
      
      midline[, y.pred := ml.pred]#add smoothed predictions
      
      #head section
      head.dat <- midline[1:(ant.per * smooth.points), ]
      head.lm <- lm(y.pred ~ x, head.dat)
      
      head.p <- summary(head.lm)$r.squared #how well does head lm fit
      
      midline$mid.pred <-
        predict(head.lm, newdata = midline)#add lm prediction to midline df
      
      midline <- midline[complete.cases(midline), ]
      
      midline[, wave.y := y.pred - mid.pred] #wave y based on midline y and straight head.lm pred points
      
      midline[, roi := r.name]
      
      
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
          width = w,
          height = h
        )
        if (image.type == "bin")
          EBImage::display(z, method = "raster")
        if (image.type == "orig")
          EBImage::display(img, method = "raster")
        
        
        if (plot.pml)
          lines(predict(lm(mid.pred ~ x, midline)),
                x = midline$x,
                col = "blue",
                lwd = 4)
        with(midline, lines(y.pred ~ x, col = "red", lwd = 4))
        if (plot.pml)
          with(midline[1:ceiling(ant.per * smooth.points), ],
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
        dim = img.dim
      )
    )
  }


#' @title  Midline tracking of free-moving ROIs over image sequences

#' @description  Automatically retrieves the midline of detected ROIs that are free to move in the spatial field of an image sequence. Does so through threshholding and segmentation; finds the midpoint coordinates from each half of the ROI sharing the same index then calculates a smoothed midline according to a chosen smoothing method (loess or smooth spline). Also outputs the midline amplitude relative to a reference line determined by an anterior section of the ROI. Supported image formats are jpeg, png, and tiff. Supports parallel processing of frames.
#' 
#' @param image.dir character, directory containing images to analyze.
#' @param frames numeric, vector indicating which images to process. Must be >1. See Details.
#' @param par logical, should the frames be processed in parallel using \code{cores.n}.
#' @param cores.n numeric, the number of CPU cores to use if \code{par=TRUE}. If \code{cores.n=NULL} (the default), the total number of cores minus 1 are used.
#' @param thr numeric or character ('otsu'); the threshold to determine binary image. See Details.
#' @param ant.per numeric; anterior percentage of ROI that establishes the reference for the midline displacement.
#' @param tips, numeric, the proportion the the midline data to use in calculation of the head and tail position.
#' @param smooth.n, numeric, the number of contour smoothing iterations. See Details.
#' @param red numeric, between 0-1 the proportion of contour coordinates to sample for midline estimates. Will speed up midline estimations. If 'NULL', the full contour retrieved from the ROI will be passed to \code{\link{free.ml}}.
#' @param ml.smooth a list of length two with unnamed components including a character string specifying the midline smoothing method, either 'loess' or "spline", and a numeric value specifying the amount of smoothing. See Details.
#' @param save logical, value indicating if images should be outputted with midline and predicted midline based on the \code{ant.per} \code{lm} overlaying original or binary images.
#' @param out.qual, numeric, a value between 0-1 representing the quality of outputted images. Ignored if \code{save=FALSE}.
#' @param out.dir character, the directory to which outputted images should be saved.
#' @param plot.pml logical, value indicating if outputted images should include an overlay of the midline, head region and theoretical midline based on \code{ant.per}.
#' @param flip logical, indicating if binary image should be flipped.
#' @param size.min numeric, indicating the minimum size of ROIs as a proportion of the pixel field to be considered in analysis. May be useful if smaller unimportant ROIs appear in the frame. Default is 0.02.
#' @param search.for character, the search parameter. See Details.
#' @param edges logical, should ROIs on image edges be evaluated. See Details.
#' @param border if \code{edges=TRUE}, size of border to add in pixels. Dee details.
#'
#' @export
#'
#' @details
#' The midline is calculated as the midpoint coodinates between the extrema of the coordinates describing each ROIs contour (i.e, the coordinates spanning the longest euclidean distance). Chooses ROIs based on relative ROI size or position.
#'
#' The position of the anterior of the ROI (that which is moving forward in the field) is determined by the displacement of the ROI between the first two frames. Thus, \code{frames} must be >1.
#' 
#' The midline is determine by first finding the tips of the ROI (i.e., the two coordinates in the outline that are farthest from one another) with \code{\link{free.ml}} and therefore assumes the ROI is elongate and moving along this long axis. Using \code{\link{free.ml}}, the function bisects the ROI contour at the tips, giving it two sides with coordinates of equal length. The midline coordinates are calculated as the midpoints defined between all pairs of coordinates with the same index value.
#'
#'Thresholding operations can be performed with an arbitrary (user defined) numeric value or with Otsu's method ('thr="otsu"'). The latter chooses a threshold value by minimizing the combined intra-class variance. See \code{\link{otsu}}.
#'
#'If 'edges=TRUE', it is best to add an artificial border so that any part of the ROI in contact with the edge can be distinguished from it.
#'
#' \code{search.for} determines how ROIs are chosen:
#' \itemize{
#' \item "offset", the ROI with a centroid that is the shortest linear distance to the center of the field
#' \item "offset.x", the ROI with a centroid x position that is closest to the x position of the center of the field
#' \item "offset.y", the ROI with a centroid y position that is closest to the y position of the center of the field
#' \item "largest", the largest ROI.
#' }
#'
#' These choices will be made on ROI sets that are not on the edge of the field if 'edges=FALSE'.
#'
#' \code{edges} Set by default to 'FALSE'. It is not advisable to include shapes that are on the edge of any frame and are therefore incomplete.	Yet, if set to 'TRUE', the \code{border} adds a black border to the image so that the intended ROI may be distinguished from the edge.
#'
#' For midline smoothing, if \code{ml.smooth} contains 'spline' (the default), \code{\link{smooth_spline}} from the \code{smoothr} package is used to interpolate points between a reduced number of vertices using piecewise cubic polynomials. The number of vertices is calculated based on the number of midline coordinates times numeric value of the list in \code{ml.smooth}. If \code{ml.smooth} contains 'loess', \code{loess} is used to fit a polynomial surface. For contours that have a complicated midline with non-unique x values, loess smoothing can produce poor results. Thus, spline smoothing is usually the advisable option. 
#' 
#' For contour smoothing \code{smooth.n} is passed to the \code{n} parameter of \code{\link{free.ml}}, which smooths coordinates using a simple moving average. Users should be wary of oversmoothing by smoothing both the contour (from which the midline is calculated) and the midline.
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
#' @seealso \code{\link{kin.search}}, \code{\link{kin.simple}},\code{\link{free.ml}}
#' @examples
#'
#' #### plot lot caudal amplitude and produce a classic midline waveform plot of swimming fish
#' ##A very long example.
#' \dontrun{
#'
#' #download example images and place in 'example' subdirectory
#' f <- "https://github.com/ckenaley/exampledata/blob/master/ropefish.avi?raw=true"
#'
#' download.file(f, paste0(tempdir(),"/ropefish.avi"))
#'
#' dir.create(paste0(tempdir(),"/images"))
#' dir.create(paste0(tempdir(),"/out"))
#' 
#' vid.to.images(paste0(tempdir(),"/ropefish.avi"), out.dir = paste0(tempdir(),"/images"))
#' 
#' kin <- kin.free(image.dir =paste0(tempdir(),"/images"),
#'       par=FALSE,
#'       save=TRUE,
#'       out.dir=paste0(tempdir(),"/out"),
#'       ml.smooth=list("spline",0.9),
#'       thr = "otsu",
#'       size.min=0.01,
#'       frame=1:20,
#'       red=0.5
#'       )
#'
#' fi <- list.files(paste0(tempdir(),"/out"),full.names=T)
#' EBImage::display(EBImage::readImage(fi[1]))
#' 
#' #plot instantaneous amplitude of tail over frames
#' p <- ggplot(dat=kin$kin.dat,aes(x=frame,y=amp))+geom_line()+geom_point()+theme_classic(15)
#' print(p)
#'
#' # midline plot
#' ml <- kin$midline
#' p <- ggplot(dat=ml,aes(x=x.sm,y=y.sm,col=frame))+geom_point()+theme_classic(15)
#' print(p)
#'
#'p <- ggplot(dat=kin$cont.sm[frame==1],aes(x=x,y=y,col=n))+geom_point()
#' print(p)
#' 
#'
#' unlink(paste0(tempdir(),"/images"),recursive=TRUE)
#' unlink(paste0(tempdir(),"/out"),recursive=TRUE)
#'}
#'

kin.free <-
  function(image.dir = NULL,frames=NULL,par=FALSE,cores.n=NULL,thr = "otsu",ant.per = 0.10,tips = 0.02,smooth.n=5,red=NULL,ml.smooth = list("spline",0.25),save = FALSE,out.qual = 1,out.dir = NULL,plot.pml = TRUE,flip = TRUE,size.min = 0.02,search.for = "largest",edges = FALSE,border = 5) {
    
    #to prevent NSE warnings and NSB notes
    size <-
      x <-
      y.pred <-
      wave.y <-
      mid.pred <-
      roi  <-
      prev.dist <-
      prev.distF <-
      nxt.dist <-
      nxt.distF <-
      x.c <-
      y.c <-
      x.sm <-
      min.x <-
      range.x <-
      max.x <-
      m <-
      b <-
      y.sm <-max.n <- n <- 
      keep <-  l <-  head.pval <-  bl <-  per.bl <-  amp <- y <- im <- or.dist <- head.dist <- closer <- n.shift <- x<- NULL
    
    
    
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
    
    if (!search.for %in% c("offset.x", "offset.y", "largest", "offset"))
      stop("'search.for' must be set to 'offset', 'offset.x','offset.y', or 'largest')")
  
    
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
    
    
    kin.res <-
      foreach(
        im = images,
        .combine = append,
        .packages = c("data.table", "Momocs", "trackter"),
        .export = c("free.ml")
      ) %fun% {
        frame <- which(im == images) - 1
        
        img <- EBImage::readImage(im, all = FALSE)
        
        img.dim <- dim(img)[1:2]
        
        
        # computes binary mask
        if (thr != "otsu" &
            !is.numeric(thr))
          stop("'thr' must be set to 'otsu' or a numeric value=0-1")
        if (thr == "otsu") {
          EBImage::colorMode(img) = EBImage::Grayscale
          thr <- EBImage::otsu(img)[1]
        }
        
        y = img > thr #contrast threshold
        
        if (flip) {
          #flip binary
          y[y == 1] <- 5
          y[y == 0] <- 1
          y[y == 5] <- 0
        }
        z = EBImage::bwlabel(y)
        
        #EBImage::display(z)
        if (edges) {
          bord <- border * 2
          
          z1 <-
            EBImage::resize(
              z,
              w = dim(z)[1],
              h = dim(z)[2],
              output.dim = c(dim(z)[1] + bord, dim(z)[2] + bord),
              bg.col = "white"
            )
          
          EBImage::translate(z1, v = c(border, border), bg.col = "black")
          
        }
        
        rois <- tabulate(z)
        
        pix <- dim(z[, , 1])[1] * dim(z[, , 1])[2]
        w <- dim(z[, , 1])[1] #width of image
        h <- dim(z[, , 1])[2] #height of image
        per <-
          rois / (w * h) #how big are rois compared to pixel field
        
        c.roi <-
          which(per >= size.min) #candidate rois, filtered by size of of pixel field
        
        names(c.roi) <-
          as.factor(letters[order(rois[c.roi], decreasing = TRUE)])
        
        z.l <- list()
        out.l <- list()
        
        for (r in c.roi) {
          r.name <- as.character(names(c.roi)[c.roi == r])
          z.r <- z
          z.r[z != r] <- 0
          z.r[z == r] <- 1
          z.m <- z.r[, , 1]
          z.m[1, 1] <- 0 #this gets a 1 when
          z.l[[r.name]] <- z.m
          
          z.c <- EBImage::ocontour(z.m)
          
          wall <-
            any(z.c[[1]][, 1] > dim(z)[1] - 2 |
                  z.c[[1]][, 1] < 2  |
                  z.c[[1]][, 2] > dim(z)[2] - 2 |
                  z.c[[1]][, 2] < 2)
          
          r.out <- Momocs::Out(EBImage::ocontour(z.m))
          if (wall)
            edge <- TRUE
          if (!wall)
            edge <- FALSE
          r.out$fac <-
            data.frame(
              shape = paste0("roi-", r.name),
              type = paste0("roi"),
              edge = edge
            )
          out.l[[r.name]] <- r.out
          rois[c.roi[r.name]]
          
        }
        
        #don't combine if only one ROI
        if (length(out.l) == 1) {
          roi.out2 <- out.l[[1]]
        } else{
          roi.out2 <- Momocs::combine(out.l)
        }
        
        #centroids
        cent <- Momocs::coo_centpos(roi.out2) #centroids
        offset.x <- abs((img.dim[1] / 2) - cent[, "x"])
        offset.y <- abs((img.dim[2] / 2) - cent[, "y"])
        
        offset <-
          apply(cent, 1, function(x)
            dist.2d(x["x"], img.dim[1] / 2, x["y"], img.dim[2] / 2))
        
        classes <-
          data.table(
            roi = gsub("roi-", "", roi.out2$fac$shape),
            edge = roi.out2$fac$edge,
            size = rois[c.roi],
            offset.x = offset.x,
            offset.y = offset.y,
            offset = offset
          )
        
        classes <-
          data.table(frame = frame, classes)
        
        if (edges == FALSE) {
          classes <- classes[edge == FALSE]
          if (nrow(classes) == 0)
            stop("no ROI found that is not on edge")
        }
        
        if (search.for == "offset") {
          z.best <- z.l[[classes[which.min(offset), ]$roi]] #best segment
          r.name <- classes[which.min(offset), ]$roi
          best.class <- classes[which.min(offset), ]
        }
        
        if (search.for == "offset.x") {
          z.best <- z.l[[classes[which.min(offset.x), ]$roi]] #best segment
          r.name <- classes[which.min(offset.x), ]$roi
          best.class <- classes[which.min(offset.x), ]
        }
        
        if (search.for == "offset.y") {
          z.best <- z.l[[classes[which.min(offset.y), ]$roi]] #best segment
          r.name <- classes[which.min(offset.y), ]$roi
          best.class <- classes[which.min(offset.y), ]
        }
        
        if (search.for == "largest") {
          z.best <- z.l[[classes[which.max(size), ]$roi]] #best segment
          r.name <- classes[which.max(size), ]$roi
          best.class <- classes[which.max(size), ]
        }
        
        
        best.cont <- data.table(EBImage::ocontour(z.best)[[1]])
        colnames(best.cont) <- c("x", "y")
        
        cont.im <- data.table(frame = frame, best.cont)
        
        
        ##using free.ml

        
        fml <- free.ml(out = as.matrix(best.cont), smooth.n = smooth.n,red=red)
        
        cont.sm <- data.table(frame = frame, fml$cont.sm)
     
        ml<- fml$ml
        
        #EBImage::display(z,"raster")
        #with(ml,points(x,y,col="blue",cex=0.1))
        
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
        ml[, roi := r.name]
        
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
          class = classes,
          dim = img.dim
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
    
   #qplot(d=cont.dat[frame==0],x,y)+geom_point(d=midline.dat2[frame==0],aes(x,y,col=n))
    #centroid in next frame
    
    keep.n <- midline.dat2[,list(max.n=round(max(n)*ant.per)),by=frame]
    midline.dat2 <- midline.dat2[keep.n, on="frame"] 
    
    head.ml <- copy(midline.dat2[n <= max.n])
      head.ml[, c("x.c","y.c"):=list(mean(x), y.c = mean(y)), by = list(frame)]
    
    ends3 <- copy(head.ml)[,keep:=n==min(n)|n==max(n),by=frame][keep==T]
    
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
    # qplot(data=midline.dat2,x=x,y=y,col=n)+facet_wrap(frame~.)
    
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
    
  
    #qplot(d=cont.dat[frame==0],x,y)+geom_point(d=midline.dat2[frame==0],aes(x,y,col=n))+geom_point(d=mid.pred2[frame==0],aes(x,mid.pred))
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


#' @title  Midline estimations for closed contours
#' @description Internal function used in \code{kin} functions for calculating a midline spanning a closed contour. 
#'
#' @param out a matrix of named x,y values describing a closed outline (contour)
#' @param smooth.n the number of smoothing iterations. See Details
#' @param red numeric, between 0-1 the proportion of contour coordinates to sample. Will speed up midline estimations. If 'NULL', the full contour in \code{out} will be used See Details.
#' 
#' 
#' @return A list with the following components:
#'
#' \code{ml}: a data table consisting of the x,y coordinates of the midline and their index value
#' 
#' \code{cont.sm}: a data table containing the smoothed outline according to \code{coo_smooth}. If 'smooth.n=0', this will be identical to the input 'out'.
#' 
#' \code{cont.sides}: a data table containing the smoothed outline factorized by side 'a' or 'b' and index.
#' 
#' 
#' @details
#' The midline is determine by first finding the tips of the contour (i.e., the two coordinates in the outline that are farthest from one another) with \code{\link{coo_truss}} and therefore assumes the contour is elongate. The function then bisects the contour across this axis using \code{\link{coo_slice}}, giving it two sides with coordinates of equal length. The midline is calculated as the midpoints defined between all pairs of coordinates with the same index value.
#' 
#' NOTE: The function does not make decisions about position, i.e, output values of 'n', although ordered along the long axis of the contour, may be differently arranged given the rotation of the contour.

#' 
#' 'smooth.n' is passed to the 'n' parameter of \code{\link{coo_smooth}}, which smooths coordinates using a simple moving average. Users should be careful not to oversmooth. If the input contour has few points (say just a 100 or so extracted from \code{kin} functions run on low resolution images), much detail will be lost. In general, \code{smooth.n} should be <5.
#' @export
#' 
#' @seealso \code{\link{coo_smooth}}, \code{\link{coo_slice}},\code{\link{coo_slide}}, \code{\link{coo_truss}}, \code{\link{kin.free}}
#' 
#' @examples
#' # a lateral midline, but a midline nonetheless
#' require(Momocs)
#' o <- Momocs::nsfishes$coo[[136]]
#' colnames(o) <- c("x","y")
#' plot(o[,1],o[,2])
#' fml <- free.ml(o,smooth.n=0)
#' points(fml$ml$x,fml$ml$y,col="red")
#' #note the difference
#' fml2 <- free.ml(o,smooth.n=10)
#' points(fml2$ml$x,fml2$ml$y,col="blue")

free.ml <- function(out = NULL,smooth.n=NULL,red=NULL) {
  
  if(!"matrix" %in% class(out)) stop("'out' must be a matrix")
  n <- side <- x <- y <- NULL
  
  if (!is.null(red) ){
    if(!is.numeric(red)) stop("'red' must be numeric and 0-1")
    if (red<0 | red>1 )
      stop("'red' must be numeric and 0-1")
  }
  
  if(!is.null(red)) red.n <- round(red*(nrow(out)/2),0)
  
  coo <- Momocs::coo_close(out)
  if(!is.null(red)) coo <- Momocs::coo_interpolate(coo,n=red.n)
  colnames(coo) <- c("x","y")
  coo <- Momocs::coo_smooth(coo,smooth.n)
  #coo <-  Momocs::coo_slide(coo, id=10)
  tr <-  Momocs::coo_truss(coo)
  tip.n <- names(tr[which.max(tr)])
  tips <-
    c(as.numeric(gsub("(\\d+)-(\\d+)", "\\1", tip.n)), as.numeric(gsub("(\\d+)-(\\d+)", "\\2", tip.n)))
  

  coo <-  Momocs::coo_slide(coo, tips[2])
  
  tr2 <-  Momocs::coo_truss(coo)
  tip.n2 <- names(tr[which.max(tr2)])
  tips2 <-
    c(as.numeric(gsub("(\\d+)-(\\d+)", "\\1", tip.n2)), as.numeric(gsub("(\\d+)-(\\d+)", "\\2", tip.n2)))
  
  
  coo.sl <-  Momocs::coo_slice(coo, ids = tips2)

  n.pts <- sapply(coo.sl,nrow)
  coo.sl <-
    lapply(coo.sl, function(x)
      Momocs::coo_sample(x,n=min(n.pts)))
  
  coo.sl <-
    lapply(coo.sl, function(x)
      Momocs::coo_smoothcurve(Opn(list(x)),n=1))
  
coo.dist <- rbind(
  data.table(coo.sl[[1]]$coo$shp1,n=1:min(n.pts),side="a"),
  data.table(coo.sl[[2]]$coo$shp1,n=min(n.pts):1,side="b")
)


colnames(coo.dist)[1:2] <- c("x","y")

#qplot(d=coo.dist,x,y,col=side)
  

  setkeyv(coo.dist,c("n"))
  
  coo.ml <- coo.dist[, list(x = sum(x) / 2, y = sum(y) / 2), by = list(n)]

  #qplot(d=coo.dist,x,y,col=n)+geom_point(d=coo.ml)
  
  return(list(ml = coo.ml,cont.sm=data.table(coo),cont.sides=coo.dist))
}


#' @title Tracking of fin-like extensions of body contours

#' @description  Estimates the amplitudes of regions along a body contour that are protruding. Useful in computing paired-fin amplitudes from contour data produced from  \link{kin.simple} and \link{kin.search}. Also computes a smoothed midline based on the body outline with the fin region removed.
#'
#' @param out a data frame or matrix with 'x' and 'y' data as columns. The first row should be a point near the head, the last near the tail
#' @param fin.pos numeric, a vector of length 2 indicating the start and end of the contour region that contains the fins of interest as a proportion of the body length.
#' @param smooth.n numeric, the number of smoothing operations undertaken by \link{coo_smooth} on the contour described by 'x'. See Details.
#' @param ml.smooth numeric (0-1), the smoothing value for the midline. See details. 
#' @param red numeric, between 0-1 the proportion of contour coordinates to sample. Will speed up fin position and midline estimations. If 'NULL', the full conour in \code{out} will be used See Details.
#' @export
#' @importFrom graphics lines
#' @importFrom stats complete.cases fitted lm loess  predict smooth.spline
#' @importFrom utils head setTxtProgressBar tail txtProgressBar
#'
#' @details
#' If \code{red} is specified, the contour \code{out} is sampled with \code{\link{coo_interpolate}} from the \code{Momocs} package. The number of points sampled (n) equals \code{red} times the number of points in \code{out}
#'
#' To establish the contour positions that are within \code{fin.pos}, a midline is estimated using \code{\link{free.ml}}. These midline points are indexed by position along the body length by calculating the cumulative distance between midline coordinates in pixels.  
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
#' \code{body} a data table consisting of x,y coordinates of the body contour
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
#' kin <- kin.free(image.dir = ti,thr=0.9,
#' ant.per = 0.25,red=0.5,smooth.n=3,ml.smooth=list("spline",0.97))
#' 
#' #fin amplitudes by frame with data.table
#' fin.pos <- c(0.25,.5)
#' fin.dat <- kin$cont.sm[, 
#' { f <- fin.kin(data.frame(x=x,y=y),fin.pos =fin.pos,smooth.n=0,red=0.75,ml.smooth=0.75);
#' list(amp=f$amp$amp2,side=f$amp$side)},by=frame]
#' 
#' p <- ggplot(dat=fin.dat,aes(x=frame,y=amp,col=side))+geom_line()+theme_classic(15)
#'print(p)
#'
#'
#' ## plot body and fin contours of frame 8
#' cont <- kin$cont.sm[frame==8,list(x,y)]
#' fins <- fin.kin(cont,fin.pos =fin.pos,red=NULL,smooth.n=0,ml.smooth=0.75)
#'
#' #plot body contour and fins
#' p <- qplot(data=fins$body,x=x,y=y)+geom_point(data=fins$fin,aes(x,y),col="red",size=3)
#' p+geom_point(data=fins$fin.pts,aes(x,y,shape=pos))+xlim(c(0,kin$dim[1]))+ylim(c(0,kin$dim[2]))
#' #plot body contour minus fins and the body midline
#' p <- qplot(data=fins$comp,x=x,y=y)+geom_point(data=fins$midline,aes(x,y),col="red",size=2)
#' p+xlim(c(0,kin$dim[1]))+ylim(c(0,kin$dim[2]))
#'
#'unlink(ti,recursive=TRUE)
#' }
#'
fin.kin <-
  function(out,fin.pos = NULL,smooth.n = 5, ml.smooth = 0.9,red=NULL) {
    y <- x <- x.sm <- y.sm <- n <- m <- b <- dist2 <- x.c <- y.c <- tip1 <- tip2 <- method <- amp2 <-pos <- y.pred <-side <- ends <- NULL # due to NSE notes in R CMD checks
   
    if (is.null(fin.pos))               
      stop("'fin.pos' not defined")
 
    if (length(fin.pos) != 2)
      stop("length of 'fin.pos' argument must be 2")
    if (!is.matrix(out))
     out<- as.matrix(out)
    if (is.null(colnames(out)))
      colnames(out) <- c("x", "y")
    
    fml <- free.ml(out=as.matrix(out),smooth.n =smooth.n,red=red)
    
    #with(fml$cont.sides,plot(x,y))
    sm.spline <- function(x, sm = 0.5) {
      s <-  round(seq(1, nrow(x), length.out = nrow(x) * (1 - sm)), 0)
      x <- as.matrix(x)
      d <- data.table(smoothr::smooth_spline(x[s, ], n = nrow(x)))
      colnames(d) <- c("x.sm", "y.sm")
      return(d)
    }
    
    #return position of a point rotated theta about another point
    point.ang.orig<- function(p,o,theta){
      xrot<-cos(theta)*(p[1]-o[1])-sin(theta)*(p[2]-o[2])+o[1]
      yrot<-sin(theta)*(p[1]-o[1])+cos(theta)*(p[2]-o[2])+o[2]
      return(c(xrot,yrot))
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
    
    ## shift points to match order in input
    
    out2 <- data.table(out,n=1:nrow(out))
    
    head <- out2[n==1]
    qplot(d=out2,x,y,col=n)
    
    
    cont.sm <- fml$cont.sm[,n:=1:.N]
    #qplot(d=cont.sm,x,y,col=n)
    
    cont.flip <- cont.sm[,head.dist:=dist.2d(x,head$x,y,head$y)][,list(close.n=n[which.min(head.dist)])]
    
  
    
    #reset n=1 to head
    cont.sm2 <- data.table(Momocs::coo_slide(as.matrix(cont.sm[,list(x,y)]),cont.flip$close.n))
    cont.sm2[,n:=1:.N]
    #qplot(d=cont.sm2,x,y,col=n)
    
    
    setkeyv(cont.sm2,c("n"))
    
    ml2 <- sm.spline(x=as.matrix(fml$ml[,list(x,y)]),sm=.8)[,n:=1:.N]
    
    qplot(d=ml2,x.sm,y.sm,col=n)
    min.dist <- with(ml2,dist.2d(x.sm,head$x,y.sm,head$y))
    
    #    #flip cont.sides  and ml2 if needed
    if(which.min(min.dist)>which.max(min.dist)){ ml2[,n:=max(n):min(n)]
      
      fml$cont.sides[,n:=max(n):min(n),by=side]
    }
    
    setkeyv(fml$cont.sides,c("side","n"))
    setkeyv(ml2,c("n"))
    
    #qplot(d=ml2,x.sm,y.sm,col=n)
    
    #qplot(d=fml$cont.sides,x,y,col=n)
    rot <- ml2[,{pa <- point.ang.orig(c(x.sm,y.sm),c(ml2$x.sm[1],ml2$y.sm[1]),pi/2);list(x=pa[1],y=pa[2])},by=n]
    rot.lm <- lm(y~x,rot)
    

    
    fml$cont.sides[,dist:=dist.2d.line(x,y,coef(rot.lm)[2],coef(rot.lm)[1]),by=list(n,side)][,bl:=dist/max(dist),by=side]
    
    fins <- fml$cont.sides[data.table::between(bl,fin.pos[1],fin.pos[2]),]
    
    #qplot(data= fins,x=x,y=y,col=side)
    
    lm.coefs <- fins[,ends:=n%in%c(min(n),max(n)),by=side][ends==TRUE,{l <- lm(y~x);list(b=coef(l)[1],m=coef(l)[2])},by=side]
    
    cents <- fins[ends==TRUE,list(x.c=mean(x),y.c=mean(y)),by=side]
    fins <- fins[lm.coefs,on="side"][cents,on="side"][,dist:=dist.2d.line(x,y,m,b),by=list(side,n)][,dist2:=dist.2d(x,x.c,y,y.c),by=side]
    
    
    fins[,tip1:=dist==dist[which.max(dist)],by=side]

    fins[tip1==T,tip1:=dist2==max(dist2)&n==max(n),by=side] #&n==max(n) for duplicated coords
    
    cp.a <- with(fins[side=='a'],features(n,dist2,smoother = "smooth.spline",spar=0.3))
    cp.an <- round(cp.a$cpts,0)
    cp.b <- with(fins[side=='b'],features(n,dist2,smoother = "smooth.spline",spar=0.3))
    cp.bn <- round(cp.b$cpts,0)
    
    
    cp.bn <- fins[side=='b'][n%in%cp.bn,][which.max(dist2),]$n
    cp.an <- fins[side=='a'][n%in%cp.an,][which.max(dist2),]$n
    
    if(length(cp.an)==0) cp.an <- fins[side=='a'][which.min(dist2)]$n
    if(length(cp.bn)==0) cp.an <- fins[side=='b'][which.min(dist2)]$n
    cp <- data.table(side=c("a","b"),n=c(cp.an,cp.bn),tip2=TRUE)
    
    
    #qplot(cont)
    #qplot(d=fml$cont.sides,x,y,col=n)+facet_grid(side~.)
    
    fins <- cp[fins,on=c("side","n")]

    
    lms <- list()
    for(i in c("a","b")){
      lms[[i]] <- lm(y~x,fins[ends==TRUE & side==i])
    }
    
    fins[, y.pred := predict(lms[[side]], newdata = data.frame(x = x)), by = list(side)]
    
    comp <- merge(fml$cont.sides, fins[,list(y.pred,n,side)], by = c("n","side"), all.x = TRUE)
    comp[!is.na(y.pred), y := y.pred]
    
  
    #p <- qplot(data= fml$cont.sides,x=x,y=y)+geom_point(dat=fins,aes(x,y.pred,col=side))
    #print(p)
  
    setkeyv(comp, c("n","side"))
    
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
    
    
  bl <- sum(copy(ml.comp.s[,dist:=dist.2d(x,dplyr::lead(x),y,dplyr::lead(y))])$dist,na.rm = TRUE)
    
    # p <- qplot(data= fml$cont.sides,x=x,y=y)+geom_point(dat=fins,aes(x,y,col=side))+geom_point(data=fins[tip1==TRUE],aes(x,y),size=3,col="black")+geom_point(d=cents,aes(x.c,y.c),col="blue")+geom_point(d= ml2,aes(x.sm,y.sm))
    # print(p)
    # 
    return(list(
      body = cont.sm2, #smoothed contour
      fin = fins[,list(side,n,x,y)],
      fin.pts = finPts,
      comp = comp[,list(n,side,x,y)],
      midline = ml.comp.s[,list(x,y)],
      bl = bl ,
      amp = amp
    ))
  }


