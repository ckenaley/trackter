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
#'Thresholding operations can be performed with an arbitrary (user defined) numeric value or w;
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


#' @title  Midline tracking of free-moving ROI over image sequences

#' @description  Automatically retrieves the midline of a detected ROI that is free to move in the spatial field of an image sequence. Does so through thresholding and segmentation; finds the y-value midpoint along the between two x values from each half of the ROI sharing the same index then calculates a smoothed midline according to a chosen smoothing method (loess or spline). Also outputs the midline amplitude relative to a reference line determined by an anterior section of the ROI. Supported image formats are jpeg, png, and tiff. Support parallel processing of frames.
#' 
#'
#'
#' @param image.dir character, directory containing images to analyze.
#' @param frames numeric, vector indicating which images to process.
#' @param par logical, should the frames be processed in parallel using \code{cores.n}.
#' @param cores.n numeric, the number of CPU cores to use if \code{par=TRUE}. If \code{cores.n=NULL} (the default), the total number of cores minus 1 are used.
#' @param thr numeric or character ('otsu') threshold to determine binary image. See Details.
#' @param ant.per numeric; left-most percentage of ROI that establishes the horizontal reference for the midline displacement.
#' @param tips, numeric, the proportion the the midline data to use in calculation of the head and tail position.
#' @param smooth.what, character of length 1 or 2, either "ml", "cont", or both.
#' @param smooth.n, numeric, the number of contour smoothing interations. See Details.
#' @param smoothing character, the midline smoothing method, either 'loess' or "spline".See Details.
#' @param smooth numeric; if \code{smoothing} is set to 'loess', smoothing parameter value for plotted midline. If \code{smooth='spline'} then this value must be less than 1. See Details.
#' @param smooth.points numeric, number of equally spaced points along the ROI midline on which the smoothed midline is computed.
#' @param flip logical, indicating if binary should be flipped.
#' @param size.min numeric, indicating the minimum size of ROIs as a proportion of the pixel field to be considered in analysis. May be useful if smaller unimportant ROIs appear in the frame. Default is 0.02.
#' @param save logical, value indicating if images should be outputted with midline and predicted midline based on the \code{ant.per} \code{lm} overlaying original or binary images.
#' @param out.qual, numeric, a value between 0-1 representing the quality of outputted images. Ignored if \code{save=FALSE}.
#' @param out.dir character, the directory to which outputted images should be saved.
#' @param plot.pml logical, value indicating if outputted images should include an overlay of the midline, head region and theoretical midline based on \code{ant.per}.
#' @param search.for character, the search parameter. See Details.
#' @param edges logical, should ROIs on image edges be evaluated. See Details.
#' @param border if \code{edges=TRUE}, size of border to add in pixels. Dee details.
#'
#' @export
#'
#' @details
#' The midline is calculated as the midpoints between the extrema of the coordinates describing the each ROIs contour (i.e, the coordinates spanning the longest euclidean distance). Chooses ROIs based on relative ROI size or position.
#'
#' The midline is determine by first finding the tips of the ROI (i.e., the two coordinates in the outline that are farthest from one another) with \code{\link{free.ml}}and therefore assumes the ROI is elongate and moving along this long axis. Using \code{\link{free.ml}}, the function bisects the ROI contour at the tips, giving it two sides with coordinates of equal length. The midline is calculated as the midpoints defined between all pairs of coordinates with the same index value.
#'
#'Thresholding operations can be performed with an arbitrary (user defined) numeric value or w;
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
#' For midline smoothing, if \code{smoothing='spline'} (the default), \code{\link{smooth_spline}} from the \code{smoothr} package is used to interpolate points between a reduced number of vertices using piecewise cubic polynomials. The number of vertices is calculated based on the number of midline coordinatate times \code{1-smooth}. If \code{smoothing='loess'}, \code{stats::loess} is used to fit a polynomial surface. For contours that have a complicated midline with non-unique x values, loess smoothing can produce poor results. Thus, spline smoothing is usually the advisable option. 
#' 
#' For contour smoothing 'smooth.n' is passed to the 'n' parameter of \code{\link{free.ml}}, which smooths coordinates using a simple moving average. Users should be wary of oversmoothing by smoothing both the contour (from which the midline is calculated) and the midline.
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
#' \item 'x' and 'y': x and y positions of the midline of the ROI
#' \item 'roi': a character indicating ROI size ('a' being the largest)
#' \item 'y.pred': the smooth predicted midline
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
#' \item 'x' and 'y.m': x and y positions of the midline of the ROI
#' \item 'x.sm' and 'y.sm': the smoothed midline positions predicted by \code{smoothing}
#' \item 'y.pred': midline points fit to a smooth spline or loess model with spar or span equal to \code{smooth} (red curve in the outputted images if \code{save=TRUE})
#' \item 'wave.y': y position of midline points (see below)
#' \item 'per.bl': the percentage of 'x.sm' along the body length calculated as the cumulative sum of distances between points
#' }
#' 
#' \code{mid.pred} the theoretical midline based on a linear model established by the anterior section of \code{midline$x.sm} in \code{ant.per}. Used to calculate \code{midline$wave.y} as the orthogonal distance between the line defined by 'x' and 'mid.pred' and each coordinate defined by '\code{midline$x.sm} and \code{midline$y.sm}. A data table that contains the following:
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
#' 
#' 
#' vid.to.images(paste0(tempdir(),"/ropefish.avi"), out.dir = paste0(tempdir(),"/images"))
#' 
#' kin <- kin.free(image.dir =paste0(tempdir(),"/images"),
#'       par=FALSE,
#'       save=FALSE,
#'       smooth.what=c("ml","cont"),
#'       thr = "otsu",
#'       smoothing="spline",
#'       smooth=0.92,
#'       size.min=0.01,
#'       frame=1:2
#'       )
#'
#' #plot instantaneous amplitude of tail (last/rightmost point) over frames
#' p <- ggplot(dat=kin$kin.dat,aes(x=frame,y=amp))+geom_line()+geom_point()+theme_classic(15)
#' print(p)
#'
#' # midline plot
#' ml <- kin$midline
#' p <- ggplot(dat=ml,aes(x=x.sm,y=y.sm,col=wave.y))+theme_classic(15)
#' p <- p+geom_point(aes(group=frame), size = 1.5)+facet_wrap(frame~.)
#' print(p)
#'
#' unlink(paste0(tempdir(),"/images"),recursive=TRUE)
#'}
#'

kin.free <-
  function(image.dir = NULL,frames=NULL,par=FALSE,cores.n=NULL,thr = "otsu",ant.per = 0.10,tips = 0.02,smooth.what="ml",smooth.n=5,smoothing = "spline",smooth = 0.25, smooth.points = 200,save = TRUE,out.qual = 1,out.dir = NULL,plot.pml = TRUE,flip = TRUE,size.min = 0.02,search.for = "largest",edges = FALSE,border = 5) {
    
    
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
      y.sm <-
      keep <-  l <-  head.pval <-  bl <-  per.bl <-  amp <- y <- im <- or.dist <- NULL
    
    
    
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
    
    if (!all(smooth.what %in% c("ml", "cont")))
      stop("'smooth.what' must be set to 'ml', 'cont', or 'c('ml','cont')'")
    
    
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
        if ("cont" %in% smooth.what)
          n <- smooth.n
        else
          n <- 0
        fml <- free.ml(out = as.matrix(best.cont), smooth.n = n)
        
        cont.sm <- data.table(frame = frame, fml$cont.sm)
        y.df <- fml$ml
        #EBImage::display(z,"raster")
        #with(y.df,points(x,y,col="blue",cex=0.1))
        
        ends <- ceiling(nrow(y.df) * tips)
        tip.y <-
          mean(tail(y.df$y[!is.na(y.df$y)], ends))#tip is mean y.m of last tips per pixels
        tip.x <-
          mean(tail(y.df$x[!is.na(y.df$y)], ends))
        
        head.y <-
          mean(head(y.df$y[!is.na(y.df$y)], ends))
        head.x <-
          mean(head(y.df$x[!is.na(y.df$y)], ends))
        
        #n midline points
        if (is.null(smooth.points))
          smooth.points <- nrow(y.df)
        midline <-
          y.df[seq(1, nrow(y.df), length.out = smooth.points), ] #two hundred points on midline
        
        midline <- midline[complete.cases(midline)]
        midline <- data.table(frame, midline)
        
        midline <- midline[complete.cases(midline), ]
        
        midline[, roi := r.name]
        
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
          mid = midline,
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
    if (length(flip.f) != 0)
      midline.dat2[frame %in% (flip.f - 1), n := max(n):1, by = list(frame)]
    
    #now flip if head is the wrong direction
    
    setkeyv(midline.dat2, c("frame", "n"))
    
    #centroid in next frame
    cent <-
      midline.dat2[n %in% 1:(ant.per * smooth.points), ][n %in% c(min(n), max(n))][, list(x.c =
                                                                                            mean(x), y.c = mean(y)), by = frame][, c("x.c", "y.c") := list(dplyr::lead(x.c), dplyr::lead(y.c))]
    cent.ref <-
      midline.dat2[n %in% 1:(ant.per * smooth.points), ][n %in% c(min(n), max(n)) &
                                                           frame == 0][, list(x.c = mean(x), y.c = mean(y)), by = frame][frame == 0, ]

    #head centroid in first frame  
    cent.ref <-
      midline.dat2[cent.ref, on = "frame"][n %in% 1:(ant.per * smooth.points), ][n %in%
                                                                                   c(min(n), max(n)) &
                                                                                   frame == 0][, or.dist := dist.2d(x.c, x, y.c, y), by = frame]
    ends3 <-
      midline.dat2[cent, on = "frame"][n %in% 1:(ant.per * smooth.points), ][n %in%
                                                                               c(min(n), max(n))]
    ends3[, dist := dist.2d(x.c, x, y.c, y), by = frame]
    ends3 <-
      ends3[cent.ref[, list(frame, n, or.dist)], on = c("frame", "n")][, diff :=
                                                                         dist - or.dist]
    direction <-
      ends3[, list(direction = ifelse(n[which.min(dist)] == max(n), "tailward", "headward"))]
    
    
    if (direction[1, ]$direction == "tailward")
      midline.dat2[, n := max(n):1, by = list(frame)]
    # qplot(data=midline.dat2,x=x,y=y,col=n)+facet_wrap(frame~.)
    
    #rerun kin, head, and ml calculations with new data
    
    ####which type of lines to be fitted, spline or loess
    if ("ml" %in% smooth.what &
        !any(c("spline", "loess") == smoothing))
      stop("'smooth.cont==F' and smoothing' must = 'loess' or 'spline'")
    
    if ("ml" %in% smooth.what & smoothing == "loess")
      midline.dat2[, y.pred :=
                     fitted(loess(y ~ x,
                                  span = smooth,
                                  degree = 0)), by = list(frame)]
    
    
    if ("ml" %in% smooth.what & smoothing == "spline") {
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
    
    if (!"ml" %in% smooth.what)
      midline.dat2[, c("x.sm", "y.sm") := list(x, y)]
    
    
    setkeyv(midline.dat2, c("frame", "n"))
    
    head.dat <-
      copy(midline.dat2[n %in% 1:(ant.per * smooth.points)])
    
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
      midline.dat2[, list(x.sm, frame)][, list(
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
          with(midline.dat2[frame == (frame2 - 1)][1:ceiling(ant.per * smooth.points),],
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
#' @param smooth.n the number of smoothing iterations. See Detail
#' 
#' 
#' @return A list with the following components:
#'
#' \code{ml}: a data table consisting of the x,y coordinates of the midline and their index value
#' 
#' \code{cont.sm}: a data table containing the smoothed outline according to \code{coo_smooth}. If 'smooth.n=0', this will be identical to the input 'out'.
#' 
#'
#' @details
#' The midline is determine by first finding the tips of the contour (i.e., the two coordinates in the outline that are farthest from one another) with \code{\link{coo_truss}} and therefore assumes the contour is elongate. The function then bisects the contour across this axis using \code{\link{coo_slice}}, giving it two sides with coordinates of equal length. The midline is calculated as the midpoints defined between all pairs of coordinates with the same index value.
#' 
#' 'smooth.n' is passed to the 'n' parameter of \code{\link{coo_smooth}}, which smooths coordinates using a simple moving average.
#' @export
#' 
#' @seealso \code{\link{coo_smooth}}, \code{\link{coo_slice}},\code{\link{coo_slide}}, \code{\link{coo_truss}}, \code{\link{kin.free}}
#' 
#' @examples
#' require(Momocs)
#' o <- Momocs::nsfishes$coo[[136]]
#' colnames(o) <- c("x","y")
#' plot(o[,1],o[,2])
#' fml <- free.ml(o,smooth.n=0)
#' points(fml$ml$x,fml$ml$y,col="red")
#' #note the difference
#' fml2 <- free.ml(o,smooth.n=10)
#' points(fml2$ml$x,fml2$ml$y,col="blue")

free.ml <- function(out = NULL,smooth.n=NULL) {
  
  if(!"matrix" %in% class(out)) stop("'out' must be a matrix")
  n <- side <- x <- y <- NULL
  
  coo <- Momocs::coo_close(out)
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
  coo.sl <-
    lapply(coo.sl, function(x)
      Momocs::coo_sample(x, n = min(sapply(coo.sl, nrow))))
  coo.sl[[2]] <- coo.sl[[2]][nrow(coo.sl[[2]]):1, ]
  coo.dist <-
    data.table(do.call(rbind, coo.sl))[, n := rep(min(sapply(coo.sl, nrow)):1, 2)][,side:=c(rep("a",min(sapply(coo.sl, nrow))),rep("b",min(sapply(coo.sl, nrow))))]
  
  setkeyv(coo.dist,c("n"))
  
  coo.ml <- coo.dist[, list(x = sum(x) / 2, y = sum(y) / 2), by = list(n)]

  
  return(list(ml = coo.ml,cont.sm=data.table(coo)))
}


#' @title Tracking of fin-like extensions of body contours

#' @description  Estimates the amplitudes of regions along a body contour that are protruding. Useful in computing paired-fin amplitudes from contour data produced from  \link{kin.simple} and \link{kin.search}. Also computes a smoothed midline based on the body outline with the fin region removed.
#'
#' @param x a data frame or matrix with 'x' and 'y' data as columns.
#' @param fin.pos numeric, a vector of length 2 indicating the start and end of the contour region that contains the fins of interest as a proportion of the body length.
#' @param smooth.n numeric, the number of smoothing operations undertaken by \link{coo_smooth} on the contour described by 'x'.
#' @param tip.ang the minimum angle, in degrees, that defines tip of each fin. See Details.
#' @param x.bins numeric, when less than or equal to 1, the proportion of contour coordinates to sample for midline estimation. If greater than 1, the absolute number of equally spaced x values from which to compute the midline. See Details.
#' @param smoothing character, the midline smoothing method, either 'loess' or 'spline'.
#' @param ml.smooth numeric the smoothing value for the midline. If \code{smoothing} is set to 'loess', passed to 'span' value for \code{\link{loess}}. If \code{smoothing} is set to 'spline', passed to 'spar' value for \code{\link{smooth.spline}}
#'
#' @export
#' @importFrom graphics lines
#' @importFrom stats complete.cases fitted lm loess  predict smooth.spline
#' @importFrom utils head setTxtProgressBar tail txtProgressBar
#'
#' @details
#'The algorithm assumes a left-right orientation, i.e., the head of the contour is left. If otherwise oriented, contour can be flipped with \code{\link{coo_flipx}} and \code{\link{coo_flipy}} after converting contour to class \code{coo}.
#'
#'  \code{tip.angle} is used to define the tip of the fin, assuming that the tip of the fin is pointed and, for a sufficiently smoothed fin contour, will have contour edges that form the highest angles within the fin region defined by \code{fin.pos}. Low values of \code{smooth.n} (<5) should be avoided if the contour is jagged, perhaps due to digitization.
#'
#'In addition to fin amplitude and contour extraction, this function also produces a composite contour of the body minus the fin area described by \code{fin.pos}. Fin contours are replaced by a simple linear prediction constructed from the coordinates of the first and last values covered by \code{fin.pos}. The result is a straight line between the start and end of each fin. From this composite body contour, a midline prediction is made based on the method indicated by \code{smoothing} and number of points indicated by \code{x.bins}.
#'
#'  \code{x.bins} controls the bin size of x values used to estimate the midline. From these bins, mean x and the range of y is calculated. The midpoint at each mean x is then calculated from the mid point of y. When less then 1, \code{x.bins} values approaching 1 may result in poor a midline as x values on one side of the contour may not have corresponding identical values on the other. Values closer to 0 will result in fewer points but a more robust midline. Higher \code{smooth.n} values will also result in a more robust midline estimation (but also a loss of contour information).
#'
#'
#' @return A list with the following components:
#'
#' \code{body} a data table consisting of x,y coordinates of the body contour
#'
#' \code{fin} a data table describing the contour of the fins consisting of the following:
#'
#' \itemize{
#' \item x,y coordinates within the range of \code{fin.pos}
#'
#' \item 'ang': the angle formed by each coordinate and its adjacent points.
#' \item 'fin': fin side, 'L' or 'R'
#' \item 'y.pred': predicted y values according to \code{lm()} from start to end of fins.
#' }
#'
#' \code{fin.pts} a data table describing fin position consisting of the following:
#'
#' \itemize{
#' \item  x,y coordinates of the fin tips, start, and end within the range of \code{fin.pos}.
#' \item 'ang': the angle formed by the coordinates and their adjacent points.
#' \item 'pos': description  of the coordinates' positions, 'start', 'end' or 'tip'.
#' }
#'
#' \code{comp} a data table describing the composite contour of the body minus the fins.
#' \itemize{
#' \item  x,y coordinates of the body except the range of x values within \code{fin.pos}. These values take on a straight line described by the prediction of \code{lm()} based on the start and end of the fin. See Details.
#' }
#'
#' \code{midline} a data table describing the estimated
#' \itemize{
#' \item 'x': the mean x position within the bin.
#' \item 'x.bin': the x bin in \code{\link{cut}} notation.
#' \item 'y.m': the y midpoint at the bind and mean x value.
#' \item 'ml.pred': the y midline value according to the smoothing parameters.
#' }
#'
#' @seealso \code{\link{kin.simple}}, \code{\link{kin.LDA}}, \code{\link{kin.search}}, \code{\link{efourier}}, \code{\link{coo_angle_edges}}, \code{\link{coo_smooth}}, \code{\link{loess}}, \code{\link{smooth.spline}}
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
#' #extract images with ffmpeg operations and reduce them to 600 px wide with a filter
#' filt.red <- " -vf scale=600:-1 " #filter
#' vid.to.images2(vid.path="sunfish.avi",filt = filt.red) #extract
#'
#' #number of frames
#' fr <- length(list.files("images"))
#' #extract contours and other data
#' kin <- kin.simple(image.dir = "images",frames=c(1:fr),thr=0.9,ant.per = 0.25)
#' #fin amplitudes by frame with data.table
#' fin.pos <- c(0.25,.5)
#' fin.dat <- kin$cont[, { f <- fin.kin(data.frame(x=x,y=y),fin.pos =fin.pos);
#' list(amp=f$amp$amp,fin=f$amp$fin,amp.bl=f$amp$amp.bl)},by=list(frame)]
#' p <- ggplot(dat=fin.dat,aes(x=frame,y=amp,col=fin))+geom_line()+theme_classic(15)
#'print(p)
#'
#'
#' ## plot body and fin contours of frame 1
#' cont <- data.frame(x=kin$cont[frame==2,list(x,y)]$x,y=kin$cont[frame==2,list(y)]$y)
#' fins <- fin.kin(cont,fin.pos =fin.pos,x.bins=100)
#'
#' #plot body contour and fins
#' p <- qplot(data=fins$body,x=x,y=y)+geom_point(data=fins$fin,aes(x,y),col="red",size=3)
#' p+geom_point(data=fins$fin.pts,aes(x,y,shape=pos))+xlim(c(0,kin$dim[1]))+ylim(c(0,kin$dim[2]))
#'
#' #plot body contour minus fins and the body midline
#' p <- qplot(data=fins$comp,x=x,y=y)+geom_point(data=fins$midline,aes(x,ml.pred),col="red",size=2)
#' p+xlim(c(0,kin$dim[1]))+ylim(c(0,kin$dim[2]))
#'
#' }
#'
#'
#'

fin.kin <-
  function(x,
           fin.pos = NULL,
           smooth.n = 50,
           tip.ang = 10,
           smoothing = "loess",
           x.bins = 0.2,
           ml.smooth = 0.25) {
    y <-
      ang <-
      pos <-
      y.pred <-
      fin <-
      y.m <- ml.pred <- x.bin <- NULL # due to NSE notes in R CMD check
    
    if (is.null(fin.pos))
      stop("'fin.pos' not defined")
    if (length(fin.pos) != 2)
      stop("length of 'fin.pos' argument must be 2")
    if (!is.matrix(x))
      x.m <- as.matrix(x)
    if (is.null(colnames(x.m)))
      colnames(x.m) <- c("x", "y")
    #
    
    x.m[, 2] <-
      x.m[, 2] - max(x.m[, 2])#to get positive y coords after flip
    
    #plot(x.m[,1],x.m[,2])
    x.o <- coo_flipx(Out(list(x.m)))
    
    #panel(x.o)
    
    #with(x.o,plot(x,y))
    x.s <- data.table(coo_smooth(x.o, n = smooth.n)[[1]][[1]])
    colnames(x.s) <- c("x", "y")
    

    bl <- diff(range(x$x))
    
    fin.range <- min(x) + fin.pos * bl
    
    #Left fin
    finL <- x.s[x >= fin.range[1] &
                  x <= fin.range[2] & y < y[which.min(x)]]
    finL.o <- Out(list(as.matrix(finL[, list(x, y)])))
    finL[, ang := (coo_angle_edges(finL.o)[[1]])]
    
    finL2 <- finL[seq(1, nrow(finL), 5)]
    finL2.o <- Out(list(as.matrix(finL2[, list(x, y)])))
    finL2 <- data.table(finL2.o[[1]][[1]])
    finL2[, ang := deg(coo_angle_edges(finL2.o)[[1]])]
    
    ptsL <- finL2[order(abs(ang)), ][1:15]
    ptsL[which.min(x), pos := "start"]
    ptsL[which.max(x), pos := "end"]
    ptsL[ang < tip.ang, ][which.min(y)]$pos <- "tip"
    ptsL <- ptsL[!is.na(pos)]
    
    #right fin
    
    finR <- x.s[x >= fin.range[1] &
                  x <= fin.range[2] & y > y[which.min(x)]]
    finR.o <- Out(list(as.matrix(finR[, list(x, y)])))
    finR[, ang := (coo_angle_edges(finR.o)[[1]])]
    
    #panel(x.o)
    finR2 <- finR[seq(1, nrow(finR), 5)]
    finR2.o <- Out(list(as.matrix(finR2[, list(x, y)])))
    finR2 <- data.table(finR2.o[[1]][[1]])
    finR2[, ang := deg(coo_angle_edges(finR2.o)[[1]])]
    
    ptsR <- finR2[order(abs(ang)), ][1:15]
    ptsR[which.min(x), pos := "start"]
    ptsR[which.max(x), pos := "end"]
    ptsR[ang < tip.ang, ][which.max(y)]$pos <- "tip"
    ptsR <- ptsR[!is.na(pos)]
    
    finPts <- rbind(data.table(ptsL, fin = "L"), data.table(ptsR, fin = "R"))
    
    fins <- rbind(data.table(finL, fin = "L"), data.table(finR, fin = "R"))
    
    fins[, y.pred := predict(lm(y ~ x, data.frame(
      y = c(y[which.min(x)], y[which.max(x)]), x = c(min(x), max(x))
    )), newdata = data.frame(x = x)), by = list(fin)]
    
    comp <- comp2 <- merge(x.s, fins, by = c("x", "y"), all.x = TRUE)
    setkey(comp, "x")
    comp[!is.na(y.pred), y := y.pred]
    comp <- comp[, list(x, y)]
    
    setkey(comp2, "x")
    
    if (x.bins > 1)
      comp2[, x.bin := cut(x, x.bins)]
    if (x.bins <= 1)
      comp2[, x.bin := cut(x, nrow(comp) * x.bins)]
    
    comp2 <- comp2[, list(x = mean(x), y = range(y)), by = list(x.bin)]
    
    mid <- comp2[, list(y.m = mean(y)), by = list(x, x.bin)]
    
    if (smoothing == "loess")
      mid[, ml.pred := predict(loess(y.m ~ x, span = ml.smooth))]
    if (smoothing == "spline")
      mid[, ml.pred := smooth.spline(x = x, y = y.m, spar = ml.smooth)$y]
    
    bl <-
      mid[, list(sum(dist.2d(
        x, dplyr::lead(x), y.m, dplyr::lead(y.m)
      ), na.rm = TRUE))]$V1
    
    amp <-
      finPts[, list(amp = abs(y[pos == "start"] - y[pos == "tip"]), amp.bl =
                      abs(y[pos == "start"] - y[pos == "tip"]) / bl), by = list(fin)]
    
    return(list(
      body = x.s,
      fin = fins,
      fin.pts = finPts,
      comp = comp,
      midline = mid,
      bl = bl,
      amp = amp
    ))
  }


