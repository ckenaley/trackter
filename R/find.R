#' @title  ROI identification through thresholding and segmentation
#' 
#' @description  Finds regions of interests (ROIs) within images according to search parameters by thresholding, binarization, and segmenting connected regions. Returns contour coordinates of ROI(s) and position and size descriptors. 
#' 
#' @param img character, the image file path or vector of paths. 
#' @param thr numeric or character ('otsu'); the threshold to determine binary image. See Details.
#' @param search.for character, the search parameter. See Details.
#' @param size.min numeric, indicating the minimum size of ROIs as a proportion of the pixel field to be considered in analysis. May be useful if smaller unimportant ROIs appear in the frame(s). Default is 0.01.
#' @param edges logical, should ROIs on image edges be evaluated. See Details.
#' @param border integer, if \code{edges=TRUE}, size of border to add in pixels. Otherwise ignored. Dee details.
#' @param results logical, should a binary image be printed to the graphics device. Value is forced to FALSE if the length of 'img' is >1.
#' @param bg character, the background color of the binary images (either 'white' or 'black'). Ignored if \code{results=FALSE}.
#' 
#' @export
#'
#' @details
#'Thresholding operations can be performed with an arbitrary (user defined, 0-1) numeric value or with Otsu's method ('thr="otsu"'). The latter chooses a threshold value by minimizing the combined intra-class variance. See \code{\link{otsu}}.
#'
#'
#' \code{search.for} determines how ROIs are chosen:
#' \itemize{
#' \item "offset", the ROI with a centroid that is the shortest linear distance to the center of the field
#' \item "offset.x", the ROI with a centroid x position that is closest to the x position of the center of the field
#' \item "offset.y", the ROI with a centroid y position that is closest to the y position of the center of the field
#' \item "largest", the largest ROI and the default.
#' }
#'
#' These choices will be made on ROI sets that are not on the edge of the field if 'edges=FALSE'.
#'
#' \code{edges} is set by default to 'FALSE'. It is not advisable to include shapes that are on the edge of any frame and are therefore incomplete.	Yet, if set to 'TRUE', the \code{border} adds a border of the background color to the image so that the intended ROI may be distinguished from the edge.
#' 
#' @return A named list or list of lists containing the following:
#'
#' \code{best} a data table consisting x,y coordinates of the best contour according to 'search.for'
#' \code{best.class} the name the best ROI found in the image(s)
#' 
#' \code{classes} a data table containing the following for all ROIs detected:
#' \itemize{
#' \item 'roi': the name (as letters) of each ROI found in the image
#' \item 'edge': indicating whether ROI was on the edge of the image field
#' \item 'size': area size of the ROI in pixels^2
#' \item 'offset.x': ROI distance from horizontal center
#' \item 'offset.y': ROI distance from vertical center
#' \item 'offset': euclidean distance of ROI's centroid to image center
#' }
#' 
#' \code{dim} the x and y dimensions of the image(s) analyzed. 
#'
#' @export
#'
#' 
#' @seealso \code{\link{kin.search}}, \code{\link{kin.simple}},\code{\link{kin.free}}, \code{\link{bwlabel}}
#' 
#' @examples
#'
#'library(data.table)
#' #acces example video of a sunfish swimming and establish directories
#' 
#' v <-system.file("extdata/vid", "sunfish_BCF.avi", package = "trackter")
#'
#' dir.create(paste0(tempdir(),"/images"))
#' dir.create(paste0(tempdir(),"/out"))
#' 
#' #extract images
#' vid.to.images(v, out.dir = paste0(tempdir(),"/images"))
#' 
#' #find ROIs and extract contours
#' 
#' f <- list.files(paste0(tempdir(),"/images"),full.names=TRUE)
#' rois <- find.roi(img=f,search.for="largest",size.min=0.10,thr=0.7)
#' 
#' #see results of a few
#' 
#' conts <- do.call(rbind,lapply(rois[1:10], function(x) data.table(img=x$img,x$best)))
#' 
#' library(ggplot2)
#' qplot(data=conts,x,y,col=img)
#' 
#' #clean up
#' unlink(paste0(tempdir(),"/images"),recursive=TRUE)
#' unlink(paste0(tempdir(),"/out"),recursive=TRUE)
#' 

find.roi <-  function(img=NULL,thr="otsu",search.for="largest",size.min=0.01,bg="white",edges=FALSE,border = 5,results=FALSE){
  
  size <- NULL
   
  if (length(img) > 1) results <-  FALSE

  if(!all(file.exists(img))) stop("at least one file in 'img' does not exist")
  
    if (!search.for %in% c("offset.x", "offset.y", "largest", "offset"))
      stop("'search.for' must be set to 'offset', 'offset.x','offset.y', or 'largest')")
    
    if (size.min <= 0 |
        size.min > 1)
      stop("'size.min' must be >0 and <=1.")
    
    r <- lapply(img,function(x) {
      im <- EBImage::readImage(x, all = FALSE)
      
      img.dim <- dim(im)[1:2]
      img.n <- basename(x)
      
      EBImage::colorMode(im) = EBImage::Grayscale
      
      # computes binary mask
      if (thr != "otsu" &
          !is.numeric(thr))
        stop("'thr' must be set to 'otsu' or a numeric value=0-1")
      if (thr == "otsu") {
        thr <- EBImage::otsu(im)[1]
      }
      
      
      x <-  im[, , 1]
      
      #display(x,"raster")
      
      ## computes binary mask
      y <- x < thr #contrast threshold
      y <-  EBImage::closing(y, EBImage::makeBrush(1, shape = "disc"))
      #display(y, title='Cell nuclei binary mask')
      
      ## bwlabel
      z <- EBImage::bwlabel(y)
      # display(normalize(z), title='Cell nuclei')
      # nbnuclei = apply(z, 1, max)
      # cat('Number of nuclei=', paste(nbnuclei, collapse=','),'\n')
      
      # ## paint nuclei in color
      # cols = c('black', sample(rainbow(max(z))))
      # zrainbow = Image(cols[1+z], dim=dim(z))
      # display(zrainbow, title='Cell nuclei (recolored)')
      
      
      #EBImage::display(z,"raster")
      if (edges) {
        bord <- border * 2
        
        z1 <-
          EBImage::resize(
            z,
            w = dim(z)[1],
            h = dim(z)[2],
            output.dim = c(dim(z)[1] + bord, dim(z)[2] + bord),
            bg.col = "black"
          )
        
        z <-
          EBImage::translate(z1, v = c(border, border), bg.col = "black")
        
      }
      
      rois <- tabulate(z)
      
      pix <- dim(z)[1] * dim(z)[2]
      w <- dim(z)[1] #width of image
      h <- dim(z)[2] #height of image
      per <-
        rois / (w * h) #how big are rois compared to pixel field
      
      c.roi <-
        which(per >= size.min) #candidate rois, filtered by size of of pixel field
      
      if(length(c.roi)==0) stop("no ROIs larger than 'size.min'")
      
      names(c.roi) <-
        as.factor(letters[order(rois[c.roi], decreasing = TRUE)])
      
      z.l <- list()
      out.l <- list()
      
      for (r in c.roi) {
        r.name <- as.character(names(c.roi)[c.roi == r])
        z.r <- z
        z.r[z != r] <- 0
        z.r[z == r] <- 1
        z.m <- z.r
        z.m[1, 1] <- 0 #this gets a 1 when
        z.l[[r.name]] <- z.m
        
        z.c <- EBImage::ocontour(z.m)
        
        
        #EBImage::display(z.m,"raster")
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
      
      if (edges == FALSE) {
        classes <- classes[edge == FALSE]
        if (nrow(classes) == 0)
          stop("no ROI found that is not on edge")
      }
      
      if (search.for == "offset") {
        z.best <- z.l[[classes[which.min(offset),]$roi]] #best segment
        r.name <- classes[which.min(offset),]$roi
        best.class <- classes[which.min(offset),]
      }
      
      if (search.for == "offset.x") {
        z.best <- z.l[[classes[which.min(offset.x),]$roi]] #best segment
        r.name <- classes[which.min(offset.x),]$roi
        best.class <- classes[which.min(offset.x),]
      }
      
      if (search.for == "offset.y") {
        z.best <- z.l[[classes[which.min(offset.y),]$roi]] #best segment
        r.name <- classes[which.min(offset.y),]$roi
        best.class <- classes[which.min(offset.y),]
      }
      
      if (search.for == "largest") {
        z.best <- z.l[[classes[which.max(size),]$roi]] #best segment
        r.name <- classes[which.max(size),]$roi
        best.class <- classes[which.max(size),]
      }
      
    
      best.cont <- data.table(EBImage::ocontour(z.best)[[1]])
      colnames(best.cont) <- c("x", "y")
      
      
    
      if (bg == "white") {
        #flip binary
        z.best[z.best == 0] <- 5
        z.best[z.best == 1] <- 0
        z.best[z.best == 5] <- 1
      }
      
      if (results & length(img) < 2)
        EBImage::display(z.best, "raster")
      
      list(
        img = img.n,
        best = best.cont,
        best.class=r.name,
        classes = classes,
        dim = img.dim
      )
    }
    )
    
  
    return(r)
  }



