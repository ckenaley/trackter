#' @title  Midline detection from video data
#'
#' @description A wrapper, decomposing videos with \code{vid.to.images2} to frames for automatical retrieval of ROI midlines using \code{kin.img2}.
#'
#' @param vid.args list; arguments to be passed to \code{vid.to.images2}
#' @param kin.args list; arguments to be passed to \code{kin.img2}
#'
#'
#' @details By default, images are outputted to an 'images' subdirectory in the working directory and processed images to a 'processed_images' subdirectory.
#
#' @return Returns the components of \code{kin.img2}
#'
#' @seealso \code{\link{kin.img2}}, \code{\link{vid.to.images2}}
#' @export
#' @import data.table
#' @examples produce a classic midline waveform plot of swimming fish
#'\dontrun{
#' require(wesanderson)
#' require(ggplot2)
#' require(dplyr)
#'
#' #download an example video (7.5 MB) and place in working directory
#' f <- "https://github.com/ckenaley/exampledata/blob/master/trout1_63_test.avi?raw=true"
#' download.file(f, file.path(getwd(), "trout1_63_test.avi"), method = "libcurl")
#' kin <- kin.vid2(vid.args=list(vid.path ="trout1_63_test.avi",silent=TRUE),kin.args=list(thr=0.7,frames=1:10,frame.rate=10,sequenced=TRUE,rem.file=F))
#'
#' ml <- kin$midline
#' #normalize x (y is normalized to midline by "kin.img/kin.vid")
#' ml <- ddply(ml,.(frame),transform,x2=x-x[1])
#'
#' #compute instantaneous amplitude of tail (last/rightmost point) and wave crest x position  by frame
#' ml2 <- ddply(ml,.(frame),summarize,amp.i=last(wave.y))
#'
#' ml <- merge(ml,ml2,by="frame") #merge these
#'
#' pal <- wes_palette("Zissou1", 100, type = "continuous") #"Zissou" color palette
#' p <- ggplot(dat=ml,aes(x=x2,y=wave.y))+theme_classic(15)+scale_color_gradientn(colours = pal)
#' p+geom_line(aes(group=frame,color=amp.i),stat="smooth",method = "loess", size = 1.5,alpha = 0.5)
#'}


kin.vid2 <-function(kin.args,vid.args){

  if(!file.exists(vid.args$vid.path)) stop(paste0(vid.args$vid.path," does not exist in ", "getwd()"))

  oldwd <- getwd()
  do.call(vid.to.images2,vid.args)
  kin.args$image.dir <- paste0(getwd(),"/images")
  do.call(kin.img2,kin.args)
}

######### kin.img2

#' @title  Midline tracking over image seqences

#' @description  Automatically retrieves the midline of a detected ROI in each image of a sequence; finds the y-value midpoint along the x-value array of the ROI and fits a midline according to a chosen smoothing method (loess or spline). Also outputs the midline amplitude relative to a reference line determined by an anterior section of the ROI. Supported image formats are jpeg, png, and tiff.
#'
#' @param image.dir Directory containing images to analyze.
#' @param frames numeric, vector indicating which images to process.
#' @param sequenced logical, do the names of the image sequence end in a number sequence. If set to 'TRUE', output video will be named the file base name of the first image minus the image number ('gsub("(.+)\\_\\d*$", "\\1", image[1])')
#' @param thr numeric, threshold to determine binary image. May require some tweaking through iteration.
#' @param plot.midline logical, value indicating if outputted images should include plotted midline and reference line based on anterior section of the ROI.
#' @param smoothing character, the midline smoothing method, either 'loess' or "spline".
#' @param smooth numeric; if \code{smoothing} is set to 'loess', smoothing parameter value for plotted midline.
#' @param smooth.points numeric, number of equally spaced points along the ROI midline on which the smoothed midline is computed.
#' @param image.type character; the type of image to be outputted.
#' @param flip logical, indicating if binary should be flipped.
#' @param show.prog logical value indicating if outputted image should be displayed during analysis.
#' @param size.min numeric, indicating the minimum size of ROIs as a proportion of the pixel field to be considered in analysis. May be useful if smaller unimportant ROIs appear in the frame. Default is 0.02.
#' @param n.blob numeric, indicating which nth largest ROI is the ROI to be analyzed. May require tweaking through iteration. Perhaps best to let the function choose by using "search.for".
#' @param make.video logical value indicating if a video should be saved of midline position overlaying origina frames.
#' @param video.name character, the name given to the outputted video.
#' @param qual numeric; quality of the outputted video from 1-100\%. Defaults to 50\%.
#' @param ant.per numeric; left-most percentage of ROI that establishes the vertical reference for the midline displacement.
#' @param frame.rate numeric; outputted video frame rate in fps.
#' @param rem.file logical value indicating if the outputted images, both from the original video and images with midline overlay, should be deleted.
#' @param search.for character; the search parameter: "offset" the ROI closest to the midpoint of the field or "largest" the largest roi (equivalent to n.blob=1). Will produce a warning if the ROI indicated by either of these settings is one that has no positive pixels along 90\% of its midline (see details).
#' @param burn numeric, how many frames (a burn-in) should be used to determine the ROI. If all ROIs in the the burn-in frames are the same, the algorithm will choose that ROI in all subsequent frames. This should increase the speed of midline detection in post burn-in frames.
#'  @param silent logical, should output of ffmpeg system calls be hidden.
#'
#' @export
#'
#' @details
#'The algorithm assumes a left-right orientation, i.e., the head of the ROI is positioned left, the tail right. The \code{ant.per} value therefor establishes the reference line (theoretical straight midline) based on that portion of the head.  By default, images are outputted to the \code{image.dir} subdirectory in the working directory. Chooses ROIs based on relative ROI size or position.
#'
#'\code{image.type} Can be set as "orig" or "bin". "orig" plots midline and reference lines over the original video frames, "bin" over binary images.
#'\code{n.blob} May be useful if there are other highly contrasted ROIs in the frame and the user expects and knows their relative size
#'\code{search.for} The algorithm attempts to resolve an ROI that has positive (i.e., dark) pixels along more than 90% of the ROI's midline. This should be useful if the ROI of interest is surrounded by irregular dark object (e.g., walls).
#'
#'\code{make.video} If "TRUE" a video of the same names as \code{video.name}, if given, is created by \code{images.to.video2} and outputted in the working directory. If \code{video.name} is "NULL", a base name minus the image sequence number will be determined and "_kin" appended to the file name.
#'
#' \code{rem.file} If "TRUE" and \code{make.video} is also "TRUE", a video of processed images is still produced.
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
#' \item 'mid.pred': the predicted linear midline based on the points/pixels defined by \code{head.per} (green points in the outputted images/video)
#' \item 'y.pred': midline points fit to a smooth spline or loess model with spar or span equal to \code{smooth} (red curve in the outputted images/video)
#' \item 'wave.y': midline points 'y.pred' normalized to 'mid.pred'
#' \item 'roi': a character indicating ROI size (a being the largest)
#' \item 'cent.x': x centroid of ROI
#' \item 'cent.y': y centroid of ROI
#' \item 'offset.x': ROI distance from horizontal center
#' \item 'offset.y': ROI distance from vertical center
#' \item 'offset.total': sum of ROI offset.x and offset.y
#' \item 'ar': aspect ration of the ROI
#' \item 'size': size of ROI in pixels
#' }

#'  \code{head.lms}  'lm' objects, one for each frame described by \code{frames} of the linear model fit to the \code{ant.per} section of the ROI
#' @seealso \code{\link{kin.vid}}
#' @export
#' @examples
#' #produce a classic midline waveform plot of swimming fish
#' \dontrun{
#' require(wesanderson)
#' require(ggplot2)
#' require(data.table)
#' require(dplyr)
#' require(EBImage)
#'
#' #download example images and place in 'example' subdirectory
#' f <- "https://github.com/ckenaley/exampledata/blob/master/example.zip?raw=true"
#' download.file(f,"temp.zip")
#' unzip("temp.zip")
#' unlink("temp.zip")
#'
#'kin <- kin.img2(image.dir ="example",search.for = "largest",smoothing = "spline",frames=1:15,show.prog = T,thr = 0.6,make.video = TRUE,silent=FALSE)

#' ml <- kin$midline
#' #normalize x (y is normalized to midline by "kin.img/kin.vid")
#' ml <- ddply(ml,.(frame),transform,x2=x-x[1])
#'
#' #compute instantaneous amplitude of tail (last/rightmost point) and wave crest x position  by frame
#' ml2 <- ddply(ml,.(frame),summarize,amp.i=abs(last(wave.y)))
#'
#' ml <- merge(ml,ml2,by="frame") #merge these
#'
#' pal <- wes_palette("Zissou1", 100, type = "continuous") #"Zissou" color palette
#' p <- ggplot(dat=ml,aes(x=x2,y=wave.y))+theme_classic(15)+scale_color_gradientn(colours = pal)
#' p <- p+geom_line(aes(group=frame,color=amp.i),stat="smooth",method = "loess", size = 1.5,alpha = 0.5)
#' print(p)
#'}
#'
kin.img2 <-function(image.dir=NULL,sequenced=TRUE,frames=NULL,thr=0.7,plot.midline=TRUE, show.prog=FALSE,ant.per=0.10,smoothing="loess",smooth=.2, smooth.points=200, image.type="orig",flip=TRUE,rem.file=FALSE,make.video=TRUE,video.name=NULL,qual=50,frame.rate=10,size.min=0.05,n.blob=NULL,search.for="largest",burn=10,silent=TRUE){

  unlink("processed_images",recursive = T)
  dir.create("processed_images")

  proc.dir <- "processed_images"

  images <- paste0(image.dir,"/",list.files(image.dir)[!grepl("Icon\r",list.files(image.dir))]) #remove pesky Icon\r
print(images)
  if(any(frames>length(images))) stop("variable 'frames' out of range of image sequence")
  if(!is.null(frames)) images <- images[frames]

  pb = txtProgressBar(min = 0, max = length(images), initial = 0,style=3)


  trial <- gsub("\\.[^.]*$", "", basename(images[1]))
  if(sequenced==TRUE) trial <- gsub("(.+)\\_\\d*$", "\\1", trial)

  kin.dat <- list()
  midline.dat <- list()
  lms <- list()

  #assume largest ROI on first frame
  prob.roi <- "a"

  for(im in images){


    frame <- which(im==images)-1 #could cause problems
    img <- EBImage::readImage(im,all=F) #if don't add package, others use "display"
    ## computes binary mask
    y = img >thr #contrast threshold
    if(flip){#flip binary
      y[y==1] <- 5
      y[y==0] <- 1
      y[y==5] <- 0
    }
    z = bwlabel(y)

    rois <- tabulate(z)

    ###idealized fish
   #  fish <- EBImage::readImage("/Users/Chris/Documents/R.scripts/trackter/fish.jpg",all=F) #if don't add package, others use "display"
   #  ## computes binary mask
   #  fish.y = fish >thr #contrast threshold
   #  if(flip){#flip binary
   #    fish.y[fish.y==1] <- 5
   #    fish.y[fish.y==0] <- 1
   #    fish.y[fish.y==5] <- 0
   #  }
   #  fish.z = bwlabel(fish.y)
   #
   #  fish.m <- fish.z[,,1]
   # fish.m[1,1] <- 0

   #fish.c <- ocontour(fish.m)
    #fish.c.df <- data.frame(fish.c)

   # fish.co <- efourier(ocontour(fish.m),nb.h = 40)
   # efourier_norm(ef=fish.co)

#plot(fish.c.df[,1],fish.c.df[,2])

    roi <- which.max(rois)#find the largest roi, can tell it to find nth largest blob

    if(search.for=="size") roi <- which(rois==rois[order(rois,decreasing =T)[n.blob]])

    pix <- dim(z[,,1])[1]*dim(z[,,1])[2]
    w <- dim(z[,,1])[1] #width of image
    h <- dim(z[,,1])[2] #height of image
    per <- rois/(w*h) #how big are rois compared to pixel field

    c.roi <-  which(per>=size.min) #candidate rois, filtered by size of of pixel field

    names(c.roi) <- as.factor(letters[order(rois[c.roi],decreasing = T)])

    setTxtProgressBar(pb,which(images==im))


    if(!any(c("largest","offset")==search.for)) stop("'search.for' must = 'offset' or 'largest'")



    #c.roi <- c.roi[which(names(c.roi)==prob.roi)]
    cand.kin <- list() #store each roi kin data
    cand.mid <- list() #store each roi midline data
    z.l <- list()
    void.l <- list() #store eval of centroids, empty or not


    if(!is.null(burn) && which(im==images)>burn && all(prob.rois==prob.rois[1])) c.roi <- c.roi[which(names(c.roi)==prob.rois[1])]

    for(r in c.roi){
      r.name <- as.character(names(c.roi)[c.roi==r])
      z.r <- z
      z.r[z!=r] <- 0
      z.r[z==r] <- 1
      z.m <- z.r[,,1]
      z.m[1,1] <- 0 #this gets a 1 when



      z.l[[r.name]] <- z.m
      length <- diff(range(unlist(apply(z.m,2,function(x) which(x==1)))))
      width <- diff(range(unlist(apply(z.m,1,function(x) which(x==1)))))


      ar <- length/width
      #centroids
      cent.x <- mean(unlist(apply(z.m,2,function(x) which(x==1))))
      cent.y <- mean(unlist(apply(z.m,1,function(x) which(x==1))))
      offset.x <- abs(cent.x-w/2)# distance
      offset.y <- abs(cent.y-h/2)
      offset.tot <- offset.x+offset.y
      points(cent.x,cent.y,col="white",pch=16)


      y.max <- apply(z.m,1,function(x) ifelse(any(x==1),max(which(x==1)),NA))

      y.min <- apply(z.m,1,function(x) ifelse(any(x==1),min(which(x==1)),NA))

      y.df <- data.table(x=1:nrow(z.m),y.max,y.min)

      y.df <- y.df[,y.m :=(y.max-y.min)/2+y.min, by=.(x,y.min,y.max)]

      # zc <- data.frame(ocontour(z.m))
      # colnames(zc) <- c("x","y")
      # with(zc,plot(x,y))



      # ef.co <- efourier(ocontour(z.m),nb.h = 40)
      # efourier_norm(ef=ef.co)

      tip.y <- mean(tail(y.df$y.m[!is.na(y.df$y.m)],30))#tip is mean y.m of last 30 pixels
      tip.x <- mean(tail(y.df$x[!is.na(y.df$y.m)],30))#tip is mean y.m of last 30 pixels

      head.y <- mean(head(y.df$y.m[!is.na(y.df$y.m)],30))#tip is mean y.m of first 30 pixels
      head.x <- mean(head(y.df$x[!is.na(y.df$y.m)],30))#tip is mean y.m of first 30 pixels

      #n midline points
      if(is.null(smooth.points)) smooth.points <- nrow(y.df)
      midline <- y.df[seq(1,nrow(y.df),length.out = smooth.points),] #two hundred points on midline

      midline <- midline[complete.cases(midline)]

      #head section
      head.dat <- midline[1:(ant.per*smooth.points),]
      head.lm <- lm(y.m~x,head.dat)

      head.p <- summary(head.lm)$r.squared #how well does head lm fit
      midline$mid.pred <- predict(head.lm,newdata=midline)#add lm prediction to midline df
      midline <- midline[complete.cases(midline),]


      #midline pixels
      z.dat <- imageData(z.m)#extract image data
      m.pix <- apply(midline[,c("x","y.m")],1,function(x) z.dat[x[1],x[2]])

      ####which type of lines to be fitted, spline or loess
      if(!any(c("spline","loess")==smoothing)) stop("'smoothing' must = 'loess' or 'spline'")

      if(smoothing=="loess")  ml.pred <- fitted(loess(midline$y.m~midline$x,span=smooth,degree=1))
      if(smoothing=="spline") ml.pred <- smooth.spline(x = midline$x,y=midline$y.m,spar=smooth)$y

      midline[,y.pred:=ml.pred]#add smoothed predictions

      midline[,wave.y:=dist.2d(x,x,y.pred,mid.pred)] #wave y based on  pred points
      midline[y.pred<mid.pred,wave.y:=wave.y*-1] #neg or pos amp
      midline[,roi:=r.name]
      n.roi <- paste0(basename(im),"-",r)
      cand.kin[[r.name]] <- data.frame(frame,x=tip.x,y=tip.y,head.x,head.y,amp=last(midline$wave.y),head.pval=head.p,roi=r.name,cent.x,cent.y,offset.x,offset.y,offset.tot,ar,size=rois[c.roi[r.name]])
      cand.mid[[r.name]] <- midline


      void.l[[r.name]] <- sum(m.pix)/length(m.pix)>0.95 #is the midline full of positive body pixels (more than 90%)

    }

    cand.kin.i <- data.table(do.call(rbind,cand.kin)) #store all candidate roi kin data
    cand.mid.i <- data.table(do.call(rbind,cand.mid))

    #store raw and best kin data

    if(search.for=="largest") prob.roi <- "a"

    off.min <- cand.kin.i[last(order(offset.tot,decreasing = F)),]$roi
    if(search.for=="offset") prob.roi <- off.min
    if(!prob.roi%in%names(void.l)[which(void.l==T)]) (warning(paste0("entire midline of probable ROI for image ", im, " does not have positive pixels. Try other 'seach.for' options")))


    if(!is.null(burn) & which(im==images)<=burn &length(prob.roi)>1) prob.rois <- combine(prob.roi)

    kin.i <- cand.kin.i[roi==prob.roi]
    mid.i <- cand.mid.i[roi==prob.roi]
    if(show.prog) {
      EBImage::display(z.l[[prob.roi]],method = "raster")
    }
    midline.dat[[basename(im)]] <- data.frame(frame,mid.i)
    kin.dat[[basename(im)]] <- kin.i
    jpeg(paste0(proc.dir,"/",trial,"_",sprintf("%03d",frame),".jpg"),quality = 0.5)
    if(image.type=="bin") EBImage::display(z,method = "raster")
    if(image.type=="orig") EBImage:: display(img,method = "raster")
    if(plot.midline) {

      lines(predict(lm(mid.pred~x,mid.i)),x=mid.i$x,col="blue",lwd=4)
      lines(y.pred~x,mid.i,col="red",lwd=4)
      with(mid.i[1:(ant.per*smooth.points),],points(x,y.m,col="green",pch=16,cex=0.75))
    }
    dev.off()
  }


  kin.dat <- do.call(rbind,kin.dat)
  midline.dat <- data.table(do.call(rbind,midline.dat))
  rownames(midline.dat) <- NULL
  
  if(make.video){vid.name <- paste0(trial,"_kin")}
  if(make.video) {images.to.video2(image.dir = proc.dir, vid.name = vid.name, qual=qual,frame.rate = frame.rate,silent=silent)
message(paste0(vid.name," outputted to ",getwd()))
  }
  #clean up
  if(rem.file){
    if(is.null(image.dir)) stop("'image.dir' not specified and 'rem.file=TRUE'. Won't delete working directory!")
    unlink(proc.dir,recursive = T)
    unlink(image.dir,recursive = T)
  }


  return(list(kin.dat=kin.dat,midline=midline.dat))
}

######### kin.LDA

#' @title  Midline tracking over image sequences with ROI search using LDA

#' @description  Experimental. Automatically retrieves the midline of a fish-like ROI class detected through linear discrimate analysis (LDA) of PCA on elliptical Fourier described shapes. Initial training of ROI is user defined or with the 'train.dat' data set loaded with \code{trackter} (see details). For each detected ROI, \code{kin.LDA} finds the y-value midpoint along the x-value array of the ROI and fits a midline according to a chosen smoothing method (loess or spline). Also outputs the midline amplitude relative to a reference line determined by an anterior section of the ROI. Supported image formats are jpeg, png, and tiff.
#'
#' @param image.dir character, directory containing images to analyze.
#' @param frames numeric, vector indicating which images to process.
#' @param thr numeric, threshold to determine binary image. May require some tweaking through iteration.
#' @param enorm logical, should the EFA coeffecients from \code{efourier} operations be normalized or not. See \code{details} and \code{?Momocs::efourier} 
#' @param harms numeric, the number of harmonics to use. If missing, \code{Momocs} sets 'nh.b' to 12. Will produce messages.
#'
#'@param rescale logical, should all shapes in PCA be rescaled. Performs best as 'FALSE'.
#'@param train.dat Classified \code{Out} and \code{Coo} outlines that are produced from \code{Momocs}. See details. 
#'@param train numeric, the number of frames on which to retrain the LDA data set. See details.
#'@param after.train character, if set to 'size', LDA will be skipped after \code{retrain} and the ROI with a size closest to the ROI found by the LDA $>=$ will be chosen. This peeds calculations considerably. If 'LDA', the default, 'LDA' will continue usining the retraining classifications from frames $<=$ 'train'.
#'
#'@param ties character, how to chose ROI's in any one frame that appear fish-like. See details.
#' @param smoothing character, the midline smoothing method, either 'loess' or 'spline'.
#' @param smooth numeric; if \code{smoothing} is set to 'loess', 'span' parameter value for \code{\link{loess}}. If \code{smoothing} is set to 'spline' 'spar' parameter value for \code{\link{smooth.spline}}
#' @param smooth.points numeric, number of equally spaced points along the ROI midline on which the smoothed midline is computed.
#' @param show.prog logical value indicating if outputted image should be displayed during analysis.
#' @param ant.per numeric; left-most percentage of ROI that establishes the vertical reference for the midline displacement.
#' #' @param tips, numeric, the proportion the the midline data to use in calculation of the head and tail position. 
#' @param size.min numeric, indicating the minimum size of ROIs as a proportion of the pixel field to be considered in analysis. May be useful if smaller unimportant ROIs appear in the frame. Default is 0.05.
#' @param save logical, value indicating if images should be outputted with midline and predicted midline based on the \code{ant.per} \code{lm()} overlaying original or binary images.
#' @param image.type character; the type of image to be outputted, either 'orig' or 'bin' representing the original or binary images, respectively. Ignored if 'save==FALSE'.
#' #' @param plot.pml logical, value indicating if outputted images should include the predicted midline (in blue) and the points according to \code{ant.per} used to construct the predicted midline (in green).
#' @param flip logical, indicating if binary should be flipped.
#' 
#' @export
#'
#' @details
#'The algorithm assumes a left-right orientation, i.e., the head of the ROI is positioned left, the tail right. ffmpeg operations or even imageJ can rotate images not in this orientation. The \code{ant.per} value therefor establishes the reference line (theoretical straight midline) based on that portion of the head.  If 'save=TRUE', images are saved as binary or the original with a body mideline overlay and, if chosen, with the theoretical midline (based on \code{ant.per}). 
#'
#'Before \code{train}, ROIs are chosen according to LDA of a PCA object constructed from \code{efourier} analysis. LDA is trained by a user define 'train.dat' when the frame $<=$ \code{retrain}. LDA will procede after \code{retrain} if \code{after.train}='LDA', but the LDA will be trained by the contours classified as 'fish' and 'not.fish' found during the chosen training period. 
#'
#'\code{enorm} Normalization of EFA coefficents is often perilous, especially for symetrical shapes, a conditional met for undulating, bilaterally sysmetical organisms at least some of the time and perhaps for many of the frames included in any analysis. Thus, 'enorm' by default is set to 'FALSE'. 'enorm=TRUE' may produce odd ROI choices and should be used cautiously.
#'
#'\code{train.dat} This should be a \code{Coo} and \code{Out} object produced by \code{efourier} analysis of predefined shapes from \code{Momocs::efourier}. A user defined dataset or the 'shapes' dataset in \code{trackter} must be used for training. 'shapes' includes several arbitrary shapes (circrles, squares, U-shapes, etc.) as well as several fish shapes (sunfish (genus Lepomis), eel (genus Anguilla), and trout (genus Onchorhynchus) swimming over one tail-beat cycle). A user-defined dataset must have shapes classified with factors identical to the 'fisshapes' contours, that is by shape, type, and edge. Shape levels should indicate what type of shape is described by the contour (e.g., 'circle', 'L-shape', 'trout', 'eel', etc). The type levels must describe the shape as 'fish' or 'not.fish'. The edge levels must be 'FALSE'. See also \link{fishshapes} and \link{Momocs}.
#'
#'\code{edges} Set by default to 'FALSE'. It is not advisable to include shapes that are on the edge of any frame and are therfore incomplete.
#'\code{retrain} After this value, the LDA analysis will use the ROIs determined as 'fish' and 'not.fish' in the frames $>=$ \code{retrain} to discrimate fish from non-fish shapes. This speeds up analysis considerably.
#'\code{ties} Determiens how to chose ROIs if more than one fish-like ROI is found in any frame. 'fish' will result in chosing the ROI with shape types in which the best *and* second-best fish-like shape (according to posterior probabilities) match a fish-like shape in the trainin and/or retraining datasets.'post' will chose the best fish-like shape according the the highest posterior probability from LDA.
#'
#' @return A list with the following components:
#'
#' \code{kin.dat} a data table consisting of frame-by-frame position parameters for the ROI determined by LDA analysis.
#' \itemize{
#' \item the frame number
#'
#' \item 'x' and ''y': the position of the tail (rightmost or posteriormost)
#' \item 'head.x' and 'head.y': the x and y position of the head (leftmost or anteriormost)
#' \item 'amp': the amplitude (\code{amp}) of the tail relative to thr theoretical midline determined by the \code{lm()} predictions from \code{ant.per}
#' \item 'roi': a character indicating the ROI ranked by size ('a' being the largest)
#' \item 'head.pval': p values of the \code{lm()} fit that describes the position of the head as determined by \code{ant.per} (green points in the outputted images/video)}
#'
#' \code{midline} A data table containing, for each frame described by \code{frames}, the following: \itemize{
#' \item 'x' and 'y.m': x and y positions of the midline of the ROI
#' #' \item 'y.min' and 'y.max': min and max y positions ROI's countour used in y.m calculation
#' \item 'mid.pred': the predicted linear midline based on the points/pixels defined by \code{head.per} (green points in the outputted images/video)
#' \item 'y.pred': midline points fit to a smooth spline or loess model with spar or span equal to \code{smooth} (red curve in the outputted images/video)
#' \item 'wave.y': midline points 'y.pred' normalized to 'mid.pred'
#' \item 'roi': a character indicating ROI size (a being the largest)
#' }
#' 
#' \code{cont} A data table containing x and y positions of the contours used to calculate the data in 'kin.dat'. Contains the following: 
#' \itemize{
#' \item 'frame': the frame
#' #' \item 'x' and 'y': the x and y positions of the countours
#' }
#' 
#' #' \code{all.classes} A data table containing the following: 
#' \itemize{
#' \item 'frame': the frame
#' \item 'roi': the name of each ROI found in a frame.
#' \item 'size': the size of each ROI
#' }
#' @seealso \code{\link{kin.simple}}, \code{\link{efourier}} \code{\link{LDA}}
#' @export
#' @import data.table
#' @import Momocs
#' @examples
#' #produce a classic midline waveform plot of swimming fish searching a image field with a two fish-like ROIs
#' \dontrun{
#' require(wesanderson)
#' require(ggplot2)
#' require(data.table)
#' require(dplyr)
#' require(EBImage)
#'
#' #download example images and place in 'example' subdirectory
#' f <- "https://github.com/ckenaley/exampledata/blob/master/example.zip?raw=true"
#' download.file(f,"temp.zip")
#' unzip("temp.zip")
#' unlink("temp.zip")
#'kin <- kin.simple(image.dir = "example",frames=1:20,thr=0.7,ant.per=.25,enorm=F,show.prog = F,retrain=2,train.dat = train,after.train="LDA",edges=F)

#' ml <- kin$midline
#' #normalize x (y is normalized to midline by "kin.img/kin.vid")
#' ml <- ddply(ml,.(frame),transform,x2=x-x[1])
#'
#' #compute instantaneous amplitude of tail (last/rightmost point) and wave crest x position  by frame
#' ml2 <- ddply(ml,.(frame),summarize,amp.i=abs(last(wave.y)))
#'
#' ml <- merge(ml,ml2,by="frame") #merge these
#'
#' pal <- wes_palette("Zissou1", 100, type = "continuous") #"Zissou" color palette
#' p <- ggplot(dat=ml,aes(x=x2,y=wave.y))+theme_classic(15)+scale_color_gradientn(colours = pal)
#' p <- p+geom_line(aes(group=frame,color=amp.i),stat="smooth",method = "loess", size = 1.5,alpha = 0.5)
#' print(p)
#'}
#'

kin.LDA <-function(image.dir=NULL,frames=NULL,thr=0.7,train.dat=NULL,rescale=F,harms=15,enorm=T,retrain=5,after.train="LDA",ties="fish",edges=F,size.min=0.05,show.prog=FALSE,ant.per=0.20,tips=0.2,smoothing="loess",smooth=.3,smooth.points=200,save=T,image.type="orig",plot.pml=TRUE,flip=TRUE){
  
  
  if(is.null(train.dat)) train.dat<- readRDS("shapes.RDS") ###load training data 
  
  if(save){unlink("processed_images",recursive = T)
  dir.create("processed_images")
  
  proc.dir <- "processed_images"
  }
  
  images <- paste0(image.dir,"/",list.files(image.dir)[!grepl("Icon\r",list.files(image.dir))]) #remove pesky Icon\r
  
  if(any(frames>length(images))) stop("variable 'frames' out of range of image sequence")
  if(!is.null(frames)) images <- images[frames]
  
  trial <- gsub("\\.[^.]*$", "", basename(images[1]))
  if(sequenced==TRUE) trial <- gsub("(.+)\\_\\d*$", "\\1", trial)
  
  kin.l <- list()
  midline.l<- list()
  classes.l <- list()
  
  lms <- list()
  conts <- list()
  pb = txtProgressBar(min = 0, max = length(images), initial = 0,style=3)
  
  roi.outs <- list() #store the rois for each image
  for(im in images){
    
    frame <- which(im==images)-1
    
    img <- EBImage::readImage(im,all=F) #if don't add package, others use "display"
  
    ## computes binary mask
    y = img >thr #contrast threshold
    if(flip){#flip binary
      y[y==1] <- 5
      y[y==0] <- 1
      y[y==5] <- 0
    }
    z = bwlabel(y)
    rois <- tabulate(z)
    pix <- dim(z[,,1])[1]*dim(z[,,1])[2]
    w <- dim(z[,,1])[1] #width of image
    h <- dim(z[,,1])[2] #height of image
    per <- rois/(w*h) #how big are rois compared to pixel field
    
    c.roi <-  which(per>=size.min) #candidate rois, filtered by size of of pixel field
    
    names(c.roi) <- as.factor(letters[order(rois[c.roi],decreasing = T)])
    
    z.l <- list()
    out.l <- list()
    size.l <- list()
    
    if(frame==retrain){
      kin.train <- do.call(rbind,kin.l)
      out.train <- roi.outs
      size.train <- do.call(rbind,size.l)
      roi.train <-combine(roi.outs)
      
    }
    
    if(is.null(after.train)) stop("'after.train' not set to 'size' or 'LDA'.")
    
    if(!after.train %in% c("LDA","size")) stop("after.train not set to 'size' or 'LDA'.")
    
    if(after.train=="size" & frame>retrain){
      sz.m <- mean(kin.train$size)
      size.diff <- rois[c.roi]-sz.m
      
      for(r in c.roi){
        r.name <- as.character(names(c.roi)[c.roi==r])
        z.r <- z
        z.r[z!=r] <- 0
        z.r[z==r] <- 1
        z.m <- z.r[,,1]
        z.m[1,1] <- 0 #this gets a 1 when
        z.l[[r.name]] <- z.m
        
        z.c <- ocontour(z.m)
        
        wall <- any(z.c[[1]][,1]>dim(z)[1]-2 | z.c[[1]][,1]<2  |z.c[[1]][,2]>dim(z)[2]-2 | z.c[[1]][,2]<2)
        
        
        r.out <- Out(ocontour(z.m))
        if(wall ) edge <- T
        if(!wall) edge <- F
        r.out$fac <- data.frame(shape=paste0("roi-",r.name),type=paste0("roi"),edge=edge)
        out.l[[r.name]] <- r.out
        rois[c.roi[r.name]]
      }
      
      roi.out2 <- combine(out.l)
      
      classes <- data.table(roi=gsub("roi-","",roi.out2$fac$shape),type=NA,shape=NA,post=NA,type2=NA,shape2=NA,post2= NA,retrained=frame>retrain,edge=roi.out2$fac$edge,size.diff=size.diff)
      
      classes.l[[paste0(frame)]] <- data.table(frame=frame,classes)
      
      if(edges==F) classes <- classes[edge==F]
      
      z.best <- z.l[[classes[which.min(size.diff)]$roi]]
      r.name <-classes[which.min(size.diff)]$roi
      best.class <- classes[which.min(size.diff)]
    }
    else{
      for(r in c.roi){
        r.name <- as.character(names(c.roi)[c.roi==r])
        z.r <- z
        z.r[z!=r] <- 0
        z.r[z==r] <- 1
        z.m <- z.r[,,1]
        z.m[1,1] <- 0 #this gets a 1 when
        z.l[[r.name]] <- z.m
        
        z.c <- ocontour(z.m)
        
        wall <- any(z.c[[1]][,1]>dim(z)[1]-2 | z.c[[1]][,1]<2  |z.c[[1]][,2]>dim(z)[2]-2 | z.c[[1]][,2]<2)
        
        r.out <- Out(ocontour(z.m))
        if(wall ) edge <- T
        if(!wall) edge <- F
        r.out$fac <- data.frame(shape=paste0("roi-",r.name),type=paste0("roi"),edge=edge)
        out.l[[r.name]] <- r.out
        
        rois[c.roi[r.name]]
      }
      
      ##### efourier and LDA analysis ####
      #train.data from training
      roi.out2 <- combine(out.l)
  
      if(frame>retrain){all.outs <- combine(roi.train,roi.out2)}else{all.outs <- combine(train.dat,roi.out2)}
      #rescale?
      if(rescale) all.outs <- coo_scale(coo_center(all.outs))
      # panel(all.outs)
      if(length(roi.out2)==1) all.outs <- combine(all.outs,roi.out2) #must have two roi rows for rePCA to work
      
      shapes.out <- all.outs %>% filter(type!="roi")

      roi.out <- filter(all.outs,type=="roi") #shape that is roi
      roi.f <- suppressMessages(efourier(roi.out,nb.h = harms,norm=enorm,start=F) )#ef analysis
      
      shapes.f <- suppressMessages(efourier(shapes.out,nb.h = harms,norm=enorm,start = F))
    
      shapes.p <- suppressWarnings(PCA(shapes.f))#PCA of non roi
      #why LDA wouldn't drop levels to prevent warnings is a stumper
      shapes.p$fac$shape <- droplevels(shapes.p$fac$shape)
      
      #if few shapes classified as same, keep a tone of variance
      if(grepl("-",shapes.f$fac$shape[1])){ suppressMessages(shape.l <- LDA(shapes.p,"shape",retain=0.9999999))}else{ 
        suppressMessages(shape.l <- LDA(shapes.p,"shape"))
      }#LDA of non roi
      
      #redo the same PCA and LDA  with roi
      roi.pca <- suppressMessages(rePCA(shapes.p, roi.f))
      
      roi.lda <- suppressMessages(reLDA(roi.pca,shape.l))
      
      roi.shape<- roi.lda$class
      roi.shape2 <- apply(roi.lda$posterior,1, function(x) names(x[order(x,decreasing = T)[2]]))
      
      #retrieve roi type
      roi.type <- sapply(roi.shape,function(x) unique(filter(shapes.p$fac,shape==x)$type))
      
      roi.post <- apply(roi.lda$posterior,1,max)
      roi.post2 <- apply(roi.lda$posterior,1, function(x) x[order(x,decreasing = T)[2]])
      roi.type2 <- sapply(roi.shape2,function(x) unique(filter(shapes.p$fac,shape==x)$type))
      
      classes <- data.table(roi=gsub("roi-","",roi.out2$fac$shape),type=roi.type,shape=roi.shape,post=roi.post,type2=roi.type2,shape2=roi.shape2,post2= roi.post2,retrained=frame>retrain,edge=roi.out2$fac$edge,size.diff=NA)
      
      classes.l[[paste0(frame)]] <- data.table(frame=frame,classes)
      
      if(edges==F) classes <- classes[edge==F]
      
      if(nrow(classes[type=="fish"])<1){stop("No shape type found matching 'fish' or 'edge==F'")}
      
      if(nrow(classes[type=="fish"])>1){
        if(ties=="post"){
          
          z.best <- z.l[[classes[type=="fish"][which.max(post)]$roi]]
          r.name <-classes[type=="fish"][which.max(post)]$roi
          best.class <- classes[type=="fish"][which.max(post)]
          warning("More than one 'fish' found in LDA. Tie broken with highest post. prob.")}
        if(ties=="fish"){
          ###is there just one roi with type and type2 fish?
          if(nrow(classes[type=="fish" & type2=="fish"])==1){
            
            z.best <- z.l[[classes[type=="fish" & type2=="fish"]$roi]]
            r.name <-classes[type=="fish" & type2=="fish"]$roi
            best.class <- classes[type=="fish" & type2=="fish"]
          }
          
          ###are there two rois with type and type2 fish?
          if(nrow(classes[type=="fish" & type2=="fish"])>1){
            
            z.best <- z.l[[classes[type=="fish" & type2=="fish"][which.max(post)]$roi]]
            r.name <-classes[type=="fish" & type2=="fish"][which.max(post)]$roi
            best.class <- classes[type=="fish" & type2=="fish"][which.max(post)]
            warning("'ties=fish' and >1 ROIs have best and second best type as 'fish'; breaking tie with sum of post. prob.")}
          if(nrow(classes[type=="fish" & type2=="fish"])==0){
            classes.l[[paste0(frame)]] <- data.table(frame=frame,classes)
            z.best <- z.l[[classes[type=="fish"][which.max(post)]$roi]]
            r.name <-classes[type=="fish"][which.max(post)]$roi
            best.class <- classes[type=="fish"][which.max(post)]
            warning("'fish' chosen to break tie, but no second best shapes='fish'. Chose 'fish' with highest post. prob")
          }
        }
      }
      else{ 
        r.name <- classes$roi
        z.best <- z.l[[r.name]]
        best.class <- classes
      }
    }  
    if(show.prog) {
      EBImage::display(z.best,method = "raster")
    }
    
    best.cont <- data.table(ocontour(z.best)[[1]])
    colnames(best.cont) <- c("x","y")
    
    conts[[paste0(frame)]] <- data.table(frame=frame,best.cont)
    
    y.df <- best.cont[,.(y.min=min(y),y.max=max(y),y.m=mean(y)),by=.(x)]
    setkey(y.df,"x")
    
    ends <- ceiling(nrow(y.df)*tips)
    tip.y <- mean(tail(y.df$y.m[!is.na(y.df$y.m)],ends))#tip is mean y.m of last 30 pixels
    tip.x <- mean(tail(y.df$x[!is.na(y.df$y.m)],ends))#tip is mean y.m of last 30 pixels
    
    head.y <- mean(head(y.df$y.m[!is.na(y.df$y.m)],ends))#tip is mean y.m of first 30 pixels
    head.x <- mean(head(y.df$x[!is.na(y.df$y.m)],ends))#tip is mean y.m of first 30 pixels
    
    #n midline points
    if(is.null(smooth.points)) smooth.points <- nrow(y.df)
    midline <- y.df[seq(1,nrow(y.df),length.out = smooth.points),] #two hundred points on midline
    
    midline <- midline[complete.cases(midline)]
    midline <- data.table(frame,midline)
    
    
    ####which type of lines to be fitted, spline or loess
    if(!any(c("spline","loess")==smoothing)) stop("'smoothing' must = 'loess' or 'spline'")
    
    if(smoothing=="loess")  ml.pred <- fitted(loess(midline$y.m~midline$x,span=smooth,degree=1))
    if(smoothing=="spline") ml.pred <- smooth.spline(x = midline$x,y=midline$y.m,spar=smooth)$y
    
    midline[,y.pred:=ml.pred]#add smoothed predictions
    
    #head section
    head.dat <- midline[1:(ant.per*smooth.points),]
    head.lm <- lm(y.pred~x,head.dat)
    
    head.p <- summary(head.lm)$r.squared #how well does head lm fit
    
    midline$mid.pred <- predict(head.lm,newdata=midline)#add lm prediction to midline df
    
    midline <- midline[complete.cases(midline),]
    
    midline[,wave.y:=y.pred-mid.pred] #wave y based on midline y and straight head.lm pred points
    
    midline[,roi:=r.name]
    n.roi <- paste0(basename(im),"-",r)
    
    kin.l[[paste(frame)]] <- data.table(frame,x=tip.x,y=tip.y,head.x,head.y,amp=last(midline$wave.y),head.pval=head.p,size=rois[c.roi[r.name]],best.class)
    midline.l[[paste(frame)]] <- midline
    
    if(frame<=retrain){
      roi.out2$fac <- data.frame(shape=roi.shape,type=roi.type,edge=roi.out2$fac$edge)
      #are all shapes the same, if so give them unique levels
      if(all(roi.shape==roi.shape[1]))   roi.out2$fac <-  data.frame(shape=paste0(roi.shape,"-",index(roi.shape)),type=roi.type,edge=roi.out2$fac$edge)
      
      roi.outs[[paste(frame)]] <- roi.out2
    } #save for retraining
    
    
    if(save){
      
      jpeg(paste0(proc.dir,"/",trial,"_",sprintf("%03d",frame),".jpg"),quality = 0.5)
      if(image.type=="bin")EBImage::display(z,method = "raster")
      if(image.type=="orig")EBImage:: display(img,method = "raster")
      
        
      if(plot.pml) lines(predict(lm(mid.pred~x,midline)),x=midline$x,col="blue",lwd=4)
        with(midline,lines(y.pred~x,col="red",lwd=4))
        if(plot.pml) with(midline[1:ceiling(ant.per*smooth.points),],points(x,y.pred,col="green",pch=16,cex=0.75))
      
      dev.off()
    }
    setTxtProgressBar(pb,which(images==im))
  }
  
  classes.dat <- do.call(rbind,classes.l)
  kin.dat <- do.call(rbind,kin.l)
  midline.dat <- data.table(do.call(rbind,midline.l))
  cont.dat <- do.call(rbind,conts)
  
  return(list(kin.dat=kin.dat,midline=midline.dat,cont=cont.dat,all.classes=classes.dat))
}

#### kin.simple

#' @title  Simplified midline tracking over image sequences 

#' @description  Automatically retrieves the midline of an ROI detected based on size. Assumes the ROI of interest is the largest detected and not interesecting the edges of the image frame, conditions often met in kinematic studies. For each ROI of interest, finds the y-value midpoint along the x-value array of the ROI and fits a midline according to a chosen smoothing method (loess or spline). Also outputs the midline amplitude relative to a reference line determined by an anterior section of the ROI and outputs contours ROIs in each frame for subsequent analysis. Supported image formats are jpeg, png, and tiff.
#'
#'
#' @param image.dir character, directory containing images to analyze.
#' @param frames numeric, vector indicating which images to process.
#' @param thr numeric, threshold to determine binary image. May require some tweaking through iteration.
#' @param size.min numeric, indicating the minimum size of ROIs as a proportion of the pixel field to be considered in analysis. May be useful if smaller unimportant ROIs appear in the frame. Default is 0.05.
#' @param ant.per numeric; left-most proportion of ROI that establishes the vertical reference for the midline displacement.
#' @param tips, numeric, the proportion the the midline data to use in calculation of the head and tail position. 
#' @param smoothing character, the midline smoothing method, either 'loess' or 'spline'.
#' @param smooth numeric; if \code{smoothing} is set to 'loess', passed to 'span' parameter of \code{\link{loess}}. If \code{smoothing} is set to 'spline', passed to 'spar' parameter of \code{\link{smooth.spline}}
#' @param smooth.points numeric, number of equally spaced points along the ROI midline on which the smoothed midline is computed.
#' @param save logical, value indicating if images should be outputted with midline and predicted midline based on the \code{ant.per} \code{lm()} overlaying original or binary images.
#' @param plot.pml logical, value indicating if outputted images should include the predicted midline (in blue) and the points according to \code{ant.per} used to construct the predicted midline (in green).
#' @param image.type character; the type of image to be outputted, either 'orig' or 'bin' representing the original or binary images, respectively. Ignored if 'save=FALSE'.
#' @param flip logical, indicating if binary image should be flipped.
#' @param show.prog logical, indicating if outputted image should be displayed during analysis.
#' 
#' @export
#'
#' @details
#'The algorithm assumes a left-right orientation, i.e., the head of the ROI is positioned left, the tail right. ffmpeg operations or even imageJ can rotate images not in this orientation. The \code{ant.per} value therefor establishes the reference line (theoretical straight midline) based on that portion of the head.  If 'save=TRUE', images are saved as binary or the original with a body mideline overlay and, if chosen, with the theoretical midline (based on \code{ant.per}). 
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
#' #' \item 'y.min' and 'y.max': min and max y positions ROI's countour used in y.m calculation
#' \item 'mid.pred': the predicted linear midline based on the points/pixels defined by \code{ant.per} (green points in the outputted images/video if 'plot.pml=TRUE')
#' \item 'y.pred': midline points fit to a smooth spline or loess model with spar or span equal to \code{smooth} (red curve in the outputted images/video)
#' \item 'wave.y': midline points 'y.pred' normalized to 'mid.pred'
#' \item 'roi': a character indicating ROI size (a being the largest)
#' }
#' 
#' \code{cont} A data table containing x and y positions of the contours used to calculate the data in 'kin.dat'. Contains the following: 
#' \itemize{
#' \item 'frame': the frame
#' \item 'x' and 'y': the x and y positions of the countours
#' }
#' 
#' \code{all.classes} A data table containing the following: 
#' \itemize{
#' \item 'frame': the frame
#' \item 'roi': the name of each ROI found in a frame.
#' \item 'size': the size of each ROI
#' }
#' @seealso \code{\link{kin.simple}}, \code{\link{efourier}} \code{\link{LDA}}
#' @export
#' @import data.table
#' @import Momocs
#' @examples
#' #produce a classic midline waveform plot of swimming fish searching a image field with a two fish-like ROIs
#' \dontrun{
#' require(wesanderson)
#' require(ggplot2)
#' require(data.table)
#' require(dplyr)
#' require(EBImage)
#'
#' #download example avi video
#' f <- "https://github.com/ckenaley/exampledata/blob/master/trout1_63_test.avi?raw=true"
#' download.file(f,"trout1_63_test.avi")
#' 
#' #extract images and reduce them to 600 px wide with a filter
#' filt.red <- " -vf scale=600:-1 " #filter
#' vid.to.images2(vid.path="trout1_63_test.avi",filt = filt.red) #extract
#' 
#' #number of frames
#' fr <- length(list.files("images"))
#' #extract midline and other data
#' kin <- kin.simple(image.dir = "images",frames=c(1:fr),thr=0.6,ant.per = 0.2,show.prog=T)
#' ml <- kin$midline
#' #normalize x (y is normalized to midline by \code{kin.simple})
#' ml <- ddply(ml,.(frame),transform,x2=x-x[1])
#'
#' #compute instantaneous amplitude of tail (last/rightmost points) and wave crest x position  by frame
#' ml2 <- ddply(ml,.(frame),summarize,amp.i=abs(last(wave.y)))
#'
#' ml <- merge(ml,ml2,by="frame") #merge these
#'
#' pal <- wes_palette("Zissou1", 100, type = "continuous") #"Zissou" color palette
#' p <- ggplot(dat=ml,aes(x=x2,y=wave.y))+theme_classic(15)+scale_color_gradientn(colours = pal)
#' p <- p+geom_line(aes(group=frame,color=amp.i),stat="smooth",method = "loess", size = 1.5,alpha = 0.5)
#' print(p)
#' 
#' #delete 'images' and 'processed_images' folders
#' unlink("processed_images",recursive = T)
#' unlink("images",,recursive = T)
#'}
#'


kin.simple <-function(image.dir=NULL,frames=NULL,thr=0.7,size.min=0.05,ant.per=0.20,tips=0.02,smoothing="loess",smooth=0.2,smooth.points=200,save=TRUE,plot.pml=TRUE,image.type="orig",flip=TRUE,show.prog=FALSE){
  
  unlink("processed_images",recursive = T)
  dir.create("processed_images")
  
  proc.dir <- "processed_images"
  
  images <- paste0(image.dir,"/",list.files(image.dir)[!grepl("Icon\r",list.files(image.dir))]) #remove pesky Icon\r
  
  if(any(frames>length(images))) stop("variable 'frames' out of range of image sequence")
  if(!is.null(frames)) images <- images[frames]
  
  
  trial <- gsub("\\.[^.]*$", "", basename(images[1]))
  if(sequenced==TRUE) trial <- gsub("(.+)\\_\\d*$", "\\1", trial)
  
  kin.l <- list()
  midline.l<- list()
  classes.l <- list()
  lms <- list()
  conts <- list()
  pb = txtProgressBar(min = 0, max = length(images), initial = 0,style=3)
  
  
  roi.outs <- list() #store the rois for each image
  for(im in images){
    
    frame <- which(im==images)-1
    
    img <- EBImage::readImage(im,all=F) #if don't add package, others use "display"
    
    #EBImage::display(img,method="raster")
    
    ## computes binary mask
    y = img >thr #contrast threshold
    
    
    
    if(flip){#flip binary
      y[y==1] <- 5
      y[y==0] <- 1
      y[y==5] <- 0
    }
    z = bwlabel(y)
    
    #EBImage::display(y,method="raster")
    
    # EBImage::display(z,method="raster")
    rois <- tabulate(z)
    
    pix <- dim(z[,,1])[1]*dim(z[,,1])[2]
    w <- dim(z[,,1])[1] #width of image
    h <- dim(z[,,1])[2] #height of image
    per <- rois/(w*h) #how big are rois compared to pixel field
    
    c.roi <-  which(per>=size.min) #candidate rois, filtered by size of of pixel field
    
    names(c.roi) <- as.factor(letters[order(rois[c.roi],decreasing = T)])
    
    z.l <- list()
    out.l <- list()
    
    for(r in c.roi){
      r.name <- as.character(names(c.roi)[c.roi==r])
      z.r <- z
      z.r[z!=r] <- 0
      z.r[z==r] <- 1
      z.m <- z.r[,,1]
      z.m[1,1] <- 0 #this gets a 1 when
      z.l[[r.name]] <- z.m
      
      z.c <- ocontour(z.m)
      #EBImage::display(z.m,method="raster")
      
      wall <- any(z.c[[1]][,1]>dim(z)[1]-2 | z.c[[1]][,1]<2  |z.c[[1]][,2]>dim(z)[2]-2 | z.c[[1]][,2]<2)
      
      
      r.out <- Out(ocontour(z.m))
      if(wall ) edge <- T
      if(!wall) edge <- F
      r.out$fac <- data.frame(shape=paste0("roi-",r.name),type=paste0("roi"),edge=edge)
      out.l[[r.name]] <- r.out
      rois[c.roi[r.name]]

    }
    #don't combine if only one ROI
    if(length(out.l)==1){roi.out2 <- out.l[[1]]}else{roi.out2 <- combine(out.l)}
    
    classes <- data.table(roi=gsub("roi-","",roi.out2$fac$shape),edge=roi.out2$fac$edge,size=rois[c.roi])
    
    classes.l[[paste0(frame)]] <- data.table(frame=frame,classes)
    z.best <- z.l[[classes[edge==F,][which.max(size)]$roi]]
    r.name <-classes[edge==F,][which.max(size)]$roi
    best.class <- classes[edge==F,][which.max(size)]
    
    
    if(show.prog) {
      EBImage::display(z.best,method = "raster")
    }
    
    best.cont <- data.table(ocontour(z.best)[[1]])
    colnames(best.cont) <- c("x","y")
    
    conts[[paste0(frame)]] <- data.table(frame=frame,best.cont)
    
    y.df <- best.cont[,.(y.min=min(y),y.max=max(y),y.m=mean(y)),by=.(x)]
    setkey(y.df,"x")
    
    ends <- ceiling(nrow(y.df)*tips)
    tip.y <- mean(tail(y.df$y.m[!is.na(y.df$y.m)],ends))#tip is mean y.m of last 30 pixels
    tip.x <- mean(tail(y.df$x[!is.na(y.df$y.m)],ends))#tip is mean y.m of last 30 pixels
    
    head.y <- mean(head(y.df$y.m[!is.na(y.df$y.m)],ends))#tip is mean y.m of first 30 pixels
    head.x <- mean(head(y.df$x[!is.na(y.df$y.m)],ends))#tip is mean y.m of first 30 pixels
    
    #n midline points
    if(is.null(smooth.points)) smooth.points <- nrow(y.df)
    midline <- y.df[seq(1,nrow(y.df),length.out = smooth.points),] #two hundred points on midline
    
    midline <- midline[complete.cases(midline)]
    midline <- data.table(frame,midline)
    
    
    ####which type of lines to be fitted, spline or loess
    if(!any(c("spline","loess")==smoothing)) stop("'smoothing' must = 'loess' or 'spline'")
    
    if(smoothing=="loess")  ml.pred <- fitted(loess(midline$y.m~midline$x,span=smooth,degree=1))
    if(smoothing=="spline") ml.pred <- smooth.spline(x = midline$x,y=midline$y.m,spar=smooth)$y
    
    midline[,y.pred:=ml.pred]#add smoothed predictions
    
    #head section
    head.dat <- midline[1:(ant.per*smooth.points),]
    head.lm <- lm(y.pred~x,head.dat)
    
    head.p <- summary(head.lm)$r.squared #how well does head lm fit
    
    midline$mid.pred <- predict(head.lm,newdata=midline)#add lm prediction to midline df
    
    midline <- midline[complete.cases(midline),]
    
    midline[,wave.y:=y.pred-mid.pred] #wave y based on midline y and straight head.lm pred points
    
    midline[,roi:=r.name]
    n.roi <- paste0(basename(im),"-",r)
    
    kin.l[[paste(frame)]] <- data.table(frame,x=tip.x,y=tip.y,head.x,head.y,amp=last(midline$wave.y),head.pval=head.p,best.class)
    midline.l[[paste(frame)]] <- midline
    
    
    if(save){
      
      jpeg(paste0(proc.dir,"/",trial,"_",sprintf("%03d",frame),".jpg"),quality = 0.5)
      if(image.type=="bin")EBImage::display(z,method = "raster")
      if(image.type=="orig")EBImage:: display(img,method = "raster")
    
  
        if(plot.pml) lines(predict(lm(mid.pred~x,midline)),x=midline$x,col="blue",lwd=4)
        with(midline,lines(y.pred~x,col="red",lwd=4))
        if(plot.pml) with(midline[1:ceiling(ant.per*smooth.points),],points(x,y.pred,col="green",pch=16,cex=0.75))
      
      dev.off()
    }
    setTxtProgressBar(pb,which(images==im))
  }
  
  classes.dat <- do.call(rbind,classes.l)
  kin.dat <- do.call(rbind,kin.l)
  midline.dat <- data.table(do.call(rbind,midline.l))
  cont.dat <- do.call(rbind,conts)
  
  return(list(kin.dat=kin.dat,midline=midline.dat,cont=cont.dat,all.classes=classes.dat))
}


######### fin.kin

#' @title Tracking of fin-like extensions of body contours 

#' @description  Estimates the amplitudes of regions along a body contour that are protruding. Useful in computing paired-fin amplitudes from contour data produce from \link{kin.LDA} and \link{kin.simple}. Also computes a smoothed midline based on the body outline with the fin region removed.
#'
#' @param x a data frame or matrix with 'x' and 'y' data as columns.
#' @param fin.pos numeric, a vector of length 2 indicating the start and end of the contour region that contains the fins of interest as a proportion of the the length. 
#' @param smooth.n numeric, the number of smoothing operations undertaken by \link{coo_smooth} on the contour described by 'x'.
#' @param tip.angle the minimum angle, in degrees, that defines tip of each fin. See Details.
#' @param smoothing character, the midline smoothing method, either 'loess' or 'spline'.
#' @param ml.smooth numeric the smoothing value for the midline. If \code{smoothing} is set to 'loess', passed to 'span' value for \code{\link{loess}}. If \code{smoothing} is set to 'spline', passed to 'spar' value for \code{\link{smooth.spline}}
#' 
#' @export
#'
#' @details
#'The algorithm assumes a left-right orientation, i.e., the head of the contour is left. If otherwise oriented, contour can be flipped with \code{\link{coo_flipx}} and \code{\link{coo_flipy}} after converting contour to class \code{coo}.
#'
#'  \code{tip.angle} is used to define the tip of the fin, assuming that the tip of the fin is pointed and, for a sufficently smoothed fin contour, will have countour edges that form the highest angles within the fin region defined by \code{fin.pos}. Low values of \code{smooth.n} ($<$5) should be avoided if the contour is jagged, perhaps due to digitization.
#'  
#'  In addition to fin amplitude and contour extraction, also produces a composite contour of the body minus the fin area described by \code{fin.pos}. Fin contours are replaced by a simple linear prediction constructed from the coordinates of the first and last values covered by \code{fin.pos}. The result is a straight line between the start and end of each fin. From this composite body contour, a midline prediction is made based on the method indicated by \code{smoothing}. 
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
#' #' \code{fin.pts} a data table describing fin position consisting of the following:
#' 
#' \itemize{
#' \item  x,y coordinates of the fin tips, start, and end within the range of \code{fin.pos}.
#' \item 'ang': the angle formed by the coordinataes and their adjacent points.
#' \item 'pos': descripiton of the coordinates' positions, 'start', 'end' or 'tip'.
#' }
#' 
#' \code{comp} a data table describing the composite contour of the body minus the fins.
#' \itemize{
#' \item  x,y coordinates of the body except the range of x values withing \code{fin.pos}. These values take on a straight line described by the prediction of \code{lm()} based on the start and end of the fin. See Details.
#' \item 'y.m': the predicted y position of the midline described by the \code{ml.smooth} parameter. See Details.
#' }
#'
#' @seealso \code{\link{kin.simple}}, \code{\link{kin.LDA}}, \code{\link{efourier}}, \code{\link{coo_angle_edges}}, \code{\link{coo_smooth}}, \code{\link{loess}}, \code{\link{smooth.spline}}
#' @export
#' @import data.table
#' @import Momocs
#' @examples
#' ###plot pectoral-fin amplitudes of a swimming sunfish
#' \dontrun{
#' require(ggplot2)
#' require(data.table)
#' require(dplyr)
#' require(EBImage)
#'
#' #download example avi video
#' f <- "https://github.com/ckenaley/exampledata/blob/master/sunfish_pect.avi?raw=true"
#' download.file(f,"sunfish.avi")
#' 
#' #extract images with and reduce them to 600 px wide with a filter
#' filt.red <- " -vf scale=600:-1 " #filter
#' vid.to.images2(vid.path="sunfish.avi",filt = filt.red) #extract
#' 
#' #number of frames
#' fr <- length(list.files("images"))
#' #extract contours and other data
#' kin <- kin.simple(image.dir = "images",frames=c(1:fr),thr=0.9,ant.per = 0.25)
#' #fin amplitudes by frame with data.table
#' fin.dat <- kin$cont[, { f <- fin.kin(data.frame(x=x,y=y),fin.pos =fin.pos);list(amp=f$amp$amp,fin=f$amp$fin,amp.bl=f$amp$amp.bl)},by=.(frame)]
#' p <- ggplot(dat=fin.dat,aes(x=frame,y=amp,col=fin))+geom_line()+theme_classic(15)
#'print(p)
#'
#'#delete 'images' and 'processed_images' folders
#' unlink("processed_images",recursive = T)
#' unlink("images",,recursive = T)
#' 
#' ## plot body and fin contours of frame 1
#' cont <- data.frame(x=kin$cont[frame==1,.(x,y)]$x,y=kin$cont[frame==1,.(y)]$y)
#' fins <- fin.kin(cont,fin.pos =fin.pos)
#'
#' #for plotting
#' p.dimX <- round(range(cont$x)*c(0.9,1.1))
#' p.dimY <- round(range(cont$y)*c(0.7,1))
#' 
#' #plot body contour and fins 
#' qplot(data=fins$body,x=x,y=y)+geom_point(data=fins$fin,aes(x,y),col="red",size=3)+geom_point(data=fins$fin.pts,aes(x,y,shape=pos))+xlim(p.dimX)+ylim(p.dimY)
#' 
#' #plot body contour minus fins and the body midline
#' qplot(data=fins$comp,x=x,y=y)+geom_point(aes(x,y.m),col="red",size=2)+xlim(p.dimX)+ylim(p.dimY)
#' 
#'#midlines based on body contours minus fins
#' fin.dat <- kin$cont[, { f <- fin.kin(data.frame(x=x,y=y),fin.pos =fin.pos);list(x=f$comp$x,fin=f$comp$y,yamp.bl=f$amp$amp.bl)},by=.(frame)]
#'}
#'

fin.kin <- function(x,fin.pos=NULL,smooth.n=50,tip.ang=10,smoothing="loess",ml.smooth=0.2){
  if(is.null(fin.pos)) stop("'fin.pos' not defined")
  if(!is.matrix(x)) x.m <- as.matrix(x)
  if(is.null(colnames(x.m))) colnames(x.m) <- c("x","y")
  x.o <- Out(x.m)
  x.s <- data.table(coo_smooth(x.o,n=smooth.n)[[1]][[1]])
  
  bl <- diff(range(x$x))
  
  fin.range <-min(x)+fin.pos*bl
  
  #Left fin
  finL <- x.s[x>=fin.range[1] & x<=fin.range[2] &y<y[which.min(x)]]
  finL.o <- Out(as.matrix(finL[,.(x,y)]))
  finL[,ang:= (coo_angle_edges(finL.o)[[1]])]
  
  finL2 <- finL[seq(1,nrow(finL),5)]
  finL2.o <- Out(as.matrix(finL2[,.(x,y)]))
  finL2 <- data.table(finL2.o[[1]][[1]])
  finL2[,ang:=deg(coo_angle_edges(finL2.o)[[1]])]
  
  ptsL <- finL2[order(abs(ang)),][1:15]
  ptsL[which.min(x),pos:="start"]
  ptsL[which.max(x),pos:="end"]
  ptsL[ang>tip.ang,][which.min(y)]$pos <- "tip"
  ptsL <- ptsL[!is.na(pos)]
  
  #right fin
  
  finR <- x.s[x>=fin.range[1] & x<=fin.range[2] &y>y[which.min(x)]]
  finR.o <- Out(as.matrix(finR[,.(x,y)]))
  finR[,ang:= (coo_angle_edges(finR.o)[[1]])]
  
  #panel(x.o)
  finR2 <- finR[seq(1,nrow(finR),5)]
  finR2.o <- Out(as.matrix(finR2[,.(x,y)]))
  finR2 <- data.table(finR2.o[[1]][[1]])
  finR2[,ang:=deg(coo_angle_edges(finR2.o)[[1]])]
  
  ptsR <- finR2[order(abs(ang)),][1:15]
  ptsR[which.min(x),pos:="start"]
  ptsR[which.max(x),pos:="end"]
  ptsR[ang>tip.ang,][which.max(y)]$pos <- "tip"
  ptsR <- ptsR[!is.na(pos)]
  
  finPts <- rbind(data.table(ptsL,fin="R"),data.table(ptsR,fin="L"))
  
  fins <- rbind(data.table(finL,fin="R"),data.table(finR,fin="L"))
  
  ###add predicted y values according to lm from start to end of fins
  fins[,y.pred:=predict(lm(y~x,data.frame(y=c(y[which.min(x)],y[which.max(x)]),x=c(min(x),max(x)))),newdata=data.frame(x=x)),by=.(fin)]
  
  comp <- merge(x.s,fins,by=c("x","y"),all.x = T)
  comp[!is.na(fin),y:=y.pred]
  setkey(comp,"x")
  comp <- comp[,x:=round(x)][,.(y=range(y)),by=.(x)]
  
  if(smoothing=="loess") comp[,y.m:=predict(loess(y~x,span = ml.smooth))]
  if(smoothing=="spline") comp[,y.m:=smooth.spline(x = x,y=y,spar=ml.smooth)$y]
  
  bl <- comp[,.(sum(dist.2d(x,lead(x),y.m,lead(y.m)),na.rm=T))]$V1
  
  amp <- finPts[,.(amp=y[pos=="start"]-y[pos=="tip"],amp.bl=(y[pos=="start"]-y[pos=="tip"])/bl),by=.(fin)]
  
  return(list(body=x.s,fin=fins,fin.pts=finPts,comp=comp,bl=bl,amp=amp))
}
