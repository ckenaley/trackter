#' @title  Midline detection from video data
#'
#' @description A wrapper, decomposing videos with \code{vid.to.images2} to frames for automatical retrieval of ROI midlines using \code{kin.img2}.
#'
#' @param vid.args list; arguments to be passed to \code{vid.to.images2}
#' @param kin.args list; arguments to be passed to \code{kin.img2}
#'
#'
#' @details By default, images are outputted to an "images" subdirectory in the working directory and processed images to a "processed_images" subdirectory.
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
#' @param sequenced logical, do the names of the image sequence end in a number sequence. If set to "TRUE", output video will be named the file base name of the first image minus the image number ('gsub("(.+)\\_\\d*$", "\\1", image[1])')
#' @param thr numeric, threshold to determine binary image. May require some tweaking through iteration.
#' @param plot.midline logical, value indicating if outputted images should include plotted midline and reference line based on anterior section of the ROI.
#' @param smoothing character, the midline smoothing method, either "loess" or "spline".
#' @param smooth numeric; if \code{smoothing} is set to "loess", smoothing parameter value for plotted midline.
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
#' \item "head.x" and "head.y": the x and y position of the head (leftmost or anteriormost)
#' \item "x" and ""y": the position of the tail (rightmost or posteriormost)
#' \item "amp": the amplitude (\code{amp}) of the tail
#' \item "cent.x" and "cent.y": centroid coordinate of ROI
#' \item "roi": a character indicating ROI size ('a' being the largest)
#' \item "head.pval": p values of the \code{lm()} fit that describes the position of the head as determined by \code{ant.per} (green points in the outputted images/video)}
#'
#' \code{midline} A data frame containing, for each frame described by \code{frames}, the following: \itemize{
#' \item "x" and "y.m": x and y positions of the midline of the ROI
#' \item "mid.pred": the predicted linear midline based on the points/pixels defined by \code{head.per} (green points in the outputted images/video)
#' \item "y.pred": midline points fit to a smooth spline or loess model with spar or span equal to \code{smooth} (red curve in the outputted images/video)
#' \item "wave.y": midline points "y.pred" normalized to "mid.pred"
#' \item "roi": a character indicating ROI size (a being the largest)
#' \item "cent.x": x centroid of ROI
#' \item "cent.y": y centroid of ROI
#' \item "offset.x": ROI distance from horizontal center
#' \item "offset.y": ROI distance from vertical center
#' \item "offset.total": sum of ROI offset.x and offset.y
#' \item "ar": aspect ration of the ROI
#' \item "size": size of ROI in pixels
#' }

#'  \code{head.lms}  "lm" objects, one for each frame described by \code{frames} of the linear model fit to the \code{ant.per} section of the ROI
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
#' #download example images and place in "example" subdirectory
#' f <- "https://github.com/ckenaley/exampledata/blob/master/example.zip?raw=true"
#' download.file(f,"temp.zip")
#' unzip("temp.zip")
#' unlink("temp.zip")
#'
#'
#' kin <- kin.img2(image.dir ="example",search.for = "largest",smoothing = "spline",frames=1:5,show.prog = T,thr = 0.6,make.video = T,video.name="test")
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
      #points(cent.x,cent.y,col="white",pch=16)

      y.max <- apply(z.m,1,function(x) ifelse(any(x==1),max(which(x==1)),NA))

      y.min <- apply(z.m,1,function(x) ifelse(any(x==1),min(which(x==1)),NA))

      y.df <- data.table(x=1:nrow(z.m),y.max,y.min)

      y.df <- y.df[,y.m :=(y.max-y.min)/2+y.min, by=.(x,y.min,y.max)]

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


      void.l[[r.name]] <- sum(m.pix)/length(m.pix)>0.90 #is the midline full of positive body pixels (more than 90%)

    }

    cand.kin.i <- data.table(do.call(rbind,cand.kin)) #store all candidate roi kin data
    cand.mid.i <- data.table(do.call(rbind,cand.mid))

    #store raw and best kin data

    if(search.for=="largest") prob.roi <- "a"
    off.min <- cand.kin.i[last(order(offset.tot,decreasing = F)),]$roi
    if(search.for=="offset") prob.roi <- off.min
    if(!prob.roi%in%names(void.l)[which(void.l==T)]) (warning(paste0("entire midline of probable ROI for image ", im, " does not have positive pixels. Try other 'seach.for' options")))


    if(!is.null(burn) & which(im==images)<=burn) prob.rois <- combine(prob.roi)

    kin.i <- cand.kin.i[roi==prob.roi]
    mid.i <- cand.mid.i[roi==prob.roi]
    if(show.prog) {
      EBImage::display(z.l[[prob.roi]],method = "raster")
    }
    midline.dat[[basename(im)]] <- data.frame(frame,mid.i)
    kin.dat[[basename(im)]] <- kin.i
    jpeg(paste0(proc.dir,"/",trial,"_",sprintf("%03d",frame),".jpg"),quality = 0.5)
    if(image.type=="bin")EBImage::display(z,method = "raster")
    if(image.type=="orig")EBImage:: display(img,method = "raster")
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
  if(make.video==T & !is.null(video.name)){vid.name <- video.name}else{vid.name <- paste0(trial,"_kin")}
  if(make.video) images.to.video2(proc.dir, vid.name = vid.name, qual=qual,frame.rate = frame.rate,silent=silent)

  #clean up
  if(rem.file){
    if(is.null(image.dir)) stop("'image.dir' not specified and 'rem.file=TRUE'. Won't delete working directory!")
    unlink(proc.dir,recursive = T)
    unlink(image.dir,recursive = T)
  }


  return(list(kin.dat=kin.dat,midline=midline.dat))
}

