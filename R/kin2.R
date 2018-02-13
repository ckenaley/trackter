#' Retrives the midline of a detected ROI in each frame of a video.
#'
#' @param vid.path Character; path of video to be analyze.
#' @param frames numeric vectors indicating which video frames to process
#' @param thr numeric, threshhold to determine binary image. May require some tweeking through iteration
#' @param plot.midline logical value indicating if outputted images should inclued plotted midline and reference line.
#' @param show.prog logical value indicating if outputted image should be display during analysis.
#' @param smooth numeric; smoothing parameter value for plotted midline
#' @param image.type character; the type of image to be outputted
#' @param n.blob numeric, indicating which nth largest ROI is the ROI to be analyzed. May require tweeking through interation
#' #' @param make.video logical value indicating if a video should be saved of midline position overlaying origina frames
#' @param qual numeric; quality of the outputted video from 1-100\%. Defaults to 50\%.
#' @param ant.per numeric; left-most percentage of ROI that establishes the vertical reference for the midline displacement.
#' @param frame.rate numeric; outputted video frame rate in fps.
#' @param rem.file logical value indicating if the outputted images, both from the original video and images with midline overlay, should be deleted. Default is "TRUE".
#'
#' @details
#' By default, images are outputted to an "images" subdirectory in the working directory.
#'
#'\code{image.type} Can be set as "orig" or "bin". "orig" plots midline and reference lines over the orginal video frames, "bin" over binary images.
#'\code{n.blob} May be useful if there are other highly contrasted ROIs in the frame.
#'
#'\code{make.video} If "TRUE" a video of the same names as \code{video.name} is outputted in the working directory.
#'
#'\code{rem.file} If "TRUE", \code{make.video} is also "TRUE" a video of processed images is still produced.
#
#' @return A list with the following components:
#'
#' \code{kin.dat} a data frame consisting of position paramters for the ROI indicated by \code{n.blob}:\itemize{
#' \item the frame number
#' \item "head.x" and "head.y": the x and y position of the head (leftmost or anteriormost)
#' \item "x" and ""y": the position of the tail (rightmost or posteriormost)
#' \item "amp": the amplitude (\code{amp}) of the tail
#' \item "head.pval": p values of the \code{lm()} fit that describes the position of the head as determined by \code{ant.per} (green points in the outputted images/video)}
#'
#' \code{midline} A data frame containing, for each frame described by \code{frames}, the following: \itemize{
#' \item "x" and "y.m": x and y positions of the midline of the ROI
#' \item "mid.pred": the predicted linear midline based on the points/pixels defined by \code{head.per} (green points in the outputted images/video)
#' \item "y.pred": midline points fit to a loess smoothing model with span equal to \code{smooth} (red curve in the outputted images/video)
#' \item "wave.y": midline points "y.pred" normalized to "mid.pred"}
#'
#'  \code{head.lms}  "lm" objects, one for each frame desrbied by \code{frames} of the linear model fit to the \code{ant.per} section of the ROI
#' @details Chooses ROIs that are big (>5\% of the pixel field) and identifies the one with the largest variance to the trailing edge amplitude (i.e., assumes that ROI is the one moving)
#' @seealso \code{\link{kin.img}}
#' @export
#' @import data.table
#' @examples
#' #produce a classic midline waveform plot of swimming fish
#'
#' require(wesanderson)
#' require(ggplot2)
#' require(plyr)
#'
#' #download an example video (7.5 MB) and place in working directory
#' f <- "https://github.com/ckenaley/exampledata/blob/master/trout1_63_test.avi?raw=true"
#' download.file(f, file.path(getwd(), "trout1_63_test.avi"), method = "libcurl")
#' kin <- kin.vid(vid.path ="trout1_63_test.avi",thr=0.7,frames=1:20,frame.rate=10)
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
#'

kin.vid <-function(vid.path=NULL,frames=NULL,thr=0.7,plot.midline=TRUE, show.prog=FALSE, ant.per=0.15,smooth=.2, image.type="orig",flip=T,n.blob=NULL,rem.file=TRUE,make.video=T,qual=50,frame.rate=10){

  if(!file.exists(vid.path)) stop(paste0(vid.path," does not exist in ", "getwd()"))

  oldwd <- getwd()
  trial <- gsub(".avi","",vid.path)

  for(i in 1:2) vid.to.images(vid.path=vid.path,qual = 20)
  image.dir <- paste0(getwd(),"/images")
  kin.img2(image.dir = image.dir,frames,thr,plot.midline, show.prog, ant.per,smooth, image.type,flip,n.blob,rem.file,make.video,qual,frame.rate)
}

######### kin.img

#' Retrives the midline of an ROI from an image sequence.
#'
#' @param image.dir Directory containing images to analyze.
#' @param frames numeric vectors indicating which images to process.
#' @param thr numeric, threshhold to determine binary image. May require some tweeking through iteration.
#' @param plot.midline logical value indicating if outputted images should inclued plotted midline and reference line.
#' @param smooth numeric; smoothing parameter value for plotted midline
#' @param image.type character; the type of image to be outputted.
#' @param flip logical, indicating if binary should be flipped.
#' @param show.prog logical value indicating if outputted image should be displayed during analysis.
#' @param n.blob numeric, indicating which nth largest ROI is the ROI to be analyzed. May require tweeking through interation. Perhaps best to let the function choose.
#' @param make.video logical value indicating if a video should be saved of midline position overlaying origina frames.
#' @param qual numeric; quality of the outputted video from 1-100\%. Defaults to 50\%.
#' @param ant.per numeric; left-most percentage of ROI that establishes the vertical reference for the midline displacement.
#' @param frame.rate numeric; outputted video frame rate in fps.
#' @param rem.file logical value indicating if the outputted images, both from the original video and images with midline overlay, should be deleted.
#' @export
#' @details
#' By default, images are outputted to the \code{image.dir} subdirectory in the working directory. Chooses ROIs that are big (>5\% of the pixel field) and identifies the one with the largest variance to the trailing edge amplitude (i.e., assumes that ROI is the one moving).
#'
#'\code{image.type} Can be set as "orig" or "bin". "orig" plots midline and reference lines over the orginal video frames, "bin" over binary images.
#'\code{n.blob} May be useful if there are other highly contrasted ROIs in the frame.
#'
#'\code{make.video} If "TRUE" a vidoe of the same names as \code{video.name} is outputted in a \code{image.dir} subdirectory.
#'
#'\code{rem.file} If "TRUE", \code{make.video} is also "TRUE" a video of processed images is still produced.
#
#' @return A list with the following components:
#'
#'\code{kin.dat} a data frame consisting of position paramters for the ROI indicated by \code{n.blob}:\itemize{
#' \item the frame number
#' \item "head.x" and "head.y": the x and y position of the head (leftmost or anteriormost)
#' \item "x" and ""y": the position of the tail (rightmost or posteriormost)
#' \item "amp": the amplitude (\code{amp}) of the tail
#' \item "head.pval": p values of the \code{lm()} fit that describes the position of the head as determined by \code{ant.per} (green points in the outputted images/video)}
#'
#' \code{midline} A data frame containing, for each frame described by \code{frames}, the following: \itemize{
#' \item "x" and "y.m": x and y positions of the midline of the ROI
#' \item "mid.pred": the predicted linear midline based on the points/pixels defined by \code{head.per} (green points in the outputted images/video)
#' \item "y.pred": midline points fit to a loess smoothing model with span equal to \code{smooth} (red curve in the outputted images/video)
#' \item "wave.y": midline points "y.pred" normalized to "mid.pred"}
#'
#'  \code{head.lms}  "lm" objects, one for each frame desribed by \code{frames} of the linear model fit to the \code{ant.per} section of the ROI
#' @seealso \code{\link{kin.vid}}
#' @export
#' @examples
#' #produce a classic midline waveform plot of swimming fish
#' require(wesanderson)
#' require(ggplot2)
#' require(data.table)
#' require(dplyr)
#' require(EBImage)
#'
#'
#' #download example images and place in "example" subdirectory
#' f <- "https://github.com/ckenaley/exampledata/blob/master/example.zip?raw=true"
#' download.file(f,"temp.zip")
#' unzip("temp.zip")
#' unlink("temp.zip")
#' kin <- kin.img2(image.dir=paste0(getwd(),"/example"),thr=0.7,frames=1:10,smooth=0.5)
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
#'
frames <- 1:10
kin.img2 <-function(image.dir=NULL,frames=NULL,thr=0.7,plot.midline=TRUE, show.prog=FALSE,ant.per=0.15,smooth=.2, image.type="orig",flip=TRUE,rem.file=TRUE,make.video=TRUE,qual=50,frame.rate=10,n.blob=NULL){

  unlink("processed_images",recursive = T)
  dir.create("processed_images")

  proc.dir <- "processed_images"

  images <- paste0(image.dir,"/",list.files(image.dir,pattern="jpg"))

  if(any(frames>length(images))) stop("variable 'frames' out of range of image sequence")
  if(!is.null(frames)) images <- images[frames]

  pb = txtProgressBar(min = 0, max = length(images), initial = 0,style=3)

  und <- grepl("_\\d+\\.\\w+\\b",images[1])

  trial <- gsub("(\\w)\\d+\\.\\w+\\b","\\1", basename(images[1]))
  if(und) trial<- gsub("(.+)_\\d+\\.\\w+\\b","\\1",basename(images[1]))


  kin.dat <- list()
  kin.dat.raw <- list()
  midline.dat.raw <- list()
  midline.dat <- list()
  lms <- list()
  for(im in images){

    frame <- which(im==images)-1 #could cause problems
    img <- EBImage::readImage(im,all=F) #if don't add package, others use "display"
    ## computes binary mask
    y = img >thr #contrast threshold
    if(flip){
      y[y==1] <- 5#flip binary
      y[y==0] <- 1
      y[y==5] <- 0
    }
    z = bwlabel(y)


    rois <- tabulate(z)

    roi <- which.max(rois)#find the largest roi, can tell it to find nth larger blob

    if(!is.null(n.blob)) roi <- which(rois==rois[order(rois,decreasing =T)[n.blob]])

    pix <- dim(z[,,1])[1]*dim(z[,,1])[2]
  per <- rois/pix #how big are rois compared to pixel field

 c.roi <-  which(per>=0.05) #candidate rois, filtered by 5 % of pixel field

names(c.roi) <- as.factor(letters[order(rois[c.roi],decreasing = T)])

 kin.burn <- NULL
 if(which(im==images)>2){
   kin.burn <- do.call(rbind,kin.dat.raw)
   rownames(kin.burn) <- NULL
   amp.var <- ddply(kin.burn,.(roi),summarize,amp.v=var(amp))
 }else{
   amp.var <- data.frame(roi=as.character("a"),amp.v=1) #assume largest ROI on first frame
 }

 cand.kin <- list()
 for(r in c.roi){
   r.name <- as.character(names(c.roi)[c.roi==r])
       z.r <- z
    z.r[z!=r] <- 0
    z.r[z==r] <- 1
    z.m <- z.r[,,1]
    z.m[1,1] <- 0 #dunno why, but this gets a 1
    if(show.prog) EBImage::display(z.m, method="raster")

    #centroids
    cent.x <- mean(unlist(apply(z.m,1,function(x) which(x==1))))
    cent.y <- mean(unlist(apply(z.m,2,function(x) which(x==1))))

    #points(cent.x,cent.y,col="white",pch=16)

    y.max <- apply(z.m,1,function(x) ifelse(any(x==1),max(which(x==1)),NA))
    y.min <- apply(z.m,1,function(x) ifelse(any(x==1),min(which(x==1)),NA))

    y.df <- data.frame(x=1:nrow(z.m),y.max,y.min)
    y.df <- as.data.table(y.df) #data.table operations speed this up
    y.df <- y.df[,y.m :=(y.max-y.min)/2+y.min, by=.(x,y.min,y.max)]

    tip.y <- mean(tail(y.df$y.m[!is.na(y.df$y.m)],30))#tip is mean y.m of last 30 pixels
    tip.x <- mean(tail(y.df$x[!is.na(y.df$y.m)],30))#tip is mean y.m of last 30 pixels

    head.y <- mean(head(y.df$y.m[!is.na(y.df$y.m)],30))#tip is mean y.m of last 30 pixels
    head.x <- mean(head(y.df$x[!is.na(y.df$y.m)],30))#tip is mean y.m of last 30 pixels

    #n midline points
    midline <- y.df[seq(1,nrow(y.df),length.out = 100),] #hundred points on midline
    midline <- midline[complete.cases(midline),]

    #head section
    head.dat <- midline[1:(ant.per*100),]
    head.lm <- lm(y.m~x,head.dat)
    head.p <- summary(head.lm)$r.squared #how well does head lm fit
    midline$mid.pred <- predict(head.lm,newdata=midline)#add lm prediction to midline df
    midline <- midline[complete.cases(midline),]
    midline$y.pred <- predict(loess(y.m~x,midline,span=smooth,degree=1))#add smoothed predictions
    midline$wave.y <- with(midline,dist.2d(x,x,y.pred,mid.pred)) #wave y based on  pred points
    midline$wave.y[midline$y.pred<midline$mid.pred] <- midline$wave.y[midline$y.pred<midline$mid.pred]*-1 #neg or pos amp

    n.roi <- paste0(basename(im),"-",r)
    cand.kin[[r]] <- data.frame(frame,x=tip.x,y=tip.y,head.x,head.y,amp=last
                                (midline$wave.y),head.pval=head.p,roi=r.name,cent.x,cent.y)



    m.var.i <- amp.var$roi[which.max(amp.var$amp.v)]
    if(r.name==m.var.i){
      midline.dat[[basename(im)]] <- data.frame(frame,midline,roi=r.name)
      lms[[basename(im)]] <- head.lm
    jpeg(paste0(proc.dir,"/",trial,"_",sprintf("%03d",frame),".jpg"),quality = 0.5)
    if(image.type=="bin")EBImage::display(z,method = "raster")
    if(image.type=="orig")EBImage:: display(img,method = "raster")
    if(plot.midline) {
      lines(predict(lm(mid.pred~x,midline)),x=midline$x,col="blue",lwd=4)
      lines(predict(loess(y.m~x,midline,span=smooth,degree=1)),x=midline$x,col="red",lwd=4)
      with(head.dat,points(x,y.m,col="green",pch=16,cex=0.75))
    }
    dev.off()
    }
 }
 #store raw and best kin and midline data
 max.var <- as.character(amp.var$roi[which.max(amp.var$amp.v)])
 cand.kin.i <- do.call(rbind,cand.kin)
 rownames(cand.kin.i) <- NULL
 cand.kin[[basename(im)]] <- cand.kin.i
 kin.dat.raw[[basename(im)]] <- do.call(rbind,cand.kin)
 kin.dat[[basename(im)]] <- cand.kin.i[cand.kin.i$roi==max.var,]
 setTxtProgressBar(pb,which(images==im))
  }


  kin.dat <- do.call(rbind,kin.dat)
  midline.dat <- do.call(rbind,midline.dat)
  if(make.video) images.to.video(image.dir =proc.dir, vid.name = trial, qual=qual,frame.rate = frame.rate)

  #clean up
  if(rem.file){

    unlink(proc.dir,recursive = T)
    unlink(image.dir,recursive = T)

  }
  return(list(kin.dat=kin.dat,midline=midline.dat,head.lms=lms))

}
