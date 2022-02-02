#' @title Extracts images from a video file
#'
#' @description Uses ffmpeg through the \code{av} package to extract images from a video.
#'
#' @param vid.path character; path of video file to be processed.
#' @param out.dir character; directory path in which to store images.
#' @param format character; image format such as "png" or "jpg"; must be available from \code{\link[av]{av_encoders}}.
#' @param vfilter character; a string defining an ffmpeg filter, the \code{-vf} argument in the ffmpeg command line utility. Passed to \code{\link[av]{av_encode_video}}
#' @param overwrite logical; should images in path described by 'out.dir' be overwritten if they exists.
#' @param ... other arguments to be passed to \code{\link[av]{av_encode_video}}.
#' @return Extracts all the images of the video and saves them to an "images" directory with appended number sequence.
#' @seealso \code{\link{images.to.video}}, \code{\link[av]{av_encode_video}} \code{\link[av]{av_video_images}}.
#' @export
#' @import av
#' 
#' @examples
#'
#' #access video that loads with package
#' v <-system.file("extdata/vid", "sunfish_BCF.avi", package = "trackter")
#'
#' #create directory in which to store images
#' dir.create(paste0(tempdir(),"/images"))
#' 
#' vid.to.images(vid.path=v,out.dir= paste0(tempdir(),"/images"))
#'
#' #see the images in the "images" subdirectory
#' list.files( paste0(tempdir(),"/images"))
#' 
#' #clean up
#' unlink(paste0(tempdir(),"/images"),recursive=TRUE)
#' 


vid.to.images <- function(vid.path=NULL,out.dir=NULL,format="jpg",vfilter=NULL,overwrite=FALSE,...)  {

  if(is.null(out.dir)|length(out.dir)==0) stop("'out.dir' not specified.")
  if(is.null(vid.path)) stop("'vid.path' not specified.")
  
  if(!file.exists(out.dir)) stop("'out.dir' (", paste0(out.dir),") does not exist")
  
  if(!file.exists(vid.path)) stop("'vid.path' (", paste0(vid.path),") does not exist")
  
  vbn<- gsub("\\..*$","",basename(vid.path))
  
  framerate <- av_media_info(vid.path)$video$framerate
  
  sz <- av_media_info(vid.path)$video[,c("width","height")]
  
  codec <- av_media_info(vid.path)$video$codec


  if(is.null(vfilter)) vfilter <-  "null"
  
  tmp <- tempfile()
  dir.create(tmp)
  output <- file.path(tmp, paste0(vbn,"_%6d.", format))
  
  nframe <-  av_media_info(vid.path)$video$frames
  
  av_encode_video(input = vid.path, output = output, vfilter = vfilter,codec = codec,...)
  
  imgs <- list.files(tmp,full.names = TRUE,pattern=paste0(vbn,"_\\d+.", format))
  
  if(nframe<length(imgs)) {imgs <- imgs[1:nframe] 
  message("Extra audio frame dropped from output")
  }
  imgs2 <- file.path(out.dir,basename(imgs))
  invisible(file.copy(imgs,imgs2, overwrite=overwrite))
  
}

#' @title Stitches images into a video file
#' @description Uses ffmpeg through the \code{av} package to stitch images into a video.
#'
#' @param image.dir character; directory containing images to stitch.
#' @param out.dir character; directory in which to store video.
#' @param vid.name character; file name given to video including extension. 
#' @param overwrite logical; should path described by \code{vid.name} be overwritten if it exists. 
#' @param ... other arguments to be passed to \code{av::av_encode_video}.
#' @return Outputs a video of name "video.name+vid.ext".
#' @export
#' @import av
#' @details Assumes images are appended with a numeric sequence.
#' @seealso \code{\link{vid.to.images}}, \code{\link[av]{encoding}} \code{\link[av]{av_video_images}}.
#' @examples
#'
#' #make some images
#' \donttest{
#' dir.create(paste0(tempdir(),"/images")) #make a directory to store images
#'
#' a <- 2
#' b <- 3
#' theta <- seq(0,10*pi,0.01)
#' r <- a + b*theta
#' df <- data.frame(x=r*cos(theta), y=r*sin(theta)) # Cartesian coords
#' every.i <- 30
#' for(i in seq(1,length(theta),30)) {
#'   jpeg(paste0(tempdir(),"/images/image_",sprintf("%03d",which(i==seq(1,length(theta),30))),".jpg"))
#'   with(df[1:i,],plot(x,y,xlim=range(df$x),ylim=range(df$y),col="red"))
#'   dev.off()
#'   }
#'
#'images.to.video(image.dir=paste0(tempdir(),"/images"),vid.name="spiral.mp4",out.dir=tempdir())
#'
#'file.exists(paste0(tempdir(),"/spiral.mp4"))
#'
#' #clean up
#' unlink(paste0(tempdir(),"/spiral.mp4"))
#' unlink(paste0(tempdir(),"/images"),recursive=TRUE)
#' }

images.to.video <- function(image.dir=NULL,out.dir=NULL,vid.name=NULL,overwrite=FALSE,...)  {
  
  if(!dir.exists(image.dir)) stop("'image.dir' (", paste0(image.dir),") does not exist")
 
  if(is.null(out.dir)) stop("'out.dir' not specified.")

  if(!file.exists(out.dir)) stop("'out.dir' (", paste0(out.dir),") does not exist")
  
  vid.path <- file.path(out.dir,vid.name)

  if(file.exists(vid.path) & overwrite==FALSE) stop("video with name 'vid.name'  exist in 'out.dir' directory. To save file of this name in this location, 'overwrite' must be 'TRUE'")
  
  if(file.exists(vid.path) & overwrite==TRUE) unlink(vid.path,recursive = TRUE)
  
  images <- list.files(image.dir,pattern="jpg|png|tiff|jpeg|bmp",full.names = TRUE)

  av_encode_video(input = images,output = vid.path,...)
  
  if(!file.exists(vid.path)) warning("failed to create video'")
  
  if(file.exists(vid.path)) message("successfully created video \'", paste0(vid.path),"\'")
  
}

