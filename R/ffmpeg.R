#' @title Extracts images from a video file with ffmpeg
#'
#' @description Uses ffmpeg systems calls to extract images from a video.
#'
#' @param vid.path Character; path of video file to be processed.
#' @param qual numeric; the quality of the jpeg images to be rendered from 1-100\%. Defaults to 50\%.
#' @return Extracts all the images of the video and saves them to an "images" directory with appended number sequence
#' @seealso \code{\link{images.to.video}}
#' @export
#' @examples
#'
#' #make a video with animation package
#' require(animation)
#' fun <- function(){
#' y <- sin(1:50)
#' x <- 1:50
#' for(i in 1:50) {
#'   plot(x[i],y[i],col="red",xlim=c(0,50),ylim=range(y))
#'   animation::ani.pause()
#'   }
#' }
#' animation::saveVideo(fun(),video.name="wave.mp4",interval = 0.2)
#'
#' vid.to.images(vid.path="wave.mp4",qual=100)
#'
#' #see the images in the "images" directory
#' list.files( paste0(getwd(),"/","images"))


vid.to.images <- function(vid.path=NULL,qual=50)  {
  
  version <-  try(system("ffmpeg -version", intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command 'ffmpeg' is not available in your system. Please install FFmpeg first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}
   
  qual <- round(30-30*(qual/100)+1,0)

  
  image.dir <- gsub(basename(vid.path),"images",vid.path)
  
  #delete or create the image directory
 
  if(!file.exists(vid.path)) stop("'vid.path' invalid")

  unlink(image.dir,recursive = T)
  dir.create(image.dir)
  
  vid.path <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",vid.path)
  
  image.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",image.dir)
  
  video.name<- gsub(".avi","",basename(vid.path))
  
  system(paste0("ffmpeg -i ", vid.path, " -q:v ",qual," ", image.dir,"/",video.name,"_%5d.jpg")) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4
  
}

#' @title Stitches images into a video file
#' @description Stitches images into a video file of type indicated by "vid.ext"
#'
#' @param image.dir Character; directory containing images to stich.
#' @param vid.name character; file name to be give video including extension.  mp4 currently works best.
#' @param qual numeric; the quality of the video rendered from 1-100\%. Defaults to 50\%.
#' @param frame.rate numeric; video frame rate in fps.
#' @param silent logical; should output of \code{system} call for ffmpeg operation be suppressed.
#' @return Outputs a video of name "video.name+vid.ext".
#' @export
#' @details Assumes images are appended with a numeric sequence.
#' @seealso \code{\link{vid.to.images}}
#' @examples
#'
#' #make some images
#' dir.create("images") #make a directory to store images
#'
#' a <- 2
#' b <- 3
#' theta <- seq(0,10*pi,0.01)
#' r <- a + b*theta
#' df <- data.frame(x=r*cos(theta), y=r*sin(theta)) # Cartesian coords
#' every.i <- 30
#' for(i in seq(1,length(theta),30)) {
#'   jpeg(paste0(getwd(),"/images/image_",sprintf("%03d",which(i==seq(1,length(theta),30))),".jpg"))
#'   with(df[1:i,],plot(x,y,xlim=range(df$x),ylim=range(df$y),col="red"))
#'   dev.off()
#'   }
#'
#'images.to.video(image.dir=paste0(getwd(),"/images"),vid.name="spiral",frame.rate=5,qual=100,silent=FALSE)

images.to.video <- function(image.dir=NULL,vid.name=NULL,qual=50,frame.rate=10,silent=TRUE)  {
  
  version <-  try(system(paste("ffmpeg -version"), intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command 'ffmpeg' is not available in your system. Please install FFmpeg first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}
  
  qual <- round(30-30*(qual/100)+1,0)
  
  if(!dir.exists(image.dir)) {
    stop("The directory \"", image.dir, "\" does not exist.")}
  
  
  
  vid.path <- paste0(getwd(),"/",vid.name)
  #delete video if it exists
  if(file.exists(vid.path)) unlink(vid.path,recursive = T)
  
  images <- paste0(image.dir,"/",list.files(image.dir,pattern="jpg|png|tiff|jpeg|bmp"))
  
  
  vid.path <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",vid.path) #if vid.path has spaces
  
  images <- paste0(image.dir,"/",list.files(image.dir,pattern="jpg|png|tiff|jpeg|bmp",ignore.case = T))
  
  image.name <- gsub("(.+)\\_\\d*\\.\\w+$", "\\1", basename(images[1]))
  
  ext <-    gsub(".*(\\.\\w+$)", "\\1", basename(images[1]))
  
  num <- gsub(".*_(\\d+)\\.\\w+$","\\1",basename(images[1]))
  num.l <- nchar(num)
  num.for <- paste0("_%",num.l,"d",ext)
  
  image.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",image.dir) #if image.dir has spaces
  print(paste0(image.dir,"/", image.name,num.for))
  
  system(paste0("ffmpeg -i ", image.dir,"/", image.name,num.for," -q:v ",qual," -r ", frame.rate," -vcodec mpeg4 ", vid.path),ignore.stderr = silent) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4
  
  if(!file.exists(vid.name)) warning("ffmpeg failed to create video. Retry and inspect output of system call with 'silent=F'")
  
}


#' @title Extracts images from a video file with ffmep
#' @description Extract images from video file using ffmpegs flexible video filters and codecs
#'
#' @param vid.path Character; path of video file to be processed.
#' @param filt character; video filter that should be applied to ffmpeg operation. See \url{https://ffmpeg.org/ffmpeg-filters.html}
#' @param codec character; video codec to apply in ffmpeg operation
#' @param silent logical; should output of \code{system} call for ffmpeg operation be suppressed.
#' @details Particularly useful for resizing images
#' @return Extracts all the images of the video and saves them to an "images" directory with appended number sequence
#' @seealso \code{\link{images.to.video}}
#' @export
#' @examples
#' #make a video with animation package
#' fun <- function(){
#' y <- sin(1:50)
#' x <- 1:50
#' for(i in 1:50) {
#'   plot(x[i],y[i],col="red",xlim=c(0,50),ylim=range(y))
#'   animation::ani.pause()
#'   }
#' }
#' animation::saveVideo(fun(),video.name="wave.mp4",interval = 0.2)
#'
#'#reduce the image images to 200 px wide maintaining aspect ratio
#'#notice the spaces at the beginning/end of string
#'filt.red <- " -vf scale=200:-1 "
#'c <- " -c:v libx264 "
#' vid.to.images2(vid.path="wave.mp4",filt=filt.red,codec=NULL)
#'
#' #see the images in the "images" directory
#' list.files( paste0(getwd(),"/","images"))

vid.to.images2 <- function(vid.path=NULL,filt=NULL,codec=NULL,silent=TRUE)  {
  
  version <-  try(system("ffmpeg -version", intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command 'ffmpeg' is not available in your system. Please install FFmpeg first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}
  
  if(!file.exists(vid.path)) stop("'vid.path' invalid")
  
  image.dir <- gsub(basename(vid.path),"images",vid.path)
  
  #delete or create the image directory
  
  unlink(image.dir,recursive = T)
  dir.create(image.dir)
  
  #vid.path <- gsub("Google Drive","\"Google Drive\"",vid.path) ## remove spaces from dir if google drive
  vid.path <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",vid.path)
  
  #image.dir <- paste0(gsub("Google Drive","\"Google Drive\"",image.dir),"/") #degooglize path
  image.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",image.dir)
  
  video.name<- gsub(".avi","",basename(vid.path))
  video.name <- gsub(".avi","_red",video.name)
  
  if(is.null(codec)) codec <- " "
  system(paste0("ffmpeg -i ", vid.path, codec, image.dir,"/",video.name,"_%5d.jpg -y"),ignore.stderr = silent) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4
  
  if(!is.null(filt)){
    system(paste0("ffmpeg -i ", image.dir,"/",video.name,"_%5d.jpg", filt, image.dir,"/",video.name,"_%5d.jpg  -y"),ignore.stderr = silent) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4
    
  }
}

#' @title  Stitches images from video file passing filters to ffmpeg
#'
#' @description Wrapper for ffmpeg video operations. Permits flexible filtering.
#'
#' @param image.dir character; directory containing images to stich.
#' @param vid.name character; file name to be given video.
#' @param qual numeric; the quality of the video rendered from 1-100\%. Defaults to 50\%.
#' @param vid.ext chacracter; video type to output. mp4 currently works best.
#' @param frame.rate numeric; video frame rate in fps.
#' @param raw logical; encodes a raw AVI video with the "rawvideo" codec.
#' @param filt character; video filter that should be applied to ffmpeg operation. See \url{https://ffmpeg.org/ffmpeg-filters.html}.
#' @param silent logical; should output of \code{system} call for ffmpeg operation be suppressed.
#'
#' @return Outputs a video of name "video.name+vid.ext".
#' @export
#' @details Assumes images are appended with a numeric sequence beginning with "_".
#' @seealso \code{\link{vid.to.images2}}
#' @examples
#'
#' #make some spiralled images and video
#'
#' dir.create("images") #make a directory to store images
#'
#' a <- 2
#' b <- 3
#' theta <- seq(0,10*pi,0.01)
#' r <- a + b*theta
#' df <- data.frame(x=r*cos(theta), y=r*sin(theta)) # Cartesian coords
#' every.i <- 30
#' for(i in seq(1,length(theta),30)) {
#'   jpeg(paste0(getwd(),"/images/image_",sprintf("%03d",which(i==seq(1,length(theta),30))),".jpg"))
#'   with(df[1:i,],plot(x,y,xlim=range(df$x),ylim=range(df$y),col="red"))
#'   dev.off()
#'   }
#'
#' images.to.video2(image.dir=paste0(getwd(),"/images"),vid.name="spiral",frame.rate=5,qual=100,raw=FALSE)
#'

images.to.video2 <- function(image.dir=NULL,vid.name=NULL,qual=50,vid.ext=".mp4",frame.rate=10,raw=TRUE,filt=NULL,silent=TRUE)  {
  
  if(!raw)   vid.name <- paste0(vid.name,"_red")
  
  version <-  try(system(paste("ffmpeg -version"), intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command 'ffmpeg' is not available on your system. Please install ffmpeg first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}
  
  qual <- round(30-30*(qual/100)+1,0)
  
  if(!dir.exists(image.dir)) {
    stop("The directory \"", image.dir, "\" does not exist.")}
  
  vid.ext <- ".mp4"
  if(raw) vid.ext <-".avi"
  
  vid.path <- paste0(getwd(),"/",vid.name)
  #delete video if it exists
  unlink(vid.path,recursive = T)
  
  
  vid.path <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",vid.path) #if vid.path has spaces
  
  
  images <- paste0(image.dir,"/",list.files(image.dir,pattern="jpg|png|tiff|jpeg|bmp",ignore.case = T))
  
  #[here]
  
  
  image.name <- gsub("(.+)\\_\\d*\\.\\w+$", "\\1", basename(images[1]))
  
  ext <-    gsub(".*(\\.\\w+$)", "\\1", basename(images[1]))
  
  num <- gsub(".*_(\\d+)\\.\\w+$","\\1",basename(images[1]))
  num.l <- nchar(num)
  num.for <- paste0("_%",num.l,"d",ext)
  
  image.dir <- normalizePath(dirname(images[1]))
  
  image.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",image.dir) #if image.dir has spaces
  
  if(!raw) system(paste0("ffmpeg -i ", image.dir,"/", image.name,num.for," -q:v ",qual," -r ", frame.rate," -f mp4", filt," -vcodec libx264 -pix_fmt yuv420p ", vid.path,vid.ext, " -y"),ignore.stderr = silent) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4
  
  if(raw) system(paste0("ffmpeg -i ", image.dir,"/", image.name,num.for," -q:v ",qual," -r ", frame.rate," -f avi -vcodec rawvideo ", vid.path,vid.ext, " -y"),ignore.stderr = silent)
  
  message(paste0("video saved to ", vid.path,vid.ext))
}
