#' Extracts images from a video file with ffmep
#'
#' @param vid.path Character; path of video file to be processed.
#' @param qual numeric; the quality of the jpeg images to be rendered from 1-100\%. Defaults to 50\%.
#' @return Extracts all the images of the video and saves them to an "images" directory with appended number sequence
#' @seealso \code{\link{images.to.videos}}
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
#' saveVideo(fun(),video.name="wave.mp4",interval = 0.2)
#'
#' vid.to.images(vid.path="wave.mp4",qual=100)
#'
#' #see the images in the "images" directory
#' list.files( paste0(getwd(),"/","images"))


vid.to.images <- function(vid.path=NULL,qual=50)  {

  version <-  try(system("ffmpeg -version", intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command \"", ffmpeg, "\" is not available in your system. Please install FFmpeg or avconv first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}

  qual <- round(30-30*(qual/100)+1,0)



  image.dir <- gsub(basename(vid.path),"images",vid.path)

  #delete or create the image directory

    unlink(image.dir,recursive = T)
      dir.create(image.dir)

  vid.path <- gsub("Google Drive","\"Google Drive\"",vid.path) ## remove spaces from dir if google drive
  image.dir <- paste0(gsub("Google Drive","\"Google Drive\"",image.dir),"/") #degooglize path
 video.name<- gsub(".avi","",basename(vid.path))


   system(paste0("ffmpeg -i ", vid.path, " -q:v ",qual," ", image.dir,"/",video.name,"_%5d.jpg")) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4

}

#' Stitches images into an AVI file
#'
#' @param image.dir Character; directory containing images to stich.
#' @param vid.path character; file name to be give video.
#' @param qual numeric; the quality of the video rendered from 1-100\%. Defaults to 50\%.
#' @param vid.ext chacracter; video type to output. mp4 currently works best.
#' @param frame.rate numeric; video frame rate in fps.
#' @return Outputs a video of name "video.name+vid.ext".
#' @export
#' @details Assumes images are appended with a numeric sequence.
#' @seealso \code{\link{vid.to.images}}
#' @examples
#'
#' #make some images
#' require(emojifont)
#'
#' y <- sin(1:50)
#' x <- 1:50
#' for(i in 1:50) {
#'   jpeg(paste0(getwd(),"/images/image",sprintf("%03d",i),".jpg"))
#'   plot(x[i],y[i],col="red",xlim=c(0,50),ylim=range(y),cex=0)
#'   text(x[i], y[i], labels=emoji('cow'), cex=6, col='steelblue', family='EmojiOne')
#'   dev.off()
#'   }
#'
#' images.to.video(image.dir="images",vid.name="flyingcow",frame.rate=5)
#'

images.to.video <- function(image.dir=NULL,vid.name=NULL,qual=50,vid.ext=".mp4",frame.rate=10)  {

  version <-  try(system(paste("ffmpeg -version"), intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command \"", ffmpeg, "\" is not available in your system. Please install FFmpeg or avconv first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}

  qual <- round(30-30*(qual/100)+1,0)

  if(!dir.exists(image.dir)) {
    stop("The directory \"", image.dir, "\" does not exist.")}


  vid.path <- paste0(getwd(),"/",vid.name,vid.ext)

  #delete video if it exists
  unlink(vid.path,recursive = T)

  images <- list.files(image.dir,pattern="jpg|png|tiff|jpeg|bmp")

  image.dir <- paste0(gsub("Google Drive","\"Google Drive\"",image.dir),"/") #degooglize path

  vid.path <- paste0(gsub("Google Drive","\"Google Drive\"",vid.path)) #degooglize path

  und <- grepl("_\\d+\\.\\w+\\b",images[1])


  image.name <- gsub("(\\w)\\d+\\.\\w+\\b","\\1", images[1])
  if(und) image.name <- gsub("(.+)_\\d+\\.\\w+\\b","\\1",images[1])

  num<- gsub(image.name,"",images[1])
  num <- gsub("_","",num);
  num <- unlist(strsplit(num,"[.]"))

  ext <- num[2]
  num.l <- nchar(num[1])
  num.for <- paste0("%",num.l,"d.",ext)
  if(und) num.for <- paste0("_",num.for)



    system(paste0("ffmpeg -i ", image.dir, image.name,num.for," -q:v ",qual," -r ", frame.rate," -vcodec mpeg4 ", vid.path)) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4

}

