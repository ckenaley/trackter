#' Extracts images from a video file with ffmep
#'
#' @param vid.path Character; path of video file to be processed.
#' @param qual numeric; the quality of the jpeg images to be rendered from 1-100\%. Defaults to 50\%.
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
#' saveVideo(fun(),video.name="wave.mp4",interval = 0.2)
#'
#' vid.to.images(vid.path="wave.mp4",qual=100)
#'
#' #see the images in the "images" directory
#' list.files( paste0(getwd(),"/","images"))


vid.to.images <- function(vid.path=NULL,qual=50)  {

  version <-  try(system("ffmpeg -version", intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command \"", ffmpeg, "\" is not available in your system. Please install FFmpeg first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}

  qual <- round(30-30*(qual/100)+1,0)



  image.dir <- gsub(basename(vid.path),"images",vid.path)

  #delete or create the image directory
  if(is.null(image.dir)) stop("'image.dir' not specified")

    unlink(image.dir,recursive = T)
      dir.create(image.dir)

  #vid.path <- gsub("Google Drive","\"Google Drive\"",vid.path) ## remove spaces from dir if google drive
  vid.path <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",vid.path)
  #image.dir <- paste0(gsub("Google Drive","\"Google Drive\"",image.dir),"/") #degooglize path
  image.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",image.dir)

   video.name<- gsub(".avi","",basename(vid.path))


   system(paste0("ffmpeg -i ", vid.path, " -q:v ",qual," ", image.dir,"/",video.name,"_%5d.jpg")) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4

}

#' Stitches images into an AVI file
#'
#' @param image.dir Character; directory containing images to stich.
#' @param vid.name character; file name to be give video.
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
#' dir.create("images") #make a directory to stor images
#' for(i in 1:50) {
#'   jpeg(paste0(getwd(),"/images/image",sprintf("%03d",i),".jpg"))
#'   plot(x[i],y[i],col="red",xlim=c(0,50),ylim=range(y),cex=0)
#'   text(x[i], y[i], labels=emoji('cow'), cex=6, col='steelblue', family='EmojiOne')
#'   dev.off()
#'   }
#'
#' images.to.video(image.dir=paste0(getwd(),"/images"),vid.name="flyingcow",frame.rate=5,qual=100)
#'

images.to.video <- function(image.dir=NULL,vid.name=NULL,qual=50,vid.ext=".mp4",frame.rate=10)  {

  version <-  try(system(paste("ffmpeg -version"), intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command \"", ffmpeg, "\" is not available in your system. Please install FFmpeg first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}

  qual <- round(30-30*(qual/100)+1,0)

  if(!dir.exists(image.dir)) {
    stop("The directory \"", image.dir, "\" does not exist.")}


  vid.path <- paste0(getwd(),"/",vid.name,vid.ext)

  #delete video if it exists
  unlink(vid.path,recursive = T)

  images <- paste0(image.dir,"/",list.files(image.dir,pattern="jpg|png|tiff|jpeg|bmp"))

  image.dir <- paste0(gsub("Google Drive","\"Google Drive\"",image.dir),"/") #degooglize path


  vid.path <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",vid.path)
  #vid.path <- paste0(gsub("Google Drive","\"Google Drive\"",vid.path)) #degooglize path

  und <- grepl("_\\d+\\.\\w+\\b",images[1])


  image.name <- gsub("(\\w)\\d+\\.\\w+\\b","\\1", basename(images[1]))
  if(und) image.name <- gsub("(.+)_\\d+\\.\\w+\\b","\\1",basename(images[1]))

  num<- gsub(image.name,"",basename(images[1]))
  num <- gsub("_","",num);
  num <- unlist(strsplit(num,"[.]"))

  ext <- num[2]
  num.l <- nchar(num[1])
  num.for <- paste0("%",num.l,"d.",ext)
  if(und) num.for <- paste0("_",num.for)



    system(paste0("ffmpeg -i ", image.dir, image.name,num.for," -q:v ",qual," -r ", frame.rate," -vcodec mpeg4 ", vid.path)) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4

}


#' Extracts images from a video file with ffmep using flexib video filters and codecs
#'
#' @param vid.path Character; path of video file to be processed.
#' @param filter character; video filter that should be applied to ffmpeg operation. See https://ffmpeg.org/ffmpeg-filters.html
#' @param codec character; video codec to apply in ffmpeg operation
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
#' saveVideo(fun(),video.name="wave.mp4",interval = 0.2)
#'
#'filt.red <- " -vf scale=200:-1 " #reduce the image images to 200 px wide maintaining aspect ratio
#'c <- " -c:v libx264 "
#' vid.to.images2(vid.path="wave.mp4",filt=filt.red,codec=NULL)
#'
#' #see the images in the "images" directory
#' list.files( paste0(getwd(),"/","images"))

vid.to.images2 <- function(vid.path=NULL,filt=NULL,codec=NULL)  {

  version <-  try(system("ffmpeg -version", intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command \"", ffmpeg, "\" is not available in your system. Please install FFmpeg first:",
            ifelse(.Platform$OS.type == "windows", "http://ffmpeg.arrozcru.org/autobuilds/",
                   "http://ffmpeg.org/download.html"))
    return()}



  image.dir <- gsub(basename(vid.path),"images",vid.path)

  #delete or create the image directory

    unlink(image.dir,recursive = T)
    dir.create(image.dir)

  #vid.path <- gsub("Google Drive","\"Google Drive\"",vid.path) ## remove spaces from dir if google drive
  vid.path <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",vid.path)

  #image.dir <- paste0(gsub("Google Drive","\"Google Drive\"",image.dir),"/") #degooglize path
  image.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",image.dir)
print(image.dir)
   video.name<- gsub(".avi","",basename(vid.path))
  video.name <- gsub(".avi","_red",video.name)

if(is.null(codec)) codec <- " "
  system(paste0("ffmpeg -i ", vid.path, codec, image.dir,"/",video.name,"_%5d.jpg -y")) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4

  if(!is.null(filt)){
  system(paste0("ffmpeg -i ", image.dir,"/",video.name,"_%5d.jpg", filt, image.dir,"/",video.name,"_%5d.jpg  -y")) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4

}
}

#' Stitches images from video file passing filters to ffmpeg
#'
#' @param image.dir Character; directory containing images to stich.
#' @param vid.name character; file name to be give video.
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
#' dir.create("images") #make a directory to stor images
#' for(i in 1:50) {
#'   jpeg(paste0(getwd(),"/images/image",sprintf("%03d",i),".jpg"))
#'   plot(x[i],y[i],col="red",xlim=c(0,50),ylim=range(y),cex=0)
#'   text(x[i], y[i], labels=emoji('cow'), cex=6, col='steelblue', family='EmojiOne')
#'   dev.off()
#'   }
#'
#' images.to.video(image.dir=paste0(getwd(),"/images"),vid.name="flyingcow",frame.rate=5,qual=100)
#'

images.to.video2 <- function(image.dir=NULL,vid.name=NULL,qual=50,vid.ext=".mp4",frame.rate=10,raw=T,filt=NULL)  {
  vid.name <- gsub(".avi","_red",vid.name)
  version <-  try(system(paste("ffmpeg -version"), intern = TRUE))
  if (inherits(version, "try-error")) {
    warning("The command \"", ffmpeg, "\" is not available in your system. Please install FFmpeg first:",
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

  images <- paste0(image.dir,"/",list.files(image.dir,pattern="jpg|png|tiff|jpeg|bmp"))

  #image.dir <- paste0(gsub("Google Drive","\"Google Drive\"",image.dir),"/") #degooglize path
  image.dir <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",image.dir)
  #vid.path <- paste0(gsub("Google Drive","\"Google Drive\"",vid.path)) #degooglize path
  vid.path <- gsub("\\/(\\w+ \\w+)\\/","/\"\\1\"/",vid.path)
  und <- grepl("_\\d+\\.\\w+\\b",images[1])


  image.name <- gsub("(\\w)\\d+\\.\\w+\\b","\\1", basename(images[1]))
  if(und) image.name <- gsub("(.+)_\\d+\\.\\w+\\b","\\1",basename(images[1]))

  num<- gsub(image.name,"",basename(images[1]))
  num <- gsub("_","",num);
  num <- unlist(strsplit(num,"[.]"))

  ext <- num[2]
  num.l <- nchar(num[1])
  num.for <- paste0("%",num.l,"d.",ext)
  if(und) num.for <- paste0("_",num.for)


  if(!raw) system(paste0("ffmpeg -i ", image.dir, image.name,num.for," -q:v ",qual," -r ", frame.rate," -f mp4", filt," -vcodec libx264 -pix_fmt yuv420p ", vid.path,vid.ext, " -y")) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4
  if(raw)
    if(raw) system(paste0("ffmpeg -i ", image.dir, image.name,num.for," -q:v ",qual," -r ", frame.rate," -f avi -vcodec rawvideo ", vid.path,vid.ext, " -y")) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4



}

## If ffmpeg is already installed, you need to uninstall it.
# brew uninstall ffmpeg

# you may very well want to specify other options (e.g. --with-faac)
# brew install ffmpeg --with-freetype

