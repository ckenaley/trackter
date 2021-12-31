params <-
list(EVAL = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("trackter")

## ----eval=FALSE---------------------------------------------------------------
#      require(devtools)
#      install_github("ckenaley/trackter")
#      require(trackter)

## ----eval=FALSE---------------------------------------------------------------
#    if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#     BiocManager::install("EBImage")

## -----------------------------------------------------------------------------
library(trackter)
library(ggplot2)
library(data.table)

## -----------------------------------------------------------------------------
i <- system.file("extdata/img","sunfish_BCF.jpg",package="trackter")
y <- EBImage::readImage(i)

EBImage::display(y,method="raster")

## -----------------------------------------------------------------------------

dir <- system.file("extdata","img",package="trackter")

im <- list.files(system.file("extdata/img", package = "trackter"),full.names = TRUE)
im<- im[grepl("sunfish.jpg",im)]
im.dir <-paste0(getwd(),"/images")
dir.create(im.dir)

file.copy(im,paste0(im.dir,"/",basename(im)))

kin.y <- kin.simple(image.dir = im.dir,save=TRUE,out.dir = getwd())



## -----------------------------------------------------------------------------
print(sapply(kin.y,class))

## -----------------------------------------------------------------------------
print(kin.y$kin.dat)

## -----------------------------------------------------------------------------
print(head(kin.y$midline))
ml <- melt(kin.y$midline[,.(x,y,y.sm,wave.y)],"x")
qplot(data=ml,x=x,y=value)+facet_wrap(variable~.)

## -----------------------------------------------------------------------------
y2 <- EBImage::readImage("sunfish_000.jpg")
EBImage::display(y2,method="raster")
#clean up
unlink(im.dir,recursive=TRUE)

## ----include=FALSE------------------------------------------------------------
#clean up
unlink("sunfish_BCF_000.jpg")

## -----------------------------------------------------------------------------

#construct video file path
v <- system.file("extdata/vid","sunfish_BCF_red.avi",package="trackter")

#create new directory to store images
dir.create(paste0(getwd(),"/img2"))

#extract images from video and output to new "img2" directory
vid.to.images(v,out.dir=paste0(getwd(),"/img2"))

## -----------------------------------------------------------------------------
  kin.y2 <- kin.simple(image.dir=paste0(getwd(),"/img2"),ant.per = 0.2,save = FALSE)

#clean up
unlink(paste0(getwd(),"/img2"),recursive = TRUE)

## -----------------------------------------------------------------------------
#position
qplot(data=kin.y2$kin.dat,x=frame,y=y)

## -----------------------------------------------------------------------------
#amplitude relative to a theoretical midline established by head
qplot(data=kin.y2$kin.dat,x=frame,y=amp) 

## -----------------------------------------------------------------------------
kin.y2$midline[,x2:=x-x[1],by=frame]

#waveform plot
qplot(data=kin.y2$midline,x=x2,y=wave.y,col=frame)

## -----------------------------------------------------------------------------
w <- halfwave(x=kin.y$midline$x,y=kin.y$midline$wave.y,method="zeros")
print(w)

## -----------------------------------------------------------------------------
qplot(data=w$names,x=x,y=y,col=wave)

## -----------------------------------------------------------------------------
 wave.dat <- kin.y2$midline[, { w <- halfwave(x,wave.y,method="zeros")$dat;
 list(l=as.numeric(w$l),
 	amp=as.numeric(w$amp1),
 	pos=as.numeric(w$pos1),
 	start=as.numeric(w$wave.begin),
 	end=as.numeric(w$wave.end))},
	 by=.(frame)]

qplot(data=wave.dat,x=pos,y=l)

