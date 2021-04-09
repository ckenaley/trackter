## ---- include = FALSE------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", out.width = "50%",
  message = FALSE, warning = FALSE, error = FALSE,fig.path=""
)


## ----eval=FALSE------------------------------------------------------------------------
## install.packages("trackter")


## ----eval=FALSE------------------------------------------------------------------------
##     require(devtools)
##     install_github("ckenaley/trackter")
##     require(trackter)


## ----eval=FALSE------------------------------------------------------------------------
##   if (!requireNamespace("BiocManager", quietly = TRUE))
##    install.packages("BiocManager")
##    BiocManager::install("EBImage")


## ----setup-----------------------------------------------------------------------------
library(trackter)
library(ggplot2)
library(data.table)


## ----sunfishimage,fig.height=2.5-------------------------------------------------------

i <- system.file("extdata/img","sunfish_BCF.jpg",package="trackter")
y <- EBImage::readImage(i)


EBImage::display(y,method="raster")




## ----results=FALSE---------------------------------------------------------------------

dir <- system.file("extdata","img",package="trackter")

kin.y <- kin.simple(image.dir = dir,out.dir = getwd())




## ----kiny print------------------------------------------------------------------------
print(sapply(kin.y,class))


## ----kinykindat------------------------------------------------------------------------
print(kin.y$kin.dat)


## ----midline, fig.pos="center",fig.width=5---------------------------------------------

print(kin.y$midline)
ml <- melt(kin.y$midline[,.(x,y.m,y.pred,wave.y)],"x")
qplot(data=ml,x=x,y=value)+facet_wrap(variable~.)




## ----sunfishoverlay, fig.pos="center",fig.width=5--------------------------------------

y2 <- EBImage::readImage("sunfish_BCF_000.jpg")
EBImage::display(y2,method="raster")

#clean up
unlink("sunfish_BCF_000.jpg")



## ----vidtoimages,results=FALSE---------------------------------------------------------

#construct video file path
v <- system.file("extdata/vid","sunfish_BCF_red.avi",package="trackter")

#create new directory to store images
dir.create(paste0(getwd(),"/img2"))

#extract images from video and output to new "img2" directory
vid.to.images(v,out.dir=paste0(getwd(),"/img2"))  


## ----kinsimple-------------------------------------------------------------------------
  kin.y2 <- kin.simple(image.dir =paste0(getwd(),"/img2"),thr=0.6,ant.per = 0.2,save = FALSE)

#clean up
unlink(paste0(getwd(),"/img2"),recursive = TRUE)



## ----posamp,fig.pos="center",fig.width=5-----------------------------------------------
qplot(data=kin.y2$kin.dat,x=frame,y=y) #position
qplot(data=kin.y2$kin.dat,x=frame,y=amp) #amplitude relative to a theoretical midline established by head


## ----midline2, fig.pos="center",fig.width=5--------------------------------------------
kin.y2$midline[,x2:=x-x[1],by=frame]

#wave form plot
qplot(data=kin.y2$midline,x=x2,y=wave.y,col=frame)




## ----wave------------------------------------------------------------------------------
w <- halfwave(x=kin.y$midline$x,y=kin.y$midline$wave.y,method="zeros")
print(w)


## ----wavenumber, fig.pos="center",fig.width=5------------------------------------------
qplot(data=w$names,x=x,y=y,col=wave)


## ----wavelengthplot, fig.pos="center",fig.width=5--------------------------------------
 wave.dat <- kin.y2$midline[, { w <- halfwave(x,wave.y,method="zeros")$dat;
 list(l=as.numeric(w$l),
 	amp=as.numeric(w$amp1),
 	pos=as.numeric(w$pos1),
 	start=as.numeric(w$wave.begin),
 	end=as.numeric(w$wave.end))},
	 by=.(frame)]

qplot(data=wave.dat,x=pos,y=l)

