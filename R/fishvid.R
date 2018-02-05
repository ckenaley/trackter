process.fish <- function(trial=NULL,movie=T,track=T,frames=NULL,thresh=0.7,plot.midline=T,ant.per=0.15,smooth.per=.2, image.type="orig",fps=100,flip=T,n.blob=NULL,rem.file=T){

  image.dir <- paste0(curdir,"/images/")

   trial <- gsub(".avi","",trial)

  saveVideo({
    for (i in 1:50) {
      data = data.frame(x=rnorm(1000),y=rnorm(1000))
      plot1 = ggplot(data, aes(x=x, y=y)) + geom_point()
      plot2 = ggplot(data, aes(x=y, y=x)) + geom_point()

      grid.arrange(arrangeGrob(plot1, plot2, heights=c(3/4, 1/4), ncol=1))

      ani.options(interval = 0.05, ani.dev="png", ani.height=800)
    }
  },video.name = "test_png.mp4", other.opts = "-b 1000k")




  system(paste0("ffmpeg -i ", curdir,"/",trial, ".avi -q:v 15 ", image.dir,"/",trial,"%3d.jpg")) #see https://trac.ffmpeg.org/wiki/Encode/MPEG-4

  setwd(paste0(getwd(),"/images"))
  images <- list.files()
  #images <- images[grep(paste0(fish,"_\\d+\\.jpg"),images)]
  images <- images[grep(".jpg",images)]

  if(!is.null(frames)) images <- list.files(pattern=fish)[frames]

  pb = txtProgressBar(min = 0, max = length(images), initial = 0,style=3)

  kin.dat <- list()
  for(im in images[2]){
    met.dat <- unlist(strsplit(im,"-"))
    met.dat2 <- unlist(strsplit(met.dat,"_"))
    #frame <-as.numeric(gsub("(\\d*)\\.jpg","\\1",met.dat2[7])) #fix here

    frame <- which(im==images)-1 #could cause problems
    img <- EBImage::readImage(im,all=F) #if don't add package, others use "display"
    ## computes binary mask
    y = img >thresh #contrast threshold

    if(flip){
      y[y==1] <- 5#flip binary
      y[y==0] <- 1
      y[y==5] <- 0
    }

    z = bwlabel(y)
    #EBImage::display(z,method = "raster")

    #if(length(dim(z))==3) z <- z[,,1] #if color maintained, dump it
    #may not need after filter
    roi <- which.max(tabulate(z))#find the largest roi, can tell it to find nth larger blob
    if(!is.null(n.blob)) roi <- which(tabulate(z)==tabulate(z)[order(tabulate(z),decreasing = T)[n.blob]])
    z[z!=roi] <- 0
    z[z==roi] <- 1
    EBImage::display(z,method = "raster")
    if(!track) writeImage(z,files = paste0(fish,"_",sprintf("%03d",frame),"_bin.jpg"))

    z.m <- as.matrix(imageData(z))
    z.m[1,1] <- 0 #dunno why, but this gets a 1


    y.max <- apply(z.m,1,function(x) ifelse(any(x==1),max(which(x==1)),NA))
    y.min <- apply(z.m,1,function(x) ifelse(any(x==1),min(which(x==1)),NA))
    y.df <- data.frame(x=1:nrow(z.m),y.max,y.min)
    y.df <- as.data.table(y.df) #data.table operations speed this up

    y.df <- y.df[,y.m :=(y.max-y.min)/2+y.min, by=.(x,y.min,y.max)]


    tip.y <- mean(tail(y.df$y.m[!is.na(y.df$y.m)],30))#tip is mean y.m of last 30 pixels
    tip.x <- mean(tail(y.df$x[!is.na(y.df$y.m)],30))#tip is mean y.m of last 30 pixels

    head.y <- mean(head(y.df$y.m[!is.na(y.df$y.m)],30))#tip is mean y.m of last 30 pixels
    head.x <- mean(head(y.df$x[!is.na(y.df$y.m)],30))#tip is mean y.m of last 30 pixels


    #50 midline points
    midline <- y.df[seq(1,nrow(y.df),length.out = 100),]
    midline <- midline[complete.cases(midline),]

    #head section
    head.dat <- midline[1:(ant.per*100),]
    head.lm <- lm(y.m~x,head.dat)
    head.p <- summary(head.lm)$r.squared #how well does head lm fit
    midline$mid.pred <- predict(head.lm,newdata=midline)#add lin prediction to midline df
    midline <- midline[complete.cases(midline),]
    midline$y.pred <- predict(loess(y.m~x,midline,span=smooth.per,degree=1))#add smooth predictions
    midline$amp <- with(midline,dist.2d(x,x,y.pred,mid.pred)) #amp based on last pred points
    midline$amp[midline$y.pred<midline$mid.pred] <- midline$amp[midline$y.pred<midline$mid.pred]*-1 #neg or pos amp
    kin.dat[[im]] <- data.frame(frame,x=tip.x,y=tip.y,head.x,head.y,amp=last(midline$amp),head.pval=head.p)

    plot.dat <- do.call(rbind,kin.dat)
    setTxtProgressBar(pb,which(images==im))

    jpeg(paste0(getwd(),"/processed.images/",trial,"_",sprintf("%03d",frame),".jpg"),quality = 0.5)
    if(image.type=="bin")EBImage::display(z,method = "raster")
    if(image.type=="orig")EBImage:: display(img,method = "raster")
    if(plot.midline) {
      lines(predict(lm(mid.pred~x,midline)),x=midline$x,col="blue",lwd=4)
      lines(predict(loess(y.m~x,midline,span=smooth.per,degree=1)),x=midline$x,col="red",lwd=4)
      with(head.dat,points(x,y.m,col="green",pch=16,cex=0.75))
    }
    dev.off()
  }
  if(rem.file) file.remove(im)
  kin.dat <- do.call(rbind,kin.dat)
  freqs<- fps/(diff(find_peaks(kin.dat$amp,m=4)))
  freq.m <- fps/(mean(diff(find_peaks(kin.dat$amp,m=4))))
  freq.se <- fps/(se(diff(find_peaks(kin.dat$amp,m=4))))
  return(list(kin.dat=kin.dat,head.p=head.p,freq.m=freq.m,freq.se=se(freqs),freqs=freqs))
}
