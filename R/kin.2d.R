#' Distance between two points in Cartesian space.
#'
#' @param x1 Numeric; x position of coordinate 1
#' @param y1 numeric; y position of coordinate 1
#' @param x2 numeric; x position of coordinate 2
#' @param y2 numeric; y position of coordinate 2
#' @return A single value of the distance between p[x1,y1] and p[x2,y2]
#' @export
#' @examples
#' #Find the lengths of the sides of a tringle and print to plot
#' x <- c(0,3,2)
#' y <- c(0,3,0)
#' plot(x,y)
#' lines(x,y)
#' lines(x[c(1,3)],y[c(1,3)])
#' hyp <- dist.2d(x[1],x[2],y[1],y[2])
#' s1 <- dist.2d(x[1],x[3],y[1],y[3])
#' s2 <- dist.2d(x[2],x[3],y[2],y[3])
#' text(mean(x[1:2],mean(y[2:3])),labels=round(hyp,1))
#' text(mean(x[c(1,3)]),y[1]+0.25,labels=round(s1,1))
#' text(mean(x[c(2:3)]),mean(y[2:3]),labels=round(s2,1))

dist.2d <- function(x1,x2,y1,y2){sqrt((x2-x1)^2+(y2-y1)^2)}


#' Compute wavelength from the last half wave of sine-like waveform
#'
#' @param x Numeric; x position
#' @param y numeric; y position
#' @param p logical; should a plot of the last half wavelength
#' @return a list with half wavelength "l" and position of the crest, "p"
#' @export
#' @import features
#' @examples
#'
#' #Find length of the last half wave
#' x <- seq(0,pi,0.1)
#' y <- sin(x*pi)
#'
#' wave(x,y,p=T)
wave <- function(x,y,p=F){
  #add pantom point if last y is zero
  added <- F
  if(last(y<1e-14)){x <- c(x,last(x)+diff(tail(x,2)))
  y <- c(y,last(y)+1e-10)
  added <- T
  }
  ft <- features(x=x,y=y)
  for(n in 1:length(ft$cpts)){
    zeros <- tail(which(abs(sign(y)-lead(sign(y)))>1),2)
    pt <- x[which.min(abs(x-tail(ft$cpts,1)))]
    l <- diff(x[zeros])


    if(length(l)>0 & length(zeros)>1) {
      if(added){x <- x[-length(x)]
      y <- y[-length(y)]
      }
      if(p){plot(x,y,ylim=c(1.1*min(y),1.1*max(y)),main=paste0("wavelength= ",l,", crest position= ",pt))
        points(x[zeros[1]:zeros[2]],y[zeros[1]:zeros[2]],col="red")
        abline(h=0)

      }
      return(list(l=l,pt=pt))
    }else{
      return(list(l=NULL,pt=NULL))
    }
  }
}

#' Compute amplitude(s) and wavelength(s) of a wave form, amongst other things
#'
#' @param x Numeric; x position
#' @param y numeric; y position
#' @param p logical; should a plot of the last half wavelength
#' @return a list with amplitude "a", frequence "f", amplitude returned from a smoothed sign function "a.f", signal to noise ratio "snr", half wavelength "wave.l".
#' @export
#' @import features
#' @examples
#'
#' #Compute waveform patterns
#' x <- seq(0,pi,0.1)
#' y <- sin(x^1.3*pi)
#' plot(x,y)
#'
#' amp.freq(x=x,y=y)

amp.freq <- function(x=NULL,y,sf=100){
  s <- 1/sf
  if(is.null(x)) x <- 1:length(y)
  amp.n <- features(x=x,y=y)
  x.n <- unlist(sapply(amp.n$cpts,function(z) which.min(abs(z-x))))
  amp<- abs(diff(y[x.n])/2)
  amp.f <- abs(diff(attributes(amp.n)$fits$fn[amp.n$cpts])/2)
  snr <- fget(amp.n)$f["snr"]
  names(snr) <-NULL
  freq<- 1/(s*diff(amp.n$cpts))
  wave.l <- diff(tail(amp.n$cpts,2))
  tail.dat <-  list(a=amp,f=freq,a.f=amp.f,snr=snr,wave.l=wave.l)#a.f is amp accourding to function

  return(tail.dat)
}

