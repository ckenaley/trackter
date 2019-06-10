#' @title Compute wavelength from a sine-like waveform
#' @description computes wavelength from a sine-like waveform based on either peak-to-peak or internodal distance
#'
#' @param x Numeric; x position
#' @param y numeric; y position
#' @param p logical; should the data be plotted with indicated waves
#' @param method character; how should wavelength be calculated, where it crosses zero ("zeros") or peak to peak ("peaks"). "peaks" method calculates full wavelength, "zeros" the half wavelength.
#' @param zero.begin logical; does wave begin at zero? Default is 'TRUE' and will help find wave beginning at first x,y values if y=0
#' @return a list with half wavelength "l" and position of the crest, "p"
#' @export
#' @import features
#' @import ggplot2
#' @import data.table
#' @examples
#'
#' #Find length of the last half wave
#' x <- seq(0,pi,0.1)
#' y <- sin(x*pi)
#'
#' wave(x,y,p=T)
#'
#'
wave <-
  function(x,
           y,
           p = F,
           method = "zeros",
           zero.begin=T

        ) {

     if (!method %in% c("peaks", "zeros"))  stop("method must be set to 'peaks' or 'zeros' (the default)")
    x = c(unlist(x))
    y = c(unlist(y))

    added <- F
    #is first value y=0?
    if(y[1]==0){x.0 <- x[1]-diff(x[1:2])
    y.0 <- y[1]-diff(y[1:2])
    x <- c(x.0,x)
    y <- c(y.0,y)
    added <- T
    }


    dt <- data.table(x, y)


#find peaks and then lengths


    if (method == "peaks") {
      ft <- features(x, y)
      x.pos <- sapply(ft$cpts, function(z)
        which.min(abs(x - z)))
      y.pk <- y[x.pos]
      x.pk <- x[x.pos]


      wave.dat <- data.table(
        zeros = 0,
        start = x[x.pos],
        end = x[lead(x.pos, lead.wave)],
        begin.index = x.pos,
        end.index = lead(x.pos, lead.wave),
        wave = as.character(1:length(lead(x.pos, lead.wave))),
        l = x[lead(x.pos, lead.wave)] - x[x.pos]
      )

      wave.dat[, amp :=y[x.pos]]
      wave.dat[, pos := start+ (end - start) / 2]

      if (all(is.na(wave.dat$l)))  warning("No full waves found in data. Try to method='zeros' to find half waves")

}


    if (method == "zeros") {
      y.signs <- sign(dt$y)
      if(zero.begin) y.signs[y.signs==0] <- -1

      z <- which(abs(diff(y.signs))>1)



      names(x) <- NULL
      wave.dat <- data.table(zeros = z, start = x[z])

      if (length(z) > 1) {
        wave.dat[, c("end", "begin.index", "end.index", "wave", "l") := list(lead(start),
                                                                             zeros,
                                                                             lead(zeros),
                                                                             as.character(1:.N),
                                                                             lead(start) - start
                                                                             )
                                                                             ]


        }
    }

    wave.dat <-  wave.dat[!is.na(l)]

    if (nrow(wave.dat) != 0) {
      if (all(is.na(wave.dat$l))) {
        wave.dat2 <- dt[, wave := NA]
      } else{
        wave.dat2 <- wave.dat[!is.na(l), .(x = x[begin.index:end.index]), by = .(wave)]
        wave.dat2 <- merge(dt, wave.dat2, all.x = T,by="x")
      }

      pks <-wave.dat2[!is.na(wave), .(amp = y[which.max(abs(y))], pos = as.numeric(x[which.max(abs(y))])), by = wave]

      if (method != "peaks" & nrow(pks)>0){
        wave.dat <- merge(wave.dat, pks, by = "wave")
      wave.dat[, l2 := lead(pos) - pos]
      wave.dat[, pos2 := pos + (lead(pos) - pos) / 2]

      if (p) {
        p1 <- qplot(data = wave.dat2,
              x,
              y,
              col = wave,
              alpha = 0.5) + theme_light(15)
        print(p1)
      }
      }
    } else{
      wave.dat <- data.table(
        zeros = 0,
        start = 0,
        end = 0,
        begin.index = 0,
        end.index = 0,
        wave = 0,
        l = as.numeric(0),
        pos = 0,
        amp = 0,
        l2 = 0,
        pos2 = 0
      )
      wave.dat2 <- dt[, wave := NA]

      if (p) {
        p1 <-   qplot(data = wave.dat2, x, y) + theme_light(15)
        print(p1)
      }
    }

    #zap dummy data to find wave starting at y=0

    if(added){
      wave.dat2 <- wave.dat2[-1,]
      wave.dat2$y <- wave.dat2$y-1e-17
      wave.dat[,zeros:=zeros-1]
      wave.dat[,begin.index:=begin.index-1]
      wave.dat[,end.index:=end.index-1]
    }

    return(list(
      method=method,
      names = wave.dat2,
      dat = wave.dat
    ))
  }

#' @title Computes amplitude and wavelength of wave-like data
#' @description Computes amplitude(s) and wavelength(s) of a wave form, amongst other things, based on a sampling frequency
#'
#' @param x Numeric; x position (or sample number)
#' @param y numeric; y position
#' @param sf numeric; sample frequency (i.e., how often was x and y sampled) in Hz
#' @return a list with amplitude "a", frequence "f", amplitude returned from a smoothed sign function "a.f" based on output from \code{features}, signal to noise ratio "snr".
#' @export
#' @import features
#'@seealso \code{features}
#' @examples
#'
#' #Compute waveform patterns
#' x <- seq(0,pi,0.1)
#' y <- sin(x^1.3*pi)
#' plot(x,y)
#'
#' amp.freq(x=x,y=y)

amp.freq <- function(x = NULL, y, sf = 100) {
  s <- 1 / sf
  if (is.null(x))
    x <- 1:length(y)
  amp.n <- features(x = x, y = y)
  x.n <- unlist(sapply(amp.n$cpts, function(z)
    which.min(abs(z - x))))
  amp <- abs(diff(y[x.n]) / 2)
  amp.f <- abs(diff(attributes(amp.n)$fits$fn[x.n]) / 2)
  snr <- fget(amp.n)$f["snr"]
  names(snr) <- NULL
  freq <-
    1 / (s * diff(amp.n$cpts[seq(1, length(amp.n$cpts), 2)])) #peak to peak or trough to trough freq
  tail.dat <-
    list(a = amp,
         f = freq,
         a.f = amp.f,
         snr = snr)#a.f is amp according to function

  return(tail.dat)
}
