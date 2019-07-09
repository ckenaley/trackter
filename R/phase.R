library(psd)
atan3 <- function(beta2, beta1)
{
  if (beta1 ==0 & beta2 ==0)
    return(NA)
  if (beta1 == Inf | beta1 == -Inf)
    return(NA)
  if (beta2 == Inf | beta2 == -Inf)
    return(NA)
  if (beta1 > 0)
    v <- atan(beta2/beta1);
  if(beta2 >=0 & beta1 <0)
    v <- pi + atan(beta2/beta1);
  if(beta2 <0 & beta1 <0)
    v <- -pi + atan(beta2/beta1);
  if(beta2 >0 & beta1==0)
    v <- pi/2;
  if(beta2 <0 & beta1==0)
    v <- - (pi/2);
  if (v < 0)
    v <- v + 2*pi;
  return(v)
}


fft.freq <- function(x){
  out <- pspectrum(x,verbose = F)
out$freq[which.max(out$spec)]
}

phase<- function(y,t,f){
  #if(f1!=f2) warning("f1 != f2; results may be compromised")
  # fitting procedure:
  fit <- lm(y ~ sin(2*pi*f*t)+cos(2*pi*f*t))

  a <- fit$coefficients[2]
  b<- fit$coefficients[3]
  deg(atan3(b,a))
}

Mode <- function(x) {
  ux <- na.omit(unique(x) )
  tab <- tabulate(match(x, ux)); ux[tab == max(tab) ]
}


#function to find peaks
peaks <- function(x,y){
  ft <- features(x, y,)
  x.pos <- sapply(ft$cpts, function(z)
    which.min(ceiling(abs(x - z))))
  y.pk <- y[x.pos]
  x.pk <- x[x.pos]

  return(list(y.pk=y.pk,x.pk=x.pk,curve=ft$curvature))
}

