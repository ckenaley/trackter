
#' @title  Midline estimations for closed contours
#' @description Internal functions used in \code{kin} functions for calculating the long-axis midlline spanning a closed contour. 
#' 
#' \code{free.ml.ang} estimates a midline by finding angular tips that lie in end regions of the long axis of a contour and bisects the contour according to these tips.
#' 
#' \code{free.ml.hull} estimates a midline by finding the most distant coordinates in the convex hull and then bisects the contour according to these tips.
#' 
#' \code{free.ml.del} estimates a midline through Delaunay triangulation and Voronoi tesselation.
#'
#' @param out a matrix of named x,y values describing a closed outline (contour)
#' @param smooth.n the number of smoothing iterations. See Details
#' @param red numeric, between 0-1 the proportion of contour coordinates to sample. Will speed up midline estimations. If 'NULL', the full contour in \code{out} will be used See Details.
#' @param dens integer,  the factor by which to densify the contour. Will slow down midline estimations but may increase the likelihood of finding a pointed tip. If 'NULL', the original contour in \code{out} will be used See Details.
#' 
#' @usage free.ml.ang(out = NULL, smooth.n = NULL, dens = NULL, red = NULL)
#' 
#' @return A list with the following components:
#'
#' \code{ml}: a data table consisting of the x,y coordinates of the midline and their index value
#' 
#' \code{cont.sm}: a data table containing the smoothed outline according to \code{coo_smooth}. If 'smooth.n=0', this will be identical to the input 'out'.
#' 
#' \code{cont.sides}: a data table containing the smoothed outline factorized by side 'a' or 'b' and index.
#' 
#' 
#' @details
#' 
#' \code{free.ml.ang} estimates the midline by first creating a convex hull of the contour and then finding candidate tips that are farthest from one another with \code{\link{coo_truss}}. This and \code{free.ml.hull} therefore assume the contour is elongate. The candidate tips are then refined by finding coordinates whose index are within the 5% quantiles of the extrema. From the additional coordinates, the final candidate tips are defined as those that form the greatest angle between adjacent coordinates. The function then **bisects** the contour across this axis determined by the tips using \code{\link{coo_slice}}, giving it two sides with coordinates of equal length. The midline is calculated as the midpoints defined between all pairs of coordinates with the same index value.
#' 
#'\code{free.ml.hull} creates a convex hull of the contour and then finds the two coordinates in the hull that are farthest from one another. As with \code{free.ml.ang}, this function then bisects the contour across an axis determined by the tips and calculates the midline as midpoints defined between all pairs of coordinates with the same index value
#' 
#' 
#' \code{free.ml.del} estimates the midline by first finding candidate tips of the contour. Midline coordinates are retrieved through Delaunay triangulation and Voronoi tesselation. Triangulated points that lie within the contour are used to build a distance-weighted minimum spanning tree that expands until the the path intersects with the contour. The extrema of the midline path are then indexed as the first and last according to their distance from the candidate tips.
#' 
#' The use of \code{free.ml.ang} and \code{free.ml.hull} may be more appropriate for complicated outlines (i.e., those with appendages). The use of \code{free.ml.del} produces better results when contour regions overlap (i.e, kinks or snakes back on itself), but produces less precise midlines for complicated contours and is slower for high resolution outlines. Reducing the resolution with 'red' may hasten the speed. 
#' 
#' 'smooth.n' is passed to the 'n' parameter of \code{\link{coo_smooth}}, which smooths coordinates using a simple moving average. Users should be careful not to oversmooth. If the input contour has few points (say just a 100 or so extracted from \code{kin} functions run on low resolution images), much detail will be lost. In general, 'smooth.n' should be <5. 
#' 
#' If 'red' is specified, the resolution of the contour is reduced by invoking \code{\link{coo_interpolate}} with 'n' of this function equal to 'red' times the number of coordinates in the original contour.  
#' 
#' For input contours with few points (say >100), users should consider densifying with 'dens'. For example, an input contour with 100 coordinates would be transformed to one of 200 with 'dens=1'.
#' 
#' A humble thanks to Barry Rowlingson of Lancaster University for his \href{https://stackoverflow.com/questions/9595117/identify-a-linear-feature-on-a-raster-map-and-return-a-linear-shape-object-using/9643004#9643004}{answer to a post on Stack Overflow} that inspired the logic behind \code{free.ml.del}.
#' 
#' @note These functions do not make decisions about position, i.e, output values of 'n', although ordered along the long axis of the contour, may be differently arranged given the rotation of the contour.
#' 
#' @name free.ml.ang
#' @export
#' @importFrom stats quantile fitted lm predict coef  
#' @importFrom deldir deldir
#' @importFrom sp Polygon Polygons SpatialPoints SpatialPolygons over
#' @importFrom rgeos gDistance
#' @importFrom igraph get.diameter graph.adjacency minimum.spanning.tree vcount
#' @importFrom sf as_Spatial
#' @import methods
#' 
#' @seealso \code{\link{coo_smooth}}, \code{\link{coo_slice}},\code{\link{coo_slide}}, \code{\link{coo_truss}}, \code{\link{kin.free}}, \code{\link{deldir}}, \code{\link{kin.free}}
#' 
#' @examples
#' # a lateral midline, but a midline nonetheless
#' require(Momocs)
#' o <- Momocs::nsfishes$coo[[136]]
#' colnames(o) <- c("x","y")
#' plot(o)
#' fml <- free.ml.hull(o)
#' points(fml$ml$x,fml$ml$y,col="red")
#' #note the difference
#' fml2 <- free.ml.ang(o) 
#' points(fml2$ml$x,fml2$ml$y,col="blue")
#' 
#' # with free.ml.del, a poorer choice
#' fml3 <- free.ml.del(o)
#' points(fml3$ml$x,fml3$ml$y,col="green")
#' 
#' #free.ml.del with the contour of a swimming worm, C. elegans
#' #from SI Video 3 of Pierce-Shimomura, et al., 2008.
#' # Genetic analysis of crawling and swimming locomotory patterns in C. elegans. 
#' #Proceedings of the National Academy of Sciences, 105(52), pp.20982-20987.
#' 
#' cel <-system.file("extdata", "celegans.csv", package = "trackter")
#' out <- as.matrix(read.csv(cel))
#' 
#' #smooth it
#' out.sm <- Momocs::coo_smooth(out,n=1)
#' 
#' fml.del <- free.ml.del(out.sm)
#' points(fml.del$ml$x,fml.del$ml$y,col="red")
#' 
#' 
free.ml.ang <- function(out = NULL,smooth.n=NULL,dens=NULL,red=NULL) {
  
  hull <- bl <- bl2 <- NULL
  if(!"matrix" %in% class(out)) stop("'out' must be a matrix")
  n <- side <- x <- y <- tip <- n2 <- ang <- is.tip <- dist <- dist2 <- is.tip2 <- tip1.dist <- tip2.dist <- is.tip1 <- ang.rm<- NULL
  
  densify <- function(xy,n=5){
    ## densify a 2-col matrix
    cbind(Dens(xy[,1],n=n),Dens(xy[,2],n=n))
  }
  
  Dens <- function(x,n=5){
    ## densify a vector
    out = rep(NA,1+(length(x)-1)*(n+1))
    ss = seq(1,length(out),by=(n+1))
    out[ss]=x
    for(s in 1:(length(x)-1)){
      out[(1+ss[s]):(ss[s+1]-1)]=seq(x[s],x[s+1],len=(n+2))[-c(1,n+2)]
    }
    out
  }
  

  if(!is.null(dens) & is.null(red)) out <- densify(out,n=dens)
  
  if(!is.null(dens)&!is.null(red)) stop("both 'red' and 'dens' are not NULL. Enter values for only one argument" )
  
  if (!is.null(red) & is.null(dens) ){
    if(!is.numeric(red)) stop("'red' must be numeric and 0-1")
    if (red<0 | red>1 )
      stop("'red' must be numeric and 0-1")
  }
  
  if(!is.null(red)) red.n <- round(red*(nrow(out)),0)
  
  #close coo
  coo <- Momocs::coo_close(out)
  
  #reduce, smooth
  if(!is.null(red)) coo <- Momocs::coo_interpolate(coo,n=red.n)
  
  if(!is.null(smooth.n)) if( smooth.n>0) coo <- Momocs::coo_smooth(coo,smooth.n)
  
  colnames(coo) <- c("x","y")

  
  #some outlines have duplicated points, nicht gut
  coo <- coo[!duplicated(coo),]
  
  coo1 <- data.table(coo,n=1:nrow(coo))
  # qplot(d=coo1,x,y,col=n)
  
  #limit to hull points, hastens truss calculation
  coo.h <- data.table(Momocs::coo_chull(coo))[,hull:=TRUE]
  coo.h <- coo1[coo.h,on=c("x","y")]
  tr.h <-  Momocs::coo_truss(as.matrix(coo.h[,list(x,y)]))
  tip.h <- names(tr.h[which.max(tr.h)])
  tips.h <-
    c(as.numeric(gsub("(\\d+)-(\\d+)", "\\1", tip.h)), as.numeric(gsub("(\\d+)-(\\d+)", "\\2", tip.h)))
  
  tips <- coo1[n%in%coo.h[tips.h,]$n,]
  tip1 <- tips[n==min(n),]
  tip2 <- tips[n==max(n),]
  #qplot(d=coo1,x,y,col=n)+geom_point(data=coo1[n%in%coo.h[tips.h,]$n,],aes(x,y),col="red")
  
  #slide to one point from tips to prevent bad slicing
  coo2 <-  Momocs::coo_slide(coo, max(tips$n))
  coo2 <- data.table(coo2)[,n:=1:.N]

  
  #find those tips in new contour
  coo2[,tip1.dist:= dist.2d(x,tip1$x,y,tip1$y)]
  coo2[,tip2.dist:= dist.2d(x,tip2$x,y,tip2$y)]
  coo2[,is.tip1:=tip1.dist==min(tip1.dist)]
  coo2[,is.tip2:=tip2.dist==min(tip2.dist)]
  
  #qplot(d=coo2,x,y,col=n)+geom_point(data=coo2[n%in%tips2,],aes(x,y),col="red")
  
  #pull tips from slidden coo (coo2) and reslice on matrix of coo2
  tips2 <- coo2[is.tip1==TRUE|is.tip2==TRUE]$n
  
  ##dont need
  coo.sl <-  Momocs::coo_slice(as.matrix(coo2[,list(x,y)]), ids = tips2)
  n.pts2 <- sapply(coo.sl,nrow)

  #resample so equal n on each side
  coo.sl <-
    lapply(coo.sl, function(x)
      Momocs::coo_sample(x,n=min(n.pts2)))
  
  #coo.sl becomes coo.sides with side factor
  coo.sides <- rbind(
    data.table(coo.sl[[1]],n=1:min(n.pts2),side="a"),
    data.table(coo.sl[[2]],n=min(n.pts2):1,side="b")
  )
  
  
  colnames(coo.sides)[1:2] <- c("x","y")
  
  ##end dont need
  
  setkeyv(coo.sides,c("n"))
  #qplot(d=coo.sides,x,y,col=n)


  #qplot(d=coo2,x,y,col=n)+geom_point(data=coo2[n%in%tips2],aes(x,y),col="red")
  
  #find 10% of ends
  ends.q <- round(quantile(min(tips2):max(tips2),probs=c(.1)))
  ends1.n <- c((max(coo2$n)-ends.q):max(coo2$n),min(tips2):ends.q)
  ends2.n <- c((max(tips2)-ends.q):(max(tips2)+ends.q))
  
  #qplot(d=coo2,x,y,col=n)+geom_point(data=coo2[n%in%ends1.n],aes(x,y),col="red")
  
  #pull ends from 
  ends1 <- coo2[n%in%ends1.n]
  ends1[,side:=ifelse(n<=ends.q,"a","b")]
  ends1[n<=ends.q,n2:=.N:1]
  ends1[n>ends.q,n2:=(.N:1)+ends.q]
  
  ends2 <- coo2[n%in%ends2.n]
  ends2[,side:=ifelse(n<=max(tips2),"a","b")]
  ends2[,n2:=1:.N]
  
  
  #qplot()+geom_point(data=ends2,aes(x,y,col=side))
  
  ang.3pts <- function(pt1,pt2,pt3){
    ang <-  atan2(pt3[2] - pt1[2], pt3[1] -pt1[1]) -atan2(pt2[2] - pt1[2], pt2[1] - pt1[1])
    return(ang)
  }
  
  setkeyv(ends1,"n2")
  setkeyv(ends2,"n2")
  
  #compute angles
  ends1.ang <- list()

  for(a in ends1$n2){
    pt1= ends1[n2==a-1,list(x,y)]
    pt2= ends1[n2==a,list(x,y)]
    pt3= ends1[n2==a+1,list(x,y)]
    ang <- ang.3pts(c(unlist(pt1)),c(unlist(pt2)),c(unlist(pt3)))
    ends1.ang[[paste0(a)]] <- data.table(n2=a,ang=abs(deg(ang)))
  }
  
  ends2.ang <- list()
  for(a in ends2$n2){
    pt1= ends2[n2==a-1,list(x,y)]
    pt2= ends2[n2==a,list(x,y)]
    pt3= ends2[n2==a+1,list(x,y)]
    ang <- ang.3pts(c(unlist(pt1)),c(unlist(pt2)),c(unlist(pt3)))
    ends2.ang[[paste0(a)]] <- data.table(n2=a,ang=abs(deg(ang)))
  }
  
  
  ends1.ang <- do.call(rbind,ends1.ang)
  ends2.ang <- do.call(rbind,ends2.ang)
  
  #rolling mean of angle
  ends1.ang[,ang.rm:= frollmean(ends1.ang$ang,3)]
  ends2.ang[,ang.rm:= frollmean(ends2.ang$ang,3)]
  
  ends1.newN <- ends1.ang[which.max(ang)]$n2
  ends2.newN <- ends2.ang[which.max(ang)]$n2
  
  ends1.newN2 <- ends1.ang[which.max(ang.rm)]$n2
  ends2.newN2 <- ends2.ang[which.max(ang.rm)]$n2
  
  ends1[,is.tip:=n2==ends1.newN2]
  ends2[,is.tip:=n2==ends2.newN2]
  
  
  #qplot()+geom_point(d=coo.dist,aes(x,y))+geom_point(d=ends1[n2%in% ends1.newN2],aes(x,y),col="red")+geom_point(d=ends1[n2%in% ends1.newN],aes(x,y),col="blue")
  #qplot()+geom_point(d=coo.dist,aes(x,y))+geom_point(d=ends2[n2%in% ends2.newN2],aes(x,y),col="red")+geom_point(d=ends2[n2%in% ends2.newN],aes(x,y),col="blue")
  
  #redefine tips
  tip1 <- ends1[is.tip==TRUE]
  tip2 <- ends2[is.tip==TRUE]
  
  #new tips based on angles
  tips3 <- c(tip1$n,tip2$n)
  
  #slide to one point from tips to prevent bad slicing
  coo3 <-  Momocs::coo_slide(as.matrix(coo2[,list(x,y)]), max(tips3))
  coo3 <- data.table(coo3)[,n:=1:.N]
  
  #find those tips in new contour
  coo3[,tip1.dist:= dist.2d(x,tip1$x,y,tip1$y)]
  coo3[,tip2.dist:= dist.2d(x,tip2$x,y,tip2$y)]
  coo3[,is.tip1:=tip1.dist==min(tip1.dist)]
  coo3[,is.tip2:=tip2.dist==min(tip2.dist)]
  
  tips4 <- coo3[is.tip1==TRUE|is.tip2==TRUE]
  
  #qplot(data=coo3,x,y,col= n)+geom_point(d=coo3[n%in%tips4],aes(x,y),col="red")
  
  #slice a final time
  coo.sl<-  Momocs::coo_slice(as.matrix(coo3[,list(x,y)]), ids = tips4$n)
  
   n.pts4 <- sapply(coo.sl,nrow)

  #make final slice coo have about as many points as original coo/out
  # coo.sl <- lapply(coo.sl, function(x)
  #   smoothr::smooth_spline(x,n = max(n.pts4)))
  
  # n.pts4.2 <- sapply(coo.sl,nrow)
  
  #new table with factors sides of equal length
  coo.sides <- rbind(
    data.table(coo.sl[[1]],n=1:nrow(coo.sl[[1]]),side="a"),
    data.table(coo.sl[[2]],n=nrow(coo.sl[[2]]):1,side="b")
  )
  
  colnames(coo.sides)[1:2] <- c("x","y")
  
  #qplot(d=coo.sides[n<10],x,y,col=as.factor(n))
  
  setkeyv(coo.sides,c("n"))
  
  coo.ml <- coo.sides[, list(x = sum(x) / 2, y = sum(y) / 2), by = list(n)]
  
  coo.sides <- rbind(
    data.table(coo.sl[[1]],n=1:nrow(coo.sl[[1]]),side="a"),
    data.table(coo.sl[[2]],n=nrow(coo.sl[[2]]):1,side="b")
  )
  
  rot <- coo.ml[,{pa <- point.ang.orig(c(x,y),c(coo.ml$x[1],coo.ml$y[1]),pi/2);list(x=pa[1],y=pa[2])},by=n][c(range(n)),]
  rot.lm <- lm(y~x,rot)
  #rot$pred <- predict(rot.lm)
  
  # qplot(d=fml$cont.sides,x,y,col=dist)+geom_point(d=rot,aes(x,),col="red")
  
  
  coo.sides[,dist:=dist.2d.line(x,y,coef(rot.lm)[2],coef(rot.lm)[1]),by=list(n,side)][,bl:=dist/max(dist),by=side]
  
  coo.sides2 <- copy(coo.sides[,bl2:=round(bl,2)])[,list(x=mean(x),y=mean(y)),by=list(bl2,side)]
  
  setkeyv(coo.sides2,c("side","bl2"))
  coo.sides2[,n:=1:.N,by=(side)]
  
 coo.sides3 <-  coo.sides2[,{s <- smoothr::smooth_spline(as.matrix(data.frame(x,y)),n = max(n.pts4));
    list(
      x=s[,1],
      y=s[,2]
    )
  },
  by=side
  ]
 
 coo.sides3[,n:=1:.N,by=(side)]

  
  coo.ml2 <- coo.sides3[, list(x = sum(x) / 2, y = sum(y) / 2), by = list(n)]
  
  
  #qplot(d=coo.sides,x,y,col=side)+geom_point(d=coo.ml2,aes(x,y),col="red")
  
  return(list(ml = coo.ml2,cont.sm=coo3[,list(n,x,y)],cont.sides=coo.sides))
}

NULL

#' @rdname free.ml.ang
#' @export
#' 

#out=as.matrix(kin$cont[frame==163,list(x,y)])
free.ml.hull <- function(out = NULL,smooth.n=NULL,dens=NULL,red=NULL) {
  hull <- bl <- bl2 <- NULL
  
  if(!"matrix" %in% class(out)) stop("'out' must be a matrix")
  n <- side <- x <- y <- tip <- n2 <- ang <- is.tip <- dist <- dist2 <- is.tip2 <- tip1.dist <- tip2.dist <- is.tip1 <- ang.rm<- NULL
  
  densify <- function(xy,n=5){
    ## densify a 2-col matrix
    cbind(Dens(xy[,1],n=n),Dens(xy[,2],n=n))
  }
  
  Dens <- function(x,n=5){
    ## densify a vector
    out = rep(NA,1+(length(x)-1)*(n+1))
    ss = seq(1,length(out),by=(n+1))
    out[ss]=x
    for(s in 1:(length(x)-1)){
      out[(1+ss[s]):(ss[s+1]-1)]=seq(x[s],x[s+1],len=(n+2))[-c(1,n+2)]
    }
    out
  }
  
  
  if(!is.null(dens) & is.null(red)) out <- densify(out,n=dens)
  
  if(!is.null(dens)&!is.null(red)) stop("both 'red' and 'dens' are not NULL. Enter values for only one argument" )
  
  if (!is.null(red) & is.null(dens) ){
    if(!is.numeric(red)) stop("'red' must be numeric and 0-1")
    if (red<0 | red>1 )
      stop("'red' must be numeric and 0-1")
  }
  
  if(!is.null(red)) red.n <- round(red*(nrow(out)),0)
  
  #close coo
  coo <- Momocs::coo_close(out)
  
  #reduce, smooth
  if(!is.null(red)) coo <- Momocs::coo_interpolate(coo,n=red.n)

  if(!is.null(smooth.n)) if( smooth.n>0) coo <- Momocs::coo_smooth(coo,smooth.n)
  
  colnames(coo) <- c("x","y")
  
  #some outlines have duplicated points, nicht gut
  coo <- coo[!duplicated(coo),]
  
  coo1 <- data.table(coo,n=1:nrow(coo))
  # qplot(d=coo1,x,y,col=n)
  
  #limit to hull points, hastens truss calculation
  coo.h <- data.table(Momocs::coo_chull(coo))[,hull:=TRUE]
  coo.h <- coo1[coo.h,on=c("x","y")]
  tr.h <-  Momocs::coo_truss(as.matrix(coo.h[,list(x,y)]))
  tip.h <- names(tr.h[which.max(tr.h)])
  tips.h <-
    c(as.numeric(gsub("(\\d+)-(\\d+)", "\\1", tip.h)), as.numeric(gsub("(\\d+)-(\\d+)", "\\2", tip.h)))
  
  tips <- coo1[n%in%coo.h[tips.h,]$n,]
  tip1 <- tips[n==min(n),]
  tip2 <- tips[n==max(n),]
  #qplot(d=coo1,x,y,col=n)+geom_point(data=coo1[n%in%coo.h[tips.h,]$n,],aes(x,y),col="red")
  
  #slide to one point from tips to prevent bad slicing
  coo2 <-  Momocs::coo_slide(coo, max(tips$n))
  coo2 <- data.table(coo2)[,n:=1:.N]
  
  
  #find those tips in new contour
  coo2[,tip1.dist:= dist.2d(x,tip1$x,y,tip1$y)]
  coo2[,tip2.dist:= dist.2d(x,tip2$x,y,tip2$y)]
  coo2[,is.tip1:=tip1.dist==min(tip1.dist)]
  coo2[,is.tip2:=tip2.dist==min(tip2.dist)]
  
  #qplot(d=coo2,x,y,col=n)+geom_point(data=coo2[n%in%tips2,],aes(x,y),col="red")
  
  #pull tips from slidden coo (coo2) and reslice on matrix of coo2
  tips2 <- coo2[is.tip1==TRUE|is.tip2==TRUE]$n
  
  coo.sl <-  Momocs::coo_slice(as.matrix(coo2[,list(x,y)]), ids = tips2)
  n.pts2 <- sapply(coo.sl,nrow)
  
  #resample so equal n on each side
  coo.sl <-
    lapply(coo.sl, function(x)
      Momocs::coo_sample(x,n=min(n.pts2)))
  
  #coo.sl becomes coo.sides with side factor
  coo.sides <- rbind(
    data.table(coo.sl[[1]],n=1:nrow(coo.sl[[1]]),side="a"),
    data.table(coo.sl[[2]],n=nrow(coo.sl[[2]]):1,side="b")
  )
  
  colnames(coo.sides)[1:2] <- c("x","y")
  
  
  coo.ml <- coo.sides[, list(x = sum(x) / 2, y = sum(y) / 2), by = list(n)]
  
  #qplot(d=coo.sides,x,y,col=n)+geom_point(d=coo.ml,aes(x,y),col="red")
  
  rot <- coo.ml[,{pa <- point.ang.orig(c(x,y),c(coo.ml$x[1],coo.ml$y[1]),pi/2);list(x=pa[1],y=pa[2])},by=n][c(range(n)),]
  rot.lm <- lm(y~x,rot)
  

  coo.sides[,dist:=dist.2d.line(x,y,coef(rot.lm)[2],coef(rot.lm)[1]),by=list(n,side)][,bl:=dist/max(dist),by=side]
  
  coo.sides2 <- copy(coo.sides[,bl2:=round(bl,2)])[,list(x=mean(x),y=mean(y)),by=list(bl2,side)]
  
  setkeyv(coo.sides2,c("side","bl2"))
  coo.sides2[,n:=1:.N,by=side]
  
  coo.sides3<-  coo.sides2[,{s <- smoothr::smooth_spline(as.matrix(data.frame(x,y)),n = max(n.pts2));
  list(
    x=s[,1],
    y=s[,2]
  )
  },
  by=side
  ]
  
  coo.sides3[,n:=1:.N,by=side]
  
  
  coo.ml2 <- coo.sides3[, list(x = sum(x) / 2, y = sum(y) / 2), by = list(n)]
  
  #qplot(d=coo.sides,x,y,col=n)+geom_point(d=coo.ml2,aes(x,y),col="red")
  
  return(list(ml = coo.ml2,cont.sm=coo2[,list(n,x,y)],cont.sides=coo.sides))
}


NULL


#' @rdname free.ml.ang
#' @export
#' 
#out <- as.matrix(kin.del$cont[frame==14,list(x,y)])
free.ml.del <- function(out = NULL,smooth.n=NULL,red=NULL,dens=NULL) {

  hull <- bl <- bl2 <- NULL
  
  if(!"matrix" %in% class(out)) stop("'out' must be a matrix")
  n <- side <- x <- y <- tip <- n2 <- ang <- is.tip <- dist <- dist2 <- is.tip2 <- tip1.dist <- tip2.dist <- is.tip1 <- ang.rm <- NULL
  densify <- function(xy,n=5){
    ## densify a 2-col matrix
    cbind(Dens(xy[,1],n=n),Dens(xy[,2],n=n))
  }
  
  Dens <- function(x,n=5){
    ## densify a vector
    out = rep(NA,1+(length(x)-1)*(n+1))
    ss = seq(1,length(out),by=(n+1))
    out[ss]=x
    for(s in 1:(length(x)-1)){
      out[(1+ss[s]):(ss[s+1]-1)]=seq(x[s],x[s+1],len=(n+2))[-c(1,n+2)]
    }
    out
  }
  
  
  if(!is.null(dens) & is.null(red)) out <- densify(out,n=dens)
  
  if(!is.null(dens)&!is.null(red)) stop("both 'red' and 'dens' are not NULL. Enter values for only one argument" )
  
  if (!is.null(red) & is.null(dens) ){
    if(!is.numeric(red)) stop("'red' must be numeric and 0-1")
    if (red<0 | red>1 )
      stop("'red' must be numeric and 0-1")
  }
  
  if(!is.null(red)) red.n <- round(red*(nrow(out)),0)
  
  coo <- Momocs::coo_close(out)
  if(!is.null(red)) coo <- Momocs::coo_interpolate(coo,n=red.n)
  if(!is.null(smooth.n)) if( smooth.n>0) coo <- Momocs::coo_smooth(coo,smooth.n)
  colnames(coo) <- c("x","y")

  
  tr <-  Momocs::coo_truss(coo)
  tip.n <- names(tr[which.max(tr)])
  
  tips <-
    c(as.numeric(gsub("(\\d+)-(\\d+)", "\\1", tip.n)), as.numeric(gsub("(\\d+)-(\\d+)", "\\2", tip.n)))
  
  coo <-  Momocs::coo_slide(coo, tips[2])
  
  
  
  ### compute triangulation
  d <- deldir(coo[,1],coo[,2],)
  
  #plot(d)
  ### find midpoints of triangle sides
  ml <- cbind((d$delsgs[,'x1']+d$delsgs[,'x2'])/2,
              (d$delsgs[,'y1']+d$delsgs[,'y2'])/2)
  
  ### get points that are inside the polygon 
  sr <- SpatialPolygons(list(Polygons(list(Polygon(coo)),ID=1)))
  ins <-  over(SpatialPoints(ml),sr)
  
  ### select the points
  ml.pts <-  ml[!is.na(ins),]
  
  G <-  graph.adjacency(as.matrix(dist(ml.pts)),weighted=TRUE,mode="upper")
  Tr <-  minimum.spanning.tree(G,weighted=TRUE)
  vct <- vcount(Tr)
  ### get a diameter
  path <-  get.diameter(Tr)
  path.l <- length(path)
  
  #plot(ml.pts)
  div <- 30
  
  while(path.l!=vct){
    dPoly <- gDistance(as(sr,"SpatialLines"),SpatialPoints(ml.pts),byid=TRUE)
    ml.pts <-  ml.pts[dPoly > max(dPoly/div),]
    
    
    ### now build a minimum spanning tree weighted on the distance
    G <-  graph.adjacency(as.matrix(dist(ml.pts)),weighted=TRUE,mode="upper")
    Tr <-  minimum.spanning.tree(G,weighted=TRUE)
    vct <- vcount(Tr)
    ### get a diameter
    path <-  get.diameter(Tr)
    path.l <- length(path)
    if(div<1){
      stop("midline path intersects outline---try a denser outline")
    }
    
    
    ### path should be the sequence of points in order
    #list(pts=pts[path+1],tree=T)
    
    
    
    div <- div-1
    
  }
  
  
  #reorder path
  ml2 <- data.table(ml.pts[as.integer(path),])
  colnames(ml2) <- c("x","y")
  ml2[,n:=1:nrow(ml.pts)]
  coo.dt <- data.table(coo)[,n:=1:.N]
  #qplot(d=coo.dt,x,y,col=n)+geom_point(d=ml2,aes(x,y),col="red")
  
  
  
  ml.end1 <- ml2[n==1]
  ml.end2 <- ml2[n==max(n)]
  end1 <- coo.dt[,d:=dist.2d(x,ml.end1$x,y,ml.end1$y),by=n][which.min(d),]
  end2 <- coo.dt[,d:=dist.2d(x,ml.end2$x,y,ml.end2$y),by=n][which.min(d),]
  
  tips2 <- c(end1$n,end2$n)
  
  coo.slide <-  Momocs::coo_slide(coo,  min(tips2))
  
  coo.slide.dt <-data.table(coo.slide)[,n:=1:.N]
  #qplot(d=coo.slide.dt,x,y,col=n)+geom_point(d=ml.end1,aes(x,y),col="red")+geom_point(d=ml.end2,aes(x,y),col="red")
  
  coo.sl <-  Momocs::coo_slice(coo.slide, ids = tips2-min(tips2))
  
  n.pts <- sapply(coo.sl,nrow)
  
  coo.sl2 <-
    lapply(coo.sl, function(x)
      Momocs::coo_sample(x,n=min(n.pts)))
  
  
  n.pts2 <- sapply(coo.sl2,nrow)
  
  coo.sl2 <- lapply(coo.sl2, function(x)
    smoothr::smooth_spline(x,n =max(n.pts2)))
  
  coo.side <- data.table(do.call(rbind,coo.sl2))[,side:=c(rep("a",max(n.pts2)),rep("b",max(n.pts2)))]
  coo.side[,n:=1:.N,by=side][side=="b",n:=rev(n)]
  colnames(coo.side)[1:2] <- c("x","y")
  setkeyv(coo.side,c("side","n"))
  
  
  
  return(list(ml=ml2[,list(x,y,n)],cont.sm=coo.dt[,list(n,x,y)],cont.sides=coo.side))
  
}

