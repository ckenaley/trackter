
context("image functions")

test_that("thr.check works fine", {
  
  y <-system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")
  
  tc <- thr.check(y)
  expect_is(tc,"numeric")
  expect_true(length(tc)==1)
  expect_message(invisible(capture.output(thr.check(y,min=0.2,max=0.3))),"Otsu value is outside defined threshold range")
  expect_error(invisible(capture.output(thr.check(y,min=0.4,max=0.3))),"'min' must be < 'max'")
  expect_error(invisible(capture.output(thr.check(y,min=0.4,max=0.4))),"'min' must be value different from 'max'")
  
  expect_error(invisible(capture.output(thr.check(y,min=NULL,max=0.3))),"both 'min' and 'max' must have value=NULL or numeric 0-1")
  
  expect_true(length(thr.check(y,otsu=FALSE))==0)
  
  
  png("rplot.png") 
  thr.check(y,otsu=FALSE)
  dev.off() 
  
  expect_true(file.exists("rplot.png"))
  
  expect_true(file.size("rplot.png")>100)
  
  unlink("rplot.png")
  
  
})


test_that("crop.img works fine", {
  
  y <-system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")
  
  od <- paste0(tempdir(),"/cropimg")
  dir.create(od)

  crop.img(img=y,ul=c(5,30),br=c(100,200),out.dir=od)

  expect_true(list.files(od)=="sunfish_BCF.jpg")
  
  crop.img(img=y,ul=c(5,30),br=c(200,200),out.dir=od,type="png")
  
  z <- system.file("extdata/vid", "sunfish_BCF.avi", package = "trackter")
  
  expect_true("sunfish_BCF.png" %in% list.files(od))
  
  
  
  expect_error(crop.img(img=z,ul=c(5,30),br=c(200,200),out.dir=od),"file in file path doesn't appear to be an image" )
  
  expect_error(crop.img(img=y,ul=c(5,30),out.dir=od),"both 'ul' and 'br' must be specified or 'locate=TRUE'")
  
  expect_error(crop.img(img=y,ul=c(5,30),br=c(200,200),out.dir=od,type="foo"),"invalid" )
  
  dimA <- dim(EBImage::readImage(y))
  dimB <- dim(EBImage::readImage(list.files(od,full.names = T)[1]))
  dimC <- dim(EBImage::readImage(list.files(od,full.names = T)[2]))
  
  expect_true(dimA[1]>dimB[1]&& dimA[1]>dimC[1] &&  dimB[1]<dimC[1])

  unlink(od,recursive = TRUE)
})


test_that("contrast.img works fine", {
  
  y <-system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")
  
  od <- paste0(tempdir(),"/cropimg")
  dir.create(od)
  
  contrast.img(img=y,c=0.5,out.dir=od)
  
  expect_true(list.files(od)=="sunfish_BCF.jpg")
  
  contrast.img(img=y,c=0.5,out.dir=od,type="png")
  
  expect_true("sunfish_BCF.png" %in% list.files(od))
  
  z <- system.file("extdata/vid", "sunfish_BCF.avi", package = "trackter")
  
  
  expect_error(contrast.img(img=z,c=0.5,out.dir=od),"file in file path doesn't appear to be an image" )
  
  expect_error(contrast.img(img=y,c=0.5,out.dir=od,type="foo"),"invalid" )
  
  unlink(od,recursive = TRUE)
  
})


test_that("data.overlay works fine", {
  
f <-system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter")

d <- dim(EBImage::readImage(f))
x <- runif(10,1,d[1])
y <- runif(10,1,d[2])
pts <- cbind(x,y)

recPlot <- function(expr) {
  pdf(NULL)
  on.exit(dev.off())
  dev.control(displaylist="enable")
  expr
  recordPlot()
}


expect_identical(
  recPlot(data.overlay(img=f,over=pts,col="red",type="p")),
  recPlot(data.overlay(img=f,over=pts,col="red",type="p"))
  ) 

expect_error(expect_identical(
  recPlot(data.overlay(img=f,over=pts,col="red",type="p")),
  recPlot(data.overlay(img=f,over=pts,col="blue",type="p"))
),"not identical to") 



  
})

test_that("gg.overlay works fine", {
  
  f <-system.file("extdata", "sunfish_kin.RDS", package = "trackter")
  
  kin <- readRDS(f)

  p <- gg.overlay(kin=kin,under="cont.sm", over="midline")
  
  expect_length(p$layers,2) 
  
  library(transformr)
  gg.overlay(kin=kin,under="cont.sm", over="midline", size=1,animate=TRUE, frames=0:1,col="red",fps=10,save=TRUE,out.dir = tempdir(),filename = "foo.gif")
  
  expect_true("foo.gif" %in% list.files(tempdir()))
  
})

