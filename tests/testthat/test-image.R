
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
  
  
  expect_error(crop.img(img=z,ul=c(5,30),br=c(200,200),out.dir=od),"file in file path doesn't appear to be an image" )
  
  expect_error(crop.img(img=y,ul=c(5,30),out.dir=od),"both 'ul' and 'br' must be specified or 'locate=TRUE'")
  
  expect_error(crop.img(img=y,ul=c(5,30),br=c(200,200),out.dir=od,type="foo"),"invalid" )
  
  dimA <- dim(EBImage::readImage(y))
  dimB <- dim(EBImage::readImage(list.files(od,full.names = T)[1]))
  dimC <- dim(EBImage::readImage(list.files(od,full.names = T)[2]))
  
  expect_true(dimA[1]>dimB[1]&& dimA[1]>dimC[1] &&  dimB[1]<dimC[1])

})
