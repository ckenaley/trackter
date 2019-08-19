
context("kin. fin.kin functions")

test_that("kin.simple works fine", {
  y <- EBImage::readImage(system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter"))
  t <-tempdir()
  dir.create(paste0(t,"/images"))
  EBImage::writeImage(y,paste0(t,"/images/sunfish001.jpg"),type = "jpeg")
  invisible(capture.output( kin.y <- kin.simple(image.dir = paste0(t,"/images"),save = T,out.dir = t)))
  
  tp <- paste0(tempdir(),"/processed_images")
  expect_length(list.files(tp),1)
  
  unlink(tp,recursive = T)
  
  invisible(capture.output( kin.y <- kin.simple(image.dir = paste0(t,"/images"),save = F)))
  
  expect_false(dir.exists(tp))
  
  expect_is(kin.y,"list")
  expect_named(kin.y,c("kin.dat", "midline","cont","all.classes","dim"))
  expect_true(kin.y$kin.dat$size>0)
  expect_type(kin.y$midline$roi,type = "character")
  expect_type(kin.y$cont$x,type = "integer")
  expect_true(kin.y$all.classes$size>0)
  
  
  expect_error(invisible(capture.output( kin.simple(image.dir = "foo",save=F))))
  expect_error(invisible(capture.output( kin.simple(image.dir = paste0(t,"/images"), frames=2,save=F))),"out of range")
  expect_error(invisible(capture.output( kin.simple(image.dir =paste0(t,"/images") , thr="foo",save=F))),"must be set to")
  expect_error(invisible(capture.output( kin.simple(image.dir =paste0(t,"/images") ,smoothing="foo",save=F))),"must =")
  unlink(paste0(t,"/images"),recursive = T)
  unlink(tp,recursive = T)
})

test_that("kin.search works fine", {
  
  y <- EBImage::readImage(system.file("extdata/img", "sunfish_BCF.jpg", package = "trackter"))
  t <-tempdir()
  dir.create(paste0(t,"/images"))
  EBImage::writeImage(y,paste0(t,"/images/sunfish001.jpg"),type = "jpeg")
  invisible(capture.output( kin.y <- kin.search(image.dir = paste0(t,"/images"),save = T,out.dir = t)))
  
  tp <- paste0(tempdir(),"/processed_images")
  expect_length(list.files(tp),1)
  
  unlink(tp,recursive = T)
  
  invisible(capture.output( kin.y <- kin.search(image.dir = paste0(t,"/images"),save = F)))
  
  expect_false(dir.exists(tp))
  
  expect_is(kin.y,"list")
  expect_named(kin.y,c("kin.dat", "midline","cont","all.classes","dim"))
  expect_true(kin.y$kin.dat$size>0)
  expect_type(kin.y$midline$roi,type = "character")
  expect_type(kin.y$cont$x,type = "integer")
  expect_true(kin.y$all.classes$size>0)
  
  
  expect_error(invisible(capture.output( kin.search(image.dir = "foo",save=F))))
  expect_error(invisible(capture.output( kin.search(image.dir = paste0(t,"/images"), frames=2,save=F))),"out of range")
  expect_error(invisible(capture.output( kin.search(image.dir =paste0(t,"/images") , thr="foo",save=F))),"must be set to")
  expect_error(invisible(capture.output( kin.search(image.dir =paste0(t,"/images") ,smoothing="foo",save=F))),"must =")
  
  unlink(paste0(t,"/images"),recursive = T)
  unlink(tp,recursive = T)
  
})


test_that("fin.kin works fine", {

  cont <- read.csv(system.file("extdata", "cont.csv", package = "trackter"))[,3:4]
  fin.y <- fin.kin(cont,fin.pos  = c(0.25,0.5))
  
  expect_is(fin.y,"list")
  expect_named(fin.y,c("body","fin","fin.pts","comp","midline","bl","amp"))
  expect_type(fin.y$body$y,type = "double")
  expect_type(fin.y$comp$y,type = "double")
  expect_type(fin.y$fin.pts$y,type = "double")
  expect_type(fin.y$fin$y,type = "double")
  expect_type(fin.y$bl[1],type = "double")
  expect_type(fin.y$amp$amp[1],type = "double")
  
  expect_error(fin.kin(kin.y$cont))
  expect_error(fin.kin(data.frame(x=kin.y$cont$x,y=kin.y$cont$y),fin.pos=0.1))
  expect_error(fin.kin(data.frame(x=kin.y$cont$x,y=kin.y$cont$y),fin.pos=NULL))
})


